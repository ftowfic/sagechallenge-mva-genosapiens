"""One bounded case projection of computed submission dependencies.

Portable roots bind the complete context parents and numeric member/singleton
maps. The operation-owned SQLite projection copies only compact pair codes for
random ranked access; it is rebuilt from verified context owners on restore.
No saved projection or serialized validation marker grants source authority.
"""
from array import array
from collections import OrderedDict
from dataclasses import dataclass
import hashlib
from genosapiens.wgs_context_domain import next_context_ordinal, require_complete_domain
import json
import sys

from genosapiens.evidence_serialization import _canonical
from genosapiens.wgs_evidence_store import EvidenceStoreError
from genosapiens.wgs_evidence_replay import EvidenceReplay
from genosapiens.wgs_locator_tape import LocatorTape, _decode
from genosapiens.wgs_packed_outcomes import outcome_chunk_values
from genosapiens.wgs_pair_outcomes import CHUNK_PAIRS
from genosapiens.wgs_pair_domain import CompactPairDomainPositions
from genosapiens.wgs_selection_evidence import SelectionOccurrence

MAX_BUILD_MEMBER_BYTES = 64 * 1024**2
MAX_CONTEXT_CACHE_BYTES = 512 * 1024**2
MAX_PAIR_CACHE_BYTES = 64 * 1024**2
MAX_CLASS_CACHE_BYTES = 16 * 1024**2
MAX_RENDER_BYTES = 16 * 1024**2


@dataclass(slots=True)
class _Context:
    rendered: array
    positions: CompactPairDomainPositions | None
    offset: int
    size: int


class _SelectionContextRows:
    """Live context projection before case-wide rendered IDs are allocated."""
    def __init__(self, scores, representations, eligibility):
        self.universe,self.members = scores.universe,scores.universe.members
        self.representations,self.eligibility = representations,eligibility
        self.keys,self.size = {},0

    def members_rows(self):
        positions = None if self.universe.pairs is None else self.universe.pairs.domain.positions()
        for number,member in enumerate(self.members.iter_members(allele_order=True),1):
            allele = member.allele
            self.size += sys.getsizeof(allele.allele_id)+112
            if self.size > MAX_BUILD_MEMBER_BYTES:
                raise EvidenceStoreError("selection member keys exceed capacity",capacity=True)
            self.keys[allele.allele_id] = number
            representation = self.representations.get(allele.allele_id)
            if representation is None or representation.native_identity_id != allele.native_identity_id:
                raise EvidenceStoreError("selection representation member changed")
            if positions is not None and positions.keys[number] != allele.allele_id:
                raise EvidenceStoreError("selection pair member order changed")
            content = None if representation.status != "eligible" else _canonical(representation.trace.rendered)
            if content is not None and len(content) > MAX_RENDER_BYTES:
                raise EvidenceStoreError("rendered submission identity exceeds capacity",capacity=True)
            yield content,number-1 if positions is None else positions.events[number],\
                0 if positions is None else int(positions.indices[number] == 2**32-1)
        if len(self.keys) != self.members.member_count:
            raise EvidenceStoreError("selection member coverage changed")

    def singleton_rows(self):
        count = 0
        for number,item in enumerate(self.universe.singletons(),1):
            row = self.eligibility.connection.execute(f"SELECT id,decision FROM {self.eligibility.prefix}_singletons WHERE ordinal=?",(number,)).fetchone()
            if row is None or row[0] != item.hypothesis.hypothesis_id or len(item.hypothesis.allele_refs) > 1:
                raise EvidenceStoreError("selection singleton coverage changed")
            member = self.keys[item.hypothesis.allele_refs[0].allele_id] if item.hypothesis.allele_refs else 0
            count += 1
            yield number,member,row[1]
        if count != self.universe.singleton_count:
            raise EvidenceStoreError("selection singleton count changed")
        self.keys.clear()

    def class_rows(self):
        for row in self.eligibility.catalog():
            yield row["ordinal"],(row["reasons"],row["complete"],row["inputs"]["incomplete"])


def _selection_context_header(ordinal, scores, representations, eligibility):
    from genosapiens.wgs_population_work_filter import work_domain
    universe,members = scores.universe,scores.universe.members
    if eligibility.buffer is not None:
        eligibility.buffer.flush()
    return {"ordinal":ordinal,"context_id":members.context_id,"member_count":members.member_count,
        **work_domain(members),
        "hypothesis_count":universe.hypothesis_count,"singleton_count":universe.singleton_count,
        "pair_count":0 if eligibility.buffer is None else eligibility.buffer.count,"pair_offset":eligibility.offset,
        "has_pair_domain":universe.pairs is not None,
        "scores":{"name":scores.name,"sha256":scores.sha256},
        "members":{"name":members.name,"sha256":members.sha256},
        "hypotheses":{"name":universe.name,"sha256":universe.sha256},
        "representations":{"name":representations.name,"sha256":representations.sha256},
        "eligibility":{"name":eligibility.name,"sha256":eligibility.sha256}}


def _selection_pair_rows(eligibility):
    if eligibility.buffer is not None:
        yield from eligibility.connection.execute(f"SELECT ordinal,count,content FROM {eligibility.buffer.table} ORDER BY ordinal")


class CompleteCaseSelectionInputs:
    def __init__(self, writer, *, name, context_domain=None):
        self.context_domain = context_domain
        self.writer,self.connection,self.name = writer,writer.connection,name
        self.prefix = "selection_inputs_" + hashlib.sha256(name.encode()).hexdigest()
        self.tables,self.closed,self.published = [],False,False
        self.context_count,self.hypothesis_count,self.rendered_count = 0,0,0
        self.context_cache,self.pair_cache,self.class_cache = OrderedDict(),OrderedDict(),OrderedDict()
        self.cache_bytes = {"context":0,"pair":0,"class":0}
        try:
            for suffix,definition in (
                ("contexts","id INTEGER PRIMARY KEY,context_id TEXT UNIQUE,content TEXT"),
                ("rendered","id INTEGER PRIMARY KEY,content TEXT UNIQUE"),
                ("members","context INTEGER,ordinal INTEGER,count INTEGER,content BLOB,PRIMARY KEY(context,ordinal)"),
                ("singletons","context INTEGER,ordinal INTEGER,member INTEGER,code INTEGER,PRIMARY KEY(context,ordinal)"),
                ("pairs","context INTEGER,ordinal INTEGER,count INTEGER,content BLOB,PRIMARY KEY(context,ordinal)"),
                ("classes","context INTEGER,id INTEGER,content TEXT,PRIMARY KEY(context,id)"),
            ):
                table = self.prefix + "_" + suffix
                self.connection.execute(f"CREATE TABLE {table} ({definition})")
                self.tables.append(table)
        except BaseException:
            writer.failed = True
            self.close()
            raise

    def _render_content(self, content):
        if content is None:
            return 0
        if len(content) > MAX_RENDER_BYTES:
            raise EvidenceStoreError("rendered submission identity exceeds capacity",capacity=True)
        row = self.connection.execute(f"SELECT id FROM {self.prefix}_rendered WHERE content=?",(content,)).fetchone()
        if row is not None:
            return row[0]
        self.rendered_count += 1
        if self.rendered_count >= 2**32:
            raise EvidenceStoreError("case rendered identities exceed capacity",capacity=True)
        self.connection.execute(f"INSERT INTO {self.prefix}_rendered VALUES (?,?)",(self.rendered_count,content))
        return self.rendered_count

    def add_context(self, ordinal, scores, representations, eligibility):
        self._require_open()
        universe,members = scores.universe,scores.universe.members
        try:
            if (self.published or type(ordinal) is not int or ordinal != next_context_ordinal(self) or ordinal >= 2**32
                    or any(o.closed or o.connection is not self.connection for o in (scores,representations,eligibility))
                    or getattr(scores,"mode",None) is not None or members.ordinal != ordinal
                    or representations.members is not members or eligibility.universe is not universe
                    or eligibility.representations is not representations or eligibility.population is not scores.population):
                raise EvidenceStoreError("foreign or incomplete case selection context")
            projection = _SelectionContextRows(scores,representations,eligibility)
            self._append_context(_selection_context_header(ordinal,scores,representations,eligibility),
                projection.members_rows(),projection.singleton_rows(),projection.class_rows(),_selection_pair_rows(eligibility))
        except BaseException:
            self.writer.failed = True
            raise

    def _append_context(self, header, members, singletons, classes, pairs):
        self._require_open()
        ordinal = header["ordinal"]
        if self.published or type(ordinal) is not int or ordinal != next_context_ordinal(self) or ordinal >= 2**32:
            raise EvidenceStoreError("foreign or incomplete case selection context")
        try:
            name = self.name + ".context-" + str(ordinal)
            member_count = singleton_count = pair_count = 0
            def member_rows():
                nonlocal member_count
                for content,event,index in members:
                    member_count += 1
                    yield self._render_content(content),event,index
            with LocatorTape(self.writer,member_rows(),width=3,name=name + ".members") as tape:
                self.connection.execute(f"INSERT INTO {self.prefix}_members SELECT ?,ordinal,count,content FROM {tape.table}",(ordinal,))
                member_reference = {"name":tape.name,"sha256":tape.sha256}
            def singleton_rows():
                nonlocal singleton_count
                for number,member,code in singletons:
                    singleton_count += 1
                    if number != singleton_count:
                        raise EvidenceStoreError("selection singleton order changed")
                    self.connection.execute(f"INSERT INTO {self.prefix}_singletons VALUES (?,?,?,?)",(ordinal,number,member,code))
                    yield member,code
            with LocatorTape(self.writer,singleton_rows(),width=2,name=name + ".singletons") as tape:
                singleton_reference = {"name":tape.name,"sha256":tape.sha256}
            for number,value in classes:
                self.connection.execute(f"INSERT INTO {self.prefix}_classes VALUES (?,?,?)",(ordinal,number,_canonical(value)))
            for row in pairs:
                pair_count += row[1]
                self.connection.execute(f"INSERT INTO {self.prefix}_pairs VALUES (?,?,?,?)",(ordinal,*row))
            if (member_count != header["member_count"] or singleton_count != header["singleton_count"]
                    or pair_count != header["pair_count"]):
                raise EvidenceStoreError("selection context coverage changed")
            value = {**header,"member_map":member_reference,"singleton_map":singleton_reference}
            sha = self.writer.add(name,"context-submission-inputs-v1",value)
            self.connection.execute(f"INSERT INTO {self.prefix}_contexts VALUES (?,?,?)",
                (ordinal,header["context_id"],_canonical({**value,"root":{"name":name,"sha256":sha}})))
            self.context_count += 1
            self.hypothesis_count += header["hypothesis_count"]
        except BaseException:
            self.writer.failed = True
            raise

    def publish(self):
        require_complete_domain(self)
        self._require_open()
        if self.published:
            raise EvidenceStoreError("case submission inputs are already published")
        try:
            rendered_sha = self.writer.add_sequence(self.name + ".rendered","case-rendered-identities-v1",
                ((i,json.loads(content)) for i,content in self.connection.execute(f"SELECT * FROM {self.prefix}_rendered ORDER BY id")))
            contexts_sha = self.writer.add_sequence(self.name + ".contexts","case-submission-contexts-v1",
                (json.loads(content) for content, in self.connection.execute(f"SELECT content FROM {self.prefix}_contexts ORDER BY id")))
            self.sha256 = self.writer.add(self.name,"case-submission-inputs-v1",{
                "schema":"genosapiens.wgs-case-submission-inputs","version":1,"context_count":self.context_count,
                "hypothesis_count":self.hypothesis_count,"rendered_count":self.rendered_count,
                "contexts":{"name":self.name + ".contexts","sha256":contexts_sha},
                "rendered":{"name":self.name + ".rendered","sha256":rendered_sha}})
            self.published = True
        except BaseException:
            self.writer.failed = True
            raise

    @property
    def reference(self):
        self._require_open()
        if not self.published:
            raise EvidenceStoreError("case submission input closure is incomplete")
        return {"name":self.name,"kind":"case-submission-inputs-v1","sha256":self.sha256}

    def _remember(self, cache, kind, key, value, size, limit):
        if size > limit:
            raise EvidenceStoreError("case submission window exceeds capacity",capacity=True)
        while cache and self.cache_bytes[kind]+size > limit:
            _,(_,old_size) = cache.popitem(last=False)
            self.cache_bytes[kind] -= old_size
        cache[key] = value,size
        self.cache_bytes[kind] += size
        return value

    def _context(self, context):
        cached = self.context_cache.get(context)
        if cached is not None:
            self.context_cache.move_to_end(context)
            return cached[0]
        row = self.connection.execute(f"SELECT content FROM {self.prefix}_contexts WHERE id=?",(context,)).fetchone()
        if row is None:
            raise EvidenceStoreError("selection context is missing")
        metadata = json.loads(row[0])
        rendered = array("I",[0])
        events,parents = array("I"),bytearray()
        expected_chunk = 0
        for ordinal,count,content in self.connection.execute(
                f"SELECT ordinal,count,content FROM {self.prefix}_members WHERE context=? ORDER BY ordinal",(context,)):
            expected_chunk += 1
            if ordinal != expected_chunk:
                raise EvidenceStoreError("selection member chunk coverage changed")
            _,values = _decode(content,count,3)
            for offset in range(0,len(values),3):
                key,event,parent = values[offset:offset+3]
                if parent not in (0,1) or len(rendered) > metadata["member_count"]:
                    raise EvidenceStoreError("selection member metadata changed")
                rendered.append(key)
                events.append(event)
                parents.append(parent)
        if len(rendered) != metadata["member_count"]+1:
            raise EvidenceStoreError("selection context membership changed")
        positions = None if not metadata["has_pair_domain"] else CompactPairDomainPositions(
            ((i,event,bool(parent)) for i,(event,parent) in enumerate(zip(events,parents,strict=True),1)),
            member_count=metadata["member_count"],pair_count=metadata["pair_count"],require_open=self._require_open)
        size = sys.getsizeof(rendered)+(0 if positions is None else positions.metadata_bytes)+256
        value = _Context(rendered,positions,metadata["pair_offset"],size)
        return self._remember(self.context_cache,"context",context,value,size,MAX_CONTEXT_CACHE_BYTES)

    def _pair_code(self, context, position):
        key = context,position//CHUNK_PAIRS+1
        cached = self.pair_cache.get(key)
        if cached is None:
            row = self.connection.execute(f"SELECT count,content FROM {self.prefix}_pairs WHERE context=? AND ordinal=?",key).fetchone()
            if row is None:
                raise EvidenceStoreError("selection pair evidence is missing")
            values = outcome_chunk_values(*row)
            self._remember(self.pair_cache,"pair",key,values,sys.getsizeof(values)+192,MAX_PAIR_CACHE_BYTES)
        else:
            self.pair_cache.move_to_end(key)
            values = cached[0]
        if position%CHUNK_PAIRS >= len(values):
            raise EvidenceStoreError("selection pair position is missing")
        return values[position%CHUNK_PAIRS]

    def _class(self, context, code):
        key = context,code
        cached = self.class_cache.get(key)
        if cached is not None:
            self.class_cache.move_to_end(key)
            return cached[0]
        row = self.connection.execute(f"SELECT content FROM {self.prefix}_classes WHERE context=? AND id=?",key).fetchone()
        if row is None:
            raise EvidenceStoreError("selection eligibility class is missing")
        reasons,complete,partial = json.loads(row[0])
        value = tuple(reasons),complete,partial
        size = len(row[0])+sum(sys.getsizeof(r) for r in reasons)+sys.getsizeof(value[0])+256
        return self._remember(self.class_cache,"class",key,value,size,MAX_CLASS_CACHE_BYTES)

    def rows_for_order(self, order):
        self._require_open()
        if (not self.published or order.closed or not order.published or order.connection is not self.connection
                or order.context_count != self.context_count or order.hypothesis_count != self.hypothesis_count):
            raise EvidenceStoreError("selection order differs from case inputs")
        for ordinal,content in self.connection.execute(f"SELECT id,content FROM {self.prefix}_contexts ORDER BY id"):
            source,ordered = json.loads(content),order.context(ordinal)
            if any(source[key] != ordered[key] for key in ("context_id","hypothesis_count","scores")):
                raise EvidenceStoreError("selection order context parent changed")
        for locator in order.tape.rows():
            context,left,right = locator
            value = self._context(context)
            if right == 0:
                row = self.connection.execute(f"SELECT member,code FROM {self.prefix}_singletons WHERE context=? AND ordinal=?",(context,left)).fetchone()
                if row is None:
                    raise EvidenceStoreError("selection singleton is missing")
                member,code = row
                keys = () if member == 0 else (member,)
            else:
                if value.positions is None:
                    raise EvidenceStoreError("selection pair domain is missing")
                position = value.positions.ordinal(left,right)-1
                code = self._pair_code(context,position)-value.offset
                if code < 0:
                    raise EvidenceStoreError("selection order names an unconstructed pair")
                keys = left,right
            reasons,complete,partial = self._class(context,code)
            rendered = tuple(sorted(value.rendered[k] for k in keys))
            if not reasons and (not rendered or any(k == 0 for k in rendered)):
                raise EvidenceStoreError("eligible selection identity is missing")
            # Failed representations have no rendered identity. All original
            # ordered exclusion reasons still survive in the complete audit.
            yield SelectionOccurrence(locator,keys,() if reasons else rendered,complete,partial,reasons)

    def _require_open(self):
        if self.closed or self.connection is None or self.writer.connection is not self.connection:
            raise EvidenceStoreError("case submission input owner is closed")

    def close(self):
        for table in reversed(self.tables):
            self.connection.execute("DROP TABLE " + table)
        self.tables.clear()
        self.context_cache.clear()
        self.pair_cache.clear()
        self.class_cache.clear()
        self.closed = True

    def __enter__(self):
        return self

    def __exit__(self,*_):
        self.close()


class CompletedCaseSelectionInputs(CompleteCaseSelectionInputs):
    """Read saved representation and eligibility maps without re-evaluation."""
    def __init__(self, reader, reference, *, context_at, progress=None, context_domain=None):
        from genosapiens.wgs_completed_science import child, fields, header, integer, require, outcome_chunks, observed_rows
        from genosapiens.wgs_locator_tape import VerifiedLocatorTape
        value = header(reader, reference, kind='case-submission-inputs-v1',
            schema='genosapiens.wgs-case-submission-inputs', keys=('context_count', 'hypothesis_count',
                'rendered_count', 'contexts', 'rendered'))
        super().__init__(EvidenceReplay(reader), name=reference.name, context_domain=context_domain)
        try:
            for row in reader.iter_array(child(reader, value['rendered'], 'case-rendered-identities-v1'),
                    max_item_bytes=MAX_RENDER_BYTES):
                require(type(row) is list and len(row) == 2
                    and integer(row[0], minimum=1, maximum=2**32-1) == self.rendered_count + 1)
                self.connection.execute(f'INSERT INTO {self.prefix}_rendered VALUES (?,?)', (row[0], _canonical(row[1])))
                self.rendered_count += 1
            for row in observed_rows(reader.iter_array(child(reader, value['contexts'], 'case-submission-contexts-v1')),
                    progress, 'saved_selection_inputs', value['context_count']):
                fields(row, ('ordinal', 'context_id', 'member_count', 'hypothesis_count', 'singleton_count',
                    'pair_count', 'pair_offset', 'has_pair_domain', 'scores', 'members', 'hypotheses',
                    'representations', 'eligibility', 'member_map', 'singleton_map', 'root',
                    *(('population_work',) if 'population_work' in row else ())))
                ordinal = integer(row['ordinal'], minimum=1, maximum=2**32-1)
                require(ordinal == next_context_ordinal(self))
                context = context_at(ordinal)
                require(all(row[key] == context[key] for key in ('context_id', 'member_count', 'hypothesis_count', 'pair_count')))
                root = child(reader, row['root'], 'context-submission-inputs-v1')
                require(reader.reconstruct(root) == {k: v for k, v in row.items() if k != 'root'})
                from genosapiens.wgs_case_resources import bound_reference
                context_header = reader.reconstruct(bound_reference(reader, context['context'], 'scientific-case-context-v1'))
                from genosapiens.wgs_population_work_filter import verify_saved_work_domain
                verify_saved_work_domain(reader, context_header, row)
                for key, kind in (('scores', 'complete-context-scores-v1'), ('members', 'context-evidence-v1'),
                        ('hypotheses', 'complete-hypothesis-universe-v1'), ('representations', 'complete-submission-representations-v1')):
                    found = child(reader, row[key], kind)
                    require(context_header['universe' if key == 'hypotheses' else key] ==
                        dict(name=found.name, kind=found.kind, sha256=found.sha256))
                found = child(reader, row['eligibility'], 'complete-submission-eligibility-v1')
                require(context_header['eligibility'] == dict(name=found.name, kind=found.kind, sha256=found.sha256))
                require(type(row['has_pair_domain']) is bool and integer(row['pair_offset'], maximum=1) in (0, 1))
                with VerifiedLocatorTape(reader, child(reader, row['member_map'], 'locator-tape-v1'), width=3) as tape:
                    require(tape.count == integer(row['member_count']))
                    for rendered, _, parent in tape.rows():
                        require(rendered <= self.rendered_count and parent in (0, 1))
                    self.connection.execute(f'INSERT INTO {self.prefix}_members SELECT ?,ordinal,count,content FROM {tape.table}', (ordinal,))
                eligibility = reader.reconstruct(child(reader, row['eligibility'], 'complete-submission-eligibility-v1'))
                require(eligibility['hypotheses'] == row['hypotheses'] and eligibility['representations'] == row['representations']
                    and all(eligibility[k] == row[k] for k in ('singleton_count', 'pair_count', 'hypothesis_count')))
                classes = 0
                for item in reader.iter_array(child(reader, eligibility['catalog'], 'submission-eligibility-catalog-v1')):
                    fields(item, ('ordinal', 'inputs', 'complete', 'reasons'))
                    require(integer(item['ordinal']) == classes and type(item['complete']) is bool
                        and type(item['inputs']) is dict and type(item['inputs'].get('incomplete')) is bool
                        and type(item['reasons']) is list and all(type(reason) is str and reason for reason in item['reasons']))
                    self.connection.execute(f'INSERT INTO {self.prefix}_classes VALUES (?,?,?)',
                        (ordinal, classes, _canonical((item['reasons'], item['complete'], item['inputs']['incomplete']))))
                    classes += 1
                with VerifiedLocatorTape(reader, child(reader, row['singleton_map'], 'locator-tape-v1'), width=2) as tape:
                    require(tape.count == integer(row['singleton_count']))
                    for number, (member, code) in enumerate(tape.rows(), 1):
                        require(member <= row['member_count'] and code < classes)
                        self.connection.execute(f'INSERT INTO {self.prefix}_singletons VALUES (?,?,?,?)', (ordinal, number, member, code))
                if eligibility['tape'] is None:
                    require(row['pair_count'] == 0)
                else:
                    universe = reader.reconstruct(child(reader, row['hypotheses'], 'complete-hypothesis-universe-v1'))
                    pair = reader.reconstruct(child(reader, universe['pair_evidence'], 'complete-pair-evidence-v1'))
                    require(row['has_pair_domain'] and pair['pair_count'] == row['pair_count'])
                    codes = []
                    for item in reader.iter_array(child(reader, eligibility['pair_codes'], 'pair-eligibility-codes-v1')):
                        fields(item, ('ordinal', 'decision'))
                        require(integer(item['ordinal']) == len(codes))
                        decision = item['decision']
                        if decision is None:
                            require(row['pair_offset'] == 1)
                            codes.append(0)
                        else:
                            require(integer(decision) < classes)
                            codes.append(decision + row['pair_offset'])
                    self.connection.executemany(f'INSERT INTO {self.prefix}_pairs VALUES (?,?,?,?)',
                        ((ordinal, *chunk) for chunk in outcome_chunks(reader, eligibility['tape_name'], eligibility['tape'],
                            catalog_sha256=eligibility['pair_codes']['sha256'], catalog_count=len(codes),
                            pair_count=row['pair_count'], domain_sha256=pair['domain_summary']['membership_sha256'], mapping=codes)))
                self.connection.execute(f'INSERT INTO {self.prefix}_contexts VALUES (?,?,?)', (ordinal, row['context_id'], _canonical(row)))
                self.context_count += 1
                self.hypothesis_count += integer(row['hypothesis_count'])
            require(all(getattr(self, key) == integer(value[key]) for key in ('context_count', 'hypothesis_count', 'rendered_count')))
            self.sha256, self.published = reference.sha256, True
        except BaseException:
            self.close()
            raise


class VerifiedCaseSelectionInputs(CompleteCaseSelectionInputs):
    def __init__(self, reader, reference, contexts):
        if reference.kind != "case-submission-inputs-v1":
            raise EvidenceStoreError("foreign case submission input root")
        reader._require_reference(reference)
        super().__init__(EvidenceReplay(reader),name=reference.name)
        try:
            for ordinal,scores,representations,eligibility in contexts:
                self.add_context(ordinal,scores,representations,eligibility)
            self.publish()
        except BaseException:
            self.close()
            raise
