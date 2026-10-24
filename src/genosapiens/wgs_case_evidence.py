"""Complete primary science, order, selection and selected evidence.

Only bounded scalar projections survive each context. Optional analyses consume
this immutable primary through their own output owner.
"""
from contextlib import ExitStack, nullcontext
import hashlib
import json

from genosapiens.evidence_serialization import _canonical
from genosapiens.evidence_model import _exact
from genosapiens.wgs_case_context import VerifiedScientificCaseContext
from genosapiens.wgs_case_identities import CompleteCaseIdentities
from genosapiens.wgs_case_order import CompleteCaseOrder
from genosapiens.wgs_case_resources import bound_reference, component_reference
from genosapiens.wgs_case_science_facts import CompleteCaseScienceFacts
from genosapiens.wgs_evidence_replay import EvidenceReplay
from genosapiens.wgs_evidence_store import EvidenceStoreError, EvidenceWriter
from genosapiens.wgs_selection_evidence import CompleteSelectionEvidence
from genosapiens.wgs_selection_inputs import CompleteCaseSelectionInputs
from genosapiens.wgs_source_closure import CompleteSourceClosure
from genosapiens.wgs_context_domain import ContextDomain, next_context_ordinal, require_complete_domain

MAX_CASE_CONTEXTS = 1_000_000
MAX_CONTEXT_DIRECTORY_ROW_BYTES = 64 * 1024
MAX_SELECTED_GROUP_BYTES = 8 * 1024**2



def _replay_selected_groups(entries, restore, identity):
    """Replay selected contexts through the same producer and recovery kernel."""
    selected_by_context, groups, group_bytes = {}, [], 0
    for rank, locator, total in entries:
        selected_by_context.setdefault(locator[0], []).append((rank, locator[1:], total))
    for ordinal in sorted(selected_by_context):
        with restore(ordinal) as context:
            from genosapiens.evidence_model import Completeness
            from genosapiens.track1_submission import SubmissionGroup, _allele_trace_key
            positions = None if context.universe.pairs is None else context.universe.pairs.domain.positions()
            for rank, locator, total in selected_by_context[ordinal]:
                constructed, _, ledger, _ = context.scores.record_at(locator, positions=positions)
                if (context.eligibility.reasons_at(locator, positions=positions)
                        or ledger.total != total
                        or constructed.hypothesis.hypothesis_id != identity((context.ordinal, *locator))):
                    raise EvidenceStoreError("selected scientific group does not replay")
                traces = tuple(sorted((context.representations.get(ref.allele_id).trace
                    for ref in constructed.hypothesis.allele_refs), key=_allele_trace_key))
                group = SubmissionGroup(constructed.hypothesis.hypothesis_id, rank, total, traces, (),
                    constructed.hypothesis.completeness is Completeness.INCOMPLETE)
                group_bytes += len(_canonical(group))
                if group_bytes > MAX_SELECTED_GROUP_BYTES:
                    raise EvidenceStoreError("selected scientific groups exceed capacity", capacity=True)
                groups.append(group)
    return tuple(sorted(groups, key=lambda group: group.source_rank))


class CompleteScientificCase:
    def __init__(self, writer, resources, *, name, context_domain=None):
        if context_domain is not None and type(context_domain) is not ContextDomain:
            raise EvidenceStoreError('invalid explicit scientific domain')
        self.context_domain = context_domain
        self.writer,self.connection,self.resources,self.name = writer,writer.connection,resources,name
        self.stack,self.closed,self.created = ExitStack(),False,False
        self.primary_published,self.published = False,False
        self.context_count,self.member_count,self.pair_count,self.hypothesis_count = 0,0,0,0
        self.positive_context_count = 0
        self.table = "scientific_case_" + hashlib.sha256(name.encode()).hexdigest()
        self.orders = {}
        try:
            resources.require(writer)
            self.connection.execute(f"CREATE TABLE {self.table} (ordinal INTEGER PRIMARY KEY,id TEXT UNIQUE,content TEXT,views TEXT)")
            self.created = True
            self.orders["primary"] = self._order("primary")
            self.inputs = self.stack.enter_context(CompleteCaseSelectionInputs(writer,name=name + ".submission-inputs",
                context_domain=context_domain))
            self.scalar_facts = self.stack.enter_context(CompleteCaseScienceFacts(writer,self.inputs,name=name+".scalar-facts",
                assessment_domain_sha256=resources.assessment_domain_sha256, context_domain=context_domain))
            self.identities = self.stack.enter_context(CompleteCaseIdentities(writer,name=name + ".identities",
                context_domain=context_domain))
            self.closure = self.stack.enter_context(CompleteSourceClosure(writer,resources.genes,resources.phenotypes,
                name=name + ".source-closure", context_domain=context_domain,
                resource_prefix=getattr(resources, '_source_prefix', None)))
        except BaseException:
            writer.failed = True
            self.close()
            raise

    def _order(self, mode):
        return self.stack.enter_context(CompleteCaseOrder(self.writer,self.resources.fusion_policy,
            policy_sha256=self.resources.fusion_policy_sha256,name=self.name + ".order." + mode,
            context_domain=self.context_domain))

    def _require_open(self):
        if self.closed or self.connection is None:
            raise EvidenceStoreError("scientific case is closed")
        self.resources.require(self.writer)

    def add_context(self, context):
        self._require_open()
        try:
            if (self.primary_published or context.closed or context.connection is not self.connection
                    or context.resources is not self.resources or context.ordinal != next_context_ordinal(self)):
                raise EvidenceStoreError("foreign or incomplete scientific case context")
            if self.context_count >= MAX_CASE_CONTEXTS:
                raise EvidenceStoreError("scientific case context count exceeds capacity",capacity=True)
            self.closure.add_context(context.members,context.universe,context.scores)
            self.closure.add_reference_queries(context.representations.queries)
            self.identities.add_context(context)
            self.inputs.add_context(context.ordinal,context.scores,context.representations,context.eligibility)
            self.scalar_facts.add_context(context)
            self.orders["primary"].add_context(context.ordinal,context.order)
            row = {"ordinal":context.ordinal,"context_id":context.context_id,"context":context.reference,
                "member_count":context.members.member_count,"pair_count":context.pair_count,
                "hypothesis_count":context.hypothesis_count}
            self._append_context_directory(row)
        except BaseException:
            self.writer.failed = True
            raise

    def _append_context_directory(self, row):
        """Advance complete coverage only after all context collectors succeed."""
        self._require_open()
        try:
            if (self.primary_published or row["ordinal"] != next_context_ordinal(self)
                    or self.context_count >= MAX_CASE_CONTEXTS
                    or any(owner.context_count != self.context_count + 1 for owner in
                        (self.closure,self.identities,self.inputs,self.scalar_facts,self.orders["primary"]))):
                raise EvidenceStoreError("scientific context collectors are incomplete")
            content = _canonical(row)
            if len(content) > MAX_CONTEXT_DIRECTORY_ROW_BYTES:
                raise EvidenceStoreError("scientific context directory row exceeds capacity",capacity=True)
            self.connection.execute(f"INSERT INTO {self.table} VALUES (?,?,?,NULL)",
                (row["ordinal"],row["context_id"],content))
            self.context_count += 1
            self.positive_context_count += int(row["member_count"] > 0)
            self.member_count += row["member_count"]
            self.pair_count += row["pair_count"]
            self.hypothesis_count += row["hypothesis_count"]
        except BaseException:
            self.writer.failed = True
            raise

    def contexts(self):
        self._require_open()
        # Do not suspend a live SQLite cursor across context construction:
        # temporary component tables must be dropped between contexts.
        previous = 0
        for _ in range(self.context_count):
            row = self.connection.execute(f'SELECT ordinal FROM {self.table} WHERE ordinal>? ORDER BY ordinal LIMIT 1',
                (previous,)).fetchone()
            if row is None:
                raise EvidenceStoreError('scientific context directory is incomplete')
            previous = row[0]
            yield self.context(previous)

    def context(self, ordinal):
        self._require_open()
        if type(ordinal) is not int:
            raise EvidenceStoreError("invalid scientific context ordinal")
        row = self.connection.execute(f"SELECT content FROM {self.table} WHERE ordinal=?",(ordinal,)).fetchone()
        if row is None:
            raise EvidenceStoreError("scientific context is unavailable")
        return json.loads(row[0])

    def ordered(self, mode="primary"):
        """All original IDs and exact ordered totals; no rich pair expansion."""
        self._require_open()
        if not self.published or mode not in self.orders:
            raise EvidenceStoreError("complete scientific order is unavailable")
        order = self.orders[mode]
        count = 0
        for family in order.families():
            if count != family["start"]:
                raise EvidenceStoreError("scientific order family coverage changed")
            for locator in order.tape.rows(start=family["start"],count=family["count"]):
                count += 1
                yield {"rank":count,"locator":locator,"hypothesis_id":self.identities.identity(locator),
                    "total":family["total"],"ties":tuple(family["ties"]),"mode":mode}
        if count != self.hypothesis_count:
            raise EvidenceStoreError("scientific order count changed")
        self._require_open()

    def ordered_at(self, rank, mode="primary"):
        self._require_open()
        if not self.published or mode not in self.orders:
            raise EvidenceStoreError("complete scientific order is unavailable")
        locator,total,ties = self.orders[mode].at_rank(rank)
        return {"rank":rank,"locator":locator,"hypothesis_id":self.identities.identity(locator),
            "total":total,"ties":ties,"mode":mode}

    def rank_of(self, locator, *, total, ties, mode="primary"):
        """Find one computed occurrence without a second case-wide H index.

        Each exact total/tie family is already sorted by hypothesis ID. Binary
        search reads bounded tape windows and verifies the final locator.
        """
        self._require_open()
        if not self.published or mode not in self.orders:
            raise EvidenceStoreError("complete scientific order is unavailable")
        order = self.orders[mode]
        family = self.connection.execute(f"SELECT start,count FROM {order.prefix}_families WHERE content=?",
            (_canonical((total, tuple(ties))),)).fetchone()
        if family is None:
            raise EvidenceStoreError("scientific rank family is missing")
        identifier = self.identities.identity(locator)
        left, right = family[0], family[0] + family[1]
        while left < right:
            middle = (left + right) // 2
            candidate, = order.tape.rows(start=middle, count=1)
            current = self.identities.identity(candidate)
            if current < identifier:
                left = middle + 1
            else:
                right = middle
        if left == family[0] + family[1]:
            raise EvidenceStoreError("scientific ranked occurrence is missing")
        candidate, = order.tape.rows(start=left, count=1)
        if candidate != locator:
            raise EvidenceStoreError("scientific ranked occurrence differs")
        return left + 1

    def scalar_records(self, mode="primary"):
        for ordered in self.ordered(mode):
            yield {"ordering":ordered,**self.scalar_facts.get(ordered["locator"])}

    def credible_complete_count(self, threshold):
        self._require_open()
        if not self.published:
            raise EvidenceStoreError("complete scientific credibility facts are unavailable")
        cutoff = _exact(threshold)
        count = 0
        order = self.orders["primary"]
        for family in order.families():
            if _exact(family["total"]) < cutoff:
                continue
            for locator in order.tape.rows(start=family["start"],count=family["count"]):
                facts = self.scalar_facts.get(locator)
                if (facts["completeness"] == "complete" and facts["model_state"] == "supported"
                        and all(state == "qualifying" for state in facts["member_states"])):
                    count += 1
        return count

    def selection_audit(self):
        self._require_open()
        if not self.published:
            raise EvidenceStoreError("complete scientific selection is unavailable")
        for row,audit in zip(self.ordered(),self.selection.audit(),strict=True):
            rank,category,reasons,primary_rank = audit
            if rank != row["rank"]:
                raise EvidenceStoreError("scientific selection rank changed")
            primary = self.ordered_at(primary_rank)["hypothesis_id"] if primary_rank else None
            yield {"rank":rank,"locator":row["locator"],"hypothesis_id":row["hypothesis_id"],
                "source_score":row["total"],"category":category,"reasons":reasons,
                "primary_hypothesis_id":primary}

    def record_at_rank(self, reader, rank, reference, reference_pin, renderer_cache, *, mode="primary"):
        """Reconstruct one rich occurrence and close its context immediately."""
        ordered = self.ordered_at(rank,mode)
        if reader.connection is not self.connection:
            raise EvidenceStoreError("foreign scientific case read owner")
        row = self.context(ordered["locator"][0])
        with self._restore(reader,row,reference,reference_pin,renderer_cache) as context,ExitStack() as stack:
            scores = self.context_scores(reader,context,mode,stack)
            locator = ordered["locator"][1:]
            positions = None if context.universe.pairs is None else context.universe.pairs.domain.positions()
            record = self.logical_record(context,locator,scores,positions=positions)
            total = record["ledger"].total
            ordering_total = total if mode == "primary" else scores.ordering_total(total)
            if (record["constructed"].hypothesis.hypothesis_id != ordered["hypothesis_id"] or ordering_total != ordered["total"]):
                raise EvidenceStoreError("scientific ordered occurrence does not replay")
            return {**record,"ordering":ordered}

    def context_scores(self, reader, context, mode, stack):
        self._require_open()
        if mode != "primary" or context.resources is not self.resources or reader.connection is not self.connection:
            raise EvidenceStoreError("foreign scientific score view")
        scores = context.scores
        if mode == "primary":
            return scores
    def _owns_context_resources(self, resources):
        return resources is self.resources

    def logical_record(self, context, locator, scores, *, positions):
        self._require_open()
        if (context.closed or not self._owns_context_resources(context.resources) or scores.universe is not context.universe
                or scores is not context.scores and getattr(scores,"primary",None) is not context.scores):
            raise EvidenceStoreError("foreign logical scientific record")
        constructed,evidence,ledger,population = scores.record_at(locator,positions=positions)
        reasons = context.eligibility.reasons_at(locator,positions=positions)
        from genosapiens.track1_submission import _allele_trace_key
        traces = (() if reasons else tuple(sorted((context.representations.get(ref.allele_id).trace
            for ref in constructed.hypothesis.allele_refs),key=_allele_trace_key)))
        record = {"context_ordinal":context.ordinal,"context_id":context.context_id,"constructed":constructed,
            "evidence":evidence,"ledger":ledger,"population_assessment":population,
            "population_policy":self.resources.population_reference,"eligibility_reasons":reasons,
            "alleles":tuple(context.members.member(ref.allele_id).allele for ref in constructed.hypothesis.allele_refs),
            "allele_traces":traces}
        # Charge the aggregate, not merely each individually bounded field.
        size = 2
        for key,value in record.items():
            size += len(_canonical(key))+len(_canonical(value))+2
            if size > 128 * 1024**2:
                raise EvidenceStoreError("logical scientific record exceeds capacity",capacity=True)
        return record

    def publish_primary(self):
        self._require_open()
        try:
            require_complete_domain(self)
            if self.primary_published:
                raise EvidenceStoreError("scientific primary is already published")
            self.orders["primary"].publish()
            self.inputs.publish()
            self.scalar_facts.publish()
            self.identities.publish()
            self.closure.publish()
            self.selection = self.stack.enter_context(CompleteSelectionEvidence(self.writer,self.orders["primary"],
                self.inputs.reference,self.inputs.rows_for_order(self.orders["primary"]),name=self.name + ".selection"))
            if (self.hypothesis_count != self.orders["primary"].hypothesis_count
                    or self.hypothesis_count != self.selection.count
                    or self.inputs.context_count != self.context_count or self.closure.context_count != self.context_count):
                raise EvidenceStoreError("scientific primary coverage changed")
            digest = self.writer.add_sequence(self.name + ".contexts","scientific-case-contexts-v1",self.contexts())
            self.contexts_reference = {"name":self.name + ".contexts","kind":"scientific-case-contexts-v1","sha256":digest}
            domain = {}
            if self.context_domain is not None:
                digest = self.writer.add_sequence(self.name + '.domain', 'scientific-context-domain-v1', self.context_domain)
                self.domain_reference = dict(name=self.name + '.domain', kind='scientific-context-domain-v1', sha256=digest)
                domain = {'context_domain': self.domain_reference}
            self.primary_header = {"schema":"genosapiens.wgs-scientific-primary","version":2 if domain else 1,
                **domain,
                "resources":self.resources.reference,"contexts":self.contexts_reference,
                "context_count":self.context_count,"member_count":self.member_count,"pair_count":self.pair_count,
                "positive_context_count":self.positive_context_count,
                "hypothesis_count":self.hypothesis_count,"order":self.orders["primary"].reference,
                "selection_inputs":self.inputs.reference,
                "scalar_facts":component_reference(self.scalar_facts,"case-scientific-scalar-facts-v1"),
                "identities":self.identities.reference,
                "selection":component_reference(self.selection,"complete-submission-selection-v1"),
                "source_closure":component_reference(self.closure,"complete-scientific-source-closure-v1")}
            digest = self.writer.add(self.name + ".primary","scientific-primary-evidence-v1",self.primary_header)
            self.primary_reference = {"name":self.name + ".primary","kind":"scientific-primary-evidence-v1","sha256":digest}
            self.primary_published = True
        except BaseException:
            self.writer.failed = True
            raise

    def _restore(self, reader, row, reference, reference_pin, renderer_cache, *, retain_identity_index=False):
        return VerifiedScientificCaseContext(reader,
            bound_reference(reader,row["context"],"scientific-case-context-v1"),
            self.resources,reference,reference_pin,renderer_cache,retain_identity_index=retain_identity_index)

    def publish(self, reference, reference_pin, renderer_cache):
        """Assemble selected groups and seal the independently complete primary."""
        self._require_open()
        try:
            if not self.primary_published or self.published:
                raise EvidenceStoreError("scientific primary is unavailable or case is already published")
            scope = self.writer.read_view() if type(self.writer) is EvidenceWriter else nullcontext(self.writer.reader)
            with scope as reader:
                primary = self.orders["primary"]
                entries, by_context = [], {}
                for rank in self.selection.selected:
                    locator, total, ties = primary.at_rank(rank)
                    entries.append((rank, locator, total))
                    by_context.setdefault(locator[0], []).append(dict(rank=rank, locator=locator, total=total,
                        ties=ties, hypothesis_id=self.identities.identity(locator), mode='primary'))
                from genosapiens.wgs_completed_science import selected_records
                from genosapiens.evidence_model import Completeness
                from genosapiens.track1_submission import SubmissionGroup
                groups, detail_bytes = [], 0
                for ordinal in sorted(by_context):
                    header = reader.reconstruct(bound_reference(reader, self.context(ordinal)['context'], 'scientific-case-context-v1'))
                    if header['version'] == 1 and type(self.writer) is not EvidenceWriter:
                        groups.extend(_replay_selected_groups([entry for entry in entries if entry[1][0] == ordinal],
                            lambda number: self._restore(reader, self.context(number), reference, reference_pin, renderer_cache),
                            self.identities.identity))
                        continue
                    for row in selected_records(reader, self, self.resources, ordinal, by_context[ordinal],
                            reference_names=reference.public_names):
                        detail_bytes += len(_canonical(row))
                        if detail_bytes > 128 * 1024**2:
                            raise EvidenceStoreError('selected scientific detail exceeds capacity', capacity=True)
                        hypothesis = row['constructed'].hypothesis
                        groups.append(SubmissionGroup(hypothesis.hypothesis_id, row['ordering']['rank'],
                            row['ledger'].total, row['allele_traces'], (), hypothesis.completeness is Completeness.INCOMPLETE))
                self.selected_groups = tuple(sorted(groups, key=lambda group: group.source_rank))
                if len(_canonical(self.selected_groups)) > MAX_SELECTED_GROUP_BYTES:
                    raise EvidenceStoreError('selected scientific groups exceed capacity', capacity=True)
                if tuple(g.source_rank for g in self.selected_groups) != tuple(self.selection.selected):
                    raise EvidenceStoreError("selected scientific group coverage changed")
                selected_sha = self.writer.add(self.name + ".selected","scientific-selected-groups-v1",{
                    "primary":self.primary_reference,"selection":component_reference(self.selection,"complete-submission-selection-v1"),
                    "reference":self.resources.reference_binding,"groups":self.selected_groups})
                self.header = {"schema":"genosapiens.wgs-complete-scientific-case",
                    "version":3 if self.context_domain is not None else 2,
                    **({'context_domain': self.domain_reference} if self.context_domain is not None else {}),
                    "resources":self.resources.reference,"primary":self.primary_reference,
                    "selected":{"name":self.name + ".selected","kind":"scientific-selected-groups-v1","sha256":selected_sha},
                    "order":self.orders["primary"].reference,
                    "assessment_domain_sha256":self.resources.assessment_domain_sha256,
                    "optional_analysis_status":"not_run",
                    "context_count":self.context_count,"member_count":self.member_count,
                    "positive_context_count":self.positive_context_count,
                    "pair_count":self.pair_count,"hypothesis_count":self.hypothesis_count}
                self.sha256 = self.writer.add(self.name,"complete-scientific-case-v1",self.header)
            self.published = True
            return self.sha256
        except BaseException:
            self.writer.failed = True
            raise

    def publish_sensitivities(self, *, primary_reader, primary_product, output_directory):
        """Explicit optional publication from the immutable, verified primary."""
        if primary_reader.case is not self or not self.published:
            raise EvidenceStoreError("sealed primary reader required for sensitivities")
        from genosapiens.wgs_scientific_sensitivities import write_scientific_sensitivities
        return write_scientific_sensitivities(primary_reader, primary_product, output_directory=output_directory)

    @property
    def reference(self):
        self._require_open()
        if not self.published:
            raise EvidenceStoreError("scientific case is incomplete")
        return component_reference(self,"complete-scientific-case-v1")

    def close(self):
        self.stack.close()
        if self.created:
            self.connection.execute("DROP TABLE " + self.table)
            self.created = False
        self.closed = True

    def __enter__(self):
        return self

    def __exit__(self,*_):
        self.close()


class CompletedScientificCase(CompleteScientificCase):
    """Read complete producer projections; no context or primary production."""
    def __init__(self, reader, root, resources, *, progress=None):
        from genosapiens.wgs_completed_science import fields, header, integer, require, observed_rows
        from genosapiens.wgs_case_order import CompletedCaseOrder
        from genosapiens.wgs_case_identities import CompletedCaseIdentities
        from genosapiens.wgs_selection_inputs import CompletedCaseSelectionInputs
        from genosapiens.wgs_case_science_facts import CompletedCaseScienceFacts
        from genosapiens.wgs_selection_evidence import CompletedSelectionEvidence
        from genosapiens.evidence_model import _typed
        from genosapiens.track1_submission import SubmissionGroup
        from genosapiens.wgs_runtime import analysis_progress
        from genosapiens.wgs_case_context import _KINDS
        self.writer, self.connection, self.resources, self.name = EvidenceReplay(reader), reader.connection, resources, root.name
        self.stack, self.closed, self.created = ExitStack(), False, False
        self.primary_published = self.published = False
        self.context_count = self.member_count = self.pair_count = self.hypothesis_count = self.positive_context_count = 0
        self.table = 'scientific_case_' + hashlib.sha256(root.name.encode()).hexdigest()
        self.orders = {}
        keys = ('context_count', 'positive_context_count', 'member_count', 'pair_count', 'hypothesis_count')
        try:
            version = reader.reconstruct(root).get('version')
            require(type(version) is int and version in (2, 3))
            extra = ('context_domain',) if version == 3 else ()
            self.header = header(reader, root, kind='complete-scientific-case-v1', schema='genosapiens.wgs-complete-scientific-case',
                version=version, keys=('resources', 'primary', 'order', 'selected', 'assessment_domain_sha256', 'optional_analysis_status', *keys, *extra))
            value = self.header
            self.context_domain = None
            if version == 3:
                self.domain_reference = value['context_domain']
                self.context_domain = ContextDomain.from_ordinals(reader.iter_array(
                    bound_reference(reader, self.domain_reference, 'scientific-context-domain-v1')))
            require(value['resources'] == resources.reference and value['assessment_domain_sha256'] == resources.assessment_domain_sha256
                and value['optional_analysis_status'] == 'not_run')
            self.primary_reference = value['primary']
            primary_root = bound_reference(reader, value['primary'], 'scientific-primary-evidence-v1')
            primary = self.primary_header = header(reader, primary_root, kind='scientific-primary-evidence-v1',
                version=2 if version == 3 else 1,
                schema='genosapiens.wgs-scientific-primary', keys=('resources', 'contexts', 'order', 'selection_inputs',
                    'scalar_facts', 'identities', 'selection', 'source_closure', *keys, *extra))
            require(version != 3 or primary['context_domain'] == self.domain_reference)
            require(primary['resources'] == resources.reference and primary['order'] == value['order'])
            require(all(integer(primary[k]) == integer(value[k]) for k in keys))
            require(value['context_count'] <= MAX_CASE_CONTEXTS)
            self.connection.execute(f'CREATE TABLE {self.table} (ordinal INTEGER PRIMARY KEY,id TEXT UNIQUE,content TEXT,views TEXT)')
            self.created = True
            self.contexts_reference = primary['contexts']
            for row in observed_rows(reader.iter_array(bound_reference(reader, primary['contexts'], 'scientific-case-contexts-v1'),
                    max_item_bytes=MAX_CONTEXT_DIRECTORY_ROW_BYTES), progress, 'saved_context_directory', value['context_count']):
                fields(row, ('ordinal', 'context_id', 'context', 'member_count', 'pair_count', 'hypothesis_count'))
                require(integer(row['ordinal'], minimum=1) == next_context_ordinal(self) and type(row['context_id']) is str and row['context_id'])
                context_root = bound_reference(reader, row['context'], 'scientific-case-context-v1')
                context_version = reader.reconstruct(context_root).get('version')
                require(type(context_version) is int and context_version in (1, 2))
                additional = ('selected_reconstruction',) if context_version == 2 else ()
                context = header(reader, context_root, kind='scientific-case-context-v1',
                    schema='genosapiens.wgs-scientific-case-context', keys=('resources', 'ordinal', 'context_id',
                        'member_count', 'pair_count', 'hypothesis_count', 'assessment_domain_sha256', *_KINDS, *additional),
                    version=context_version)
                if context_version == 2:
                    from genosapiens.wgs_completed_science import SelectedReconstructionIndex
                    SelectedReconstructionIndex(reader, context, resources)
                require(context['resources'] == resources.reference and context['assessment_domain_sha256'] == resources.assessment_domain_sha256
                    and all(context[k] == row[k] for k in ('ordinal', 'context_id', 'member_count', 'pair_count', 'hypothesis_count')))
                for key, kind in _KINDS.items():
                    bound_reference(reader, context[key], kind)
                self.connection.execute(f'INSERT INTO {self.table} VALUES (?,?,?,NULL)', (row['ordinal'], row['context_id'], _canonical(row)))
                self.context_count += 1
                self.positive_context_count += integer(row['member_count']) > 0
                for key in ('member_count', 'pair_count', 'hypothesis_count'):
                    setattr(self, key, getattr(self, key) + integer(row[key]))
            require(all(getattr(self, key) == value[key] for key in keys))
            require_complete_domain(self)
            analysis_progress(progress, 'publication', 'saved_scientific_projections', unit='contexts', completed=0, total=self.context_count)
            self.orders['primary'] = self.stack.enter_context(CompletedCaseOrder(reader,
                bound_reference(reader, primary['order'], 'complete-case-order-v1'), policy=resources.fusion_policy,
                policy_sha256=resources.fusion_policy_sha256, progress=progress, context_domain=self.context_domain))
            self.identities = self.stack.enter_context(CompletedCaseIdentities(reader,
                bound_reference(reader, primary['identities'], 'case-scientific-identities-v1'), context_at=self.context, progress=progress,
                context_domain=self.context_domain))
            self.inputs = self.stack.enter_context(CompletedCaseSelectionInputs(reader,
                bound_reference(reader, primary['selection_inputs'], 'case-submission-inputs-v1'), context_at=self.context, progress=progress,
                context_domain=self.context_domain))
            self.scalar_facts = self.stack.enter_context(CompletedCaseScienceFacts(reader,
                bound_reference(reader, primary['scalar_facts'], 'case-scientific-scalar-facts-v1'), self.inputs,
                assessment_domain_sha256=resources.assessment_domain_sha256, context_at=self.context, progress=progress,
                context_domain=self.context_domain))
            analysis_progress(progress, 'publication', 'saved_selection', force=True)
            self.selection = self.stack.enter_context(CompletedSelectionEvidence(reader,
                bound_reference(reader, primary['selection'], 'complete-submission-selection-v1'), self.orders['primary'], self.inputs.reference))
            require(all(owner.context_count == self.context_count and owner.hypothesis_count == self.hypothesis_count
                for owner in (self.orders['primary'], self.identities, self.inputs, self.scalar_facts)))
            closure_ref = bound_reference(reader, primary['source_closure'], 'complete-scientific-source-closure-v1')
            closure_version = reader.reconstruct(closure_ref).get('version')
            require(type(closure_version) is int and closure_version in (1, 2))
            closure = header(reader, closure_ref,
                kind='complete-scientific-source-closure-v1', schema='genosapiens.wgs-scientific-source-closure',
                version=closure_version, keys=('context_count', 'genes', 'phenotypes', 'owners', 'records', 'alleles', 'queries', 'contexts',
                    *(('original_members',) if closure_version == 2 else ())))
            require(integer(closure['context_count']) == self.context_count and closure['genes'] == resources.header['genes']
                and closure['phenotypes'] == resources.header['phenotypes'])
            from genosapiens.wgs_completed_science import child
            if closure_version == 2:
                child(reader, closure['original_members'], 'scientific-original-members-v1')
            for key, kind in (('owners', 'scientific-source-owners-v1'), ('records', 'scientific-source-records-v1'),
                    ('alleles', 'scientific-source-alleles-v1'), ('queries', 'scientific-source-reference-queries-v1'),
                    ('contexts', 'scientific-source-contexts-v1')):
                child(reader, closure[key], kind)
            selected = fields(reader.reconstruct(bound_reference(reader, value['selected'], 'scientific-selected-groups-v1')),
                ('primary', 'selection', 'reference', 'groups'))
            require(selected['primary'] == self.primary_reference)
            require(selected['selection'] == primary['selection'] and _canonical(selected['reference']) == _canonical(resources.reference_binding))
            self.selected_groups = _typed(selected['groups'], tuple[SubmissionGroup, ...], decode=True)
            require(len(_canonical(self.selected_groups)) <= MAX_SELECTED_GROUP_BYTES
                and tuple(g.source_rank for g in self.selected_groups) == tuple(self.selection.selected))
            self.primary_published = self.published = True
            self.sha256 = root.sha256
            for group in self.selected_groups:
                ordered = self.ordered_at(group.source_rank)
                require(group.primary_hypothesis_id == ordered['hypothesis_id'] and group.source_score == ordered['total'])
            analysis_progress(progress, 'publication', 'saved_scientific_projections', unit='contexts',
                completed=self.context_count, total=self.context_count)
        except BaseException:
            self.close()
            raise

    def _restore(self, reader, row, reference, reference_pin, renderer_cache, *, retain_identity_index=False):
        return VerifiedScientificCaseContext(reader, bound_reference(reader, row['context'], 'scientific-case-context-v1'),
            self.resources.restore_for_selected(), reference, reference_pin, renderer_cache,
            retain_identity_index=retain_identity_index)

    def _owns_context_resources(self, resources):
        return resources is self.resources._loaded


class VerifiedScientificCase(CompleteScientificCase):
    def __init__(self, reader, root, resources, reference, reference_pin, renderer_cache, *, progress=None):
        self.progress = progress
        try:
            with reader.read_scope():
                self._restore_case(reader,root,resources,reference,reference_pin,renderer_cache)
        except BaseException:
            if hasattr(self,"stack"):
                self.close()
            raise

    def _restore_case(self, reader, root, resources, reference, reference_pin, renderer_cache):
        if root.kind != "complete-scientific-case-v1":
            raise EvidenceStoreError("foreign complete scientific case")
        value = reader.reconstruct(root)
        version = value.get('version')
        extra = {'context_domain'} if version == 3 else set()
        if (set(value) != {"schema","version","resources","primary","order","selected","assessment_domain_sha256","optional_analysis_status",
                "context_count","positive_context_count","member_count","pair_count","hypothesis_count", *extra}
                or value["schema"] != "genosapiens.wgs-complete-scientific-case"
                or type(value["version"]) is not int or value["version"] not in (2, 3)
                or value["resources"] != resources.reference
                or value["assessment_domain_sha256"] != resources.assessment_domain_sha256
                or value["optional_analysis_status"] != "not_run"):
            raise EvidenceStoreError("invalid complete scientific case closure")
        primary = reader.reconstruct(bound_reference(reader,value["primary"],"scientific-primary-evidence-v1"))
        contexts = bound_reference(reader,primary["contexts"],"scientific-case-contexts-v1")
        domain = None if version == 2 else ContextDomain.from_ordinals(reader.iter_array(
            bound_reference(reader, value['context_domain'], 'scientific-context-domain-v1')))
        super().__init__(EvidenceReplay(reader),resources,name=root.name, context_domain=domain)
        try:
            from genosapiens.wgs_runtime import analysis_progress
            analysis_progress(self.progress, 'publication', 'scientific_replay',
                unit='contexts', completed=0, total=value['context_count'])
            for row in reader.iter_array(contexts,max_item_bytes=MAX_CONTEXT_DIRECTORY_ROW_BYTES):
                with self._restore(reader,row,reference,reference_pin,renderer_cache) as context:
                    self.add_context(context)
                analysis_progress(self.progress, 'publication', 'scientific_replay',
                    unit='contexts', completed=self.context_count, total=value['context_count'])
            self.publish_primary()
            self.publish(reference,reference_pin,renderer_cache)
        except BaseException:
            self.close()
            raise
