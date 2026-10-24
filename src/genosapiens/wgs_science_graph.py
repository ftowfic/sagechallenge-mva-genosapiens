"""Scientific-only graph preparation with indexed categorical support.

The upstream candidate admission producer remains unchanged. Preparation and
worker wrappers preserve its checks and transport. Tests compare the optimized
support kernel and complete graph against that original producer.
"""
from __future__ import annotations
from genosapiens.wgs_validation_scope import member_validation_scope
from collections import deque
from types import MappingProxyType
from genosapiens.candidate_admission import (
    _relation_edges, _PreparedAdmissionPhenotype, PhenotypeAdmission,
    _reach as _ordinary_reach,
)
from genosapiens.evidence_model import EvidenceReference
from genosapiens.candidate_admission import _record_bytes, _unique as _ordinary_unique

from genosapiens.candidate_admission import (
    AssociationTarget, _association_reaches_context, _bytes, _requirements_disagree,
    CandidateAdmissionError, DiseaseSupportEntry, HypothesisContext, LocusBinding, PreparedShardedResources, _PREPARED_ADMISSION_SEAL, _PREPARED_SHARDED_SEAL, _PreparedAdmissionGraph, _PreparedShardedAdmission, _freeze_policy, _hash, _prepare_admission_evidence, _prepare_admission_phenotype, _ref, _requirement_key, _target, _unique, assertion_native_disease_id, closing, defaultdict, hashlib, json, prepare_sharded_resources, re, replace, sharded_admission_policy_bytes
)
from genosapiens.wgs_admission_graph import (
    GraphSupportError, MAX_JOB_BYTES, MAX_JOB_CONTEXTS, MAX_RESULT_BYTES, _Pool, _job_packet, _jobs, _pack, hashlib, live_worker_count, os, pickle, run_cleanup, static_failure, sys, time, worker_failure_message
)


def _unique(items):
    """Deduplicate immutable reference values before canonical serialization.

    Three exact string fields fully define this record. Other inputs retain
    the ordinary path. Canonical byte keys still determine the final order.
    No input or serialization result survives this one operation.
    """
    if type(items) in (tuple, list) and all(type(item) is EvidenceReference
            and all(type(value) is str for value in
                (item.evidence_id, item.record_id, item.public_manifest_sha256)) for item in items):
        values = {(item.evidence_id, item.record_id, item.public_manifest_sha256): item for item in items}
        return tuple(sorted(values.values(), key=_record_bytes))
    return _ordinary_unique(items)


class _ReachIndex:
    """One case-local encoding of relationship edges, with exact proof values."""
    def __init__(self, edges):
        self.edges = edges
        self.prepared = {}
        self.identities = {}
        self.keys = []

    def get(self, node):
        if node not in self.prepared:
            rows = []
            for target, step, proof in sorted(self.edges[node], key=_bytes):
                records = []
                for record in proof:
                    key = _record_bytes(record)
                    identity = self.identities.get(key)
                    if identity is None:
                        identity = len(self.keys)
                        self.identities[key] = identity
                        self.keys.append(key)
                    records.append((identity, record))
                rows.append((target, step, tuple(records), frozenset(i for i, _ in records)))
            self.prepared[node] = tuple(rows)
        return self.prepared[node]


def _reach(anchor, edges):
    """Preserve ordinary queue order and shortest-path proof convergence."""
    if type(edges) is not _ReachIndex:
        return _ordinary_reach(anchor, edges)
    records, visited = {}, set()
    reached, queue = {anchor: (0, frozenset())}, deque([anchor])
    while queue:
        node = queue.popleft()
        distance, refs = reached[node]
        rows = edges.get(node)
        if node not in visited:
            for _, _, values, _ in rows:
                records.update(values)
            visited.add(node)
        for target, step, _, proof in rows:
            candidate = distance + step, refs | proof
            if target not in reached or candidate[0] < reached[target][0]:
                reached[target] = candidate
                queue.append(target)
            elif candidate[0] == reached[target][0]:
                combined = reached[target][1] | candidate[1]
                if combined != reached[target][1]:
                    reached[target] = candidate[0], combined
                    queue.append(target)
    return {node: (distance, tuple(records[i] for i in sorted(refs, key=edges.keys.__getitem__)))
            for node, (distance, refs) in reached.items()}


def _producer_record(cls, values):
    """Publish locally computed fields without repeating model verification."""
    result = object.__new__(cls)
    for name, value in zip(cls.__dataclass_fields__, values, strict=True):
        object.__setattr__(result, name, value)
    return result


def _producer_replace(value, **changes):
    return _producer_record(type(value), tuple(changes.get(name, getattr(value, name))
        for name in value.__dataclass_fields__))


def _prepare_live_admission_evidence(*, view, profiles, genes, ic_store,
        repeat_catalog, authenticated_resources=None):
    """Apply ordinary source selection while omitting resource verification."""
    evidence = list(genes or ())
    if repeat_catalog is not None:
        catalog = repeat_catalog.catalog
        if not any(item.evidence_id == catalog.evidence_id for item in evidence):
            evidence.append(catalog)
    return view, tuple(sorted((item if not item.aggregations else
        _producer_record(type(item), tuple(() if name == 'aggregations' else getattr(item, name)
            for name in item.__dataclass_fields__)) for item in evidence),
        key=lambda item: item.evidence_id))


def _prepare_admission_phenotype(view, profiles, evidence, ic_store):
    """Case-local producer; all members are immutable and bound to source owners."""
    observations = defaultdict(list)
    limits = []
    for o in view.observations:
        if o.status != "present":
            continue
        resolution = ic_store.graph.resolve(o.source_hpo_id) if ic_store is not None else None
        if resolution is not None and resolution.canonical_id:
            observations[resolution.canonical_id].append(o)
        else:
            limits.append("unresolved-positive-observation")
    if ic_store is None:
        limits.append("ic-source-unavailable")
    if profiles is None:
        limits.append("profile-source-unavailable")
    if not observations:
        limits.append("no-resolved-positive-observation")
    reach_cache, edges = {}, _ReachIndex(_relation_edges(evidence))
    anchors = {assertion_native_disease_id(g.store, a) for g in evidence for a in g.store.assertions}
    annotations = profiles.store.annotations if profiles else ()
    anchors.update(a.source_fields["source_disease_id"] for a in annotations)
    for anchor in sorted(anchors - {None}):
        reach_cache[anchor] = _reach(anchor, edges)
    pheno, profile_status, matches = [], defaultdict(list), {}
    grouped = defaultdict(list)
    for a in annotations:
        native = a.source_fields["source_disease_id"]
        if native and a.source_fields["aspect"] == "P":
            profile_status[native].append(a)
        if not native or a.polarity != "positive" or a.source_fields["aspect"] != "P" or not a.resolution or not a.resolution.canonical_id:
            continue
        for term in sorted(observations):
            pair = (term, a.resolution.canonical_id)
            if pair not in matches:
                matches[pair] = ic_store.resnik(*pair)
            match = matches[pair]
            if match.information_content is not None and match.information_content > 0:
                grouped[native, *pair].append(a)
    ids_by_term = {term: tuple(sorted(o.observation_id for o in rows))
        for term, rows in observations.items()}
    contradictions = {term: tuple(sorted(c.contradiction_id for c in view.contradictions
        if set(c.observation_ids).intersection(ids))) for term, ids in ids_by_term.items()}
    refs = {a.annotation_id: _producer_record(EvidenceReference,
        (profiles.evidence_id, a.annotation_id, profiles.public_manifest_sha256))
        for a in annotations}
    profile_refs = {}
    for (native, term, target), members in sorted(grouped.items()):
        match = matches[term, target]
        ids = ids_by_term[term]
        reference_key = native, target
        if reference_key not in profile_refs:
            profile_refs[reference_key] = _unique(tuple(refs[a.annotation_id] for a in members))
        proof_hash = hashlib.sha256((json.dumps([native, term, target, match.mica_ids,
            match.information_content], sort_keys=True, separators=(',', ':'),
            ensure_ascii=True, allow_nan=False) + '\n').encode('ascii')).hexdigest()
        pheno.append(_producer_record(PhenotypeAdmission, ("phenotype-" + proof_hash,
            native, term, target, match.mica_ids, match.information_content, ids,
            contradictions[term],
            profile_refs[reference_key], "unresolved", (), ())))
    return _PreparedAdmissionPhenotype(
        view, profiles, evidence, ic_store,
        MappingProxyType({key: tuple(value) for key, value in observations.items()}),
        tuple(limits),
        MappingProxyType({key: MappingProxyType(value) for key, value in reach_cache.items()}),
        MappingProxyType({key: tuple(value) for key, value in profile_status.items()}),
        MappingProxyType(matches), tuple(pheno), _PREPARED_ADMISSION_SEAL,
    )


def prepare_sharded_candidate_admission(
    *, sources: ShardedAdmissionSources, view: ScoringView,
    profiles: DiseaseProfile | None, genes: tuple[DiseaseGeneEvidence, ...] | None,
    ic_store: HpoICStore | None, reached_targets: tuple[AssociationTarget, ...],
    clinical_index_manifest_sha256: str | None = None,
    execution_config=None, graph_progress_path=None, resource_preparation=None, progress=None, authenticated_resources=None,
) -> _PreparedShardedAdmission:
    """Prepare fixed policy and phenotype relationships once for this case/view.

    The caller owns this bounded-lifetime object. It holds strong references to
    its immutable resource owners and is never shared through a global cache.
    """
    from genosapiens.wgs_runtime import analysis_progress
    from genosapiens.wgs_failures import WgsFailure
    analysis_progress(progress, 'context_assembly', 'resource_prepare')
    try:
        if resource_preparation is None:
            resource_preparation=prepare_sharded_resources(sources)
        if (type(resource_preparation) is not PreparedShardedResources
                or resource_preparation.seal is not _PREPARED_SHARDED_SEAL
                or resource_preparation.sources is not sources):
            raise ValueError
        source_view = view
        policy = sharded_admission_policy_bytes(
            sources=sources, profiles=profiles, genes=genes, ic_store=ic_store,
            clinical_index_manifest_sha256=clinical_index_manifest_sha256,
        )
        view, evidence = _prepare_live_admission_evidence(
            view=view, profiles=profiles, genes=genes, ic_store=ic_store,
            repeat_catalog=sources.repeat_catalog, authenticated_resources=authenticated_resources,
        )
        phenotype = _prepare_admission_phenotype(view,profiles,evidence,ic_store)
        return _PreparedShardedAdmission(
            sources, source_view, view, profiles, genes, ic_store,
            clinical_index_manifest_sha256, policy,
            _freeze_policy(json.loads(policy)), hashlib.sha256(policy).hexdigest(),
            _hash(view), evidence,
            phenotype,_prepare_admission_graph(phenotype, reached_targets=reached_targets,execution_config=execution_config,
                progress_path=graph_progress_path,progress=progress),resource_preparation,
            frozenset(_bytes(value) for value in reached_targets),
            _PREPARED_SHARDED_SEAL,
        )
    except WgsFailure:
        raise
    except (ValueError, TypeError, RuntimeError, OSError, KeyError, AttributeError):
        raise CandidateAdmissionError("invalid sharded admission preparation") from None


def _prepare_admission_graph(prepared, *, reached_targets, execution_config=None, progress_path=None, progress=None):
    """Close only contexts reached by positive assessed source associations."""
    if type(reached_targets) is not tuple or any(type(v) is not AssociationTarget for v in reached_targets):
        raise CandidateAdmissionError('invalid graph target closure')
    from genosapiens.wgs_runtime import analysis_progress
    processed = 0
    analysis_progress(progress, 'context_assembly', 'candidate_construction', counts={'candidates_processed':0})
    observations,limits,reach_cache = prepared.observations,prepared.limits,prepared.reach_cache
    profile_status,matches,pheno = prepared.profile_status,prepared.matches,prepared.pheno
    evidence,ic_store = prepared.evidence,prepared.ic_store
    # Each candidate alternative is anchored in actual source requirements. No
    # union-find over partial models: an unknown cannot bridge AD and AR.
    candidates, candidate_origins = {}, defaultdict(set)
    phenotype_anchors, phenotype_by_anchor = defaultdict(set), defaultdict(list)
    for p in pheno:
        for native in reach_cache[p.native_disease_id]:
            phenotype_anchors[native].add(p.native_disease_id)
            phenotype_by_anchor[native].append(p)
    support_index = defaultdict(lambda: defaultdict(list))
    public_by_id = {g.evidence_id: g for g in evidence}
    for public in evidence:
        for a in public.store.assertions:
            with member_validation_scope():
                native = assertion_native_disease_id(public.store, a)
                namespace, target = _target(public, a)
                if not native or not target:
                    continue
                support_index[namespace, target][public.evidence_id].append(a)
                possible = {native, *phenotype_anchors[native]}
                for anchor in sorted(possible):
                    relation_refs = reach_cache[anchor][native][1]
                    for requirement in a.requirements:
                        key = "context-" + _hash([anchor, namespace, target, _requirement_key(requirement)])
                        refs = _unique((*(_ref(public, k) for k in (a.record_id, *requirement.record_ids)), *relation_refs))
                        locus_bindings = (LocusBinding(_ref(public, a.locus_record_id), ("LocusId",)),) if a.locus_id else ()
                        context = HypothesisContext(key, anchor, anchor, a.gene_id, a.locus_id, requirement.inheritance,
                            "source-requirement-context-only", (requirement,), refs, locus_bindings)
                        if key in candidates:
                            old, _ = candidates[key]
                            context = _producer_replace(old, requirements=_unique((*old.requirements, requirement)),
                                evidence_refs=_unique((*old.evidence_refs, *refs)),
                                locus_bindings=_unique((*old.locus_bindings, *locus_bindings)))
                        candidates[key] = (context, namespace)
                        processed += 1
                        if processed % 256 == 0:
                            analysis_progress(progress, 'context_assembly', 'candidate_construction',
                                counts={'candidates_processed':processed}, unit='candidates', completed=processed)
                        if anchor == native:
                            candidate_origins[key].add("variant")
    targets_by_gene,targets_by_condition=defaultdict(list),defaultdict(list)
    for target in reached_targets:
        if target.target is not None:
            targets_by_gene[target.namespace,target.target].append(target)
        else:
            targets_by_condition[target.condition].append(target)
    reached={}
    for key,(context,namespace) in candidates.items():
        if 'variant' not in candidate_origins[key]:
            continue
        target=context.locus_id or context.gene_id
        stable=target.split('.')[0] if namespace=='Ensembl' else target
        possible=(*targets_by_gene.get((namespace,stable),()),
            *(value for native in reach_cache[context.native_disease_id]
                for value in targets_by_condition.get(native,())))
        if any(_association_reaches_context(value,context,namespace,reach_cache) for value in possible):
            reached[key]=(context,namespace)
    candidates=reached
    analysis_progress(progress, 'context_assembly', 'candidate_construction',
        counts={'candidates_processed':processed, 'candidates_retained':len(candidates)},
        unit='candidates', completed=processed)
    ledgers, supported, phenotype_state, context_pheno = {}, {}, {}, {}
    progress_callback = progress
    progress = None
    if execution_config is not None:
        from genosapiens.bounded_wgs_config import BoundedWgsExecutionConfig
        from genosapiens.wgs_admission_graph import GraphProgress
        if type(execution_config) is not BoundedWgsExecutionConfig or execution_config.public_queries.concurrent_batches != 10:
            raise CandidateAdmissionError('invalid graph execution configuration')
        progress = GraphProgress(progress_path, len(candidates), progress=progress_callback)
        with closing(parallel_support(candidates, support_index, public_by_id, reach_cache, progress=progress)) as results:
            for key, ledger in results:
                if key not in candidates or key in ledgers:
                    raise CandidateAdmissionError('invalid graph support coverage')
                ledgers[key] = ledger
        if len(ledgers) != len(candidates):
            raise CandidateAdmissionError('incomplete graph support')
        ledgers = {key: ledgers[key] for key in sorted(candidates)}
    else:
        for key, (context, namespace) in sorted(candidates.items()):
            reached = reach_cache[context.native_disease_id]
            indexed = tuple((public_by_id[public_id], assertions) for public_id, assertions in
                            support_index[namespace, context.locus_id or context.gene_id].items())
            ledgers[key] = _support(context, namespace, indexed, reached)

    # Requirements and gene/locus alternatives do not affect phenotype state.
    # Cache the exact former expression once per native disease in this graph.
    state_by_native = {}
    for key, (context, namespace) in sorted(candidates.items()):
        supported[key] = any(e.positive for e in ledgers[key])
        # A native profile match can reach this native gene context by an actual
        # directed/equivalence proof, never by sharing an optional mapped view.
        native_id = context.native_disease_id
        if native_id not in state_by_native:
            reached = reach_cache[native_id]
            proofs = tuple(phenotype_by_anchor[native_id])
            rows = [a for native in reached for a in profile_status.get(native, ())]
            informative_rows = [a for a in rows if a.polarity == "positive" and a.resolution and a.resolution.canonical_id]
            unknown_rows = any(a.polarity in ("unknown", "unresolved") or not a.resolution or not a.resolution.canonical_id for a in rows)
            unknown_comparison = any(matches[term, a.resolution.canonical_id].information_content is None
                                     for term in observations for a in informative_rows)
            state = "present" if proofs else "unknown" if (
                not informative_rows or not observations or ic_store is None or unknown_rows or unknown_comparison
                or "unresolved-positive-observation" in limits) else "absent"
            state_by_native[native_id] = (proofs, state)
        context_pheno[key], phenotype_state[key] = state_by_native[native_id]
    phenotype_contexts = defaultdict(list)
    for key in candidates:
        if supported[key]:
            for p in context_pheno[key]:
                phenotype_contexts[p.proof_id].append(key)
    pheno = tuple(_producer_replace(p, context_ids=tuple(sorted(phenotype_contexts[p.proof_id])),
        evidence_refs=p.evidence_refs if not phenotype_contexts[p.proof_id] else _unique((*p.evidence_refs, *(ref for key in phenotype_contexts[p.proof_id]
            for ref in reach_cache[p.native_disease_id][candidates[key][0].native_disease_id][1]))),
        outcome="admitted" if phenotype_contexts[p.proof_id] else "unresolved",
        reasons=("informative-positive-phenotype-association", *(() if phenotype_contexts[p.proof_id]
                  else ("no-applicable-positive-disease-support",)))) for p in pheno)
    candidate_by_target, candidate_by_condition = defaultdict(list), defaultdict(list)
    for key, (context, namespace) in candidates.items():
        if "variant" not in candidate_origins[key]:
            continue
        target = context.locus_id or context.gene_id
        stable = target.split(".")[0] if namespace == "Ensembl" else target
        candidate_by_target[namespace, stable].append(key)
        for native in reach_cache[context.native_disease_id]:
            candidate_by_condition[native].append(key)
    seed_contexts = {}
    by_context = defaultdict(list)
    for proof in pheno:
        for key in proof.context_ids:
            by_context[key].append(proof)
    for key,(context,namespace) in candidates.items():
        proofs = by_context.get(key,())
        refs = _unique((*context.evidence_refs, *(e.source_ref for e in ledgers[key]),
            *(r for e in ledgers[key] for r in e.relation_refs),
            *(r for p in proofs for r in p.evidence_refs)))
        seed_contexts[key] = (_producer_replace(context,evidence_refs=refs),namespace,
            phenotype_state[key],tuple(sorted(p.proof_id for p in proofs)),ledgers[key])
    graph = _PreparedAdmissionGraph(prepared,_freeze_policy(candidates),_freeze_policy(ledgers),
        _freeze_policy(supported),_freeze_policy(phenotype_state),pheno,
        _freeze_policy(dict(candidate_by_target)),_freeze_policy(dict(candidate_by_condition)),
        _freeze_policy(seed_contexts),_PREPARED_ADMISSION_SEAL)
    if progress is not None:
        progress.write('complete', force=True)
    return graph


@member_validation_scope()
def _support(context, namespace, evidence, reached):
    """Task 7.4 categorical applicability only: no weights, caps or penalties."""
    entries, details = [], []
    selected = context.requirements[0]
    for public, assertions in evidence:
        for assertion in assertions:
            ns, target = _target(public, assertion)
            if ns != namespace or target != (context.locus_id or context.gene_id):
                continue
            native = assertion_native_disease_id(public.store, assertion)
            distance, relation_refs = reached.get(native, (None, ()))
            in_view = assertion.record_id in public.store.views.get(context.disease_id, ())
            if distance is None and not in_view:
                continue
            row = public.store.records[assertion.record_id].source_fields
            compatible = tuple(r for r in assertion.requirements if not _requirements_disagree(r, selected))
            matching = any(r.inheritance is not None and r.inheritance == selected.inheritance for r in compatible)
            # Unknown requirements do not contribute to a different known alternative.
            exact_model = any(_requirement_key(r) == _requirement_key(selected) for r in compatible)
            applicability = "unknown" if distance is None else "exact" if native == context.native_disease_id else (
                "equivalent" if distance == 0 else "broader")
            scope = "native_locus" if assertion.locus_id is not None else "indirect" if row.get("is_ontology_propagated") is True else (
                "broad" if distance and distance > 0 or row.get("disease_specificity") == "broad" else
                "direct" if row.get("association_scope") in ("source_direct", "ot_mapped_evidence")
                and row.get("is_ontology_propagated") is False and row.get("disease_specificity") in ("subtype", "disorder") else "contextual")
            direct = scope == "direct" and distance == 0 and assertion.evidence_kind == "curated"
            reasons = list(assertion.diagnostics)
            if not compatible:
                reasons.append("different_inheritance_or_requirements")
            elif not exact_model:
                reasons.append("separate_requirement_alternative")
            if distance is None:
                reasons.append("mapped_view_not_equivalence_applicability_unknown")
            entries.append(DiseaseSupportEntry(_ref(public, assertion.record_id), assertion.polarity, applicability,
                scope, "different" if not compatible else "exact_alternative" if exact_model else "unknown_or_other_alternative",
                distance is not None and exact_model and assertion.polarity in ("positive", "contextual"),
                relation_refs, tuple(reasons), (), ""))
            details.append((assertion, row, distance, direct, matching, compatible))
    specific = [i for i, (a, _, _, direct, matching, _) in enumerate(details)
                if direct and matching and a.polarity in ("positive", "negative")]
    opposing = defaultdict(list)
    curated_families = {details[j][0].source_family for j, entry in enumerate(entries)
        if details[j][0].evidence_kind == 'curated' and entry.applicability in ('exact', 'equivalent') and entry.positive}
    for j, (assertion, _, distance, *_rest) in enumerate(details):
        if assertion.polarity in ('positive', 'negative') and distance is not None:
            opposing[assertion.polarity, distance].append(j)
    for i, entry in enumerate(entries):
        a, row, distance, _, matching, compatible = details[i]
        reasons = list(entry.reasons)
        positive = entry.positive
        if specific and distance is not None and distance > 0:
            positive = False
            reasons.append("umbrella_evidence_superseded")
            if any(details[j][0].polarity != a.polarity for j in specific):
                reasons.append("umbrella_validity_conflict")
            if not compatible:
                reasons.append("umbrella_inheritance_conflict")
        conflicts = []
        opposite = 'negative' if a.polarity == 'positive' else 'positive' if a.polarity == 'negative' else None
        for j in opposing.get((opposite, distance), ()):
            peer = entries[j]
            b, other_row, other_distance, _, other_matching, other_requirements = details[j]
            if ({a.polarity, b.polarity} == {"positive", "negative"}
                    and distance is not None and distance == other_distance):
                known_scope = all(assertion.evidence_kind == "curated"
                    and source.get("is_ontology_propagated") is False
                    and source.get("association_scope") in ("source_direct", "ot_mapped_evidence")
                    and source.get("disease_specificity") in ("subtype", "disorder", "broad")
                    for assertion, source in ((a, row), (b, other_row)))
                if not compatible or not other_requirements:
                    reasons.append("opposing_assertion_has_incompatible_requirements")
                elif not known_scope or not matching or not other_matching:
                    reasons.append("opposing_assertion_scope_or_requirements_unknown")
                elif any(not _requirements_disagree(left, right) for left in compatible for right in other_requirements):
                    conflicts.append(peer.source_ref)
                    reasons.append("equal_specificity_validity_conflict_unpenalized")
                else:
                    reasons.append("opposing_assertion_has_incompatible_requirements")
        if a.evidence_kind == "aggregate" and a.source_family in curated_families:
            reasons.append("aggregate_source_already_represented")
            # Keep an unknown alternative explicit; no numeric duplicate credit exists here.
        entries[i] = _producer_replace(entry, positive=positive, reasons=tuple(sorted(set(reasons))), conflict_refs=_unique(conflicts))
    # Unscored dependence groups retain every record. Unknown models never merge contexts.
    neighbours = defaultdict(set)
    for i, (a, row, *_rest) in enumerate(details):
        tokens = {"family:" + a.source_family}
        for publication in row.get("publications") or ():
            if publication and (match := re.fullmatch(r"(?:PMID:)?([1-9][0-9]*)", publication, re.I)):
                tokens.add("PMID:" + match[1])
        for token in tokens:
            neighbours[token].add(i)
    parents = list(range(len(entries)))
    sizes = [1] * len(entries)
    def root(index):
        while parents[index] != index:
            parents[index] = parents[parents[index]]
            index = parents[index]
        return index
    for group in neighbours.values():
        first = next(iter(group))
        for index in group:
            a, b = root(first), root(index)
            if a != b:
                if sizes[a] < sizes[b]:
                    a, b = b, a
                parents[b] = a
                sizes[a] += sizes[b]
    components, tokens = defaultdict(list), defaultdict(list)
    for index in range(len(entries)):
        components[root(index)].append(index)
    for token, group in neighbours.items():
        tokens[root(next(iter(group)))].append(token)
    for representative, members in components.items():
        group_id = "dependency-" + _hash(sorted(tokens[representative]))
        for i in members:
            entries[i] = _producer_replace(entries[i], dependency_group=group_id)
    return tuple(sorted(entries, key=lambda e: _bytes(e.source_ref)))


def _worker(index, jobs, results):
    job = digest = None
    operation = 'graph_worker_execute'
    try:
        results.put(('ready', index, os.getpid()))
        while True:
            idle_started=time.monotonic_ns()
            packet = jobs.get()
            idle=time.monotonic_ns()-idle_started
            if packet is None:
                results.put(('stopped', index))
                return
            job, payload = _job_packet(packet)
            digest = hashlib.sha256(payload).hexdigest()
            operation = 'graph_worker_execute'
            if type(payload) is not bytes or len(payload) > MAX_JOB_BYTES:
                raise GraphSupportError('graph input is invalid')
            indexed, rows = pickle.loads(payload)
            if type(rows) is not tuple or not 0 < len(rows) <= MAX_JOB_CONTEXTS:
                raise GraphSupportError('graph input is invalid')
            started = time.process_time_ns()
            wall_started=time.monotonic_ns()
            output = tuple((key, _support(context, namespace, indexed, reached))
                           for key, context, namespace, reached in rows)
            operation = 'graph_result_validate'
            message = _pack(output, MAX_RESULT_BYTES, operation=operation)
            results.put(('result', index, hashlib.sha256(payload).hexdigest(),
                         message, time.process_time_ns() - started,time.monotonic_ns()-wall_started,idle,job))
            del payload, indexed, rows, output, message
    except BaseException as error:
        results.put(worker_failure_message(error, worker=index, job=job, digest=digest, operation=operation))


def parallel_support(candidates, support_index, public_by_id, reach_cache, *, progress):
    pool = _Pool(worker_count=live_worker_count(len(candidates)), worker_target=_worker)
    pool.progress_callback = lambda: progress.write('support')
    jobs = _jobs(candidates, support_index, public_by_id, reach_cache)
    progress.write('worker_startup', force=True)
    try:
        pool.start()
        progress.metrics = pool.metrics
        for index in range(pool.worker_count):
            payload = next(jobs, None)
            if payload is None:
                break
            pool.submit(index, payload)
        progress.write('support', force=True)
        while pool.pending:
            index, rows = pool.receive()
            yield from rows
            progress.completed += len(rows)
            progress.write('support')
            payload = next(jobs, None)
            if payload is not None:
                pool.submit(index, payload)
        if progress.completed != len(candidates):
            raise GraphSupportError(category='malformed_evidence', operation=pool.result_operation)
        pool.finish()
        progress.write('assembling_graph', force=True)
    except BaseException as error:
        primary = error if isinstance(error, (KeyboardInterrupt, GeneratorExit)) else static_failure(error, 'bounded_products', operation='graph_worker_execute')
        run_cleanup([('diagnostic_publish', lambda: progress.write('failed', force=True))], primary=primary)
        raise primary from None
    finally:
        run_cleanup([('worker_shutdown', jobs.close), ('worker_shutdown', pool.close)], primary=sys.exc_info()[1])
