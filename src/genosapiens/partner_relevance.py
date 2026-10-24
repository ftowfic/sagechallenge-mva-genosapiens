"""Pure, context-scoped supplied-ALT relevance, not partner qualification.

The upstream consumer MUST replay ``validate_candidate_admission`` with the
independently held original inputs, source selections and policy pin BEFORE
calling this kernel. Internal consistency and hashes are not authentication.
The consumer must also pin this executable and its dependencies, recompute this
result, and compare its full ``as_bytes()`` before using persisted accounting.

Every inventory record and ALT slot is visited, including non-seeds, duplicate
source representations, and unsupported native sources. Positive association
proofs establish relevance only to their actual source target. Gene annotations
and links are not exclusive target sets. Missing clinical evidence, optional
identifier mappings, or a no-match in one rule cannot exclude a gene context.

The bounded exclusion route is an exact source-native repeat locus, NOT its
optional gene-wide interpretation: catalog ReferenceRegion/MainReferenceRegion
defines that locus, and a replayed repeat proof supplies its reference-coordinate
binding. Disjoint bounded supplied geometry excludes local locus membership.
The closed gene/condition-link adapters do not assert membership of a named
repeat locus. They remain visible, but cannot enlarge a catalog-defined region.
An overlapping record without a locus association stays unresolved. Unknown
source geometry, including a nominal SV without both endpoint bounds, is never
made precise. This is neither molecular contradiction nor biological absence.

No I/O, phase, carriage, genotype, qualification, population, CSV decision,
enumeration-completion flag, or no-partner certificate is produced here.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping
from dataclasses import dataclass
import hashlib
from typing import Literal

from genosapiens import admission_inventory as inventory_api
from genosapiens import admission_links as links, annotation_admission as annotation
from genosapiens import candidate_admission as candidates, repeat_catalog_admission as repeats
from genosapiens.allele_identity import (
    GenomicInterval, MitochondrialAllele, MosaicAllele, RepeatAllele, SequenceAllele,
    StructuralAllele,
)
from genosapiens.candidate_admission import CandidateAdmission
from genosapiens.evidence_model import (
    EvidenceReference, HypothesisContext, IdentityStage, SourcedRecord, _Controlled,
    _validate_context_sources,
)


SCHEMA = "genosapiens.partner-relevance"
POLICY_BYTES = candidates._bytes({
    "schema": SCHEMA,
    "membership": "all-records-and-source-ALT-slots; never-seed-alleles-only",
    "positive": "actual-source-target-association; namespace-aware; Ensembl-stable-ID",
    "gene_exclusion": "no-closed-gene-wide-bounds-or-exclusive-link-set; unresolved",
    "locus_exclusion": "source-defined-repeat-region; reference-bound; all-source-envelopes-disjoint",
    "geometry": "substitution-reference-span; indel-placement-unresolved; explicit-repeat-region; symbolic-SV-both-endpoint-CIs",
    "unknown": "unresolved-source-or-geometry-never-nominally-bounded",
    "non_variant": "replayed-source-non-variant-slot; no-carriage-inference",
    "admission_replay": "mandatory-upstream; internal-consistency-not-authentication",
})
POLICY_SHA256 = hashlib.sha256(POLICY_BYTES).hexdigest()
RelevanceState = Literal["relevant", "excluded_by_evidence", "unresolved_relevance"]


class PartnerRelevanceError(ValueError):
    """Static diagnostic; invalid/missing accounting is not scientific uncertainty."""


@dataclass(frozen=True, slots=True, repr=False)
class RelevanceProofReference(_Controlled):
    """Exact tuple locator and full-content hash within the bound admission."""

    collection: Literal["annotation_proofs", "link_proofs", "repeat_proofs"]
    ordinal: int
    sha256: str


@dataclass(frozen=True, slots=True, repr=False)
class RelevanceRoute(_Controlled):
    subject_id: str
    rule_id: str
    adapter: str
    association_state: Literal["match", "no_match", "unresolved"]
    target_namespace: str | None
    target_id: str | None
    condition_id: str | None
    matches_context_target: bool
    decision_ids: tuple[str, ...]
    proof_ref: RelevanceProofReference | None
    source_refs: tuple[SourcedRecord, ...]
    evidence_refs: tuple[EvidenceReference, ...]
    reasons: tuple[str, ...]


@dataclass(frozen=True, slots=True, repr=False)
class SourceContextRelevance(_Controlled):
    source_record_id: str
    source_allele_id: str
    derived_allele_ids: tuple[str, ...]
    source_status: Literal["derived", "unresolved", "non_variant"]
    unresolved_entry_ids: tuple[str, ...]
    state: RelevanceState
    routes: tuple[RelevanceRoute, ...]
    source_refs: tuple[SourcedRecord, ...]
    evidence_refs: tuple[EvidenceReference, ...]
    geometry_proof_refs: tuple[RelevanceProofReference, ...]
    source_bounds: tuple[GenomicInterval, ...]
    context_bounds: tuple[GenomicInterval, ...]
    reasons: tuple[str, ...]
    binding_sha256: str


@dataclass(frozen=True, slots=True, repr=False)
class ContextRelevanceAssessment(_Controlled):
    context: HypothesisContext
    entries: tuple[SourceContextRelevance, ...]
    binding: Mapping[str, object]

    def as_bytes(self) -> bytes:
        """Controlled/no-LLM audit, never a trusted replay or completion token."""
        return candidates._bytes({"schema": SCHEMA, "classification": "controlled",
                                  "llm_policy": "no-llm", "assessment": self})


_CHECKED_ASSOCIATIONS: dict[int, tuple[CandidateAdmission, tuple, str, str, str]] = {}


def _checked_associations_uncached(admission, context):
    """Reconcile existing source proofs/decisions, without pretending to replay I/O."""
    inv = admission.inventory
    if (admission.binding["inventory_sha256"] != hashlib.sha256(
            inventory_api.admission_inventory_as_json(inv).encode("ascii")).hexdigest()
            or any(admission.binding[key] != value for key, value in (
                ("input_sha256", inv.input_sha256), ("reference_role", inv.reference_role.value),
                ("reference_identity_id", inv.reference_identity_id),
                ("reference_manifest_sha256", inv.reference_manifest_sha256),
                ("reference_concordance_sha256", inv.reference_concordance_sha256)))):
        raise ValueError
    seeds = [s for s in admission.contexts if s.context.context_id == context.context_id]
    if len(seeds) != 1 or seeds[0].context != context:
        raise ValueError
    _validate_context_sources(context, admission.disease_gene_evidence)
    alleles = {a.allele_id: a for a in inv.alleles}
    slots = {s.source.source_allele_id: s for s in inv.source_alleles}
    if len(alleles) != len(inv.alleles) or len(slots) != len(inv.source_alleles):
        raise ValueError
    seen, projected = set(), defaultdict(set)
    for record in inv.records:
        qc, source = record.qc, record.source
        columns = source.source_fields["columns"]
        if (source.source_record_id != qc.source_record_id or source.source_sha256 not in inv.source_input_sha256s
                or source.source_dataset != "supplied-vcf"
                or (columns[0], int(columns[1]), columns[3], tuple(columns[4].split(",")))
                != (qc.contig, qc.position, qc.reference, qc.alternates)
                or len(record.source_allele_ids) != len(qc.alternates)):
            raise ValueError
        for number, slot_id in enumerate(record.source_allele_ids, 1):
            s = slots[slot_id].source
            if (slot_id in seen or s.source_record_id != qc.source_record_id
                    or s.input_sha256 != source.source_sha256 or s.allele_number != number
                    or s.alternate != qc.alternates[number - 1]
                    or (s.record.record_number, s.record.contig, s.record.position,
                        s.record.reference, s.record.alternates)
                    != (qc.record_number, qc.contig, qc.position, qc.reference, qc.alternates)):
                raise ValueError
            seen.add(slot_id)
            for allele_id in slots[slot_id].derived_allele_ids:
                projected[allele_id].add(slot_id)
    if seen != set(slots) or set(projected) != set(alleles):
        raise ValueError
    for key, allele in alleles.items():
        if {i.identity_id for i in allele.identities if i.stage is IdentityStage.SOURCE} != projected[key]:
            raise ValueError
        native = candidates._native(allele)
        identity = native.identity
        if type(identity) in (MosaicAllele, MitochondrialAllele):
            identity = identity.sequence
        if (identity.reference_role != inv.reference_role or identity.reference_identity_id != inv.reference_identity_id
                or native.reference_manifest_sha256 != inv.reference_manifest_sha256):
            raise ValueError
    unresolved = defaultdict(list)
    for entry in inv.unresolved_entries:
        if entry.source_record_id != slots[entry.source_allele_id].source.source_record_id:
            raise ValueError
        unresolved[entry.source_allele_id].append(entry)
    if (len({e.entry_id for e in inv.unresolved_entries}) != len(inv.unresolved_entries)
            or any(s.status == "unresolved" and not unresolved[k] for k, s in slots.items())):
        raise ValueError
    required = {(a, rule) for s in inv.source_alleles for a in s.derived_allele_ids for rule in s.required_rule_ids}
    required.update((e.source_allele_id, rule) for e in inv.unresolved_entries for rule in e.required_rule_ids)
    if len(required) != len(inv.required_rule_keys) or required != set(inv.required_rule_keys):
        raise ValueError
    expected = {(subject, rule, adapter) for subject, rule in required
                for adapter in candidates._RULE_ADAPTERS[rule]}
    registry = {EvidenceReference(p.evidence_id, key, p.public_manifest_sha256): record
                for p in admission.public_sources for key, record in p.records.items()}
    binding_sha = candidates._hash(admission.binding)
    decisions = defaultdict(list)
    for decision in admission.decisions:
        if decision.binding_sha256 != binding_sha:
            raise ValueError
        for source in decision.source_refs:
            if source.public_manifest_sha256 is not None:
                ref = EvidenceReference("admission-public-" + candidates._hash(source),
                                        source.source_record_id, source.public_manifest_sha256)
                if ref not in decision.evidence_refs or registry.get(ref) != source:
                    raise ValueError
        decisions[decision.subject_id, decision.rule_id, decision.adapter].append(decision)
    if set(decisions) != expected or len({candidates._bytes(d) for d in admission.decisions}) != len(admission.decisions):
        raise ValueError

    associations = defaultdict(list)
    used_decisions = set()
    for collection in ("annotation_proofs", "link_proofs", "repeat_proofs"):
        for ordinal, proof in enumerate(getattr(admission, collection)):
            ref = RelevanceProofReference(collection, ordinal, candidates._hash(proof))
            if collection == "repeat_proofs":
                subject = proof.allele_ref.allele_id
                native = candidates._native(alleles[subject])
                if (proof.allele_ref.identity_id != native.identity_id or proof.native_identity != native.identity
                        or proof.input_source != native.source
                        or proof.reference_manifest_sha256 != native.reference_manifest_sha256):
                    raise ValueError
                catalogs = [g for g in admission.disease_gene_evidence
                            if g.public_manifest_sha256 == proof.catalog_manifest_sha256
                            and g.store.source_binding.get("contract_sha256") == proof.catalog_contract_sha256
                            and g.store.evidence_sha256 == proof.catalog_evidence_sha256]
                if len(catalogs) != 1 or proof.adapter_version != catalogs[0].store.source_binding["adapter_version"]:
                    raise ValueError
                catalog = catalogs[0]
                ns = catalog.store.assertions[0].source_family + ":locus" if proof.context else (
                    "Ensembl" if proof.gene_id else None)
                target = proof.locus_id if proof.context else proof.gene_id
                condition = proof.context.native_disease_id if proof.context else None
                state, rule, adapter = proof.association_status, inventory_api.RULE_IDS[3], "repeat"
                sources, evidence = (), proof.evidence_refs
                for evidence_ref in evidence:
                    if (evidence_ref.evidence_id != catalog.evidence_id
                            or evidence_ref.public_manifest_sha256 != catalog.public_manifest_sha256
                            or evidence_ref.record_id not in catalog.store.records):
                        raise ValueError
                if proof.context:
                    _validate_context_sources(proof.context, admission.disease_gene_evidence)
                    if proof.context.locus_id != proof.locus_id or proof.context.gene_id != proof.gene_id:
                        raise ValueError
                if state == "match":
                    records = [catalog.store.records[r.record_id] for r in evidence
                               if catalog.store.records[r.record_id].table == "native_repeat_locus"]
                    if (not records or any(r.source_fields["LocusId"] != proof.locus_id for r in records)
                            or any(repeats._selected_region(r.source_fields)[0] != proof.selected_region for r in records)):
                        raise ValueError
            else:
                subject = proof.allele_id
                if proof.native != candidates._native(alleles[subject]):
                    raise ValueError
                state = {"admitted": "match", "evaluated_no_match": "no_match", "unresolved": "unresolved"}[proof.outcome]
                ns, target = proof.target_namespace, proof.target_id
                sources, evidence = (*proof.source_refs, proof.provenance), ()
                if collection == "annotation_proofs":
                    rule, adapter = proof.rule, "annotation"
                    condition = None
                    if (proof.policy_sha256 != annotation.POLICY_SHA256
                            or proof.provenance.source_dataset != annotation.SCHEMA
                            or proof.provenance.source_sha256 != admission.binding["policy"]["sources"]["annotation"]):
                        raise ValueError
                    raw_targets = [annotation._gene(r.source_fields["csq"]["Gene"])
                                   for r in proof.source_refs if "csq" in r.source_fields]
                    raw_targets += [annotation._gene(r.source_fields["raw_record"].get("gene_id"))
                                    for r in proof.source_refs if "raw_record" in r.source_fields]
                    if target is not None and (target not in raw_targets or ns != "Ensembl"):
                        raise ValueError
                    if state == "match" and (not target or not proof.source_refs):
                        raise ValueError
                else:
                    rule = proof.rule.replace("_", "-")
                    condition = proof.condition_id
                    adapter = {"encode-registry-v4-re2g": "encode", "open-targets-26.06-eva": "open_targets",
                               "clinvar-grch38-2026-08-22-vcf": "clinvar"}[proof.provenance.source_fields["adapter"]]
                    if proof.provenance.source_dataset != links.SCHEMA:
                        raise ValueError
                    if state == "match":
                        if adapter == "encode":
                            element, link = (r.source_fields["record"] for r in proof.source_refs[:2])
                            if (element["legacy_id"] != link["legacy_id"] or target != link["gene_id"]
                                    or ns != "Ensembl"):
                                raise ValueError
                        elif adapter == "open_targets":
                            variant, assertion = (r.source_fields["record"] for r in proof.source_refs[:2])
                            links.parse_open_targets_variant_row(variant)
                            if (links.parse_open_targets_eva_row(assertion) != (target, proof.condition_id)
                                    or assertion["variantId"] != variant["variantId"]):
                                raise ValueError
                        else:
                            raise ValueError
            key = subject, rule, adapter
            if key not in expected:
                raise ValueError
            matching = [d for d in decisions[key]
                        if (d.association_state, d.target_namespace, d.target_id, d.condition_id, d.source_refs)
                        == (state, ns, target, condition, sources) and set(proof.reasons) <= set(d.reasons)
                        and set(evidence) <= set(d.evidence_refs)]
            if not matching:
                raise ValueError
            used_decisions.update(candidates._bytes(d) for d in matching)
            associations[key].append((ref, state, ns, target, condition, sources,
                candidates._unique(r for d in matching for r in d.evidence_refs), proof.reasons,
                tuple(sorted({d.decision_id for d in matching}))))
    selected = admission.binding["policy"]["sources"]
    for key in sorted(expected):
        subject, _, adapter = key
        available = selected["repeat" if adapter == "repeat" else adapter] is not None
        if associations[key]:
            if not available:
                raise ValueError
            continue
        if subject in alleles and available:
            raise ValueError
        reason = "native-source-unresolved" if available else "source-unavailable"
        ds = decisions[key]
        if (len(ds) != 1 or (ds[0].association_state, ds[0].target_namespace, ds[0].target_id,
                            ds[0].condition_id, ds[0].source_refs) != ("unresolved", None, None, None, ())
                or reason not in ds[0].reasons):
            raise ValueError
        used_decisions.add(candidates._bytes(ds[0]))
        associations[key].append((None, "unresolved", None, None, None, (), ds[0].evidence_refs, (reason,), (ds[0].decision_id,)))
    if used_decisions != {candidates._bytes(d) for d in admission.decisions}:
        raise ValueError
    return seeds[0].target_namespace, alleles, slots, unresolved, associations


def _checked_associations(admission, context):
    """Validate one immutable admission once per process object.

    Genotype construction assesses every context from the same freshly replayed,
    frozen ``CandidateAdmission``.  The expensive inventory/proof reconciliation
    is context-independent, so retain its exact result behind a strong object
    reference.  The reference prevents Python object-id reuse, and a distinct or
    deserialized admission is always validated independently.
    """
    cached = _CHECKED_ASSOCIATIONS.get(id(admission))
    if cached is not None and cached[0] is admission:
        seeds = [
            seed
            for seed in admission.contexts
            if seed.context.context_id == context.context_id
        ]
        if len(seeds) != 1 or seeds[0].context != context:
            raise ValueError
        _validate_context_sources(context, admission.disease_gene_evidence)
        common, admission_sha256, binding_sha256, proofs_sha256 = cached[1:]
        return (
            seeds[0].target_namespace,
            *common,
            admission_sha256,
            binding_sha256,
            proofs_sha256,
        )
    namespace, alleles, slots, unresolved, associations = (
        _checked_associations_uncached(admission, context)
    )
    admission_sha256 = hashlib.sha256(admission.as_bytes()).hexdigest()
    binding_sha256 = candidates._hash(admission.binding)
    proofs_sha256 = candidates._hash(
        (
            admission.annotation_proofs,
            admission.link_proofs,
            admission.repeat_proofs,
        )
    )
    common = (alleles, slots, unresolved, associations)
    if id(admission) not in _CHECKED_ASSOCIATIONS and len(_CHECKED_ASSOCIATIONS) >= 4:
        _CHECKED_ASSOCIATIONS.pop(next(iter(_CHECKED_ASSOCIATIONS)))
    _CHECKED_ASSOCIATIONS[id(admission)] = (
        admission,
        common,
        admission_sha256,
        binding_sha256,
        proofs_sha256,
    )
    return (
        namespace,
        *common,
        admission_sha256,
        binding_sha256,
        proofs_sha256,
    )


def _locus_bounds(admission, context, namespace):
    """Use actual matched native/source coordinates, never guess contig aliases."""
    if context.locus_id is None:
        return (), (), ()
    bounds, refs, proofs = [], [], []
    stores = {g.evidence_id: g for g in admission.disease_gene_evidence}
    for binding in context.locus_bindings:
        ref = binding.source_ref
        public = stores[ref.evidence_id]
        if (public.store.source_binding.get("adapter_version") not in (repeats.ADAPTER_VERSION, repeats.SYNTHETIC_ADAPTER_VERSION)
                or binding.field_path != ("LocusId",)):
            return (), (), ()
        record = public.store.records[ref.record_id]
        if record.table != "native_repeat_locus" or record.source_fields["LocusId"] != context.locus_id:
            raise ValueError
        region, _, reason = repeats._selected_region(record.source_fields)
        if reason:
            return (), (), ()
        matched = []
        for ordinal, proof in enumerate(admission.repeat_proofs):
            if (proof.association_status != "match" or proof.context is None
                    or proof.context.locus_id != context.locus_id or ref not in proof.evidence_refs
                    or not proof.reference_binding_sha256):
                continue
            family = public.store.assertions[0].source_family + ":locus"
            if family != namespace or proof.selected_region != region:
                raise ValueError
            matched.append(proof.native_identity.interval)
            proofs.append(RelevanceProofReference("repeat_proofs", ordinal, candidates._hash(proof)))
        if not matched:
            return (), (), ()
        bounds.extend(matched)
        refs.append(ref)
    return candidates._unique(bounds), candidates._unique(refs), candidates._unique(proofs)


def _source_bounds(slot, alleles, unresolved):
    if slot.status != "derived" or unresolved:
        return ()
    bounds = []
    for key in slot.derived_allele_ids:
        native = candidates._native(alleles[key])
        identity = native.identity
        if type(identity) in (MosaicAllele, MitochondrialAllele):
            identity = identity.sequence
        if type(identity) is SequenceAllele:
            raw = slot.source.record
            # Interbase insertions and equivalent repeat-shifted indels need
            # placement evidence beyond nominal VCF reference spans.
            if (len(identity.reference) != len(identity.alternate)
                    or len(raw.reference) != len(slot.source.alternate)):
                return ()
            bounds.append(GenomicInterval(identity.contig, identity.position,
                                          identity.position + len(identity.reference) - 1))
            bounds.append(GenomicInterval(raw.contig, raw.position, raw.position + len(raw.reference) - 1))
        elif type(identity) is RepeatAllele:
            interpretation = native.source.source_fields.get("native_interpretation") or {}
            if interpretation.get("interval_semantics") != "source-format-closed-repeat-region":
                return ()
            bounds.append(identity.interval)
        elif type(identity) is StructuralAllele:
            envelope = inventory_api.supplied_structural_envelope(alleles[key])
            if envelope is None:
                return ()
            bounds.append(envelope)
        else:
            return ()
    return candidates._unique(bounds)


def assess_context_relevance(admission: CandidateAdmission, context: HypothesisContext) -> ContextRelevanceAssessment:
    """Assess every supplied ALT for this exact admission context, without I/O.

    Upstream MUST first replay candidate admission against independent original
    inputs and source pins. A forged, self-consistent source snapshot cannot be
    authenticated by this two-argument kernel. Missing/inconsistent accounting
    raises PartnerRelevanceError; accounted unsupported semantics remain explicit
    unresolved_relevance. Relevant context-only evidence does not assess molecular
    qualification. Consumers own hypothesis construction and completion decisions.
    """
    try:
        if type(admission) is not CandidateAdmission or type(context) is not HypothesisContext:
            raise TypeError
        (
            namespace,
            alleles,
            slots,
            unresolved,
            associations,
            admission_sha256,
            admission_binding_sha256,
            proofs_sha256,
        ) = _checked_associations(admission, context)
        bounds, bound_refs, bound_proofs = _locus_bounds(admission, context, namespace)
        binding = {"schema": SCHEMA, "policy_sha256": POLICY_SHA256,
            "admission_sha256": admission_sha256,
            "admission_binding_sha256": admission_binding_sha256,
            "input_sha256": admission.inventory.input_sha256,
            "inventory_sha256": admission.binding["inventory_sha256"],
            "context_sha256": candidates._hash(context), "target_namespace": namespace,
            "proofs_sha256": proofs_sha256,
            "scope": "supplied-source-ALT-context-relevance-only; upstream-admission-replay-required"}
        binding_sha = candidates._hash(binding)
        by_subject = defaultdict(list)
        for (subject, rule, adapter), values in sorted(associations.items()):
            for ref, state, ns, target, condition, sources, evidence, reasons, decisions in values:
                match = (state == "match" and target is not None and ns == namespace
                         and candidates._gene_join(namespace, target, context.locus_id or context.gene_id))
                by_subject[subject].append(RelevanceRoute(subject, rule, adapter, state, ns, target, condition,
                    match, decisions, ref, sources, evidence, reasons))
        entries = []
        for record in admission.inventory.records:
            for key in record.source_allele_ids:
                slot = slots[key]
                gaps = unresolved[key]
                routes = tuple(r for subject in (key, *slot.derived_allele_ids) for r in by_subject[subject])
                geometry = _source_bounds(slot, alleles, gaps)
                reasons = list(slot.reasons)
                reasons.extend(reason for gap in gaps for reason in gap.reasons)
                if any(r.matches_context_target for r in routes):
                    state = "relevant"
                    reasons.append("source_target_association_context_only_not_qualification")
                elif slot.status == "non_variant":
                    if slot.derived_allele_ids or slot.required_rule_ids or gaps:
                        raise ValueError
                    state = "excluded_by_evidence"
                    reasons.append("replayed_source_non_variant_ALT_slot_not_carriage_assessment")
                elif bounds and geometry and all(
                    a.contig == b.contig and (a.end < b.start or b.end < a.start)
                    for a in geometry for b in bounds
                ):
                    state = "excluded_by_evidence"
                    reasons.extend(("all_source_bounds_disjoint_from_source_defined_native_locus",
                                    "closed_gene_condition_links_do_not_assert_native_locus_membership"))
                else:
                    state = "unresolved_relevance"
                    reasons.append("no_source_proved_context_target_or_exclusion")
                    if not geometry:
                        reasons.append("source_geometry_unavailable_or_not_defensibly_bounded")
                    if not bounds:
                        reasons.append("no_source_supported_complete_context_bounds")
                if gaps:
                    reasons.append("accounted_unresolved_source_semantics")
                entries.append(SourceContextRelevance(record.qc.source_record_id, key, slot.derived_allele_ids,
                    slot.status, tuple(g.entry_id for g in gaps), state, routes,
                    candidates._unique((record.source, *(r for route in routes for r in route.source_refs))),
                    candidates._unique((*context.evidence_refs, *bound_refs,
                                        *(r for route in routes for r in route.evidence_refs))),
                    bound_proofs, geometry, bounds, tuple(sorted(set(reasons))), binding_sha))
        return ContextRelevanceAssessment(context, tuple(entries), binding)
    except (ValueError, TypeError, KeyError, IndexError, AttributeError, RuntimeError, StopIteration, OverflowError):
        raise PartnerRelevanceError("invalid or incomplete context relevance inputs") from None
