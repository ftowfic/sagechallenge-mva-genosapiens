"""Lossless factorization of evaluated phase bases by member source positions.

Every phase, basis, namespace, block and limit comes from an already evaluated
PairPhaseProperties. Expansion only substitutes exact source identities; it does
not infer phase, drop bases or replay science in place of a computed outcome.
"""
from dataclasses import dataclass
from typing import Literal

from genosapiens.evidence_model import HypothesisPhase, _Controlled, _ids, _typed
from genosapiens.evidence_serialization import _canonical
from genosapiens.genotype_properties import PairPhaseBasis, PairPhaseProperties
from genosapiens.wgs_evidence_store import EvidenceStoreError


@dataclass(frozen=True, slots=True, repr=False)
class PhaseBasisPattern(_Controlled):
    source_slots: tuple[tuple[int, int], tuple[int, int]]
    relation: HypothesisPhase
    rule: Literal["same-record-distinct-copies", "shared-phase-block", "source-genotype-conflict"]
    namespace: str | None
    block: str | None

    def __post_init__(self):
        _Controlled.__post_init__(self)
        if any(side not in (0,1) or index < 0 for side,index in self.source_slots):
            raise EvidenceStoreError("invalid phase source slot")


@dataclass(frozen=True, slots=True, repr=False)
class PairPhasePattern(_Controlled):
    phase: HypothesisPhase
    bases: tuple[PhaseBasisPattern, ...]
    limits: tuple[str, ...]
    version: Literal[1] = 1

    def __post_init__(self):
        _Controlled.__post_init__(self)
        if self.phase is HypothesisPhase.NOT_APPLICABLE:
            raise EvidenceStoreError("pair phase cannot be not applicable")


def _sources(left, right):
    for values in (left,right):
        _typed(values,tuple[str,...])
        _ids(values)
    return left,right


def phase_pattern(value, *, left_sources, right_sources):
    value = _typed(value,PairPhaseProperties)
    sources = _sources(left_sources,right_sources)
    slots = {}
    for side,values in enumerate(sources):
        for index,key in enumerate(values):
            # Shared source views retain the same literal identity. Choose its
            # first deterministic address; no biological equivalence is added.
            slots.setdefault(key,(side,index))
    try:
        result = PairPhasePattern(value.phase,tuple(PhaseBasisPattern(
            tuple(slots[key] for key in basis.source_allele_ids),basis.relation,basis.rule,basis.namespace,basis.block)
            for basis in value.bases),value.limits)
    except KeyError:
        raise EvidenceStoreError("phase basis has no member source owner") from None
    if expand_phase_pattern(result,allele_ids=value.allele_ids,
            left_sources=left_sources,right_sources=right_sources) != value:
        raise EvidenceStoreError("phase pattern did not preserve evaluated evidence")
    return result


def expand_phase_pattern(pattern, *, allele_ids, left_sources, right_sources):
    pattern = _typed(pattern,PairPhasePattern)
    sources = _sources(left_sources,right_sources)
    _typed(allele_ids,tuple[str,str])
    _ids(allele_ids)
    if allele_ids[0] >= allele_ids[1]:
        raise EvidenceStoreError("phase pattern pair order changed")
    try:
        bases = tuple(PairPhaseBasis(tuple(sources[side][index] for side,index in basis.source_slots),
            basis.relation,basis.rule,basis.namespace,basis.block) for basis in pattern.bases)
    except IndexError:
        raise EvidenceStoreError("phase pattern source slot is unresolved") from None
    if len(set(bases)) != len(bases):
        raise EvidenceStoreError("phase pattern source uniqueness changed")
    # Source IDs differ across pairs sharing a slot pattern. Their original
    # canonical order is restored after substitution, without dropping a basis.
    return PairPhaseProperties(allele_ids,pattern.phase,tuple(sorted(bases,key=_canonical)),pattern.limits)
