"""Canonical active-ingredient union across bounded candidate routes."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from genosapiens.track2_drug_identity import DrugIdentityResolution
from genosapiens.track2_policy import ROUTES


NEGATIVE_STATES = {"no_result", "truncated", "unavailable"}


class CandidateUnionError(ValueError):
    """Raised when candidate identity, rationale, or negative search loses lineage."""


@dataclass(frozen=True, slots=True)
class CandidateRouteObservation:
    identity: DrugIdentityResolution
    objective_id: str
    route: str
    rationale_id: str
    rationale: str
    causal_path: tuple[str, ...]
    evidence_ids: tuple[str, ...]
    lineage_group_ids: tuple[str, ...]

    def __post_init__(self) -> None:
        if (
            len(self.identity.active_ingredient_ids) != 1
            or not self.objective_id
            or self.route not in ROUTES
            or not self.rationale_id
            or not self.rationale
            or not self.causal_path
            or not self.evidence_ids
            or not self.lineage_group_ids
        ):
            raise CandidateUnionError("candidate route observation is incomplete")


@dataclass(frozen=True, slots=True)
class NegativeRouteSearch:
    objective_id: str
    route: str
    state: str
    query_ids: tuple[str, ...]
    sources: tuple[str, ...]
    reason: str

    def __post_init__(self) -> None:
        if (
            not self.objective_id
            or self.route not in ROUTES
            or self.state not in NEGATIVE_STATES
            or not self.query_ids
            or not self.sources
            or not self.reason
        ):
            raise CandidateUnionError("negative route search is incomplete")


@dataclass(frozen=True, slots=True)
class CanonicalCandidate:
    active_ingredient_id: str
    generation_routes: tuple[str, ...]
    rationales: tuple[CandidateRouteObservation, ...]
    independent_rationale_count: int
    identity_source_ids: tuple[str, ...]
    product_ids: tuple[str, ...]
    exact_substance_ids: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class CandidateUnion:
    candidates: tuple[CanonicalCandidate, ...]
    negative_searches: tuple[NegativeRouteSearch, ...]


def _rationale_key(record: CandidateRouteObservation) -> tuple:
    return (
        record.objective_id,
        record.route,
        record.causal_path,
        tuple(sorted(set(record.lineage_group_ids))),
    )


def build_candidate_union(
    observations: Iterable[CandidateRouteObservation],
    negative_searches: Iterable[NegativeRouteSearch] = (),
) -> CandidateUnion:
    records = tuple(observations)
    grouped: dict[str, list[CandidateRouteObservation]] = {}
    for record in records:
        active_id = record.identity.active_ingredient_ids[0]
        grouped.setdefault(active_id, []).append(record)
    candidates: list[CanonicalCandidate] = []
    for active_id, candidate_records in sorted(grouped.items()):
        rationales: dict[tuple, CandidateRouteObservation] = {}
        for record in sorted(candidate_records, key=lambda item: item.rationale_id):
            key = _rationale_key(record)
            existing = rationales.get(key)
            if existing is not None and (
                existing.rationale != record.rationale
                or existing.evidence_ids != record.evidence_ids
            ):
                raise CandidateUnionError(
                    "duplicate identity records disagree on one causal rationale"
                )
            rationales.setdefault(key, record)
        kept = tuple(sorted(rationales.values(), key=lambda item: item.rationale_id))
        candidates.append(
            CanonicalCandidate(
                active_ingredient_id=active_id,
                generation_routes=tuple(sorted({item.route for item in kept})),
                rationales=kept,
                independent_rationale_count=len({lineage for item in kept for lineage in item.lineage_group_ids}),
                identity_source_ids=tuple(
                    sorted(
                        {
                            source_id
                            for item in candidate_records
                            for source_id in item.identity.source_identifiers
                        }
                    )
                ),
                product_ids=tuple(
                    sorted(
                        {
                            product_id
                            for item in candidate_records
                            for product_id in item.identity.product_identifiers
                        }
                    )
                ),
                exact_substance_ids=tuple(
                    sorted(
                        {
                            item.identity.exact_substance_id
                            for item in candidate_records
                            if item.identity.exact_substance_id is not None
                        }
                    )
                ),
            )
        )
    negatives = tuple(
        sorted(negative_searches, key=lambda item: (item.objective_id, item.route, item.state))
    )
    return CandidateUnion(candidates=tuple(candidates), negative_searches=negatives)
