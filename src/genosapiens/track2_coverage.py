"""Coverage records that distinguish bounded searches from true negatives."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import date
import hashlib
import json
from typing import Iterable

from genosapiens.track2_policy import ROUTES


COVERAGE_STATES = {
    "analysis_incomplete_due_to_envelope",
    "identifier_mismatch",
    "invalid",
    "no_result",
    "not_assessed",
    "not_assessed_due_to_frozen_envelope",
    "screened_exclusion",
    "success",
    "truncated",
    "unavailable",
}
ENVELOPE_LIMITED_STATES = {
    "analysis_incomplete_due_to_envelope",
    "not_assessed_due_to_frozen_envelope",
}
SEARCH_STATES = {
    "identifier_mismatch",
    "invalid",
    "no_result",
    "screened_exclusion",
    "success",
    "truncated",
}


class CoverageError(ValueError):
    """Raised when coverage overstates what a bounded envelope established."""


@dataclass(frozen=True, slots=True)
class CoverageBound:
    maximum_results: int
    results_returned: int
    screened_exclusions: int
    graph_depth: int
    date_start: str
    date_end: str
    source_exhausted: bool
    configured_bound_reached: bool
    envelope_limited: bool

    def __post_init__(self) -> None:
        if any(
            not isinstance(value, int) or isinstance(value, bool) or value < minimum
            for value, minimum in (
                (self.maximum_results, 1),
                (self.results_returned, 0),
                (self.screened_exclusions, 0),
                (self.graph_depth, 1),
            )
        ):
            raise CoverageError("coverage bounds contain an invalid count or depth")
        try:
            start = date.fromisoformat(self.date_start)
            end = date.fromisoformat(self.date_end)
        except (TypeError, ValueError) as error:
            raise CoverageError("coverage bounds require ISO dates") from error
        if start > end:
            raise CoverageError("coverage date range is inverted")
        if self.results_returned > self.maximum_results:
            raise CoverageError("coverage results exceed the configured bound")
        if self.configured_bound_reached and self.results_returned != self.maximum_results:
            raise CoverageError("configured-bound status contradicts result count")
        if self.source_exhausted and self.configured_bound_reached:
            raise CoverageError("a bound-limited search cannot claim source exhaustion")


@dataclass(frozen=True, slots=True)
class CoverageRecord:
    coverage_id: str
    seed_id: str
    subtype_id: str | None
    route: str
    source: str
    domain: str
    state: str
    bound: CoverageBound | None

    def __post_init__(self) -> None:
        if self.state not in COVERAGE_STATES:
            raise CoverageError(f"undefined coverage state: {self.state}")
        if (
            not self.seed_id
            or self.route not in ROUTES
            or not self.source
            or not self.domain
            or self.subtype_id == ""
        ):
            raise CoverageError("coverage dimensions are incomplete")
        if self.state in SEARCH_STATES and self.bound is None:
            raise CoverageError("search coverage states require configured bounds")
        if self.bound is not None:
            if self.state == "no_result" and (
                self.bound.results_returned != 0
                or not self.bound.source_exhausted
                or self.bound.configured_bound_reached
                or self.bound.envelope_limited
            ):
                raise CoverageError(
                    "bounded or envelope-limited search cannot claim no_result"
                )
            if self.state == "truncated" and (
                not self.bound.configured_bound_reached
                or self.bound.source_exhausted
            ):
                raise CoverageError("truncated coverage must identify a reached bound")
            if self.bound.envelope_limited and self.state not in ENVELOPE_LIMITED_STATES:
                raise CoverageError(
                    "envelope-limited search requires an explicit incomplete state"
                )
            if self.state in ENVELOPE_LIMITED_STATES and not self.bound.envelope_limited:
                raise CoverageError("envelope incomplete state requires an envelope limit")
        elif self.state in ENVELOPE_LIMITED_STATES:
            raise CoverageError("envelope incomplete state requires bound metadata")
        expected_id = make_coverage_id(
            self.seed_id, self.subtype_id, self.route, self.source, self.domain
        )
        if self.coverage_id != expected_id:
            raise CoverageError("coverage ID is not stable for its dimensions")

    def as_track2_payload(self) -> dict:
        return {
            "bound": asdict(self.bound) if self.bound is not None else None,
            "domain": self.domain,
            "route": self.route,
            "seed_id": self.seed_id,
            "source": self.source,
            "state": self.state,
            "subtype_id": self.subtype_id,
        }


def make_coverage_id(
    seed_id: str,
    subtype_id: str | None,
    route: str,
    source: str,
    domain: str,
) -> str:
    content = json.dumps(
        [seed_id, subtype_id, route, source, domain],
        separators=(",", ":"),
    )
    return "coverage-record:" + hashlib.sha256(content.encode()).hexdigest()[:24]


def make_coverage_record(
    *,
    seed_id: str,
    subtype_id: str | None,
    route: str,
    source: str,
    domain: str,
    state: str,
    bound: CoverageBound | None,
) -> CoverageRecord:
    return CoverageRecord(
        coverage_id=make_coverage_id(seed_id, subtype_id, route, source, domain),
        seed_id=seed_id,
        subtype_id=subtype_id,
        route=route,
        source=source,
        domain=domain,
        state=state,
        bound=bound,
    )


def validate_coverage_matrix(
    records: Iterable[CoverageRecord],
    expected_dimensions: Iterable[tuple[str, str | None, str, str, str]],
) -> None:
    records_tuple = tuple(records)
    actual = {
        (item.seed_id, item.subtype_id, item.route, item.source, item.domain)
        for item in records_tuple
    }
    expected = set(expected_dimensions)
    if len(actual) != len(records_tuple):
        raise CoverageError("coverage matrix contains duplicate dimensions")
    missing = expected - actual
    extra = actual - expected
    if missing or extra:
        raise CoverageError(
            f"coverage matrix is incomplete: missing={len(missing)}, extra={len(extra)}"
        )
