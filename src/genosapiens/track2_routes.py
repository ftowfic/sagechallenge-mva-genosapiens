"""Complete, bounded assessment matrix for every candidate-generation route."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Mapping

from genosapiens.track2_policy import ROUTES
from genosapiens.track2_assessment import assess_objective, route_evidence


ATTEMPT_STATES = {
    "attempted_complete",
    "deferred",
    "no_result",
    "truncated",
    "unavailable",
    "unattempted",
    "unimplemented",
}
INCOMPLETE_STATES = {
    "deferred",
    "truncated",
    "unavailable",
    "unattempted",
    "unimplemented",
}


class RouteAssessmentError(ValueError):
    """Raised when a route is skipped, unbounded, or mislabeled not applicable."""


@dataclass(frozen=True, slots=True)
class RouteAttempt:
    state: str
    executed: bool
    sources: tuple[str, ...]
    query_ids: tuple[str, ...]
    screening_exclusions: tuple[str, ...]
    max_depth: int
    max_results: int
    lookback_years: int
    screening: str
    stopping: str

    def __post_init__(self) -> None:
        if self.state not in ATTEMPT_STATES:
            raise RouteAssessmentError("route attempt state is undefined")
        if not self.sources or any(not item for item in self.sources):
            raise RouteAssessmentError("route attempt sources must be explicit")
        if self.executed and (not self.query_ids or any(not item for item in self.query_ids)):
            raise RouteAssessmentError("executed route queries must be explicit")
        if any(not item for item in self.screening_exclusions):
            raise RouteAssessmentError("screening exclusions must be explicit")
        if self.state in {"deferred", "unattempted", "unimplemented"} and self.executed:
            raise RouteAssessmentError("unattempted, deferred, or unimplemented route cannot claim execution")
        if self.state not in {"deferred", "unattempted", "unimplemented"} and not self.executed:
            raise RouteAssessmentError("attempted route cannot be marked unexecuted")


@dataclass(frozen=True, slots=True)
class RouteAssessment:
    objective_id: str
    route: str
    applicability: str
    biological_rationale: str
    reviewer_role: str | None
    review_record_id: str | None
    evidence_ids: tuple[str, ...]
    attempt: RouteAttempt | None

    def __post_init__(self) -> None:
        if (
            not self.objective_id
            or self.route not in ROUTES
            or self.applicability not in {"applicable", "not_applicable"}
            or not self.biological_rationale
        ):
            raise RouteAssessmentError("route assessment is incomplete")
        if self.applicability == "not_applicable":
            if (
                not self.reviewer_role
                or not self.review_record_id
                or not self.evidence_ids
                or self.attempt is not None
            ):
                raise RouteAssessmentError(
                    "not_applicable requires reviewed biological rationale"
                )
        elif self.attempt is None:
            raise RouteAssessmentError(
                "applicable route requires an explicit bounded attempt state"
            )


@dataclass(frozen=True, slots=True)
class RoutePortfolio:
    assessments: tuple[RouteAssessment, ...]
    all_routes_assessed: bool
    all_applicable_routes_attempted: bool
    analysis_complete: bool
    incomplete_routes: tuple[tuple[str, str, str], ...]
    selected_routes: tuple[str, ...] = ()
    deferred_routes: tuple[tuple[str, str, str], ...] = ()
    public_scope_sha256: str | None = None


def _validate_attempt_bound(
    assessment: RouteAssessment, route_bounds: Mapping[str, object]
) -> None:
    if assessment.attempt is None:
        return
    expected = route_bounds[assessment.route]
    attempt = assessment.attempt
    if not isinstance(expected, Mapping) or (
        attempt.max_depth != expected["max_depth"]
        or attempt.max_results != expected["max_results"]
        or attempt.lookback_years != expected["lookback_years"]
        or attempt.screening != expected["screening"]
        or attempt.stopping != expected["stopping"]
    ):
        raise RouteAssessmentError("route attempt does not match configured bounds")


def validate_route_attempt(
    assessment: RouteAssessment, decision_policy: Mapping[str, object]
) -> None:
    route_bounds = decision_policy.get("route_bounds")
    if not isinstance(route_bounds, Mapping) or assessment.route not in route_bounds:
        raise RouteAssessmentError("decision policy route bounds are missing")
    _validate_attempt_bound(assessment, route_bounds)


def build_route_portfolio(
    objective_ids: Iterable[str],
    assessments: Iterable[RouteAssessment],
    decision_policy: Mapping[str, object],
    *,
    public_scope: Mapping[str, object] | None = None,
) -> RoutePortfolio:
    objectives = tuple(sorted(set(objective_ids)))
    if not objectives:
        raise RouteAssessmentError("at least one intervention objective is required")
    records = tuple(
        sorted(assessments, key=lambda item: (item.objective_id, item.route))
    )
    actual = {(item.objective_id, item.route) for item in records}
    expected = {(objective_id, route) for objective_id in objectives for route in ROUTES}
    if len(actual) != len(records) or actual != expected:
        raise RouteAssessmentError("route assessment matrix is incomplete")
    route_bounds = decision_policy.get("route_bounds")
    if not isinstance(route_bounds, Mapping):
        raise RouteAssessmentError("decision policy route bounds are missing")
    for record in records:
        _validate_attempt_bound(record, route_bounds)
    selected = set(ROUTES)
    deferred = ()
    scope_sha256 = None
    if public_scope is not None:
        from genosapiens.track2_freeze import validate_public_research_scope

        validate_public_research_scope(public_scope, decision_policy=decision_policy)
        selected = set(public_scope["selected_routes"])
        scope_sha256 = public_scope["content_sha256"]
        if public_scope["route_bounds"] != decision_policy["route_bounds"]:
            raise RouteAssessmentError("public scope has stale route bounds")
        omitted = tuple(item for item in records if item.route not in selected)
        if any(
            item.applicability != "applicable"
            or item.attempt is None
            or item.attempt.state != "deferred"
            or item.attempt.executed
            for item in omitted
        ):
            raise RouteAssessmentError("unselected routes must remain explicitly deferred")
        deferred = tuple((item.objective_id, item.route, "deferred") for item in omitted)
    applicable = tuple(
        item for item in records
        if item.applicability == "applicable" and item.route in selected
    )
    incomplete = tuple(
        (
            item.objective_id,
            item.route,
            item.attempt.state if item.attempt is not None else "invalid",
        )
        for item in applicable
        if item.attempt is None or item.attempt.state in INCOMPLETE_STATES
    )
    all_attempted = all(
        item.attempt is not None
        and item.attempt.state not in {"deferred", "unattempted", "unimplemented"}
        for item in applicable
    )
    return RoutePortfolio(
        assessments=records,
        all_routes_assessed=True,
        all_applicable_routes_attempted=all_attempted,
        analysis_complete=all(assess_objective(objective, None, True,
            route_evidence((row for row in records if row.objective_id == objective), selected))['status']
            == 'complete' for objective in objectives),
        incomplete_routes=incomplete,
        selected_routes=tuple(sorted(selected)),
        deferred_routes=deferred,
        public_scope_sha256=scope_sha256,
    )
