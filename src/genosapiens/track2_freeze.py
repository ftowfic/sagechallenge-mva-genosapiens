"""Immutable public-envelope freeze and offline synthesis gate."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import date
import hashlib
import json
from pathlib import Path
import re
from typing import Any, Callable, Iterable, Mapping

from genosapiens.track2_acquisition import publish_immutable_public_bytes
from genosapiens.track2_coverage import (
    CoverageBound,
    CoverageRecord,
    make_coverage_record,
)
from genosapiens.track2_policy import ROUTES
from genosapiens.track2_public_query import PublicSeedRegistry


SHA256 = re.compile(r"^[0-9a-f]{64}$")
SCHEMA = "genosapiens.track2-frozen-public-envelope/v3"
SCOPE_SCHEMA = "genosapiens.track2-public-research-scope/v1"


def declare_public_research_scope(
    registry: PublicSeedRegistry,
    *,
    included_entities: Iterable[tuple[str, str]],
    selected_routes: Iterable[str],
    source_domains: Iterable[str],
    rationale: str,
    declared_on: str,
    expansion_policy: Mapping[str, Any],
    decision_policy: Mapping[str, Any],
    reference_root: str | Path,
) -> dict[str, Any]:
    """Publish scope in public-only research, before acquisition and handoff access."""
    entities = [asdict(entry) for entry in sorted(
        registry.entries, key=lambda entry: (entry.entity_type, entry.entity_id, entry.token)
    )]
    inventory = {(entry.entity_type, entry.entity_id) for entry in registry.entries}
    included = set(included_entities)
    selected_routes = tuple(selected_routes)
    semantic = {
        "schema": SCOPE_SCHEMA,
        "classification": "public",
        "declared_on": declared_on,
        "selection_origin": "patient_independent_public_seeds",
        "entity_provenance": entities,
        "included_entities": [list(item) for item in sorted(included)],
        "deferred_entities": [list(item) for item in sorted(inventory - included)],
        "selected_routes": sorted(set(selected_routes)),
        "deferred_routes": sorted(ROUTES - set(selected_routes)),
        "source_domains": sorted(set(source_domains)),
        "selection_rationale": rationale,
        "expansion_policy_sha256": hashlib.sha256(_canonical_bytes(expansion_policy)).hexdigest(),
        "decision_policy_sha256": hashlib.sha256(_canonical_bytes(decision_policy)).hexdigest(),
        "route_bounds": decision_policy["route_bounds"],
    }
    scope = json.loads(_canonical_bytes(semantic))
    scope["content_sha256"] = _content_sha256(scope)
    validate_public_research_scope(scope)
    path = Path(reference_root) / "derived/track2/scopes" / (scope["content_sha256"] + ".json")
    publish_immutable_public_bytes(path, _canonical_bytes(scope))
    return scope


def validate_public_research_scope(
    scope: Mapping[str, Any], *, decision_policy: Mapping[str, Any] | None = None,
    expansion_policy: Mapping[str, Any] | None = None,
) -> None:
    required = {
        "schema", "classification", "declared_on", "selection_origin", "entity_provenance",
        "included_entities", "deferred_entities", "selected_routes", "deferred_routes",
        "source_domains", "selection_rationale", "expansion_policy_sha256",
        "decision_policy_sha256", "route_bounds", "content_sha256",
    }
    if (set(scope) != required or scope.get("schema") != SCOPE_SCHEMA
            or scope.get("classification") != "public"
            or scope.get("selection_origin") != "patient_independent_public_seeds"
            or scope.get("content_sha256") != _content_sha256(scope)):
        raise FrozenEnvelopeError("public research scope failed integrity or provenance validation")
    try:
        date.fromisoformat(scope["declared_on"])
        from genosapiens.track2_public_query import PublicSeedProvenance
        entries = tuple(PublicSeedProvenance(**entry) for entry in scope["entity_provenance"])
        for entry in entries:
            recreated = PublicSeedProvenance.create(**{
                key: value for key, value in asdict(entry).items()
                if key not in {"token", "classification"}
            })
            if (entry != recreated or any(word in entry.origin.lower()
                    for word in ("controlled", "handoff", "challenge_result"))
                    or SHA256.fullmatch(entry.manifest_sha256) is None):
                raise FrozenEnvelopeError("public research scope has controlled or invalid seed provenance")
        inventory = {(entry.entity_type, entry.entity_id) for entry in entries}
        included = {tuple(item) for item in scope["included_entities"]}
        deferred = {tuple(item) for item in scope["deferred_entities"]}
        selected_routes = set(scope["selected_routes"])
        deferred_routes = set(scope["deferred_routes"])
        if (not inventory or not included or included & deferred or included | deferred != inventory
                or not selected_routes or selected_routes & deferred_routes
                or selected_routes | deferred_routes != ROUTES
                or not scope["source_domains"] or not scope["selection_rationale"].strip()
                or set(scope["route_bounds"]) != ROUTES
                or any(SHA256.fullmatch(scope[field]) is None for field in (
                    "expansion_policy_sha256", "decision_policy_sha256"))):
            raise FrozenEnvelopeError("public research scope selection is incomplete")
        root = Path(__file__).resolve().parents[2]
        domains = json.loads((root / "config/track2/source-domain-matrix-v1.json").read_text())["domains"]
        if set(scope["source_domains"]) - set(domains):
            raise FrozenEnvelopeError("public research scope uses undefined source domains")
        for field in ("included_entities", "deferred_entities", "selected_routes", "deferred_routes", "source_domains"):
            values = scope[field]
            normalized = [tuple(item) if isinstance(item, list) else item for item in values]
            if normalized != sorted(set(normalized)):
                raise FrozenEnvelopeError("public research scope selection is noncanonical")
        for name, policy in (("decision", decision_policy), ("expansion", expansion_policy)):
            if policy is not None and scope[name + "_policy_sha256"] != hashlib.sha256(
                _canonical_bytes(policy)
            ).hexdigest():
                raise FrozenEnvelopeError("public research scope has a stale " + name + " policy")
    except (KeyError, TypeError, ValueError) as error:
        if isinstance(error, FrozenEnvelopeError):
            raise
        raise FrozenEnvelopeError("public research scope is malformed") from error


class FrozenEnvelopeError(ValueError):
    """Raised when public evidence was not frozen before offline synthesis."""


@dataclass(frozen=True, slots=True)
class FrozenEnvelopePublication:
    envelope_id: str
    path: Path
    sha256: str
    reused: bool


@dataclass(frozen=True, slots=True)
class SynthesisRequest:
    seed_id: str
    subtype_id: str | None
    route: str
    source: str
    domain: str
    required_entities: tuple[tuple[str, str], ...]

    def __post_init__(self) -> None:
        if (
            not self.seed_id
            or self.route not in ROUTES
            or not self.source
            or not self.domain
            or not self.required_entities
            or any(not entity_type or not entity_id for entity_type, entity_id in self.required_entities)
        ):
            raise FrozenEnvelopeError("offline synthesis request is incomplete")


@dataclass(frozen=True, slots=True)
class SynthesisResult:
    authorized: bool
    submission_ready: bool
    offline_synthesis_executed: bool
    missing_entities: tuple[tuple[str, str], ...]
    coverage: tuple[CoverageRecord, ...]
    output: Any = None


def _canonical_bytes(value: Mapping[str, Any]) -> bytes:
    return (
        json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
        + "\n"
    ).encode()


def _content_sha256(value: Mapping[str, Any]) -> str:
    semantic = {key: item for key, item in value.items() if key != "content_sha256"}
    return hashlib.sha256(_canonical_bytes(semantic)).hexdigest()


def freeze_public_envelope(
    registry: PublicSeedRegistry,
    *,
    source_manifest_sha256s: Mapping[str, str],
    runtime: Mapping[str, str],
    coverage_records: Iterable[CoverageRecord],
    expansion_policy: Mapping[str, Any],
    frozen_on: str,
    reference_root: str | Path,
    public_scope: Mapping[str, Any],
) -> FrozenEnvelopePublication:
    try:
        date.fromisoformat(frozen_on)
    except (TypeError, ValueError) as error:
        raise FrozenEnvelopeError("frozen_on must be an ISO date") from error
    required_bindings = {"input", "input_review", "plan", "plan_review", "code", "policy", "input_schema"}
    if (not required_bindings <= set(source_manifest_sha256s)
            or any(not isinstance(value, str) or not SHA256.fullmatch(value)
                   for value in source_manifest_sha256s.values())
            or set(runtime) != {"image", "code_revision"}
            or not runtime["image"].startswith("genosapiens/")
            or "@sha256:" not in runtime["image"] or not runtime["code_revision"]):
        raise FrozenEnvelopeError("current source, review, code and runtime bindings are required")
    if expansion_policy.get("schema") != "genosapiens.track2-public-expansion-policy/v1":
        raise FrozenEnvelopeError("validated public expansion policy is required")
    validate_public_research_scope(public_scope, expansion_policy=expansion_policy)
    if (public_scope["declared_on"] > frozen_on
            or public_scope["expansion_policy_sha256"] != hashlib.sha256(
                _canonical_bytes(expansion_policy)).hexdigest()):
        raise FrozenEnvelopeError("public scope is stale or was declared after freeze")
    entities = tuple(asdict(entry) for entry in registry.entries)
    if not entities or any(
        entry["classification"] != "public"
        or "controlled" in entry["origin"].lower()
        or "handoff" in entry["origin"].lower()
        for entry in entities
    ):
        raise FrozenEnvelopeError("frozen envelope contains non-public provenance")
    if {entry["token"] for entry in public_scope["entity_provenance"]} - {
        entry["token"] for entry in entities
    }:
        raise FrozenEnvelopeError("scope seed identities are absent from the frozen envelope")
    coverage_records = tuple(coverage_records)
    coverage = tuple(sorted(record.coverage_id for record in coverage_records))
    date_range = expansion_policy["date_range"]
    semantic = {
        "schema": SCHEMA,
        "classification": "public",
        "frozen_on": frozen_on,
        "network_acquisition_allowed": False,
        "date_range": {"start": date_range["start"], "end": date_range["end"]},
        "manifest_checksums": dict(sorted(source_manifest_sha256s.items())),
        "runtime": dict(runtime),
        "entity_provenance": entities,
        "coverage_record_ids": coverage,
        "coverage_records": [asdict(record) for record in sorted(
            coverage_records, key=lambda record: record.coverage_id)],
        "public_scope": public_scope,
    }
    identity = hashlib.sha256(_canonical_bytes(semantic)).hexdigest()
    envelope = {
        **semantic,
        "envelope_id": f"frozen-envelope:{identity[:24]}",
    }
    envelope["content_sha256"] = _content_sha256(envelope)
    content = _canonical_bytes(envelope)
    root = Path(reference_root).resolve()
    path = root / f"derived/track2/envelopes/{identity[:24]}.json"
    reused = path.exists()
    published_sha256 = publish_immutable_public_bytes(path, content)
    return FrozenEnvelopePublication(
        envelope_id=envelope["envelope_id"],
        path=path,
        sha256=published_sha256,
        reused=reused,
    )


def load_frozen_public_envelope(path: str | Path) -> dict[str, Any]:
    source = Path(path)
    if source.is_symlink() or not source.is_file():
        raise FrozenEnvelopeError("frozen public envelope is missing or invalid")
    try:
        content = source.read_bytes()
        envelope = json.loads(content)
    except (OSError, json.JSONDecodeError) as error:
        raise FrozenEnvelopeError(f"cannot load frozen public envelope: {error}") from error
    required = {
        "classification",
        "content_sha256",
        "coverage_record_ids",
        "coverage_records",
        "public_scope",
        "date_range",
        "entity_provenance",
        "envelope_id",
        "frozen_on",
        "manifest_checksums",
        "runtime",
        "network_acquisition_allowed",
        "schema",
    }
    if (
        not isinstance(envelope, dict)
        or set(envelope) != required
        or envelope.get("schema") != SCHEMA
        or envelope.get("classification") != "public"
        or envelope.get("network_acquisition_allowed") is not False
        or envelope.get("content_sha256") != _content_sha256(envelope)
        or content != _canonical_bytes(envelope)
    ):
        raise FrozenEnvelopeError("frozen public envelope failed integrity validation")
    validate_public_research_scope(envelope["public_scope"])
    try:
        records = tuple(CoverageRecord(**{
            **record, "bound": CoverageBound(**record["bound"]) if record["bound"] else None,
        }) for record in envelope["coverage_records"])
    except (KeyError, TypeError, ValueError) as error:
        raise FrozenEnvelopeError("frozen coverage records are invalid") from error
    if sorted(record.coverage_id for record in records) != envelope["coverage_record_ids"]:
        raise FrozenEnvelopeError("frozen coverage identities do not resolve")
    identity = hashlib.sha256(
        _canonical_bytes(
            {
                key: value
                for key, value in envelope.items()
                if key not in {"content_sha256", "envelope_id"}
            }
        )
    ).hexdigest()[:24]
    if envelope["envelope_id"] != f"frozen-envelope:{identity}":
        raise FrozenEnvelopeError("frozen public envelope ID is unstable")
    return envelope


def run_frozen_synthesis(
    request: SynthesisRequest,
    frozen_envelope: Mapping[str, Any],
    offline_synthesis: Callable[[SynthesisRequest], Any],
) -> SynthesisResult:
    if (
        frozen_envelope.get("schema") != SCHEMA
        or frozen_envelope.get("network_acquisition_allowed") is not False
        or frozen_envelope.get("content_sha256") != _content_sha256(frozen_envelope)
    ):
        raise FrozenEnvelopeError("offline synthesis requires a frozen offline envelope")
    available = {
        (entry.get("entity_type"), entry.get("entity_id"))
        for entry in frozen_envelope.get("entity_provenance", ())
        if isinstance(entry, Mapping) and entry.get("classification") == "public"
    }
    scope = frozen_envelope.get("public_scope")
    if not isinstance(scope, Mapping):
        raise FrozenEnvelopeError("offline synthesis requires predeclared public scope")
    validate_public_research_scope(scope)
    available -= {tuple(item) for item in scope["deferred_entities"]}
    missing = tuple(sorted(set(request.required_entities) - available))
    if missing or request.route not in scope["selected_routes"] or request.domain not in scope["source_domains"]:
        date_range = frozen_envelope["date_range"]
        coverage = make_coverage_record(
            seed_id=request.seed_id,
            subtype_id=request.subtype_id,
            route=request.route,
            source=request.source,
            domain=request.domain,
            state="analysis_incomplete_due_to_envelope",
            bound=CoverageBound(
                maximum_results=1,
                results_returned=0,
                screened_exclusions=0,
                graph_depth=1,
                date_start=date_range["start"],
                date_end=date_range["end"],
                source_exhausted=False,
                configured_bound_reached=False,
                envelope_limited=True,
            ),
        )
        return SynthesisResult(
            authorized=False,
            submission_ready=False,
            offline_synthesis_executed=False,
            missing_entities=missing,
            coverage=(coverage,),
        )
    output = offline_synthesis(request)
    return SynthesisResult(
        authorized=True,
        submission_ready=False,
        offline_synthesis_executed=True,
        missing_entities=(),
        coverage=(),
        output=output,
    )
