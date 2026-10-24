"""Provenance-aware boundary for public Track 2 network queries."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from typing import Any, Callable, Mapping

from genosapiens.track2_policy import (
    Track2PolicyError,
    validate_public_expansion_request,
)


FORBIDDEN_FIELDS = {
    "case_id",
    "challenge_result",
    "controlled_path",
    "handoff",
    "handoff_path",
    "phenotype_vector",
}
PROVENANCE_FIELDS = {
    "classification",
    "entity_id",
    "entity_type",
    "manifest_sha256",
    "origin",
    "source_reference_id",
    "token",
}
QUERY_FIELDS = {
    "adapter",
    "entity_type",
    "parameters",
    "provenance",
    "query_id",
    "seed_id",
}
AuditLogger = Callable[[str, str], None]


class PublicQueryError(ValueError):
    """Raised before I/O when a public query cannot prove public provenance."""


def _token_content(
    entity_type: str,
    entity_id: str,
    origin: str,
    manifest_sha256: str,
    source_reference_id: str,
) -> str:
    return "\0".join(
        (entity_type, entity_id, origin, manifest_sha256, source_reference_id)
    )


@dataclass(frozen=True, slots=True)
class PublicSeedProvenance:
    entity_type: str
    entity_id: str
    origin: str
    manifest_sha256: str
    source_reference_id: str
    token: str
    classification: str = "public"

    @classmethod
    def create(
        cls,
        *,
        entity_type: str,
        entity_id: str,
        origin: str,
        manifest_sha256: str,
        source_reference_id: str,
    ) -> "PublicSeedProvenance":
        if not all(
            isinstance(value, str) and value
            for value in (
                entity_type,
                entity_id,
                origin,
                manifest_sha256,
                source_reference_id,
            )
        ) or len(manifest_sha256) != 64:
            raise PublicQueryError("public seed provenance is incomplete")
        token = hashlib.sha256(
            _token_content(
                entity_type,
                entity_id,
                origin,
                manifest_sha256,
                source_reference_id,
            ).encode()
        ).hexdigest()
        return cls(
            entity_type=entity_type,
            entity_id=entity_id,
            origin=origin,
            manifest_sha256=manifest_sha256,
            source_reference_id=source_reference_id,
            token=token,
        )


@dataclass(frozen=True, slots=True)
class PublicSeedRegistry:
    entries: tuple[PublicSeedProvenance, ...]

    def for_entity(
        self, entity_type: str, entity_id: str
    ) -> tuple[PublicSeedProvenance, ...]:
        return tuple(
            entry
            for entry in self.entries
            if entry.entity_type == entity_type and entry.entity_id == entity_id
        )

    def contains(self, provenance: Mapping[str, Any]) -> bool:
        return any(asdict(entry) == dict(provenance) for entry in self.entries)


def build_public_seed_registry(bundle) -> PublicSeedRegistry:
    """Allowlist only independently public source-bound genotype entities.

    Patient labels are deliberately absent. The source loader verifies actual
    document bytes; synthetic fixtures cannot seed live requests.
    """
    if (bundle.data.get('input_origin') != 'published_patient_genotype'
            or bundle.data.get('classification') != 'public'
            or bundle.data.get('evidence_origin') != 'primary_literature'):
        raise PublicQueryError("independent_public_source_required")
    entries = {}
    for allele in bundle.data["genotype"]["alleles"]:
        for kind, value in (("gene", allele["gene"]),
                            ("variant", allele["gene"] + " " + allele["reported_coding"])):
            entry = PublicSeedProvenance.create(
                entity_type=kind, entity_id=value, origin="published_primary_literature",
                manifest_sha256=bundle.input_sha256,
                source_reference_id=allele["evidence_ids"][0],
            )
            entries[entry.token] = entry
    return PublicSeedRegistry(tuple(sorted(entries.values(), key=lambda item: item.token)))


def make_public_query_id(
    parameters: Mapping[str, Any], provenance_token: str, adapter: str
) -> str:
    content = json.dumps(
        {"adapter": adapter, "parameters": parameters, "provenance": provenance_token},
        sort_keys=True,
        separators=(",", ":"),
    )
    return "public-query:" + hashlib.sha256(content.encode()).hexdigest()[:24]


def _has_forbidden_field(value: Any) -> bool:
    if isinstance(value, Mapping):
        return bool(set(value) & FORBIDDEN_FIELDS) or any(
            _has_forbidden_field(item) for item in value.values()
        )
    if isinstance(value, (list, tuple)):
        return any(_has_forbidden_field(item) for item in value)
    if isinstance(value, str):
        return "controlled-work/" in value or value.startswith("/handoff/")
    return False


def _reject(logger: AuditLogger, reason: str) -> None:
    logger("public_query_rejected", reason)
    raise PublicQueryError(reason)


def validate_public_query(
    query: Mapping[str, Any],
    registry: PublicSeedRegistry,
    expansion_policy: Mapping[str, Any],
    logger: AuditLogger,
    *,
    public_scope: Mapping[str, Any] | None = None,
    route: str | None = None,
    domain: str | None = None,
) -> None:
    if set(query) != QUERY_FIELDS:
        _reject(logger, "query_shape_invalid")
    if _has_forbidden_field(query):
        _reject(logger, "controlled_field_present")
    provenance = query.get("provenance")
    if not isinstance(provenance, Mapping) or set(provenance) != PROVENANCE_FIELDS:
        _reject(logger, "query_shape_invalid")
    if (provenance.get("classification") != "public"
            or provenance.get("origin") not in {"published_primary_literature", "reviewed_public_target"}
            or not registry.contains(provenance)):
        _reject(logger, "unregistered_public_seed_provenance")
    parameters = query.get("parameters")
    if not isinstance(parameters, Mapping):
        _reject(logger, "query_shape_invalid")
    if (
        query.get("seed_id") != provenance.get("entity_id")
        or query.get("entity_type") != provenance.get("entity_type")
        or parameters.get("seed_id") != query.get("seed_id")
        or parameters.get("from_entity_type") != query.get("entity_type")
        or query.get("adapter") != parameters.get("source")
    ):
        _reject(logger, "query_shape_invalid")
    try:
        validate_public_expansion_request(parameters, expansion_policy)
    except Track2PolicyError:
        _reject(logger, "out_of_policy_query")
    if public_scope is not None:
        from genosapiens.track2_freeze import FrozenEnvelopeError, validate_public_research_scope

        try:
            validate_public_research_scope(public_scope, expansion_policy=expansion_policy)
        except FrozenEnvelopeError:
            _reject(logger, "public_scope_invalid")
        if (
            [provenance["entity_type"], provenance["entity_id"]] not in public_scope["included_entities"]
            or dict(provenance) not in public_scope["entity_provenance"]
            or route not in public_scope["selected_routes"]
            or domain not in public_scope["source_domains"]
        ):
            _reject(logger, "outside_declared_public_scope")
    elif route is not None or domain is not None:
        _reject(logger, "public_scope_required")
    expected_id = make_public_query_id(
        parameters, str(provenance["token"]), str(query["adapter"])
    )
    if query.get("query_id") != expected_id:
        _reject(logger, "query_shape_invalid")
    logger("public_query_accepted", "validated")


def execute_public_query(
    query: Mapping[str, Any],
    registry: PublicSeedRegistry,
    expansion_policy: Mapping[str, Any],
    http_call: Callable[[Mapping[str, Any]], Any],
    logger: AuditLogger,
    *,
    public_scope: Mapping[str, Any] | None = None,
    route: str | None = None,
    domain: str | None = None,
) -> Any:
    """Validate provenance and policy completely before invoking network I/O."""
    validate_public_query(query, registry, expansion_policy, logger,
                          public_scope=public_scope, route=route, domain=domain)
    return http_call(query["parameters"])
