"""Lossless public disease profiles over immutable rare-disease v6 evidence.

This is task 7.2, not a phenotype scorer. Source-dependency groups are scoped to
one disease view, source-native disease, resolved term, and annotation aspect.
They retain every frequency/context/assertion variant and source reference;
membership is not a claim of independent corroboration or a choice of frequency.
Neither the legacy canonical IDs nor global preference flags select evidence.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import asdict, dataclass
from datetime import date
import hashlib
import json
from pathlib import Path
import re
from types import MappingProxyType

from genosapiens.bounded_json import BoundedJSONError, load_bounded_json
from genosapiens.hpo_ic import HpoGraph, HpoICError, build_hpo_graph
from genosapiens.ontology_bundle import (
    OBSOLETE_RESOLUTION_POLICY,
    OntologyBundleContract,
    OntologyBundleError,
    OntologyResolution,
    _file_identity, _file_signature,
    _parse_obo,
    _require_within_root,
)
from genosapiens.rare_disease import (
    BUNDLE_MANIFEST_SCHEMA,
    RareDiseaseContract,
    RareDiseaseError,
    TABLE_COLUMNS,
    TABLE_SCHEMAS,
    _artifact_contract_sha256,
    _bundle_manifest_path,
    _source_manifest_path,
)


SCHEMA = "genosapiens.disease-profiles"
POLICY = MappingProxyType({
    "obsolete_resolution": OBSOLETE_RESOLUTION_POLICY,
    "rows": "all-v6-phenotype-assertions-no-preference-filter",
    "disease_views": "source-native-and-explicit-mapped-no-catalog-join",
    "dependency": "disease-view-source-family-native-disease-new-term-aspect-scope",
    "missing_native_identity": "dataset-release-source-record-no-inferred-equivalence",
    "missing_hpo_identity": "record-scoped-unresolved",
    "assertion_variants": "retain-all-frequency-context-qualifier-and-reference-variants",
    "qualifier_disagreement": "unresolved-no-positive-coercion",
    "frequency": "source-values-only-no-probability-conversion",
})
SOURCE_COLUMNS = TABLE_COLUMNS["phenotype_assertion"]
_ARRAY_COLUMNS = frozenset({"onset", "modifiers", "references"})
_BOOL_COLUMNS = frozenset({
    "is_ontology_propagated", "preferred_representation", "counts_as_independent", "qualifier_not",
})
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_HPO = re.compile(r"^HP:[0-9]{7}$")


class DiseaseProfileError(RuntimeError):
    """Malformed evidence or unbound public snapshot inputs."""


def _freeze(value: object) -> object:
    if isinstance(value, Mapping):
        if any(not isinstance(key, str) for key in value):
            raise DiseaseProfileError("profile metadata keys must be strings")
        return MappingProxyType({key: _freeze(child) for key, child in value.items()})
    if isinstance(value, (list, tuple)):
        return tuple(_freeze(child) for child in value)
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    raise DiseaseProfileError("profile evidence must contain JSON-compatible values")


def _json_default(value: object) -> object:
    """Expose immutable mappings without recursively copying the whole store."""
    if isinstance(value, Mapping):
        return dict(value)
    raise TypeError("profile evidence must contain JSON-compatible values")


def _bytes(value: object) -> bytes:
    try:
        return (json.dumps(value, default=_json_default, sort_keys=True,
                           separators=(",", ":"), allow_nan=False) + "\n").encode()
    except (TypeError, ValueError, RecursionError):
        raise DiseaseProfileError("profile evidence is not finite JSON") from None


@dataclass(frozen=True, slots=True)
class ProfileAnnotation:
    annotation_id: str
    source_fields: Mapping[str, object]
    normalized_source_hpo_id: str | None
    resolution: OntologyResolution | None
    polarity: str
    diagnostics: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class ProfileDependencyGroup:
    group_id: str
    aspect: str | None
    canonical_hpo_id: str | None
    annotation_ids: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class DiseaseProfile:
    disease_id: str
    groups: tuple[ProfileDependencyGroup, ...]

    def groups_for_aspect(self, aspect: str) -> tuple[ProfileDependencyGroup, ...]:
        """Select an explicit source aspect; never infer P from ontology ancestry."""
        return tuple(group for group in self.groups if group.aspect == aspect)


@dataclass(frozen=True, slots=True)
class DiseaseProfileStore:
    annotations: tuple[ProfileAnnotation, ...]
    profiles: Mapping[str, DiseaseProfile]
    unassigned_annotation_ids: tuple[str, ...]
    source_binding: Mapping[str, object]
    graph_sha256: str


def build_disease_profiles(
    graph: HpoGraph,
    rows: Iterable[Mapping[str, object]],
    *,
    source_binding: Mapping[str, object],
) -> DiseaseProfileStore:
    """Pure builder for public/synthetic v6-shaped rows, including nullable fields.

    A source-native view survives absent mappings; explicit mapped views do not
    imply disease equivalence. Each view groups shared source lineages afresh,
    separating P/I/C (and unknown aspects) and propagated/other scopes. Original
    lineage/dependency/assertion keys are retained only as source references.
    No representative annotation is chosen: downstream context interpretation
    must inspect the retained variants, including unresolved qualifier conflicts.
    """
    if not isinstance(source_binding, Mapping):
        raise DiseaseProfileError("profile source binding must be an object")
    annotations = []
    occurrences: dict[str, int] = {}
    groups: dict[tuple[str, str], list[ProfileAnnotation]] = {}
    unassigned = []
    for row in rows:
        if not isinstance(row, Mapping) or not set(SOURCE_COLUMNS).issubset(row):
            raise DiseaseProfileError("phenotype assertion is missing v6 source fields")
        for name in SOURCE_COLUMNS:
            value = row[name]
            valid = value is None
            if name in _ARRAY_COLUMNS:
                valid |= isinstance(value, (list, tuple)) and all(item is None or isinstance(item, str) for item in value)
            elif name in _BOOL_COLUMNS:
                valid |= type(value) is bool
            elif name == "representation_priority":
                valid |= type(value) is int
            else:
                valid |= isinstance(value, str)
            if not valid:
                raise DiseaseProfileError(f"invalid phenotype source field type: {name}")
        if any(not row[name] for name in ("source_family", "source_dataset", "source_release", "source_record_id")):
            raise DiseaseProfileError("phenotype assertion lacks source record provenance")
        fields = _freeze(row)
        identity = hashlib.sha256(_bytes(fields)).hexdigest()
        occurrences[identity] = occurrences.get(identity, 0) + 1
        annotation_id = f"{identity}:{occurrences[identity]}"
        source_hpo = row["source_hpo_id"]
        if (
            row["source_dataset"] == "open_targets_disease_phenotype"
            and isinstance(source_hpo, str) and re.fullmatch(r"HP_[0-9]{7}", source_hpo)
        ):
            source_hpo = source_hpo.replace("_", ":", 1)
        resolution = graph.resolve(source_hpo) if source_hpo else None
        diagnostics = []
        if source_hpo is None or not _HPO.fullmatch(source_hpo):
            diagnostics.append("missing_or_invalid_source_hpo_id")
        if resolution is None or resolution.canonical_id is None:
            diagnostics.append("unresolved_hpo")
        qualifier, qualifier_not = row["qualifier"], row["qualifier_not"]
        if qualifier not in (None, "", "NOT"):
            polarity = "unresolved"
            diagnostics.append("unsupported_qualifier")
        elif qualifier is not None and qualifier_not is not None and (qualifier == "NOT") != qualifier_not:
            polarity = "unresolved"
            diagnostics.append("qualifier_disagreement")
        elif qualifier == "NOT" or qualifier_not is True:
            polarity = "negative"
        elif qualifier == "" or qualifier_not is False:
            polarity = "positive"
        else:
            polarity = "unknown"
        if row["aspect"] not in ("P", "I", "C"):
            diagnostics.append("unknown_aspect")
        annotation = ProfileAnnotation(annotation_id, fields, source_hpo, resolution, polarity, tuple(diagnostics))
        annotations.append(annotation)
        disease_ids = {value for value in (row["source_disease_id"], row["mapped_disease_id"]) if value}
        if not disease_ids:
            unassigned.append(annotation_id)
        # The old lineage embeds a legacy canonical ID; never reuse it as a key.
        native_identity = (
            ("disease", row["source_disease_id"]) if row["source_disease_id"]
            else ("record", row["source_dataset"], row["source_release"], row["source_record_id"])
        )
        term_identity = (
            ("resolved", resolution.canonical_id) if resolution and resolution.canonical_id
            else ("unresolved", source_hpo)
        )
        if not source_hpo:
            term_identity = ("missing", row["source_dataset"], row["source_release"], row["source_record_id"])
        scope = row["association_scope"]
        if scope in ("source_direct", "ot_mapped_evidence"):
            scope = "source-assertion"
        for disease_id in disease_ids:
            key = hashlib.sha256(_bytes((
                SCHEMA, POLICY, disease_id, row["source_family"], native_identity,
                term_identity, row["aspect"], scope, row["is_ontology_propagated"],
            ))).hexdigest()
            groups.setdefault((disease_id, key), []).append(annotation)
    by_disease: dict[str, list[ProfileDependencyGroup]] = {}
    for (disease_id, key), members in sorted(groups.items()):
        first = members[0]
        group = ProfileDependencyGroup(
            key, first.source_fields["aspect"],
            first.resolution.canonical_id if first.resolution else None,
            tuple(sorted(item.annotation_id for item in members)),
        )
        by_disease.setdefault(disease_id, []).append(group)
    graph_sha256 = hashlib.sha256(_bytes({
        "policy": OBSOLETE_RESOLUTION_POLICY,
        "terms": [asdict(term) for _, term in sorted(graph.ontology.terms.items())],
    })).hexdigest()
    binding = _freeze(source_binding)
    _bytes(binding)
    return DiseaseProfileStore(
        tuple(sorted(annotations, key=lambda item: item.annotation_id)),
        MappingProxyType({key: DiseaseProfile(key, tuple(value)) for key, value in by_disease.items()}),
        tuple(sorted(unassigned)), binding, graph_sha256,
    )


def disease_profiles_as_bytes(store: DiseaseProfileStore) -> bytes:
    """Deterministic JSON, storing evidence once and referencing it from views."""
    return _bytes({
        "schema": SCHEMA, "policy": POLICY, "source_binding": store.source_binding,
        "graph_sha256": store.graph_sha256,
        "annotations": [{
            "annotation_id": item.annotation_id, "source_fields": item.source_fields,
            "normalized_source_hpo_id": item.normalized_source_hpo_id,
            "resolution": asdict(item.resolution) if item.resolution else None,
            "polarity": item.polarity, "diagnostics": item.diagnostics,
        } for item in store.annotations],
        "profiles": [asdict(profile) for profile in store.profiles.values()],
        "unassigned_annotation_ids": store.unassigned_annotation_ids,
    })


def read_pinned_disease_profiles(
    *,
    contract_path: str | Path,
    ontology_contract_path: str | Path,
    reference_root: str | Path,
    expected_bundle_sha256: str,
) -> DiseaseProfileStore:
    """Build from an installed, checksum-pinned public v6 snapshot, read-only.

    Reads only governing JSON, the HPO OBO, and phenotype_assertion.parquet.
    Upstream source identities are bound through the pinned completed bundle;
    unused source files, disease catalogs, and other projections are not opened.
    DuckDB uses an ephemeral in-memory connection, no persistent DB or downloads.
    No output or source path is created, updated, repaired, or published here.
    """
    try:
        return _read_pinned_disease_profiles(
            contract_path=contract_path, ontology_contract_path=ontology_contract_path,
            reference_root=reference_root, expected_bundle_sha256=expected_bundle_sha256,
        )
    except (BoundedJSONError, OntologyBundleError, HpoICError, RareDiseaseError, OSError, KeyError, TypeError, ValueError):
        raise DiseaseProfileError("invalid or mismatched public disease-profile input provenance") from None


def _read_pinned_disease_profiles(
    *, contract_path: str | Path, ontology_contract_path: str | Path,
    reference_root: str | Path, expected_bundle_sha256: str,
) -> DiseaseProfileStore:
    control_paths = (Path(contract_path), Path(ontology_contract_path))
    signatures = {path: _file_signature(path) for path in control_paths}
    root = Path(reference_root).resolve(strict=True)
    if not isinstance(expected_bundle_sha256, str) or not _SHA256.fullmatch(expected_bundle_sha256):
        raise DiseaseProfileError("a pinned bundle SHA-256 is required")
    contract_data, contract_bytes = load_bounded_json(contract_path)
    if not isinstance(contract_data, Mapping):
        raise DiseaseProfileError("v6 source contract must be an object")
    contract = RareDiseaseContract.from_mapping(contract_data)
    if contract.projection.schema_version != "v6":
        raise DiseaseProfileError("disease profiles require the immutable v6 source contract")
    ontology_data, ontology_bytes = load_bounded_json(ontology_contract_path)
    if not isinstance(ontology_data, Mapping):
        raise DiseaseProfileError("ontology source contract must be an object")
    ontology = OntologyBundleContract.from_mapping(ontology_data)
    hpo = ontology.artifact("hpo_ontology")
    hpoa = ontology.artifact("hpo_disease_annotations")
    if (
        (hpo.path, hpo.sha256, hpo.release) != (contract.hpo.ontology_path, contract.hpo.ontology_sha256, contract.hpo.release)
        or (hpoa.path, hpoa.sha256, hpoa.release) != (contract.hpo.annotations_path, contract.hpo.annotations_sha256, contract.hpo.release)
    ):
        raise DiseaseProfileError("HPO ontology/annotation pins differ from the v6 source contract")
    contract_sha256 = hashlib.sha256(contract_bytes).hexdigest()
    manifest_path = _bundle_manifest_path(root, contract, contract_sha256)
    legacy_manifest_path = _bundle_manifest_path(root, contract)
    if (
        not manifest_path.exists()
        and legacy_manifest_path.exists()
        and _file_identity(legacy_manifest_path)[0] == expected_bundle_sha256
    ):
        manifest_path = legacy_manifest_path
    _require_within_root(root, manifest_path, "disease profile bundle manifest")
    signatures[manifest_path] = _file_signature(manifest_path)
    manifest, manifest_bytes = load_bounded_json(manifest_path)
    if hashlib.sha256(manifest_bytes).hexdigest() != expected_bundle_sha256:
        raise DiseaseProfileError("disease profile bundle checksum mismatch")
    if not isinstance(manifest, Mapping):
        raise DiseaseProfileError("public v6 bundle manifest must be an object")
    identity = {
        "schema": BUNDLE_MANIFEST_SCHEMA, "classification": "public",
        "bundle_id": contract.bundle_id, "projection_id": contract.projection_id,
        "projection_schema": "v6", "contract_sha256": contract_sha256,
        "assembly": contract.assembly,
    }
    if any(manifest.get(key) != value for key, value in identity.items()):
        raise DiseaseProfileError("bundle manifest is not bound to the public v6 source contract")
    if not re.fullmatch(r"sha256:[0-9a-f]{64}", manifest["image_digest"]):
        raise DiseaseProfileError("bundle manifest image identity is invalid")
    date.fromisoformat(manifest["retrieved_on"])
    upstream = manifest["upstream"]
    expected_upstream = {
        "hpo_release": contract.hpo.release,
        "hpo_ontology_sha256": contract.hpo.ontology_sha256,
        "hpo_annotations_sha256": contract.hpo.annotations_sha256,
        "open_targets_release": contract.open_targets.release,
        "open_targets_croissant_sha256": contract.open_targets.croissant_sha256,
        "open_targets_manifest_path": contract.open_targets.manifest_path,
        "open_targets_manifest_sha256": upstream["open_targets_manifest_sha256"],
    }
    if upstream != expected_upstream or not _SHA256.fullmatch(upstream["open_targets_manifest_sha256"]):
        raise DiseaseProfileError("bundle upstream identities differ from the source pins")
    if len(manifest["sources"]) != len(contract.artifacts):
        raise DiseaseProfileError("bundle source manifest index is incomplete")
    for source, artifact in zip(manifest["sources"], contract.artifacts, strict=True):
        expected = {
            "role": artifact.role, "resource_id": artifact.resource_id,
            "release": artifact.release, "extraction_date": artifact.extraction_date,
            "source_revision": artifact.source_revision, "source_url": artifact.url,
            "source_sha256": artifact.sha256, "license": artifact.license,
            "license_url": artifact.license_url,
            "manifest_path": str(
                _source_manifest_path(
                    root, artifact, _artifact_contract_sha256(artifact)
                ).relative_to(root)
            ),
            "manifest_sha256": source["manifest_sha256"],
        }
        if source != expected or not _SHA256.fullmatch(source["manifest_sha256"]):
            raise DiseaseProfileError("bundle source provenance differs from the v6 source contract")
    outputs = {}
    for output in manifest["outputs"]:
        table = output["table"]
        if table not in contract.projection.expected_row_counts or table in outputs:
            raise DiseaseProfileError("bundle checksum table contains unknown or duplicate tables")
        if (
            set(output) != {"table", "path", "size_bytes", "sha256", "row_count", "retention"}
            or output["path"] != f"{contract.projection.output_prefix}/{table}.parquet"
            or type(output["row_count"]) is not int
            or output["row_count"] != contract.projection.expected_row_counts[table]
            or type(output["size_bytes"]) is not int or output["size_bytes"] <= 0
            or not _SHA256.fullmatch(output["sha256"]) or output["retention"] != "durable"
        ):
            raise DiseaseProfileError("bundle output checksum table differs from its contract")
        outputs[table] = output
    if set(outputs) != set(contract.projection.expected_row_counts):
        raise DiseaseProfileError("bundle output checksum table is incomplete")
    if manifest["validation"]["row_counts"] != dict(contract.projection.expected_row_counts):
        raise DiseaseProfileError("bundle validation counts differ from the source contract")
    output = outputs["phenotype_assertion"]
    files = ((root / hpo.path, hpo.sha256, hpo.size_bytes),
             (root / output["path"], output["sha256"], output["size_bytes"]))
    for path, digest, size in files:
        _require_within_root(root, path, "public disease profile input")
        if any(character in str(path) for character in "*?[]{}"):
            raise DiseaseProfileError("profile input paths cannot contain glob metacharacters")
        signatures[path] = _file_signature(path)
        if signatures[path][2] != size:
            raise DiseaseProfileError("public disease profile input size differs from completed checksum metadata")
    graph = build_hpo_graph(_parse_obo(files[0][0], ontology.hpo_ontology, retain_unresolved_obsolete=True))
    binding = {
        **identity, "bundle_manifest_path": str(manifest_path.relative_to(root)),
        "bundle_manifest_sha256": expected_bundle_sha256,
        "ontology_contract_sha256": hashlib.sha256(ontology_bytes).hexdigest(),
        "hpo_ontology": asdict(hpo), "hpo_annotations": asdict(hpoa),
        "sources": manifest["sources"], "upstream": upstream, "phenotype_assertion": output,
    }
    releases = {
        "hpo_disease_annotations": contract.hpo.release,
        "orphadata_product4": contract.artifact("orphadata_phenotypes").release,
        "open_targets_disease_phenotype": contract.open_targets.release,
    }
    import duckdb

    try:
        with duckdb.connect(":memory:", config={
            "threads": "1", "temp_directory": "",
            "autoinstall_known_extensions": "false", "autoload_known_extensions": "false",
        }) as connection:
            cursor = connection.execute("SELECT * FROM read_parquet(?, hive_partitioning = false)", [str(files[1][0])])
            expected_schema = tuple(
                tuple(part.strip().replace('"', '').split())
                for part in TABLE_SCHEMAS["phenotype_assertion"].split(",")
            )
            if tuple((item[0], str(item[1])) for item in cursor.description) != expected_schema:
                raise DiseaseProfileError("phenotype Parquet schema differs from v6")
            counts = dict.fromkeys(releases, 0)

            def rows():
                while batch := cursor.fetchmany(2048):
                    for values in batch:
                        row = dict(zip(SOURCE_COLUMNS, values, strict=True))
                        dataset = row["source_dataset"]
                        if dataset not in releases or row["source_release"] != releases[dataset]:
                            raise DiseaseProfileError("phenotype row source/release is not bound to the bundle")
                        counts[dataset] += 1
                        yield row

            store = build_disease_profiles(graph, rows(), source_binding=binding)
    except duckdb.Error:
        raise DiseaseProfileError("cannot read the pinned v6 phenotype Parquet") from None
    if (
        len(store.annotations) != output["row_count"]
        or counts["hpo_disease_annotations"] != contract.hpo.expected_annotation_records
        or counts["open_targets_disease_phenotype"] != contract.validation_expectations["open_targets_hpo_evidence"]
    ):
        raise DiseaseProfileError("phenotype Parquet row counts differ from the source contract")
    if any(_file_signature(path) != before for path, before in signatures.items()):
        raise DiseaseProfileError("disease profile input changed during the read-only build")
    return store
