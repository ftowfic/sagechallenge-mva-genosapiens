"""Public, release-bound disease-count HPO information content and Resnik similarity."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import date
from graphlib import CycleError, TopologicalSorter
import hashlib
import json
import math
import os
from pathlib import Path
from types import MappingProxyType
from typing import Iterable, Mapping
import uuid

from genosapiens.bounded_json import BoundedJSONError, load_bounded_json
from genosapiens.manifest import (
    PublicResourceManifest, ResourceOutput, RetentionPolicy, SourceIdentity,
    Transformation, publish_manifest,
)
from genosapiens.ontology_bundle import (
    HPOA_COLUMNS, OBSOLETE_RESOLUTION_POLICY, OntologyBundleContract, OntologyResolution, OntologyTerm,
    ParsedOntology, _file_identity, _parse_comment_table, _parse_obo,
    _require_within_root, _verify_artifact,
)
from genosapiens.storage import DEFAULT_RESERVE_BYTES, CapacityPlan, require_capacity


SCHEMA = "genosapiens.hpo-ic"
MAX_STORE_BYTES = 512 * 1024 * 1024
TRANSFORM = {
    "schema": SCHEMA,
    "obsolete_resolution": OBSOLETE_RESOLUTION_POLICY,
    "corpus": "hpo_disease_annotations",
    "positive_aspect": "P",
    "excluded_qualifier": "NOT",
    "disease_identity": "source-database-id-no-cross-namespace-merging",
    "counting": "distinct-diseases-inclusive-ancestor-propagation",
    "ic": "-ln(term-disease-count/positive-corpus-disease-count)",
    "unseen": "maximum-supported-ancestor-ic-or-null",
    "mica": "all-supported-common-ancestors-with-maximum-ic-sorted-by-id",
    "frequency": "retained-separately-not-used-in-counts",
}


class HpoICError(RuntimeError):
    """Invalid graph, corpus, or immutable resource."""


class HpoICBlocked(HpoICError):
    def __init__(self, store: HpoICStore):
        super().__init__("HPO IC publication blocked: unresolved corpus mappings or empty positive corpus")
        self.store = store


def _json_mapping(value):
    if isinstance(value, Mapping):
        return dict(value)
    raise TypeError("IC source binding must contain JSON values")


def _freeze_binding(value):
    if isinstance(value, Mapping):
        return MappingProxyType({key: _freeze_binding(child) for key, child in value.items()})
    if isinstance(value, (list, tuple)):
        return tuple(_freeze_binding(child) for child in value)
    return value


def _json_bytes(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False,
                       default=_json_mapping) + "\n").encode()


@dataclass(frozen=True, slots=True)
class HpoGraph:
    ontology: ParsedOntology
    parents: Mapping[str, tuple[str, ...]]
    ancestors: Mapping[str, frozenset[str]]
    resolutions: Mapping[str, OntologyResolution]

    def resolve(self, term_id: str) -> OntologyResolution:
        return self.resolutions.get(term_id) or self.ontology.resolve_term(term_id)


def build_hpo_graph(ontology: ParsedOntology) -> HpoGraph:
    ontology = ParsedOntology(
        MappingProxyType(dict(sorted(ontology.terms.items()))),
        MappingProxyType(dict(ontology.active_aliases)),
        ontology.parent_edges, ontology.obsolete_mappings,
    )
    identifiers = set(ontology.terms)
    identifiers.update(alias for term in ontology.terms.values() for alias in term.alt_ids)
    resolutions = {term: ontology.resolve_term(term) for term in sorted(identifiers)}
    parents = {}
    for term_id, term in ontology.terms.items():
        if term.obsolete:
            continue
        if resolutions[term_id].canonical_id != term_id:
            raise HpoICError("HPO active graph contains an ambiguous identifier")
        resolved = [ontology.resolve_term(parent).canonical_id for parent in term.parents]
        if None in resolved:
            raise HpoICError("HPO graph contains an unresolved parent")
        parents[term_id] = tuple(sorted(set(resolved)))
    if not parents:
        raise HpoICError("HPO graph has no active terms")
    ancestors: dict[str, frozenset[str]] = {}
    try:
        for term_id in TopologicalSorter(parents).static_order():
            closure = {term_id}
            for parent in parents[term_id]:
                closure.update(ancestors[parent])
            ancestors[term_id] = frozenset(closure)
    except (CycleError, KeyError):
        raise HpoICError("HPO parent graph is cyclic or incomplete") from None
    return HpoGraph(ontology, MappingProxyType(parents), MappingProxyType(ancestors), MappingProxyType(resolutions))


@dataclass(frozen=True, slots=True)
class CorpusAnnotation:
    source_record_id: str
    fields: tuple[tuple[str, str], ...]
    resolution: OntologyResolution
    eligible_positive: bool

    @property
    def disease_id(self) -> str:
        return dict(self.fields)["database_id"]

    @property
    def frequency(self) -> str | None:
        return dict(self.fields)["frequency"] or None


@dataclass(frozen=True, slots=True)
class ICValue:
    source_id: str
    canonical_id: str | None
    disease_count: int | None
    raw_ic: float | None
    effective_ic: float | None
    support_terms: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class ResnikMatch:
    left: OntologyResolution
    right: OntologyResolution
    mica_ids: tuple[str, ...]
    information_content: float | None

    @property
    def similarity(self) -> float:
        return self.information_content if self.information_content is not None else 0.0


@dataclass(frozen=True, slots=True)
class HpoICStore:
    graph: HpoGraph
    annotations: tuple[CorpusAnnotation, ...]
    diseases: tuple[str, ...]
    disease_counts: Mapping[str, int | None]
    ic: Mapping[str, float | None]
    unresolved_required: tuple[str, ...]
    source_binding: Mapping[str, object]

    @property
    def ready(self) -> bool:
        return bool(self.diseases) and not self.unresolved_required

    def value(self, term_id: str) -> ICValue:
        if not self.ready:
            raise HpoICBlocked(self)
        canonical = self.graph.resolve(term_id).canonical_id
        if canonical is None:
            return ICValue(term_id, None, None, None, None, ())
        supported = [term for term in self.graph.ancestors[canonical] if self.ic[term] is not None]
        maximum = max((self.ic[term] for term in supported), default=None)
        support = tuple(sorted(term for term in supported if self.ic[term] == maximum))
        return ICValue(term_id, canonical, self.disease_counts[canonical], self.ic[canonical], maximum, support)

    def resnik(self, left: str, right: str) -> ResnikMatch:
        if not self.ready:
            raise HpoICBlocked(self)
        a, b = self.graph.resolve(left), self.graph.resolve(right)
        common = (
            self.graph.ancestors[a.canonical_id] & self.graph.ancestors[b.canonical_id]
            if a.canonical_id is not None and b.canonical_id is not None else frozenset()
        )
        supported = [term for term in common if self.ic[term] is not None]
        maximum = max((self.ic[term] for term in supported), default=None)
        return ResnikMatch(a, b, tuple(sorted(term for term in supported if self.ic[term] == maximum)), maximum)


def build_hpo_ic(
    graph: HpoGraph, rows: Iterable[Mapping[str, str]], *, source_binding: Mapping[str, object]
) -> HpoICStore:
    """Keep all annotation evidence; use disease sets, never row/frequency weights."""
    annotations = []
    diseases: set[str] = set()
    term_diseases = {term: set() for term in graph.parents}
    unresolved = []
    for index, row in enumerate(rows, 1):
        if (
            set(row) != set(HPOA_COLUMNS) or any(not isinstance(value, str) for value in row.values())
            or ":" not in row["database_id"] or not row["disease_name"]
            or not all(row[key] for key in ("hpo_id", "reference", "evidence", "aspect", "biocuration"))
            or row["qualifier"] not in {"", "NOT"}
        ):
            raise HpoICError("HPO corpus row is malformed or lacks provenance")
        positive = row["aspect"] == "P" and row["qualifier"] != "NOT"
        resolution = graph.resolve(row["hpo_id"])
        annotation = CorpusAnnotation(f"hpoa:{index}", tuple((key, row[key]) for key in HPOA_COLUMNS), resolution, positive)
        annotations.append(annotation)
        if not positive:
            continue
        diseases.add(row["database_id"])
        if resolution.canonical_id is None:
            unresolved.append(annotation.source_record_id)
            continue
        for ancestor in graph.ancestors[resolution.canonical_id]:
            term_diseases[ancestor].add(row["database_id"])
    ready = bool(diseases) and not unresolved
    counts = {term: len(term_diseases[term]) if ready else None for term in graph.parents}
    ic = {
        term: (0.0 if count == len(diseases) else -math.log(count / len(diseases))) if count else None
        for term, count in counts.items()
    }
    return HpoICStore(
        graph, tuple(annotations), tuple(sorted(diseases)), MappingProxyType(counts),
        MappingProxyType(ic), tuple(unresolved), _freeze_binding(source_binding),
    )


def hpo_ic_as_bytes(store: HpoICStore) -> bytes:
    """Blocked builds can be inspected as diagnostics, but cannot be published."""
    return _json_bytes({
        "schema": SCHEMA, "classification": "public", "transform": TRANSFORM,
        "source_binding": dict(store.source_binding), "ready": store.ready,
        "positive_diseases": store.diseases, "positive_disease_count": len(store.diseases),
        "unresolved_required": store.unresolved_required,
        "terms": [asdict(term) for term in store.graph.ontology.terms.values()],
        "resolutions": [asdict(item) for item in store.graph.resolutions.values()],
        "parents": dict(store.graph.parents),
        "ancestors": {term: sorted(store.graph.ancestors[term]) for term in store.graph.parents},
        "disease_counts": dict(store.disease_counts), "ic": dict(store.ic),
        "annotations": [asdict(item) for item in store.annotations],
    })


def read_pinned_hpo_ic_inputs(
    contract_path: str | Path, reference_root: str | Path
) -> HpoICStore:
    """Read only the pinned public OBO/HPOA, with hashes checked before and after."""
    root = Path(reference_root).resolve(strict=True)
    try:
        contract_data, contract_bytes = load_bounded_json(contract_path)
    except BoundedJSONError:
        raise HpoICError("HPO IC contract is not a valid bounded regular JSON file") from None
    # Parse from the bytes that were hashed, not a second contract-file read.
    contract = OntologyBundleContract.from_mapping(contract_data)
    source_artifacts = (contract.artifact("hpo_ontology"), contract.artifact("hpo_disease_annotations"))
    for artifact in source_artifacts:
        _require_within_root(root, root / artifact.path, "public HPO input")
        _verify_artifact(root / artifact.path, artifact)
    graph = build_hpo_graph(_parse_obo(
        root / source_artifacts[0].path, contract.hpo_ontology, retain_unresolved_obsolete=True,
    ))
    metadata, reader, handle = _parse_comment_table(root / source_artifacts[1].path)
    try:
        if (
            tuple(reader.fieldnames or ()) != HPOA_COLUMNS
            or metadata.get("version") != contract.hpo_annotations.version
            or metadata.get("hpo-version") != contract.hpo_annotations.hpo_version
            or "HPO annotations for rare diseases" not in metadata.get("description", "")
            or not metadata.get("tracker", "").startswith("https://github.com/obophenotype/")
        ):
            raise HpoICError("HPO corpus header or release provenance mismatch")
        binding = {
            "bundle_id": contract.bundle_id, "contract_sha256": hashlib.sha256(contract_bytes).hexdigest(),
            "sources": [asdict(artifact) for artifact in source_artifacts],
        }
        store = build_hpo_ic(graph, reader, source_binding=binding)
    finally:
        handle.close()
    if {item.disease_id.split(":", 1)[0] for item in store.annotations} != set(contract.hpo_annotations.required_database_prefixes):
        raise HpoICError("HPO corpus source namespaces do not match the pin")
    for artifact in source_artifacts:
        _verify_artifact(root / artifact.path, artifact)
    return store


@dataclass(frozen=True, slots=True)
class HpoICPublication:
    store_path: str
    manifest_path: str
    manifest_sha256: str
    positive_diseases: int
    supported_terms: int
    reused: bool


def publish_hpo_ic(
    *, contract_path: str | Path, reference_root: str | Path, image_digest: str,
    retrieved_on: str | None = None, reserve_bytes: int = DEFAULT_RESERVE_BYTES,
) -> HpoICPublication:
    """Rebuild deterministically before immutable publication; no caller-provided store."""
    root = Path(reference_root).resolve(strict=True)
    store = read_pinned_hpo_ic_inputs(contract_path, root)
    if not store.ready:
        raise HpoICBlocked(store)
    content = hpo_ic_as_bytes(store)
    if len(content) > MAX_STORE_BYTES:
        raise HpoICError("HPO IC store exceeds the bounded resource size")
    parameters = hashlib.sha256(_json_bytes({"transform": TRANSFORM, "sources": dict(store.source_binding)})).hexdigest()
    relative = f"derived/phenotype/hpo-ic/v1/{parameters}/store.json"
    target = root / relative
    manifest_path = root / f"manifests/phenotype/hpo-ic-{parameters}.json"
    partial_root = target.parent / ".partial"
    for path in (target, manifest_path, partial_root):
        _require_within_root(root, path, "HPO IC output")
    expected = ResourceOutput(relative, len(content), hashlib.sha256(content).hexdigest(), RetentionPolicy("durable", None))
    sources = store.source_binding["sources"]
    source = sources[1]
    manifest = PublicResourceManifest(
        resource_id=f"hpo-ic-{parameters}", release=source["release"],
        source=SourceIdentity(source["url"], source["sha256"]),
        license=" AND ".join(sorted({item["license"] for item in sources})),
        assembly="not-applicable", retrieved_on=retrieved_on or date.today().isoformat(),
        transformation=Transformation(SCHEMA, parameters), image_digest=image_digest, outputs=(expected,),
    )
    if manifest_path.exists():
        try:
            existing_data, existing_bytes = load_bounded_json(manifest_path)
        except BoundedJSONError:
            raise HpoICError("HPO IC completion manifest is not a bounded regular JSON file") from None
        existing = PublicResourceManifest.from_mapping(existing_data)
        if (existing.resource_id, existing.release, existing.source, existing.license, existing.assembly,
            existing.transformation, existing.outputs) != (
            manifest.resource_id, manifest.release, manifest.source, manifest.license, manifest.assembly,
            manifest.transformation, manifest.outputs,
        ):
            raise HpoICError("completed HPO IC manifest differs from its inputs or transformation")
        if _file_identity(target) != (expected.sha256, expected.size_bytes):
            raise HpoICError("completed HPO IC output checksum mismatch")
        digest = hashlib.sha256(existing_bytes).hexdigest()
        return HpoICPublication(str(target), str(manifest_path), digest, len(store.diseases), sum(bool(n) for n in store.disease_counts.values()), True)
    if target.exists():
        raise HpoICError("HPO IC output lacks its completion manifest; inspect the incomplete publication")
    # Public-only containers mount output subtrees, not the controlled reference root.
    capacity_root = target.parent
    while not capacity_root.exists():
        capacity_root = capacity_root.parent
    require_capacity(capacity_root, CapacityPlan(operation="publish HPO IC", output_bytes=len(content), scratch_bytes=len(content), reserve_bytes=reserve_bytes))
    partial_root.mkdir(parents=True, mode=0o700, exist_ok=True)
    partial = partial_root / f"{uuid.uuid4().hex}.partial"
    try:
        with partial.open("xb") as output:
            output.write(content)
            output.flush()
            os.fsync(output.fileno())
        os.chmod(partial, 0o444)
        os.link(partial, target)
        directory_fd = os.open(target.parent, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
        publication = publish_manifest(manifest, manifest_path, root)
    finally:
        partial.unlink(missing_ok=True)
    return HpoICPublication(str(target), str(manifest_path), publication.manifest_sha256, len(store.diseases), sum(bool(n) for n in store.disease_counts.values()), False)


def load_hpo_ic(store_path: str | Path, *, expected_sha256: str) -> HpoICStore:
    """Verify the bound file and recompute graph/counts instead of trusting stored IC."""
    data, content = load_bounded_json(store_path, maximum_bytes=MAX_STORE_BYTES)
    if hashlib.sha256(content).hexdigest() != expected_sha256:
        raise HpoICError("HPO IC checksum mismatch")
    try:
        if data["schema"] != SCHEMA or data["transform"] != TRANSFORM or not data["ready"]:
            raise HpoICError("HPO IC schema, policy, or publication status is invalid")
        terms = {}
        for item in data["terms"]:
            term = OntologyTerm(**{key: tuple(value) if isinstance(value, list) else value for key, value in item.items()})
            if term.term_id in terms:
                raise HpoICError("HPO IC contains duplicate terms")
            terms[term.term_id] = term
        aliases = {alias: term.term_id for term in terms.values() if not term.obsolete for alias in term.alt_ids}
        graph = build_hpo_graph(ParsedOntology(terms, aliases, 0, 0))
        store = build_hpo_ic(graph, (dict(item["fields"]) for item in data["annotations"]), source_binding=data["source_binding"])
        if not store.ready or hpo_ic_as_bytes(store) != content:
            raise HpoICError("HPO IC content is noncanonical or fails recomputation")
        return store
    except (KeyError, TypeError, ValueError):
        raise HpoICError("HPO IC content is malformed") from None


def load_completed_hpo_ic(*, reference_root: str | Path, manifest_path: str | Path,
                          expected_manifest_sha256: str, expected_store_sha256: str) -> HpoICStore:
    """Load the original completed resource without rebuilding its graph/corpus.

    The separately pinned producer manifest owns scientific validation. Decode
    and check the bytes needed to materialize this immutable scoring resource;
    do not parse original OBO/HPOA or repeat ancestor/count/IC computations.
    """
    from genosapiens.evidence_model import _typed
    root, manifest_path = Path(reference_root), Path(manifest_path)
    _require_within_root(root,manifest_path,'completed HPO manifest')
    metadata, manifest_bytes = load_bounded_json(manifest_path)
    if hashlib.sha256(manifest_bytes).hexdigest() != expected_manifest_sha256:
        raise HpoICError('completed HPO manifest identity changed')
    manifest = PublicResourceManifest.from_mapping(metadata)
    if (manifest.transformation.schema != SCHEMA or len(manifest.outputs) != 1
            or manifest.outputs[0].sha256 != expected_store_sha256
            or manifest.outputs[0].size_bytes > MAX_STORE_BYTES):
        raise HpoICError('completed HPO resource binding changed')
    output = manifest.outputs[0]
    path = root / output.path
    _require_within_root(root,path,'completed HPO store')
    data, content = load_bounded_json(path,maximum_bytes=MAX_STORE_BYTES)
    if len(content) != output.size_bytes or hashlib.sha256(content).hexdigest() != output.sha256:
        raise HpoICError('consumed completed HPO bytes changed')
    try:
        if (set(data) != {'schema','classification','transform','source_binding','ready',
                'positive_diseases','positive_disease_count','unresolved_required','terms',
                'resolutions','parents','ancestors','disease_counts','ic','annotations'}
                or data['schema'] != SCHEMA or data['classification'] != 'public'
                or data['transform'] != TRANSFORM or data['ready'] is not True):
            raise HpoICError('completed HPO schema or readiness changed')
        terms = _typed(data['terms'],tuple[OntologyTerm,...],decode=True)
        term_map = {term.term_id:term for term in terms}
        resolutions = _typed(data['resolutions'],tuple[OntologyResolution,...],decode=True)
        resolution_map = {value.source_id:value for value in resolutions}
        if len(term_map) != len(terms) or len(resolution_map) != len(resolutions):
            raise HpoICError('completed HPO identifiers are duplicated')
        aliases = {alias:term.term_id for term in terms if not term.obsolete for alias in term.alt_ids}
        ontology = ParsedOntology(MappingProxyType(term_map),MappingProxyType(aliases),0,0)
        parents = _typed(data['parents'],Mapping[str,tuple[str,...]],decode=True)
        ancestors = _typed(data['ancestors'],Mapping[str,tuple[str,...]],decode=True)
        graph = HpoGraph(ontology,parents,MappingProxyType({k:frozenset(v) for k,v in ancestors.items()}),
                         MappingProxyType(resolution_map))
        store = HpoICStore(graph,_typed(data['annotations'],tuple[CorpusAnnotation,...],decode=True),
            _typed(data['positive_diseases'],tuple[str,...],decode=True),
            _typed(data['disease_counts'],Mapping[str,int | None],decode=True),
            _typed(data['ic'],Mapping[str,float | None],decode=True),
            _typed(data['unresolved_required'],tuple[str,...],decode=True),_freeze_binding(data['source_binding']))
        parameters = hashlib.sha256(_json_bytes({'transform':TRANSFORM,'sources':dict(store.source_binding)})).hexdigest()
        if (not store.ready or type(data['positive_disease_count']) is not int
                or len(store.diseases) != data['positive_disease_count']
                or parameters != manifest.transformation.parameters_sha256
                or set(parents) != set(ancestors) or set(parents) != set(store.ic)
                or set(parents) != set(store.disease_counts)):
            raise HpoICError('completed HPO metadata coverage changed')
        return store
    except (KeyError, TypeError, ValueError):
        raise HpoICError('completed HPO content is malformed') from None
