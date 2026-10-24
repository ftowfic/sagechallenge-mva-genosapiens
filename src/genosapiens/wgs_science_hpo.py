"""Science-only materialization of authenticated completed HPO resources.

Keep the original loader's byte, pin, schema, and coverage checks. Specialize
only the structural decoding of its large annotation collection.
"""
import hashlib
import json
import time
from collections.abc import Mapping
from pathlib import Path
from types import MappingProxyType
from contextvars import ContextVar

from genosapiens.bounded_json import BoundedJSONError, read_bounded_bytes
from genosapiens.manifest import PublicResourceManifest
from genosapiens.ontology_bundle import (
    OntologyResolution, OntologyTerm, ParsedOntology, _require_within_root,
)
from genosapiens.hpo_ic import (
    CorpusAnnotation, HpoGraph, HpoICError, HpoICStore, MAX_STORE_BYTES,
    SCHEMA, TRANSFORM, _freeze_binding, _json_bytes,
)

_METRICS = ContextVar('science_hpo_metrics', default=None)


def _record(name, seconds):
    metrics = _METRICS.get()
    if metrics is not None:
        metrics[name] = metrics.get(name, 0.0) + seconds


def load_bounded_json(path, **kwargs):
    """Keep the ordinary bounded JSON contract and separate read/decode cost."""
    started = time.monotonic()
    content = read_bounded_bytes(path, **kwargs)
    _record('file_read_seconds', time.monotonic()-started)
    def unique(pairs):
        value = {}
        for key, child in pairs:
            if key in value:
                raise BoundedJSONError('JSON input contains duplicate fields')
            value[key] = child
        return value
    started = time.monotonic()
    try:
        return json.loads(content, object_pairs_hook=unique), content
    except (UnicodeDecodeError, json.JSONDecodeError, RecursionError) as error:
        raise BoundedJSONError('JSON input is invalid') from error
    finally:
        _record('json_decode_seconds', time.monotonic()-started)


def _resource_fingerprint(data, aliases):
    """Hash canonical JSON with one resource record resident in the encoder."""
    digest = hashlib.sha256()
    def encoded(value):
        digest.update(json.dumps(value, sort_keys=True, separators=(',', ':'),
            ensure_ascii=True, allow_nan=False).encode('ascii'))
    def mapping(value, *, records=False):
        digest.update(b'{')
        for index, key in enumerate(sorted(value)):
            if index:
                digest.update(b',')
            encoded(key)
            digest.update(b':')
            child = value[key]
            if records and type(child) is list:
                digest.update(b'[')
                for number, record in enumerate(child):
                    if number:
                        digest.update(b',')
                    encoded(record)
                digest.update(b']')
            elif records and type(child) is dict:
                mapping(child)
            else:
                encoded(child)
        digest.update(b'}')
    digest.update(b'{"active_aliases":')
    mapping(aliases)
    digest.update(b',"store":')
    mapping(data, records=True)
    digest.update(b'}\n')
    return digest.hexdigest()


def _text(value):
    if type(value) is not str:
        raise ValueError('invalid completed HPO text')
    value.encode('utf-8', errors='strict')
    return value


def _sequence(value):
    if type(value) not in (list, tuple):
        raise ValueError('invalid completed HPO sequence')
    return value


def _strings(value):
    return tuple(_text(item) for item in _sequence(value))


def _fields(value, names):
    if not isinstance(value, Mapping) or set(value) != names:
        raise ValueError('invalid completed HPO fields')


def _resolution(value):
    _fields(value, {'source_id', 'canonical_id', 'status', 'chain', 'replaced_by', 'consider'})
    return OntologyResolution(_text(value['source_id']),
        None if value['canonical_id'] is None else _text(value['canonical_id']),
        _text(value['status']), _strings(value['chain']),
        _strings(value['replaced_by']), _strings(value['consider']))


def decode_annotations(value):
    """Validate every field and produce the ordinary immutable record types."""
    started = time.monotonic()
    result = []
    for row in _sequence(value):
        _fields(row, {'source_record_id', 'fields', 'resolution', 'eligible_positive'})
        if type(row['eligible_positive']) is not bool:
            raise ValueError('invalid completed HPO eligibility')
        pairs = []
        for pair in _sequence(row['fields']):
            if len(_sequence(pair)) != 2:
                raise ValueError('invalid completed HPO field pair')
            pairs.append((_text(pair[0]), _text(pair[1])))
        result.append(CorpusAnnotation(_text(row['source_record_id']), tuple(pairs),
            _resolution(row['resolution']), row['eligible_positive']))
    _record('annotation_decode_seconds', time.monotonic()-started)
    return tuple(result)


def load_completed_hpo_ic(*, reference_root: str | Path, manifest_path: str | Path,
                          expected_manifest_sha256: str, expected_store_sha256: str,
                          _loaded=None) -> HpoICStore:
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
        store = HpoICStore(graph,decode_annotations(data['annotations']),
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
        if _loaded is not None:
            _loaded(store, data, expected_store_sha256)
        return store
    except (KeyError, TypeError, ValueError):
        raise HpoICError('completed HPO content is malformed') from None


class ScienceHpoLoader:
    """Retain the loaded store binding without serializing its typed tree again.

    Omit the independent typed-store checksum replay under user authorization.
    Compute the ordinary scientific fingerprint from the already decoded data.
    """
    def __init__(self):
        self.binding = None
        self.metrics = {}

    def __call__(self, **kwargs):
        started = time.monotonic()
        token = _METRICS.set(self.metrics)
        try:
            result = load_completed_hpo_ic(**kwargs, _loaded=self._bind)
        finally:
            _METRICS.reset(token)
        self.metrics['total_seconds'] = time.monotonic()-started
        self.metrics['other_graph_index_and_binding_seconds'] = max(0,
            self.metrics['total_seconds'] - sum(v for k, v in self.metrics.items() if k != 'total_seconds'))
        print(json.dumps(dict(schema='genosapiens.science-hpo-initialization/v1',
            **self.metrics)), flush=True)
        return result

    def _bind(self, store, data, store_sha256):
        from genosapiens.track1_sensitivity import (
            _AuthenticatedSensitivityResources, _AUTHENTICATED_RESOURCES_SEAL,
        )
        graph, ontology = store.graph, store.graph.ontology
        started = time.monotonic()
        resource_hash = _resource_fingerprint(data, ontology.active_aliases)
        self.metrics['resource_fingerprint_seconds'] = time.monotonic()-started
        result = object.__new__(_AuthenticatedSensitivityResources)
        binding = (store, graph, ontology, ontology.terms, ontology.active_aliases,
            ontology._alias_targets, graph.parents, graph.ancestors, graph.resolutions,
            store.annotations, store.diseases, store.disease_counts, store.ic,
            store.unresolved_required, store.source_binding, resource_hash, store_sha256)
        for name, value in (('ic_store', store), ('ic_sha256', resource_hash),
                ('store_sha256', store_sha256), ('_seal', _AUTHENTICATED_RESOURCES_SEAL),
                ('_binding', binding)):
            object.__setattr__(result, name, value)
        self.binding = result

    def bind_sensitivity(self, store, *, expected_store_sha256):
        result = self.binding
        if (result is None or result.ic_store is not store
                or result.store_sha256 != expected_store_sha256):
            raise HpoICError('science HPO loader binding unavailable')
        self.binding = None
        return result
