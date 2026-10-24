"""Portable case HPO source and computed IC evidence, retained once per case.

The original HPO graph/corpus builders replay saved source terms and annotation
fields. Every computed graph resolution, closure, count and IC value is also
persisted and compared; an executable identity is not a computed result.
"""
from genosapiens.evidence_model import _typed
from genosapiens.evidence_serialization import _canonical
from genosapiens.hpo_ic import build_hpo_graph, build_hpo_ic
from genosapiens.ontology_bundle import OntologyTerm, ParsedOntology
from genosapiens.wgs_evidence_store import EvidenceStoreError
from genosapiens.wgs_evidence_replay import EvidenceReplay, ScopedSequenceReplay
from genosapiens.wgs_sensitivity import _checked_ic_hash
from genosapiens.wgs_resource_rows import (
    MAX_CHUNK_JSON_BYTES, pack_resource_rows, unpack_resource_rows,
)

MAX_IC_LOGICAL_BYTES = 1024**3
MAX_IC_ROW_BYTES = 8 * 1024**2
MAX_IC_MAP_BYTES = 64 * 1024**2
MAX_IC_TERMS = 100_000
MAX_IC_ANNOTATIONS = 5_000_000


class CaseICEvidence:
    def __init__(self, writer, store, *, expected_ic_sha256, name, authenticated_resources=None):
        self.writer,self.connection,self.name,self.store = writer,writer.connection,name,store
        self.logical_bytes = 0
        try:
            self.ic_sha256 = _checked_ic_hash(store,authenticated_resources)
            if self.ic_sha256 != expected_ic_sha256:
                raise EvidenceStoreError("case IC source identity changed")
            ontology,graph = store.graph.ontology,store.graph
            if len(ontology.terms) > MAX_IC_TERMS or len(store.annotations) > MAX_IC_ANNOTATIONS:
                raise EvidenceStoreError("case IC source counts exceed capacity",capacity=True)
            def observe_size(size):
                self.logical_bytes += size
                if size > MAX_IC_ROW_BYTES or self.logical_bytes > MAX_IC_LOGICAL_BYTES:
                    raise EvidenceStoreError("case IC evidence exceeds capacity",capacity=True)
            roots = {}
            for suffix,kind,values in (
                ("terms","case-ic-terms-v2",iter(ontology.terms.values())),
                ("annotations","case-ic-annotations-v2",iter(store.annotations)),
                ("graph","case-ic-graph-v2",(
                    {"term":key,"parents":parents,"ancestors":tuple(sorted(graph.ancestors[key])),
                     "disease_count":store.disease_counts[key],"ic":store.ic[key]}
                    for key,parents in graph.parents.items())),
                ("resolutions","case-ic-resolutions-v2",iter(graph.resolutions.values())),
            ):
                values = pack_resource_rows(values, max_row_bytes=MAX_IC_ROW_BYTES,
                                            observe_size=observe_size,
                                            _skip_verification=getattr(writer, '_skip_verification',
                                                getattr(getattr(writer, 'reader', None), '_skip_verification', False)))
                digest = writer.add_sequence(name + "." + suffix,kind,values)
                roots[suffix] = {"name":name + "." + suffix,"kind":kind,"sha256":digest}
            aliases = ontology.active_aliases
            aliases_size = len(_canonical(aliases))
            metadata = {"source_binding":store.source_binding,"diseases":store.diseases,"unresolved_required":store.unresolved_required}
            metadata_size = len(_canonical(metadata))
            self.logical_bytes += aliases_size+metadata_size
            if max(aliases_size,metadata_size) > MAX_IC_MAP_BYTES or self.logical_bytes > MAX_IC_LOGICAL_BYTES:
                raise EvidenceStoreError("case IC metadata exceeds capacity",capacity=True)
            aliases_sha = writer.add(name + ".aliases","case-ic-aliases-v1",aliases)
            metadata_sha = writer.add(name + ".metadata","case-ic-metadata-v1",metadata)
            self.sha256 = writer.add(name,"case-ic-evidence-v1",{
                "schema":"genosapiens.wgs-case-ic-evidence","version":2,"ic_sha256":self.ic_sha256,
                "term_count":len(ontology.terms),"annotation_count":len(store.annotations),"logical_bytes":self.logical_bytes,
                **roots,"aliases":{"name":name + ".aliases","kind":"case-ic-aliases-v1","sha256":aliases_sha},
                "metadata":{"name":name + ".metadata","kind":"case-ic-metadata-v1","sha256":metadata_sha}})
        except BaseException:
            writer.failed = True
            raise

    @property
    def reference(self):
        return {"name":self.name,"kind":"case-ic-evidence-v1","sha256":self.sha256,"ic_sha256":self.ic_sha256}

    def require(self, owner, store):
        if self.connection is None or owner.connection is not self.connection or store is not self.store:
            raise EvidenceStoreError("foreign case IC owner")
        completed = getattr(self,'authenticated_resources',None)
        if completed is not None and _checked_ic_hash(store,completed) != self.ic_sha256:
            raise EvidenceStoreError('completed IC owner changed')


class VerifiedCaseICEvidence(CaseICEvidence):
    def __init__(self, reader, reference, *, expected_ic_sha256, authenticated_resources=None):
        if reference.kind != "case-ic-evidence-v1":
            raise EvidenceStoreError("foreign case IC root")
        header = reader.reconstruct(reference,max_bytes=MAX_IC_ROW_BYTES)
        if (set(header) != {"schema","version","ic_sha256","term_count","annotation_count","logical_bytes",
                            "terms","annotations","graph","resolutions","aliases","metadata"}
                or header["schema"] != "genosapiens.wgs-case-ic-evidence" or type(header["version"]) is not int
                or header["version"] != 2 or header["ic_sha256"] != expected_ic_sha256
                or type(header["term_count"]) is not int or not 0 <= header["term_count"] <= MAX_IC_TERMS
                or type(header["annotation_count"]) is not int or not 0 <= header["annotation_count"] <= MAX_IC_ANNOTATIONS
                or type(header["logical_bytes"]) is not int or not 0 <= header["logical_bytes"] <= MAX_IC_LOGICAL_BYTES):
            raise EvidenceStoreError("case IC binding changed")
        def child(key):
            value = header[key]
            kind = "case-ic-" + key + ('-v1' if key in ('aliases', 'metadata') else '-v2')
            if (type(value) is not dict or set(value) != {"name","kind","sha256"}
                    or value["name"] != reference.name + "." + key or value["kind"] != kind):
                raise EvidenceStoreError("case IC child binding changed")
            ref = reader.reference(value["name"],kind=kind)
            if ref.sha256 != value["sha256"]:
                raise EvidenceStoreError("case IC child identity changed")
            return ref
        if authenticated_resources is not None:
            # D15: the original completed resource owns validation of its
            # corpus and computed graph. Bind that exact live immutable owner
            # instead of rebuilding it from its portable case copy.
            from genosapiens.track1_sensitivity import _authenticated_ic_hash
            store = authenticated_resources.ic_store
            if (_authenticated_ic_hash(authenticated_resources,store) != expected_ic_sha256
                    or len(store.graph.ontology.terms) != header['term_count']
                    or len(store.annotations) != header['annotation_count']):
                raise EvidenceStoreError('completed IC resource differs from case')
            references = {key:child(key) for key in
                          ('terms','annotations','graph','resolutions','aliases','metadata')}
            # These roots are covered by completed-resource identity reuse,
            # not by a fresh source-corpus semantic replay.
            for ref in references.values():
                reader._require_reference(ref)
            aliases = reader.reconstruct(references['aliases'],max_bytes=MAX_IC_MAP_BYTES)
            metadata = reader.reconstruct(references['metadata'],max_bytes=MAX_IC_MAP_BYTES)
            if (_canonical(aliases) != _canonical(store.graph.ontology.active_aliases)
                    or _canonical(metadata) != _canonical(dict(source_binding=store.source_binding,
                        diseases=store.diseases,unresolved_required=store.unresolved_required))):
                raise EvidenceStoreError('completed IC source binding differs from case')
            self.writer,self.connection,self.name = EvidenceReplay(reader),reader.connection,reference.name
            self.store,self.ic_sha256,self.sha256 = store,expected_ic_sha256,reference.sha256
            self.logical_bytes = header['logical_bytes']
            self.authenticated_resources = authenticated_resources
            self.validation_origin = 'completed_resource_reuse'
            return
        replay = (ScopedSequenceReplay(reader, max_sequences=2, max_item_bytes=MAX_IC_ROW_BYTES)
                  if reader._read_scope_depth or reader._active_derivative is not None else None)
        def source_rows(key):
            ref = child(key)
            rows = reader.iter_array(ref,max_item_bytes=MAX_CHUNK_JSON_BYTES)
            if replay is not None:
                rows = replay.observe(ref,rows)
            return unpack_resource_rows(rows, max_row_bytes=MAX_IC_ROW_BYTES,
                _skip_verification=reader._skip_verification)
        terms,size = {},0
        for value in source_rows("terms"):
            term = _typed(value,OntologyTerm,decode=True)
            size += len(_canonical(value))
            if term.term_id in terms or len(terms) >= header["term_count"] or size > MAX_IC_LOGICAL_BYTES:
                raise EvidenceStoreError("case IC term coverage changed")
            terms[term.term_id] = term
        if len(terms) != header["term_count"]:
            raise EvidenceStoreError("case IC term count changed")
        aliases = reader.reconstruct(child("aliases"),max_bytes=MAX_IC_MAP_BYTES)
        if type(aliases) is not dict or any(type(k) is not str or type(v) is not str for k,v in aliases.items()):
            raise EvidenceStoreError("invalid case IC aliases")
        graph = build_hpo_graph(ParsedOntology(terms,aliases,0,0))
        metadata = reader.reconstruct(child("metadata"),max_bytes=MAX_IC_MAP_BYTES)
        if set(metadata) != {"source_binding","diseases","unresolved_required"}:
            raise EvidenceStoreError("invalid case IC source metadata")
        def annotations():
            nonlocal size
            count = 0
            for value in source_rows("annotations"):
                count += 1
                size += len(_canonical(value))
                if count > header["annotation_count"] or size > MAX_IC_LOGICAL_BYTES:
                    raise EvidenceStoreError("case IC annotation coverage changed")
                yield dict(value["fields"])
            if count != header["annotation_count"]:
                raise EvidenceStoreError("case IC annotation count changed")
        store = build_hpo_ic(graph,annotations(),source_binding=metadata["source_binding"])
        # Replaying all roots also closes computed values, source annotation
        # identities, aliases and original sequence order against the builders.
        try:
            super().__init__(replay or EvidenceReplay(reader),store,
                             expected_ic_sha256=expected_ic_sha256,name=reference.name)
        finally:
            if replay is not None:
                replay.rows.clear()
