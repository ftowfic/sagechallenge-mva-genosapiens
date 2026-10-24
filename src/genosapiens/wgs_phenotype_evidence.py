"""Case-wide observation/profile ownership and source-bound phenotype results.

The existing phenotype scorer still consumes only its closed observation view,
resource owners and disease identifier. No variant, pair or diagnosis input is
introduced into primary scoring. Source/profile and IC replay occur once per
case owner, before context results can be published or restored.
"""
import hashlib
import json

from genosapiens.evidence_model import (
    DiseaseProfile, EvidenceDocument, EvidenceReference, HypothesisContext, ScoringView, _typed,
)
from genosapiens.evidence_serialization import _canonical, _owned_canonical_bytes, _owned_json_default
from genosapiens.phenotype_scoring import prepare_phenotype_scores
from genosapiens.wgs_evidence_replay import EvidenceReplay, ScopedSequenceReplay
from genosapiens.wgs_evidence_store import EvidenceStoreError, MAX_LOGICAL_BYTES
from genosapiens.wgs_sensitivity import _profiles
from genosapiens.wgs_ic_evidence import CaseICEvidence, VerifiedCaseICEvidence
from genosapiens.wgs_resource_rows import (
    MAX_CHUNK_JSON_BYTES, pack_resource_rows, unpack_resource_rows,
)

MAX_PHENOTYPE_RESOURCE_BYTES = 128 * 1024**2
MAX_CONTEXT_PHENOTYPE_BYTES = 64 * 1024**2


def _hash(value):
    if type(value) is DiseaseProfile:
        # Hash the same canonical bytes with bounded record encodings. The C
        # encoder handles each record without a Python loop for every token.
        digest = hashlib.sha256()
        for part in _profile_hash_parts(value):
            digest.update(part.encode('ascii'))
        digest.update(b'\n')
        return digest.hexdigest()
    return hashlib.sha256((_canonical(value)+"\n").encode("ascii")).hexdigest()


def _profile_hash_parts(profile):
    """Canonical profile stream. Retain at most one encoded resource record."""
    encoder = json.JSONEncoder(sort_keys=True, separators=(',', ':'), ensure_ascii=True,
                               allow_nan=False, default=_owned_json_default)

    def mapping(values, encode):
        yield '{'
        for ordinal, key in enumerate(sorted(values)):
            if ordinal:
                yield ','
            yield encoder.encode(key)
            yield ':'
            yield from encode(key, values[key])
        yield '}'

    def leaf(_key, value):
        yield encoder.encode(value)

    def store_field(key, value):
        if key == 'profiles':
            yield from mapping(value, leaf)
        elif key in ('annotations', 'unassigned_annotation_ids'):
            yield '['
            for ordinal, row in enumerate(value):
                if ordinal:
                    yield ','
                yield encoder.encode(row)
            yield ']'
        else:
            yield encoder.encode(value)

    def field(key, value):
        if key == 'store':
            yield from mapping(_owned_json_default(value), store_field)
        else:
            yield encoder.encode(value)

    yield from mapping(_owned_json_default(profile), field)


_PROFILE_FIELDS = ('annotations', 'profiles', 'unassigned_annotation_ids')
_PROFILE_KIND = 'case-phenotype-profile-v1'
_PROFILE_ROWS_KIND = 'case-phenotype-profile-records-v2'


class _ProfileReplay(ScopedSequenceReplay):
    """Three live row commitments confined to the enclosing authenticated read."""
    def __init__(self, reader):
        super().__init__(reader, max_sequences=3, max_item_bytes=MAX_PHENOTYPE_RESOURCE_BYTES)


def _write_profile(writer, profile, name):
    store = dict(source_binding=profile.store.source_binding, graph_sha256=profile.store.graph_sha256)
    for field in _PROFILE_FIELDS:
        values = getattr(profile.store, field)
        rows = ((key, values[key]) for key in sorted(values)) if field == 'profiles' else iter(values)
        child = name + '.' + field
        store[field] = dict(name=child, count=len(values),
            sha256=writer.add_sequence(child, _PROFILE_ROWS_KIND,
                # The original profile item limit includes its final newline.
                pack_resource_rows(rows, max_row_bytes=MAX_PHENOTYPE_RESOURCE_BYTES - 1,
                    _skip_verification=getattr(writer, '_skip_verification',
                        getattr(getattr(writer, 'reader', None), '_skip_verification', False)))))
    metadata = dict(schema='genosapiens.wgs-case-phenotype-profile', version=2,
                    evidence_id=profile.evidence_id, public_manifest_sha256=profile.public_manifest_sha256,
                    store=store)
    if len(_owned_canonical_bytes(metadata)) > MAX_PHENOTYPE_RESOURCE_BYTES:
        raise EvidenceStoreError('phenotype profile metadata exceeds capacity', capacity=True)
    return dict(name=name, kind=_PROFILE_KIND, sha256=writer.add(name, _PROFILE_KIND, metadata))


def _read_profile_rows(reader, descriptor, name, replay=None):
    if (type(descriptor) is not dict or set(descriptor) != {'name', 'sha256', 'count'}
            or descriptor['name'] != name or type(descriptor['count']) is not int
            or not 0 <= descriptor['count'] <= MAX_LOGICAL_BYTES):
        raise EvidenceStoreError('invalid phenotype resource sequence descriptor')
    reference = reader.reference(name, kind=_PROFILE_ROWS_KIND)
    if reference.sha256 != descriptor['sha256']:
        raise EvidenceStoreError('phenotype resource child changed')
    count = 0
    rows = reader.iter_array(reference, max_item_bytes=MAX_CHUNK_JSON_BYTES)
    if replay is not None:
        rows = replay.observe(reference, rows)
    rows = unpack_resource_rows(rows, max_row_bytes=MAX_PHENOTYPE_RESOURCE_BYTES - 1,
        _skip_verification=reader._skip_verification)
    for row in rows:
        count += 1
        if count > descriptor['count']:
            raise EvidenceStoreError('phenotype resource sequence count changed')
        yield row
    if count != descriptor['count']:
        raise EvidenceStoreError('phenotype resource sequence count changed')


def _read_profile(reader, descriptor, name, replay=None):
    if (type(descriptor) is not dict or set(descriptor) != {'name', 'kind', 'sha256'}
            or descriptor['name'] != name or descriptor['kind'] != _PROFILE_KIND):
        raise EvidenceStoreError('invalid phenotype profile reference')
    reference = reader.reference(name, kind=_PROFILE_KIND)
    if reference.sha256 != descriptor['sha256']:
        raise EvidenceStoreError('phenotype profile root changed')
    metadata = reader.reconstruct(reference, max_bytes=MAX_PHENOTYPE_RESOURCE_BYTES)
    if (type(metadata) is not dict or set(metadata) != {'schema', 'version', 'evidence_id', 'public_manifest_sha256', 'store'}
            or metadata['schema'] != 'genosapiens.wgs-case-phenotype-profile'
            or type(metadata['version']) is not int or metadata['version'] != 2
            or type(metadata['store']) is not dict
            or set(metadata['store']) != {*_PROFILE_FIELDS, 'source_binding', 'graph_sha256'}):
        raise EvidenceStoreError('invalid phenotype profile metadata')
    store = metadata['store']
    for field in _PROFILE_FIELDS:
        rows = _read_profile_rows(reader, store[field], name + '.' + field, replay)
        if field == 'profiles':
            values = {}
            previous = None
            for row in rows:
                if (type(row) is not list or len(row) != 2 or type(row[0]) is not str
                        or (previous is not None and row[0] <= previous)):
                    raise EvidenceStoreError('phenotype profile mapping order changed')
                previous = row[0]
                values[row[0]] = row[1]
            store[field] = values
        else:
            store[field] = tuple(rows)
    return _typed({key:metadata[key] for key in ('evidence_id', 'public_manifest_sha256', 'store')},
                  DiseaseProfile, decode=True)


class CasePhenotypeEvidence:
    def __init__(self, writer, *, run_manifest_sha256, observations, view, profile,
                 ic_store, expected_ic_sha256, name, authenticated_resources=None, ic_resource=None):
        self.writer,self.connection,self.name = writer,writer.connection,name
        self.run_manifest_sha256 = run_manifest_sha256
        skip_verification = getattr(writer, "_skip_verification",
            getattr(getattr(writer, "reader", None), "_skip_verification", False))
        try:
            if skip_verification:
                from types import SimpleNamespace
                document = SimpleNamespace(observations=observations,
                    disease_profiles=() if profile is None else (profile,))
            else:
                view = _typed(view,ScoringView)
                document = EvidenceDocument(run_manifest_sha256,observations=observations,contradictions=view.contradictions,
                                            disease_profiles=() if profile is None else (profile,))
                if document.scoring_view() != view:
                    raise EvidenceStoreError("case phenotype observation view changed")
            self.view,self.observations = view,document.observations
            self.profile = None if not document.disease_profiles else document.disease_profiles[0]
            self.ic_store,self.authenticated_resources = ic_store,authenticated_resources
            self.ic_resource = ic_resource or CaseICEvidence(writer,ic_store,expected_ic_sha256=expected_ic_sha256,
                name=name + ".ic",authenticated_resources=authenticated_resources)
            self.ic_resource.require(writer,ic_store)
            self.ic_sha256 = self.ic_resource.ic_sha256
            if self.ic_sha256 != expected_ic_sha256:
                raise EvidenceStoreError("case phenotype IC identity changed")
            profile_store = ((None if self.profile is None else self.profile.store) if skip_verification else
                _profiles(self.profile,profile_evidence_id=None if self.profile is None else self.profile.evidence_id,
                    ic_store=ic_store))
            self.prepared = prepare_phenotype_scores(view.observations,ic_store=ic_store,profile_store=profile_store,
                prior_identity="diagnosis-blind-primary:" + _hash(view),
                _skip_verification=skip_verification,
                _ic_store_sha256=getattr(authenticated_resources, 'store_sha256', None))
            profile_value = (None if self.profile is None else
                             _write_profile(writer, self.profile, name + '.profile'))
            self.view_sha256, self.profile_sha256 = _hash(view), _hash(self.profile)
            self._logical_owners = self._logical_resource_owners()
            header = {"schema":"genosapiens.wgs-case-phenotype-evidence","version":3,
                "run_manifest_sha256":run_manifest_sha256,"ic_sha256":self.ic_sha256,
                "ic_evidence":self.ic_resource.reference,
                "view":view,"observations":document.observations,"profile":profile_value,
                "view_sha256":self.view_sha256,"profile_sha256":self.profile_sha256}
            if len(_canonical(header)) > MAX_PHENOTYPE_RESOURCE_BYTES:
                raise EvidenceStoreError("case phenotype resources exceed capacity",capacity=True)
            self.sha256 = writer.add(name,"case-phenotype-evidence-v1",header)
        except BaseException:
            writer.failed = True
            raise

    @property
    def reference(self):
        return {"name":self.name,"kind":"case-phenotype-evidence-v1","sha256":self.sha256}

    def _logical_resource_owners(self):
        from dataclasses import fields
        objects = (self.view,) if self.profile is None else (self.view, self.profile, self.profile.store)
        return (self.view_sha256, self.profile_sha256, *objects,
            *(getattr(owner, field.name) for owner in objects for field in fields(owner)))

    def logical_hashes(self, writer):
        """Reuse hashes of this producer's private immutable resource copies."""
        owners = self._logical_resource_owners()
        if (writer.connection is not self.connection or self.connection is None
                or len(owners) != len(self._logical_owners)
                or any(a is not b for a, b in zip(owners, self._logical_owners, strict=True))):
            raise EvidenceStoreError('case phenotype logical resource owner changed')
        return self.view_sha256, self.profile_sha256

    def score(self, context):
        if self.connection is None or self.writer.connection is not self.connection:
            raise EvidenceStoreError("case phenotype owner is closed")
        context = _typed(context,HypothesisContext)
        score = self.prepared.score(context.disease_id)
        if len(score.as_bytes()) > MAX_CONTEXT_PHENOTYPE_BYTES:
            raise EvidenceStoreError("context phenotype evidence exceeds capacity",capacity=True)
        return score

    def profile_refs(self, score):
        return () if self.profile is None else tuple(EvidenceReference(
            self.profile.evidence_id,item["annotation_id"],self.profile.public_manifest_sha256)
            for item in score.ledger["profile_annotations"])


class VerifiedCasePhenotypeEvidence(CasePhenotypeEvidence):
    def __init__(self, reader, reference, *, run_manifest_sha256, expected_view_sha256,
                 expected_profile_sha256, expected_ic_sha256, ic_store=None, authenticated_resources=None, ic_resource=None):
        from genosapiens.evidence_model import Observation
        if reference.kind != "case-phenotype-evidence-v1":
            raise EvidenceStoreError("foreign case phenotype root")
        header = reader.reconstruct(reference,max_bytes=MAX_PHENOTYPE_RESOURCE_BYTES)
        if (set(header) != {"schema","version","run_manifest_sha256","ic_sha256","ic_evidence","view","observations",
                            "profile","view_sha256","profile_sha256"}
                or header["schema"] != "genosapiens.wgs-case-phenotype-evidence"
                or type(header["version"]) is not int or header["version"] != 3
                or header["run_manifest_sha256"] != run_manifest_sha256
                or header["view_sha256"] != expected_view_sha256 or header["profile_sha256"] != expected_profile_sha256
                or header["ic_sha256"] != expected_ic_sha256):
            raise EvidenceStoreError("case phenotype source binding changed")
        view = _typed(header["view"],ScoringView,decode=True)
        replay = (_ProfileReplay(reader) if reader._read_scope_depth or reader._active_derivative is not None
                  else None)
        profile = (None if header['profile'] is None else
                   _read_profile(reader, header['profile'], reference.name + '.profile', replay))
        if _hash(view) != expected_view_sha256 or _hash(profile) != expected_profile_sha256:
            raise EvidenceStoreError("case phenotype logical identity changed")
        dependency = header["ic_evidence"]
        if (type(dependency) is not dict or set(dependency) != {"name","kind","sha256","ic_sha256"}
                or dependency["kind"] != "case-ic-evidence-v1" or dependency["ic_sha256"] != expected_ic_sha256):
            raise EvidenceStoreError("case phenotype IC dependency changed")
        ic_ref = reader.reference(dependency["name"],kind="case-ic-evidence-v1")
        if ic_ref.sha256 != dependency["sha256"]:
            raise EvidenceStoreError("case phenotype IC evidence changed")
        if ic_resource is None:
            if authenticated_resources is not None:
                if ic_store is not authenticated_resources.ic_store:
                    raise EvidenceStoreError('foreign completed phenotype IC owner')
                ic_resource = VerifiedCaseICEvidence(reader,ic_ref,expected_ic_sha256=expected_ic_sha256,
                                                    authenticated_resources=authenticated_resources)
            else:
                ic_resource = (VerifiedCaseICEvidence(reader,ic_ref,expected_ic_sha256=expected_ic_sha256)
                if ic_store is None else CaseICEvidence(EvidenceReplay(reader),ic_store,
                    expected_ic_sha256=expected_ic_sha256,name=ic_ref.name,authenticated_resources=authenticated_resources))
        if ic_resource.reference != dependency:
            raise EvidenceStoreError("foreign phenotype IC resource owner")
        try:
            super().__init__(replay or EvidenceReplay(reader),run_manifest_sha256=run_manifest_sha256,
                observations=_typed(header["observations"],tuple[Observation,...],decode=True),view=view,profile=profile,
                ic_store=ic_resource.store,expected_ic_sha256=expected_ic_sha256,name=reference.name,
                authenticated_resources=authenticated_resources,ic_resource=ic_resource)
        finally:
            if replay is not None:
                replay.rows.clear()


class ContextPhenotypeEvidence:
    def __init__(self, writer, resources, context, *, name):
        if resources.connection is not writer.connection:
            raise EvidenceStoreError("foreign phenotype resource owner")
        self.name,self.connection,self.resources = name,writer.connection,resources
        try:
            self.context_sha256 = _hash(context)
            self.score = resources.score(context)
            self.profile_refs = resources.profile_refs(self.score)
            self.sha256 = writer.add(name,"context-phenotype-evidence-v1",{
                "schema":"genosapiens.wgs-context-phenotype-evidence","version":1,
                "resources":resources.reference,"context_sha256":self.context_sha256,
                "score":self.score.ledger,"profile_refs":self.profile_refs})
        except BaseException:
            writer.failed = True
            raise

    @property
    def reference(self):
        return {"name":self.name,"kind":"context-phenotype-evidence-v1","sha256":self.sha256}

    def require(self, writer, context, score, refs):
        if (writer.connection is not self.connection or self.connection is None
                or _hash(context) != self.context_sha256 or score.as_bytes() != self.score.as_bytes()
                or refs != self.profile_refs):
            raise EvidenceStoreError("context phenotype dependency changed")


class VerifiedContextPhenotypeEvidence(ContextPhenotypeEvidence):
    def __init__(self, reader, reference, resources, context):
        if reference.kind != "context-phenotype-evidence-v1":
            raise EvidenceStoreError("foreign context phenotype root")
        reader._require_reference(reference)
        super().__init__(EvidenceReplay(reader),resources,context,name=reference.name)
