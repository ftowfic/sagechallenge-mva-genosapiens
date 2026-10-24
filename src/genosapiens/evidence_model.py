"""Task 9.1 controlled evidence, not a ranker, genotype caller, or handoff.

Public profile/assertion stores and native allele identities retain their own
vocabularies. Construction freezes containers and checks structure, not the
truth of source assertions. No files or services are accessed by this module.
"""

from __future__ import annotations

from collections.abc import Mapping
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass, fields, is_dataclass
from enum import Enum, StrEnum
from fractions import Fraction
from functools import lru_cache
from graphlib import CycleError, TopologicalSorter
import math
import re
from types import MappingProxyType, UnionType
from typing import TYPE_CHECKING, ClassVar, Literal, Union, get_args, get_origin, get_type_hints

from genosapiens.allele_identity import (
    AlleleIdentity, identity_from_dict, identity_to_dict,
)
from genosapiens.disease_gene_evidence import (
    DiseaseGeneAggregation, DiseaseGeneContext, DiseaseGeneEvidenceError,
    DiseaseGeneEvidenceStore, GeneAssertion, GeneRequirements,
    NATIVE_LOCUS_TABLES, assertion_native_disease_id,
)
from genosapiens.disease_profiles import DiseaseProfileStore
from genosapiens.governance import ArtifactClassification
from genosapiens.llm_guard import LLMPolicy
from genosapiens.ontology_bundle import (
    OBSOLETE_RESOLUTION_POLICY, OntologyResolution,
)
from genosapiens.phenotype_scoring import ObservationProvenance, ScoringObservation
from genosapiens.rare_disease import TABLE_SCHEMAS
from genosapiens.vcf_qc import SampleCallQc

if TYPE_CHECKING:
    from genosapiens.phenotype_ingestion import PhenotypeIngestionResult


SCHEMA = "genosapiens.track1-evidence-model"


class EvidenceModelError(ValueError):
    """Static diagnostics only; serialization and explicit field access are sensitive."""


# Validation owns frozen copies, never trusts an incoming object's frozen flag.
# Entries hold both input and result strongly and are published only after the
# entire constructor (including subclass semantic checks) returns successfully.
# External codecs isolate each row; no identity survives that operation.
_VALIDATION_WALK: ContextVar[
    tuple[set[tuple], dict[tuple, tuple[object, object]]] | None
] = ContextVar("evidence_validation_walk", default=None)

# The decoder can hand raw field values to a known generated constructor so
# that its ordinary base/subclass checks perform the single strict field walk.
# The base consumes this identity-bound handoff before validating descendants.
_DECODE_CONSTRUCTION: ContextVar[tuple[type, Mapping] | None] = ContextVar(
    "evidence_decode_construction", default=None)


@contextmanager
def _validation_session(*, isolated=False):
    """Reuse exact immutable-child validation within one synchronous operation."""
    state = _VALIDATION_WALK.get()
    owner = isolated or state is None
    token = None
    if owner:
        state = (set(), {})
        token = _VALIDATION_WALK.set(state)
    try:
        yield
    finally:
        if owner:
            state[0].clear()
            state[1].clear()
            _VALIDATION_WALK.reset(token)


def _fail() -> None:
    raise EvidenceModelError("invalid controlled evidence structure")


@lru_cache(maxsize=None)
def _hints(cls):
    return get_type_hints(cls)


@lru_cache(maxsize=512)
def _field_contract(cls):
    hints=_hints(cls)
    return tuple((item,hints[item.name]) for item in fields(cls))


@lru_cache(maxsize=512)
def _single_walk_constructor(cls):
    # These current records use generated assignment-only __init__ methods and
    # call _Controlled first in every subclass post-init. Unknown/custom record
    # constructors retain the ordinary typed-argument construction path.
    return (isinstance(cls,type) and issubclass(cls,_Controlled)
        and cls.__module__ in ('genosapiens.evidence_model',
            'genosapiens.annotation_admission','genosapiens.wgs_admission')
        and getattr(getattr(cls.__init__,'__code__',None),'co_filename',None)=='<string>'
        and all(item.init for item,_ in _field_contract(cls)))


def _type_shape(hint):
    return _cached_type_shape(id(hint), hint)


@lru_cache(maxsize=512)
def _cached_type_shape(_identity, hint):
    # Type metadata only; never retain a record or its validation result here.
    # Equal unions can retain different argument orders. Identity in this key
    # preserves the caller's exact branch order; hint holds it strongly.
    origin,args=get_origin(hint),get_args(hint)
    scalar=None
    if origin in (UnionType,Union) and len(args)==2 and type(None) in args:
        other=next(item for item in args if item is not type(None))
        if other in (str,int,float,bool):scalar=other
    return origin,args,is_dataclass(hint),scalar


@lru_cache(maxsize=512)
def _constructor_contract(cls):
    """Resolve scalar/null dispatch once per type, retaining no field values.

    Every value still crosses the ordinary scalar or structural validator.
    Ordered unions and containers keep the existing owner/copy/rollback path.
    """
    result = []
    for item, hint in _field_contract(cls):
        origin, args, _, optional_scalar = _type_shape(hint)
        scalar = hint if hint in (str, int, float, bool, type(None)) else optional_scalar
        nullable = origin in (UnionType, Union) and type(None) in args
        result.append((item.name, hint, scalar, nullable))
    return tuple(result)


def _json_value(value):
    # Scalars dominate source fields. Keep the same exact-type/UTF-8/finite
    # checks without asking the Mapping ABC about every scalar leaf.
    if value is None or type(value) in (bool, int):
        return value
    if type(value) is str:
        value.encode("utf-8", errors="strict")
        return value
    if type(value) is float:
        if not math.isfinite(value):
            _fail()
        return value
    if isinstance(value, Mapping):
        if any(type(key) is not str for key in value):
            _fail()
        return MappingProxyType({_json_value(key): _typed_value(child, object)
                                 for key, child in value.items()})
    if type(value) in (list, tuple):
        return tuple(_typed_value(child, object) for child in value)
    _fail()


def _typed_scalar(value, hint):
    """The ordinary scalar contract, including exact numeric representation."""
    if hint is float:
        if type(value) not in (int, float) or not math.isfinite(value):
            _fail()
        if type(value) is int and not -(2**63) <= value < 2**63:
            _fail()
        return value
    if type(value) is not hint:
        _fail()
    if hint is int and not -(2**63) <= value < 2**63:
        _fail()
    if hint is str:
        value.encode("utf-8", errors="strict")
    return value


def _typed(value, hint, *, decode=False, _reuse_existing=False):
    """One strict type walker for constructors and codecs; no field inference."""
    # Scalars cannot own structural descendants or validation-cache entries.
    # Use the same exact validator without opening an empty ownership scope.
    if hint is str or hint is int or hint is bool or hint is float or hint is type(None):
        return _typed_scalar(value, hint)
    if decode:
        # A decoded row owns its validation walk. An enclosing execution may
        # reuse constructed immutable children, but must not retain every row
        # subsequently read from disk until the entire execution finishes.
        with _validation_session(isolated=True):
            return _typed_value(value, hint, decode=True,
                                _reuse_existing=_reuse_existing)
    with _validation_session():
        return _typed_value(value, hint, _reuse_existing=_reuse_existing)


def _construct_owned(cls, /, *args, **kwargs):
    """Construct a private immutable record, publishing only after all checks.

    The ordinary constructor copies every field and runs subclass semantics.
    Builders already inside a bounded operation can reuse that result without
    immediately cloning it through the same constructor again. An ordinary call
    owns and clears its own scope; public codecs always revalidate independently.
    """
    if not isinstance(cls,type) or not issubclass(cls,_Controlled):_fail()
    with _validation_session():
        _,complete=_VALIDATION_WALK.get()
        checkpoint=len(complete)
        try:
            result=cls(*args,**kwargs)
            if type(result) is not cls:_fail()
            complete[(id(result),cls,False)]=(result,result)
            return result
        except BaseException:
            while len(complete)>checkpoint:complete.popitem()
            raise


def _typed_value(value, hint, *, decode=False, _reuse_existing=False):
    # Only supported structural branches participate. Primitive validation and
    # unsupported types never acquire cache entries. The clone's entry avoids
    # walking its descendants again in its parent's validating constructor.
    if hint is str or hint is int or hint is bool or hint is float or hint is type(None):
        return _typed_scalar(value, hint)
    if hint is object and (value is None or type(value) in (str, int, float, bool)):
        return _json_value(value)
    # Raw JSON containers obey the same copied-child ownership as typed
    # containers. A caller's proxy is never trusted merely for being read-only;
    # only this operation's deep immutable copy is reused and then discarded.
    origin,args,record_type,optional_scalar = _type_shape(hint)
    if optional_scalar is not None:
        if value is None:return None
        try:return _typed_scalar(value,optional_scalar)
        except (ValueError,TypeError,RuntimeError,AttributeError):_fail()
    if value is None and origin in (UnionType, Union) and type(None) in args:
        return None
    structural = hint is object or record_type or origin in (tuple, Mapping)
    state = _VALIDATION_WALK.get()
    if state is None or not structural:
        return _typed_uncached(value, hint, decode=decode)
    active, complete = state
    key = (id(value), hint, decode)
    entry = complete.get(key)
    if entry is not None:
        if entry[0] is not value:
            _fail()
        return entry[1]
    if key in active:
        _fail()
    checkpoint = len(complete)
    active.add(key)
    try:
        result = _typed_uncached(value, hint, decode=decode)
        complete[key] = (value, result)
        complete[(id(result), hint, False)] = (result, result)
        return result
    except BaseException:
        # Discard descendants first encountered in this failed branch. Keep
        # previously validated siblings: a union commonly rejects one option
        # before accepting another, including the ordinary optional None case.
        # Insertion order makes rollback proportional to the failed branch.
        while len(complete) > checkpoint:
            complete.popitem()
        raise
    finally:
        active.discard(key)


def _typed_uncached(value, hint, *, decode=False):
    _reuse_existing = False
    origin, args, record_type, _ = _type_shape(hint)
    if hint == AlleleIdentity:
        identity = identity_from_dict(value if decode else identity_to_dict(value))
        return _typed_value(identity, type(identity))
    if origin in (UnionType, Union):
        for option in args:
            try:
                return _typed_value(
                    value,
                    option,
                    decode=decode,
                    _reuse_existing=_reuse_existing,
                )
            except (ValueError, TypeError, RuntimeError, AttributeError):
                pass
        _fail()
    if hint is type(None):
        if value is not None:
            _fail()
        return None
    if hint is object:
        return _json_value(value)
    if origin is Literal:
        if not any(type(value) is type(item) and value == item for item in args):
            _fail()
        return value
    if origin is tuple:
        if type(value) is not tuple and not (decode and type(value) is list):
            _fail()
        if len(args) == 2 and args[1] is Ellipsis:
            return tuple(
                _typed_value(
                    item,
                    args[0],
                    decode=decode,
                    _reuse_existing=_reuse_existing,
                )
                for item in value
            )
        if len(value) != len(args):
            _fail()
        return tuple(
            _typed_value(
                item,
                subtype,
                decode=decode,
                _reuse_existing=_reuse_existing,
            )
            for item, subtype in zip(value, args, strict=True)
        )
    if origin is Mapping:
        if not isinstance(value, Mapping):
            _fail()
        return MappingProxyType({
            _typed_value(
                key,
                args[0],
                decode=decode,
                _reuse_existing=_reuse_existing,
            ): _typed_value(
                child,
                args[1],
                decode=decode,
                _reuse_existing=_reuse_existing,
            )
            for key, child in value.items()
        })
    if isinstance(hint, type) and issubclass(hint, Enum):
        if decode:
            if type(value) is not str:
                _fail()
            return hint(value)
        if type(value) is not hint:
            _fail()
        return value
    if record_type:
        contract=_field_contract(hint)
        if decode:
            if not isinstance(value, Mapping) or set(value) != {item.name for item,_ in contract}:
                _fail()
            data = value
            if _single_walk_constructor(hint):
                token=_DECODE_CONSTRUCTION.set((hint,data))
                try:return hint(**data)
                finally:_DECODE_CONSTRUCTION.reset(token)
        else:
            if type(value) is not hint:
                _fail()
            data = {item.name: getattr(value, item.name) for item,_ in contract}
            if _single_walk_constructor(hint):
                # The ordinary constructor strictly copies every incoming field
                # itself; no pre-walk or unchecked instance allocation is needed.
                return hint(**data)
        arguments = {}
        for item,field_hint in contract:
            child = _typed_value(
                data[item.name],
                field_hint,
                decode=decode,
                _reuse_existing=_reuse_existing,
            )
            if item.init:
                arguments[item.name] = child
            elif child != item.default:
                _fail()
        return hint(**arguments)
    return _typed_scalar(value, hint)


class _Controlled:
    __slots__ = ()
    classification: ClassVar[ArtifactClassification] = ArtifactClassification.CONTROLLED
    llm_policy: ClassVar[LLMPolicy] = LLMPolicy.NO_LLM

    def __repr__(self):
        return "<controlled evidence; no-llm>"

    def __post_init__(self):
        construction=_DECODE_CONSTRUCTION.get()
        decode=construction is not None and construction[0] is type(self)
        contract=_field_contract(type(self))
        if decode:
            # Generated __init__ must have assigned the exact provided values.
            # Consume the handoff before any nested constructor can see it.
            _DECODE_CONSTRUCTION.set(None)
            if any(getattr(self,item.name) is not construction[1][item.name]
                    for item,_ in contract):_fail()
        with _validation_session():
            active, complete = _VALIDATION_WALK.get()
            key = (id(self), type(self), "construction")
            checkpoint = len(complete)
            try:
                if key in active:
                    _fail()
                active.add(key)
                for name, field_hint, scalar, nullable in _constructor_contract(type(self)):
                    value = getattr(self, name)
                    if nullable and value is None:
                        continue
                    copied = (_typed_scalar(value, scalar) if scalar is not None else
                              _typed_value(value, field_hint, decode=decode))
                    object.__setattr__(self, name, copied)
                # Do not register self here: subclass semantic checks have not
                # run yet. Only the enclosing _typed_value may publish it.
            except (ValueError, TypeError, RuntimeError, AttributeError, OverflowError, RecursionError):
                while len(complete) > checkpoint:
                    complete.popitem()
                raise EvidenceModelError("invalid controlled evidence structure") from None
            finally:
                active.discard(key)


def _sha(value):
    if type(value) is not str or re.fullmatch(r"[0-9a-f]{64}", value) is None:
        _fail()


def _ids(values):
    if any(not value or not value.strip() for value in values) or len(set(values)) != len(values):
        _fail()


def _exact(value: str) -> Fraction:
    if type(value) is not str or len(value) > 1024 or re.fullmatch(
        r"-?(?:0|[1-9][0-9]*)(?:\.[0-9]+|/[1-9][0-9]*)?", value,
    ) is None:
        _fail()
    return Fraction(value)


@dataclass(frozen=True, slots=True, repr=False)
class EvidenceReference(_Controlled):
    """Resolve a native record within an immutable public store, not a display join."""

    evidence_id: str
    record_id: str
    public_manifest_sha256: str

    def __post_init__(self):
        _Controlled.__post_init__(self)
        _ids((self.evidence_id,))
        _ids((self.record_id,))
        _sha(self.public_manifest_sha256)


@dataclass(frozen=True, slots=True, repr=False)
class SourcedRecord(_Controlled):
    """Uninterpreted annotation/input record, including all publisher-native fields."""

    source_dataset: str
    source_release: str
    source_record_id: str
    source_sha256: str
    public_manifest_sha256: str | None
    source_fields: Mapping[str, object] | None

    def __post_init__(self):
        _Controlled.__post_init__(self)
        for value in (self.source_dataset, self.source_release, self.source_record_id):
            _ids((value,))
        _sha(self.source_sha256)
        if self.public_manifest_sha256 is not None:
            _sha(self.public_manifest_sha256)


@dataclass(frozen=True, slots=True, repr=False)
class PublicSourcedEvidence(_Controlled):
    """Addressable raw public records, not authenticated merely by construction.

    Source adapters must verify actual artifact bytes and manifest pins before
    using these records. Keys are exact source_record_id values, not display IDs.
    """

    evidence_id: str
    public_manifest_sha256: str
    records: Mapping[str, SourcedRecord]

    def __post_init__(self):
        _Controlled.__post_init__(self)
        _ids((self.evidence_id,))
        _sha(self.public_manifest_sha256)
        if any(key != record.source_record_id or record.public_manifest_sha256 != self.public_manifest_sha256
               for key, record in self.records.items()):
            _fail()


@dataclass(frozen=True, slots=True, repr=False)
class ObservationMeasurement(_Controlled):
    name: str
    value: str
    unit: str | None
    provenance: ObservationProvenance

    def __post_init__(self):
        _Controlled.__post_init__(self)
        if any(
            not value or any(ord(character) < 32 for character in value)
            for value in (self.name, self.value)
        ) or self.unit is not None and (
            not self.unit or any(ord(character) < 32 for character in self.unit)
        ):
            _fail()


@dataclass(frozen=True, slots=True, repr=False)
class ObservationNote(_Controlled):
    text: str
    provenance: ObservationProvenance

    def __post_init__(self):
        _Controlled.__post_init__(self)
        if not self.text or any(ord(character) < 32 for character in self.text):
            _fail()


@dataclass(frozen=True, slots=True, repr=False)
class Observation(_Controlled):
    """Canonical observation whose scoring field is the only scoring projection."""

    scoring: ScoringObservation
    resolution: OntologyResolution
    ontology_manifest_sha256: str
    resolution_policy: str
    input_convention: Literal[
        "listed-positive-default-proband", "typed-structured"
    ]
    measurements: tuple[ObservationMeasurement, ...] = ()
    notes: tuple[ObservationNote, ...] = ()
    context_provenance: tuple[ObservationProvenance, ...] = ()
    source_fields: Mapping[str, object] | None = None

    def __post_init__(self):
        _Controlled.__post_init__(self)
        view = self.scoring
        successful = {"active", "alias", "replaced_by"}
        failed = {
            "unknown", "ambiguous_identifier", "ambiguous_replacement",
            "consider_only", "no_replacement", "missing_target", "cycle",
        }
        identifiers = (
            *self.resolution.chain,
            *self.resolution.replaced_by,
            *self.resolution.consider,
        )
        if (
            self.resolution.source_id != view.source_hpo_id
            or not self.resolution.chain
            or self.resolution.chain[0] != view.source_hpo_id
            or any(re.fullmatch(r"HP:\d{7}", value) is None for value in identifiers)
            or self.resolution.canonical_id is not None
            and re.fullmatch(r"HP:\d{7}", self.resolution.canonical_id) is None
            or (self.resolution.status in successful)
            != (self.resolution.canonical_id is not None)
            or self.resolution.status not in successful | failed
            or self.resolution.canonical_id is not None
            and self.resolution.chain[-1] != self.resolution.canonical_id
            or self.resolution_policy != OBSOLETE_RESOLUTION_POLICY
        ):
            _fail()
        _sha(self.ontology_manifest_sha256)
        if any(getattr(view, key) is not None for key in (
            "subject_id", "onset", "age_years", "sex", "domain", "dependency_cluster", "quality",
        )) and not self.context_provenance:
            _fail()


@dataclass(frozen=True, slots=True, repr=False)
class ObservationContradiction(_Controlled):
    """Caller-supplied unresolved provenance, including references across scopes."""

    contradiction_id: str
    observation_ids: tuple[str, ...]
    provenance: tuple[ObservationProvenance, ...]
    status: Literal["unresolved"] = "unresolved"

    def __post_init__(self):
        _Controlled.__post_init__(self)
        _ids((self.contradiction_id,))
        _ids(self.observation_ids)
        if len(self.observation_ids) < 2 or not self.provenance:
            _fail()


@dataclass(frozen=True, slots=True, repr=False)
class ScoringView(_Controlled):
    observations: tuple[ScoringObservation, ...]
    contradictions: tuple[ObservationContradiction, ...]
    context_provenance: Mapping[str, tuple[ObservationProvenance, ...]]

    def __post_init__(self):
        _Controlled.__post_init__(self)
        identifiers = tuple(item.observation_id for item in self.observations)
        _ids(identifiers)
        _ids(tuple(item.contradiction_id for item in self.contradictions))
        if set(self.context_provenance) != set(identifiers):
            _fail()


def scoring_view(
    observations: tuple[Observation, ...],
    contradictions: tuple[ObservationContradiction, ...] = (),
) -> ScoringView:
    """Closed projection only. No diagnosis extraction, reconciliation, or ablation."""
    try:
        records = _typed(observations, tuple[Observation, ...])
        _ids(tuple(item.scoring.observation_id for item in records))
        return ScoringView(
            tuple(item.scoring for item in records), contradictions,
            {item.scoring.observation_id: item.context_provenance for item in records},
        )
    except (ValueError, TypeError, RuntimeError, AttributeError):
        raise EvidenceModelError("invalid scoring projection") from None


def _public_binding(binding, manifest_sha256):
    _sha(manifest_sha256)
    if not binding or binding.get("classification", "public") != "public":
        _fail()
    if "bundle_manifest_sha256" in binding and binding["bundle_manifest_sha256"] != manifest_sha256:
        _fail()


@lru_cache(maxsize=128)
def _source_field_contract(schema_text):
    """Compile only static schema text; no case values or validation verdicts."""
    schema = dict(part.strip().replace('"', '').split() for part in schema_text.split(','))
    hints = {"VARCHAR": str, "VARCHAR[]": tuple[str | None, ...], "BOOLEAN": bool,
             "INTEGER": int, "BIGINT": int, "DOUBLE": float}
    return frozenset(schema), tuple((name, hints[dtype]) for name, dtype in schema.items())


def _source_fields(table, row):
    """Validate the existing v6 contract without rebuilding its interpretations."""
    if table == 'track2_concept':
        from genosapiens.track2_concept_supplement import validate_concept_record
        _json_value(row)
        validate_concept_record(row)
        return
    if table in NATIVE_LOCUS_TABLES:
        _json_value(row)
        if table == "native_repeat_locus" and not isinstance(row.get("LocusId"), str):
            _fail()
        return
    if table not in TABLE_SCHEMAS:
        _fail()
    names, contract = _source_field_contract(TABLE_SCHEMAS[table])
    if not names <= set(row):
        _fail()
    for name, hint in contract:
        value = row[name]
        if value is None:
            continue
        _typed(value, hint)


@dataclass(frozen=True, slots=True, repr=False)
class DiseaseProfile(_Controlled):
    """Composition of the complete task 7.2 store, not a lossy selected profile."""

    evidence_id: str
    public_manifest_sha256: str
    store: DiseaseProfileStore

    def __post_init__(self):
        _Controlled.__post_init__(self)
        _ids((self.evidence_id,))
        _public_binding(self.store.source_binding, self.public_manifest_sha256)
        _sha(self.store.graph_sha256)
        for annotation in self.store.annotations:
            _source_fields("phenotype_assertion", annotation.source_fields)
        identifiers = tuple(item.annotation_id for item in self.store.annotations)
        _ids(identifiers)
        identifier_set = set(identifiers)
        if not set(self.store.unassigned_annotation_ids) <= identifier_set:
            _fail()
        for key, profile in self.store.profiles.items():
            if key != profile.disease_id:
                _fail()
            _ids(tuple(group.group_id for group in profile.groups))
            if any(not set(group.annotation_ids) <= identifier_set for group in profile.groups):
                _fail()


@dataclass(frozen=True, slots=True, repr=False)
class DiseaseGeneEvidence(_Controlled):
    """Complete task 7.4 store and optional existing context/conflict ledgers."""

    evidence_id: str
    public_manifest_sha256: str
    store: DiseaseGeneEvidenceStore
    aggregations: tuple[DiseaseGeneAggregation, ...] = ()

    def __post_init__(self):
        _Controlled.__post_init__(self)
        _ids((self.evidence_id,))
        _public_binding(self.store.source_binding, self.public_manifest_sha256)
        _sha(self.store.evidence_sha256)
        records = self.store.records
        for record in records.values():
            _source_fields(record.table, record.source_fields)
        if any(key != record.record_id for key, record in records.items()):
            _fail()
        refs = [*self.store.unassigned_record_ids]
        for assertion in self.store.assertions:
            refs.append(assertion.record_id)
            if assertion.locus_id is not None:
                if (not assertion.locus_id.strip() or not assertion.native_disease_id
                        or not assertion.locus_record_id
                        or assertion.record_id not in records or assertion.locus_record_id not in records
                        or records[assertion.record_id].table != "native_repeat_disease"
                        or records[assertion.locus_record_id].table != "native_repeat_locus"
                        or records[assertion.locus_record_id].source_fields["LocusId"] != assertion.locus_id):
                    _fail()
                refs.append(assertion.locus_record_id)
            elif assertion.native_disease_id is not None or assertion.locus_record_id is not None:
                _fail()
            refs.extend(key for requirement in assertion.requirements for key in requirement.record_ids)
        refs.extend(key for ids in self.store.views.values() for key in ids)
        refs.extend(key for relation in self.store.relations for key in relation.record_ids)
        if not set(refs) <= set(records):
            _fail()
        for result in self.aggregations:
            if result.evidence_sha256 != self.store.evidence_sha256 or result.source_binding != self.store.source_binding:
                _fail()
            references = {item.record_id for item in result.ledger}
            references.update(key for item in result.ledger for key in (*item.relation_record_ids, *item.discounted_by_record_ids))
            references.update(key for item in result.penalties for key in (*item.positive_record_ids, *item.negative_record_ids))
            if not references <= set(records):
                _fail()


@dataclass(frozen=True, slots=True, repr=False)
class TranscriptConsequence(_Controlled):
    consequence_id: str
    source: SourcedRecord
    ensembl_gene_id: str | None = None
    ensembl_transcript_id: str | None = None
    ensembl_protein_id: str | None = None
    hgnc_symbol: str | None = None
    hgnc_id: str | None = None
    mane_refseq: str | None = None
    uniprot_accessions: tuple[str | None, ...] | None = None
    uniprot_isoforms: tuple[str | None, ...] | None = None
    hgvs_genomic: str | None = None
    hgvs_coding: str | None = None
    hgvs_protein: str | None = None
    exon: str | None = None
    intron: str | None = None
    protein_position: str | None = None
    protein_start: int | None = None
    protein_end: int | None = None
    amino_acids: str | None = None
    domain_overlaps: tuple[SourcedRecord, ...] | None = None
    transcript_status: tuple[str | None, ...] | None = None
    canonical: bool | None = None
    variant_class: str | None = None
    consequences: tuple[str | None, ...] | None = None
    biotype: str | None = None
    conservation: tuple[SourcedRecord, ...] | None = None
    constraint: tuple[SourcedRecord, ...] | None = None
    predictions: tuple[SourcedRecord, ...] | None = None

    def __post_init__(self):
        _Controlled.__post_init__(self)
        _ids((self.consequence_id,))
        for prefix, value in (("ENSG", self.ensembl_gene_id), ("ENST", self.ensembl_transcript_id),
                              ("ENSP", self.ensembl_protein_id)):
            if value is not None and re.fullmatch(prefix + r"[0-9]{11}\.[1-9][0-9]*", value) is None:
                _fail()
        if any(value is not None and value < 1 for value in (self.protein_start, self.protein_end)):
            _fail()
        if self.protein_start is not None and self.protein_end is not None and self.protein_end < self.protein_start:
            _fail()


@dataclass(frozen=True, slots=True, repr=False)
class MechanismAssertion(_Controlled):
    """Curated source assertion; direction/alternatives live in GeneRequirements.

    Publications, native/mapped disease IDs and full original fields resolve via
    source_refs. Predicted transcript severity is deliberately not an input.
    """

    mechanism_id: str
    assertion: GeneAssertion
    source_refs: tuple[EvidenceReference, ...]
    conflict_refs: tuple[EvidenceReference, ...] = ()

    def __post_init__(self):
        _Controlled.__post_init__(self)
        _ids((self.mechanism_id,))
        required = {self.assertion.record_id}
        required.update(key for item in self.assertion.requirements for key in item.record_ids)
        if not required <= {item.record_id for item in self.source_refs}:
            _fail()


class IdentityStage(StrEnum):
    SOURCE = "source"
    NATIVE = "native"
    NORMALIZED = "normalized"
    RENDERED = "rendered"


@dataclass(frozen=True, slots=True, repr=False)
class IdentityMapping(_Controlled):
    identity_id: str
    stage: IdentityStage
    identity: AlleleIdentity
    reference_manifest_sha256: str
    source: SourcedRecord
    derived_from_identity_ids: tuple[str, ...] = ()
    transformation_id: str | None = None

    def __post_init__(self):
        _Controlled.__post_init__(self)
        _ids((self.identity_id,))
        _sha(self.reference_manifest_sha256)
        _ids(self.derived_from_identity_ids)
        if self.identity_id in self.derived_from_identity_ids:
            _fail()


@dataclass(frozen=True, slots=True, repr=False)
class CallEvidence(_Controlled):
    source: SourcedRecord
    sample: SampleCallQc | None
    filters: tuple[str, ...] | None
    quality: float | None


@dataclass(frozen=True, slots=True, repr=False)
class PopulationEvidence(_Controlled):
    """One source stratum, not an aggregation or a hypothesis frequency decision."""

    source: SourcedRecord
    subset: str | None = None
    ancestry: str | None = None
    allele_count: int | None = None
    allele_number: int | None = None
    allele_frequency: float | None = None
    filtering_allele_frequency: float | None = None
    homozygote_count: int | None = None
    hemizygote_count: int | None = None
    callable: bool | None = None
    filters: tuple[str | None, ...] | None = None
    discordant: bool | None = None
    filtering_allele_frequency_99: float | None = None
    stratum: Literal["adjusted", "ancestry", "grpmax", "fafmax95", "fafmax99", "all_sites"] | None = None

    def __post_init__(self):
        _Controlled.__post_init__(self)
        if any(value is not None and value < 0 for value in (
            self.allele_count, self.allele_number, self.homozygote_count, self.hemizygote_count,
        )) or any(value is not None and not 0 <= value <= 1 for value in (
            self.allele_frequency, self.filtering_allele_frequency, self.filtering_allele_frequency_99,
        )):
            _fail()


@dataclass(frozen=True, slots=True, repr=False)
class ClinicalEvidence(_Controlled):
    """Aggregate records remain aggregates; condition links must be supplied."""

    source: SourcedRecord
    assertion_id: str | None = None
    condition_ids: tuple[str | None, ...] | None = None
    significance: tuple[str | None, ...] | None = None
    review_status: str | None = None
    explicit_condition_link: bool | None = None
    target_id: str | None = None
    confidence: str | None = None
    evaluated_on: str | None = None
    conflict_refs: tuple[EvidenceReference, ...] = ()

    def __post_init__(self):
        _Controlled.__post_init__(self)
        if self.explicit_condition_link is True and (
            not self.assertion_id or not self.assertion_id.strip()
            or not any(value is not None and value.strip() for value in self.condition_ids or ())
        ):
            _fail()


@dataclass(frozen=True, slots=True, repr=False)
class AlleleEvidence(_Controlled):
    allele_id: str
    native_identity_id: str
    identities: tuple[IdentityMapping, ...]
    calls: tuple[CallEvidence, ...] | None = None
    population: tuple[PopulationEvidence, ...] | None = None
    clinical: tuple[ClinicalEvidence, ...] | None = None
    transcript_consequences: tuple[TranscriptConsequence, ...] | None = None
    mechanism_ids: tuple[str, ...] = ()
    evidence_refs: tuple[EvidenceReference, ...] = ()

    def __post_init__(self):
        _Controlled.__post_init__(self)
        _ids((self.allele_id,))
        identifiers = tuple(item.identity_id for item in self.identities)
        _ids(identifiers)
        if not any(
            item.identity_id == self.native_identity_id and item.stage is IdentityStage.NATIVE
            for item in self.identities
        ):
            _fail()
        if any(not set(item.derived_from_identity_ids) <= set(identifiers) for item in self.identities):
            _fail()
        try:
            TopologicalSorter({item.identity_id: item.derived_from_identity_ids for item in self.identities}).prepare()
        except CycleError:
            raise EvidenceModelError("cyclic identity mapping") from None


class Completeness(StrEnum):
    COMPLETE = "complete"
    INCOMPLETE = "incomplete"
    UNKNOWN = "unknown"


class HypothesisPhase(StrEnum):
    UNKNOWN = "unknown"
    IN_TRANS = "in_trans"
    IN_CIS = "in_cis"
    CONFLICTING = "conflicting"
    NOT_APPLICABLE = "not_applicable"


class EnumerationStatus(StrEnum):
    COMPLETE = "complete"
    TRUNCATED = "truncated"
    FAILED = "failed"
    NOT_ASSESSED = "not_assessed"


@dataclass(frozen=True, slots=True, repr=False)
class LocusBinding(_Controlled):
    """Explicit native-locus field of the source record owning a requirement.

    The source adapter declares the dictionary-key path; the model does not guess
    publisher field names, reinterpret a gene ID as a locus, or acquire evidence.
    """

    source_ref: EvidenceReference
    field_path: tuple[str, ...]

    def __post_init__(self):
        _Controlled.__post_init__(self)
        if not self.field_path or any(not key.strip() for key in self.field_path):
            _fail()


@dataclass(frozen=True, slots=True, repr=False)
class HypothesisContext(_Controlled):
    context_id: str
    disease_id: str
    native_disease_id: str | None
    gene_id: str | None
    locus_id: str | None
    inheritance: str | None
    model: str
    requirements: tuple[GeneRequirements, ...]
    evidence_refs: tuple[EvidenceReference, ...]
    locus_bindings: tuple[LocusBinding, ...] = ()

    def __post_init__(self):
        _Controlled.__post_init__(self)
        for value in (self.context_id, self.disease_id, self.model):
            _ids((value,))
        if not self.gene_id and not self.locus_id or any(
            value is not None and not value.strip() for value in (self.gene_id, self.locus_id, self.native_disease_id)
        ):
            _fail()
        try:
            inheritance = DiseaseGeneContext(
                self.disease_id, self.gene_id, self.inheritance, self.native_disease_id, locus_id=self.locus_id,
            ).inheritance
        except DiseaseGeneEvidenceError:
            raise EvidenceModelError("invalid hypothesis inheritance context") from None
        if any(item.inheritance is not None and item.inheritance != inheritance for item in self.requirements):
            _fail()
        if any(item.required_allele_count is not None and item.required_allele_count <= 0 for item in self.requirements):
            _fail()
        if not {key for item in self.requirements for key in item.record_ids} <= {ref.record_id for ref in self.evidence_refs}:
            _fail()
        if self.locus_bindings and self.locus_id is None or any(
            binding.source_ref not in self.evidence_refs for binding in self.locus_bindings
        ):
            _fail()
        if self.locus_id is not None and self.requirements and not self.locus_bindings:
            _fail()


def _validate_context_sources(context, evidence, *, prepared_gene_evidence=None):
    """Non-scoring adapter of task 7.4's native anchor and relation semantics.

    Only exact source alternatives may supply requirements. Mapped-only records
    remain evidence, not established applicable requirements. Unlike aggregation,
    this does not select/discount sources or consult numerical support/penalties.
    Store-local applicability is unchanged; cross-store paths require every
    relation record to be referenced under its owning store's manifest.
    """
    if prepared_gene_evidence is not None:
        from genosapiens.genotype_evaluation import PreparedGeneEvidence
        if type(prepared_gene_evidence) is not PreparedGeneEvidence:
            _fail()
        prepared_gene_evidence.require(evidence)
    remaining = list(context.requirements)
    resolved_locus_bindings = set()
    if any(not item.record_ids for item in remaining):
        _fail()
    referenced_relations = []
    for public in evidence:
        refs = {ref.record_id for ref in context.evidence_refs
                if ref.evidence_id == public.evidence_id and ref.public_manifest_sha256 == public.public_manifest_sha256
                and ref.record_id in public.store.records}
        if prepared_gene_evidence is None:
            referenced_relations.extend(relation for relation in public.store.relations
                                        if relation.record_ids and set(relation.record_ids) <= refs)
        else:
            referenced_relations.extend(prepared_gene_evidence.stores[public.evidence_id].referenced_relations(refs))
    referenced_edges = {}
    for relation in referenced_relations:
        if not relation.subject_id or not relation.object_id:
            continue
        if relation.kind in ("equivalent", "narrower"):
            referenced_edges.setdefault(relation.subject_id, set()).add(relation.object_id)
        if relation.kind in ("equivalent", "broader"):
            referenced_edges.setdefault(relation.object_id, set()).add(relation.subject_id)
    for public in evidence:
        store = public.store
        prepared = None if prepared_gene_evidence is None else prepared_gene_evidence.stores[public.evidence_id]
        if prepared is not None:
            prepared.require(store)
        refs = {ref.record_id for ref in context.evidence_refs
                if ref.evidence_id == public.evidence_id and ref.public_manifest_sha256 == public.public_manifest_sha256}
        locus_records = set()
        for binding in context.locus_bindings:
            if binding.source_ref.evidence_id != public.evidence_id:
                continue
            value = store.records[binding.source_ref.record_id].source_fields
            for key in binding.field_path:
                if not isinstance(value, Mapping) or key not in value:
                    _fail()
                value = value[key]
            if type(value) is not str or value != context.locus_id:
                _fail()
            locus_records.add(binding.source_ref.record_id)
            resolved_locus_bindings.add(binding)
        target_assertions = tuple(
            assertion for assertion in (store.assertions if prepared is None else prepared.source_assertions(context))
            if (context.gene_id is None or assertion.gene_id == context.gene_id)
            and (assertion.locus_id is None or assertion.locus_id == context.locus_id)
            and (context.locus_id is None or (assertion.locus_record_id or assertion.record_id) in locus_records)
        )
        anchor = context.native_disease_id or context.disease_id
        # Same anchor witness as aggregate_disease_gene_evidence: a mapped view
        # must explicitly contain this native disease AND this target.
        if anchor != context.disease_id and not any(
            item.record_id in store.views.get(context.disease_id, ())
            and assertion_native_disease_id(store, item) == anchor
            for item in target_assertions
        ):
            continue
        edges = {} if prepared is None else prepared.source_edges
        for relation in store.relations if prepared is None else ():
            if not relation.subject_id or not relation.object_id:
                continue
            if relation.kind in ("equivalent", "narrower"):
                edges.setdefault(relation.subject_id, set()).add(relation.object_id)
            if relation.kind in ("equivalent", "broader"):
                edges.setdefault(relation.object_id, set()).add(relation.subject_id)
        reached, pending = {anchor}, [anchor]
        while pending:
            current = pending.pop()
            for target in (edges.get(current, frozenset()) | referenced_edges.get(current, frozenset())) - reached:
                reached.add(target)
                pending.append(target)
        for assertion in target_assertions:
            if (assertion.record_id not in refs
                    or assertion_native_disease_id(store, assertion) not in reached):
                continue
            remaining = [item for item in remaining
                         if item not in assertion.requirements or not set(item.record_ids) <= refs]
    if remaining or set(context.locus_bindings) != resolved_locus_bindings:
        _fail()


@dataclass(frozen=True, slots=True, repr=False)
class AlleleReference(_Controlled):
    allele_id: str
    identity_id: str


@dataclass(frozen=True, slots=True, repr=False)
class GenotypeHypothesis(_Controlled):
    hypothesis_id: str
    context: HypothesisContext
    allele_refs: tuple[AlleleReference, ...]
    completeness: Completeness
    phase: HypothesisPhase
    enumeration_status: EnumerationStatus
    evidence_refs: tuple[EvidenceReference, ...]
    observation_ids: tuple[str, ...] = ()
    mechanism_ids: tuple[str, ...] = ()
    completeness_evidence: tuple[SourcedRecord, ...] = ()
    phase_evidence: tuple[SourcedRecord, ...] = ()
    enumeration_evidence: tuple[SourcedRecord, ...] = ()
    # Explicit references to authoritative repeat calls, not individual alleles.
    # Both use the existing identity/evidence graph; serialization never projects
    # a locus call into allele_refs merely to satisfy an allele-only consumer.
    locus_call_refs: tuple[AlleleReference, ...] = ()

    @property
    def allele_ids(self) -> tuple[str, ...]:
        return tuple(item.allele_id for item in self.allele_refs)

    @property
    def evidence_subject_refs(self) -> tuple[AlleleReference, ...]:
        return self.allele_refs + self.locus_call_refs

    @property
    def evidence_subject_ids(self) -> tuple[str, ...]:
        return tuple(item.allele_id for item in self.evidence_subject_refs)

    def __post_init__(self):
        _Controlled.__post_init__(self)
        _ids((self.hypothesis_id,))
        _ids(self.evidence_subject_ids)
        if not self.evidence_subject_refs or self.allele_refs and self.locus_call_refs:
            _fail()
        if self.locus_call_refs and (self.context.model != "repeat" or len(self.locus_call_refs) != 1
                                     or self.phase is not HypothesisPhase.NOT_APPLICABLE):
            _fail()
        counts = {item.required_allele_count for item in self.context.requirements}
        if self.completeness is not Completeness.UNKNOWN and (
            not counts or None in counts or len(counts) != 1 or not self.completeness_evidence
        ):
            _fail()
        if self.phase in (HypothesisPhase.IN_TRANS, HypothesisPhase.IN_CIS) and not self.phase_evidence:
            _fail()
        # This checks claimed states only. Allele counts, GT, and phase are never inferred.
        if self.enumeration_status is not EnumerationStatus.NOT_ASSESSED and not self.enumeration_evidence:
            _fail()


ComponentName = Literal[
    "phenotype", "disease_gene", "molecular", "inheritance", "population", "clinical",
    "technical_and_conflict_penalties",
]


@dataclass(frozen=True, slots=True, repr=False)
class ScoreAdjustment(_Controlled):
    """Audit of an adjustment already included in the declared normalized value."""

    adjustment_id: str
    amount: str
    reason: str
    evidence_refs: tuple[EvidenceReference, ...]

    def __post_init__(self):
        _Controlled.__post_init__(self)
        _exact(self.amount)


@dataclass(frozen=True, slots=True, repr=False)
class ScoreComponent(_Controlled):
    name: ComponentName
    raw: str | None
    normalized: str | None
    neutral_value: str
    weight: str
    transform_id: str
    transform_sha256: str
    transform_parameters: Mapping[str, object]
    evidence_refs: tuple[EvidenceReference, ...] = ()
    source_values: tuple[SourcedRecord, ...] | None = None
    adjustments: tuple[ScoreAdjustment, ...] = ()

    def __post_init__(self):
        _Controlled.__post_init__(self)
        for value in (self.raw, self.normalized, self.neutral_value, self.weight):
            if value is not None:
                _exact(value)
        if self.raw is None and self.normalized is not None or _exact(self.weight) < 0:
            _fail()
        _ids((self.transform_id,))
        _sha(self.transform_sha256)


@dataclass(frozen=True, slots=True, repr=False)
class TieKey(_Controlled):
    name: str
    value: str | None
    direction: Literal["ascending", "descending"]


@dataclass(frozen=True, slots=True, repr=False)
class ScoreLedger(_Controlled):
    hypothesis_id: str
    components: tuple[ScoreComponent, ...]
    total: str
    tie_keys: tuple[TieKey, ...]
    configuration_id: str
    configuration_sha256: str
    configuration_provenance: Mapping[str, object]

    def recomputed_total(self) -> Fraction:
        """Check the approved additive equation, never execute a transformation."""
        return sum((
            (-1 if item.name == "technical_and_conflict_penalties" else 1)
            * _exact(item.weight)
            * _exact(item.normalized if item.normalized is not None else item.neutral_value)
            for item in self.components
        ), Fraction(0))

    def __post_init__(self):
        _Controlled.__post_init__(self)
        names = tuple(item.name for item in self.components)
        _ids(names)
        if set(names) != set(get_args(ComponentName)) or _exact(self.total) != self.recomputed_total():
            _fail()
        _ids((self.hypothesis_id,))
        _ids((self.configuration_id,))
        _sha(self.configuration_sha256)
        _ids(tuple(item.name for item in self.tie_keys))
        _ids(tuple(item.adjustment_id for component in self.components for item in component.adjustments))


@dataclass(frozen=True, slots=True, repr=False)
class EvidenceDocument(_Controlled):
    """Case-scoped evidence once, referenced by hypotheses; no readiness claim."""

    run_manifest_sha256: str
    observations: tuple[Observation, ...] = ()
    contradictions: tuple[ObservationContradiction, ...] = ()
    disease_profiles: tuple[DiseaseProfile, ...] = ()
    disease_gene_evidence: tuple[DiseaseGeneEvidence, ...] = ()
    mechanisms: tuple[MechanismAssertion, ...] = ()
    alleles: tuple[AlleleEvidence, ...] = ()
    hypotheses: tuple[GenotypeHypothesis, ...] = ()
    score_ledgers: tuple[ScoreLedger, ...] = ()
    public_sourced_evidence: tuple[PublicSourcedEvidence, ...] = ()

    def __post_init__(self):
        _Controlled.__post_init__(self)
        _sha(self.run_manifest_sha256)
        collections = (
            tuple(item.scoring.observation_id for item in self.observations),
            tuple(item.contradiction_id for item in self.contradictions),
            tuple(item.evidence_id for item in (*self.disease_profiles, *self.disease_gene_evidence,
                                               *self.public_sourced_evidence)),
            tuple(item.mechanism_id for item in self.mechanisms),
            tuple(item.allele_id for item in self.alleles),
            tuple(item.hypothesis_id for item in self.hypotheses),
            tuple(item.hypothesis_id for item in self.score_ledgers),
        )
        for identifiers in collections:
            _ids(identifiers)
        observation_ids, _, _, mechanism_ids, _, hypothesis_ids, ledger_ids = map(set, collections)
        if not ledger_ids <= hypothesis_ids:
            _fail()
        for contradiction in self.contradictions:
            if not set(contradiction.observation_ids) <= observation_ids:
                _fail()
        identities = {item.allele_id: item.native_identity_id for item in self.alleles}
        native_mappings = {a.allele_id: next(i for i in a.identities if i.identity_id == a.native_identity_id)
                           for a in self.alleles}
        contexts = {}
        for hypothesis in self.hypotheses:
            if (
                any(identities.get(ref.allele_id) != ref.identity_id for ref in hypothesis.evidence_subject_refs)
                or not set(hypothesis.observation_ids) <= observation_ids
                or not set(hypothesis.mechanism_ids) <= mechanism_ids
            ):
                _fail()
            for ref in hypothesis.evidence_subject_refs:
                native = native_mappings[ref.allele_id]
                # Structural validation only: the supplied adapter must replay the
                # original bytes before its interpretation can be consumed.
                interpretation = (native.source.source_fields or {}).get("native_interpretation") or {}
                if not isinstance(interpretation, Mapping):
                    _fail()
                locus_call = native.identity.kind.value == "repeat" and interpretation.get("repeat_representation") == "locus_call"
                if (ref in hypothesis.locus_call_refs) != locus_call:
                    _fail()
            context = hypothesis.context
            if context.context_id in contexts and context != contexts[context.context_id]:
                _fail()
            contexts[context.context_id] = context
        if any(not set(item.mechanism_ids) <= mechanism_ids for item in self.alleles):
            _fail()
        sources = {
            item.evidence_id: (item.public_manifest_sha256, {row.annotation_id for row in item.store.annotations})
            for item in self.disease_profiles
        }
        sources.update({
            item.evidence_id: (item.public_manifest_sha256, set(item.store.records))
            for item in self.disease_gene_evidence
        })
        sources.update({
            item.evidence_id: (item.public_manifest_sha256, set(item.records))
            for item in self.public_sourced_evidence
        })
        assertions = {
            (item.evidence_id, assertion.record_id): assertion
            for item in self.disease_gene_evidence for assertion in item.store.assertions
        }
        for mechanism in self.mechanisms:
            if not any(
                assertions.get((ref.evidence_id, ref.record_id)) == mechanism.assertion
                for ref in mechanism.source_refs
            ):
                _fail()

        def check_refs(value):
            if type(value) is EvidenceReference:
                pin, records = sources.get(value.evidence_id, (None, ()))
                if pin != value.public_manifest_sha256 or value.record_id not in records:
                    _fail()
            elif is_dataclass(value):
                for item in fields(value):
                    check_refs(getattr(value, item.name))
            elif isinstance(value, Mapping):
                for child in value.values():
                    check_refs(child)
            elif type(value) is tuple:
                for child in value:
                    check_refs(child)

        check_refs(self)
        for context in contexts.values():
            _validate_context_sources(context, self.disease_gene_evidence)

    @property
    def proband_vector(self) -> tuple[str, ...]:
        """Resolved positive proband HPO terms only; family context stays separate."""
        return tuple(
            item.resolution.canonical_id
            for item in self.observations
            if item.scoring.scope == "proband"
            and item.scoring.status == "present"
            and item.resolution.canonical_id is not None
        )

    def scoring_view(self) -> ScoringView:
        return scoring_view(self.observations, self.contradictions)


def evidence_document_from_phenotype_ingestion(
    result: PhenotypeIngestionResult, *, run_manifest_sha256: str,
) -> EvidenceDocument:
    """Package canonical list observations without inferring structured context."""
    try:
        from genosapiens.phenotype_ingestion import PhenotypeIngestionResult

        result = _typed(result, PhenotypeIngestionResult)
        return EvidenceDocument(
            run_manifest_sha256,
            observations=result.observations,
        )
    except (ValueError, TypeError, RuntimeError, AttributeError):
        raise EvidenceModelError("invalid phenotype ingestion evidence") from None


def resolve_evidence_reference(document: EvidenceDocument, ref: EvidenceReference) -> Mapping[str, object] | None:
    """Resolve exact evidence ID + manifest + record to raw fields, not a display join.

    Resolution is structural, not source authentication. Adapters must separately
    hash the original source bytes and verify their governing resource pins.
    """
    try:
        ref = _typed(ref, EvidenceReference)
        for public in (*document.disease_profiles, *document.disease_gene_evidence, *document.public_sourced_evidence):
            if public.evidence_id != ref.evidence_id:
                continue
            if public.public_manifest_sha256 != ref.public_manifest_sha256:
                _fail()
            if type(public) is DiseaseProfile:
                return next(row.source_fields for row in public.store.annotations if row.annotation_id == ref.record_id)
            records = public.records if type(public) is PublicSourcedEvidence else public.store.records
            return records[ref.record_id].source_fields
    except (KeyError, StopIteration, TypeError, AttributeError):
        pass
    raise EvidenceModelError("unresolved public evidence reference") from None


ENTITY_TYPES = (
    Observation, DiseaseProfile, DiseaseGeneEvidence, TranscriptConsequence,
    MechanismAssertion, AlleleEvidence, GenotypeHypothesis, ScoreLedger,
    EvidenceDocument, ScoringView, PublicSourcedEvidence,
)
