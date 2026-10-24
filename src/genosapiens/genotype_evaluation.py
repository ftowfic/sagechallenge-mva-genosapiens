"""Private, operation-scoped full-context inputs for genotype kernels.

Dependency direction is deliberately downward: evidence/property types only,
never model, qualification, requirement or hypothesis kernels. Consumers own
scientific decisions and put their invariant results in this exact-owner cache.
Nothing here is a published evidence type or a cross-context trust flag.
"""
from __future__ import annotations

from collections import ChainMap, defaultdict
from collections.abc import Mapping
from dataclasses import fields, is_dataclass
import hashlib
import json
from types import MappingProxyType

from genosapiens.disease_gene_evidence import PreparedDiseaseGeneStore, _requirements_disagree
from genosapiens.evidence_model import (
    AlleleEvidence, DiseaseGeneEvidence, EvidenceReference, HypothesisContext,
    _Controlled, _ids, _sha, _typed, _validate_context_sources, _validation_session,
)
from genosapiens.evidence_serialization import _canonical, _owned_json_default
from genosapiens.genotype_properties import (
    POLICY_SHA256, SuppliedGenotypeProperties, selected_source_calls,
)

_SEAL = object()
_GENE_SEAL = object()
_POLICY_SEAL = object()


def _hash(value):
    return hashlib.sha256((_canonical(value) + '\n').encode('ascii')).hexdigest()


def _validate_policy(value):
    if type(value) is tuple:
        for item in value:
            _validate_policy(item)
    elif is_dataclass(value):
        _typed(value, type(value))
    else:
        raise ValueError('invalid prepared genotype policy')


class PreparedGenotypePolicies:
    """Private immutable parent and exact child policies for one case."""
    __slots__ = ('_seal', '_ready', 'parent', 'policies', 'hashes', 'sizes', 'loci')

    def __init__(self, *, _seal=None):
        if _seal is not _POLICY_SEAL:
            raise ValueError('genotype policies require preparation')
        self._seal = _seal

    def __setattr__(self, name, value):
        if getattr(self, '_ready', False):
            raise AttributeError('prepared genotype policies are immutable')
        object.__setattr__(self, name, value)

    def require(self, policies):
        if (self._seal is not _POLICY_SEAL or set(policies) != set(self.policies)
                or any(policies[name] is not value for name, value in self.policies.items())):
            raise ValueError('prepared genotype policy owner changed')


def prepare_genotype_policies(parent, names):
    """Copy and validate a controlled parent; never trust a caller's mutable tree."""
    if (not isinstance(parent, _Controlled) or not is_dataclass(parent)
            or not type(parent).__dataclass_params__.frozen
            or type(names) is not tuple or not names
            or any(type(name) is not str for name in names) or len(set(names)) != len(names)
            or any(name not in {field.name for field in fields(parent)} for name in names)):
        raise ValueError('invalid prepared genotype policy parent')
    owner = PreparedGenotypePolicies(_seal=_POLICY_SEAL)
    with _validation_session(isolated=True):
        owner.parent = _typed(parent, type(parent))
    policies = {name:getattr(owner.parent, name) for name in names}
    for value in policies.values():
        _validate_policy(value)
    owner.policies = MappingProxyType(policies)
    owner.sizes = MappingProxyType({name:len(_canonical(value)) for name, value in policies.items()})
    hashes = {name:_hash(value) for name, value in policies.items()}
    owner.loci = tuple(sorted(policies.get('loci', ()), key=_canonical))
    if 'loci' in policies:
        hashes['loci'] = _hash(owner.loci)
    owner.hashes = MappingProxyType(hashes)
    owner._ready = True
    return owner


class PreparedGeneEvidence:
    """One immutable case resource owner; no serialized authentication flag."""
    __slots__ = ("_seal","_ready","evidence","sha256","canonical_size","registry",
                 "gene_stores","assertions_by_record","stores","relation_edges")

    def __init__(self, *, _seal=None):
        if _seal is not _GENE_SEAL:
            raise ValueError('gene evidence requires preparation')
        self._seal = _seal

    def __setattr__(self, name, value):
        if getattr(self,"_ready",False):
            raise AttributeError("prepared case gene evidence is immutable")
        object.__setattr__(self,name,value)

    def require(self, evidence):
        if self._seal is not _GENE_SEAL or evidence is not self.evidence:
            raise ValueError('prepared case gene evidence owner changed')

    def __repr__(self):
        return '<PreparedGeneEvidence>'


def prepare_gene_evidence(evidence):
    owner = PreparedGeneEvidence(_seal=_GENE_SEAL)
    # Retain the canonical immutable copy, not the unvalidated caller's object.
    owner.evidence = _typed(evidence, tuple[DiseaseGeneEvidence, ...])
    _ids(tuple(e.evidence_id for e in owner.evidence))
    # Hash the same canonical logical document without a second aggregate JSON
    # tree/string/byte buffer. Publication and replay own bounded record streams.
    digest = hashlib.sha256()
    size = 0
    encoder = json.JSONEncoder(sort_keys=True, separators=(',', ':'), ensure_ascii=True,
                               allow_nan=False, default=_owned_json_default)
    for part in encoder.iterencode(owner.evidence):
        content = part.encode('ascii')
        digest.update(content)
        size += len(content)
    digest.update(b'\n')
    owner.sha256 = digest.hexdigest()
    owner.canonical_size = size + 1
    owner.registry = MappingProxyType({EvidenceReference(e.evidence_id,key,e.public_manifest_sha256):record
        for e in owner.evidence for key,record in e.store.records.items()})
    owner.gene_stores = frozenset(e.evidence_id for e in owner.evidence)
    assertions = {}
    for public in owner.evidence:
        for assertion in public.store.assertions:
            assertions.setdefault((public.evidence_id,assertion.record_id),[]).append(assertion)
    owner.assertions_by_record = MappingProxyType({key:tuple(rows) for key,rows in assertions.items()})
    owner.stores = MappingProxyType({public.evidence_id:PreparedDiseaseGeneStore(public.store)
                                    for public in owner.evidence})
    # The source adapter defines direction and proof ordering. Build this
    # resource index once, alongside the other immutable case indexes.
    from genosapiens.candidate_admission import _relation_edges
    owner.relation_edges = MappingProxyType({key:tuple(rows)
        for key,rows in _relation_edges(owner.evidence).items()})
    owner._ready = True
    return owner


class PreparedGenotypeEvaluation:
    """Strong exact owners, full hashes and bounded context-local indexes.

    Created only by prepare_genotype_evaluation. The owning complete-context
    operation drops this value and all caches together. Source/model variants
    are accepted only when their full source context matches the prepared base.
    """
    def __init__(self, *, _seal=None):
        if _seal is not _SEAL:
            raise ValueError('prepared genotype evaluation requires preparation')
        self._seal = _seal

    def require(self, *, context, alleles=(), genotypes=None, gene_evidence=None, policies=None):
        if self._seal is not _SEAL:
            raise ValueError('invalid prepared genotype owner')
        if genotypes is not None and genotypes is not self.genotypes:
            raise ValueError('prepared genotype source owner changed')
        if gene_evidence is not None and gene_evidence is not self.gene_evidence:
            raise ValueError('prepared gene evidence owner changed')
        for name, policy in (policies or {}).items():
            if self.policies.get(name) is not policy:
                raise ValueError('prepared genotype policy owner changed')
        for allele in alleles:
            if self.alleles.get(allele.allele_id) is not allele:
                raise ValueError('prepared allele owner changed')
        self.context_key(context)
        return self

    def context_key(self, context):
        found = self.contexts.get(id(context))
        if found is not None:
            if found[0] is not context:
                raise ValueError('prepared context identity changed')
            return found[1]
        _typed(context, HypothesisContext)
        if any(getattr(context, field.name) != getattr(self.context, field.name)
               for field in fields(HypothesisContext) if field.name not in ('context_id', 'model')):
            raise ValueError('prepared context source scope changed')
        content = _canonical(context)
        digest = hashlib.sha256((content + '\n').encode('ascii')).hexdigest()
        self.contexts[id(context)] = (context, digest, content)
        self.context_source_membership[digest] = frozenset(context.evidence_refs)
        return digest

    def context_allele_hash(self, context, allele):
        self.require(context=context, alleles=(allele,))
        content = ('{"allele":' + self.allele_bytes[allele.allele_id] + ',"context":'
                   + self.contexts[id(context)][2] + '}\n')
        return hashlib.sha256(content.encode('ascii')).hexdigest()

    def alleles_hash(self, alleles):
        content = '[' + ','.join(self.allele_bytes[a.allele_id] for a in alleles) + ']\n'
        return hashlib.sha256(content.encode('ascii')).hexdigest()

    def genotype_binding_fields(self):
        return {'genotype_properties_sha256': self.genotypes_sha256}

    def validate_sources(self, context):
        key = self.context_key(context)
        # context_key already proves equality of every field except the model
        # and context ID. Neither excluded field is read by this source check.
        if not self.validated_sources:
            if any(ref.evidence_id in self.gene_stores and ref not in self.registry
                   for ref in context.evidence_refs):
                raise ValueError('prepared source reference is missing')
            _validate_context_sources(context, self.gene_evidence,
                prepared_gene_evidence=getattr(self, 'gene_evidence_owner', None))
            if any(_requirements_disagree(a, b) or (a.association_requirement is not None
                    and b.association_requirement is not None
                    and a.association_requirement != b.association_requirement)
                   for a in context.requirements for b in context.requirements):
                raise ValueError('prepared context requirements disagree')
            self.validated_sources.add(key)

    def calls_for(self, allele):
        key = allele.allele_id
        if key not in self.calls:
            self.calls[key] = selected_source_calls(allele, self.genotypes, _sources=self.source_properties)
        return self.calls[key]

    def clinical_reach(self, context):
        self.context_key(context)
        if self._clinical_reached is None:
            from genosapiens.candidate_admission import _relation_edges, _reach
            resource = getattr(self, 'gene_evidence_owner', None)
            if resource is not None:
                resource.require(self.gene_evidence)
            edges = _relation_edges(self.gene_evidence) if resource is None else resource.relation_edges
            # _reach expects a default-empty adjacency map. Missing anchors
            # belong to this call, never to the immutable shared resource.
            self._clinical_reached = MappingProxyType(_reach(
                context.native_disease_id or context.disease_id, ChainMap(edges, defaultdict(tuple))))
        return self._clinical_reached


def prepare_genotype_evaluation(*, context, alleles, genotypes, gene_evidence, policies, pairs,
                               prepared_gene_evidence=None, prepared_policies=None):
    """Validate the complete supplied objects once; retain their exact owners.

    Ordinary requirement-only consumers supply genotypes=None. Full genotype
    consumers always bind the original full properties, never a pair subset.
    Policies are a role->exact object mapping supplied by their typed kernel.
    """
    owner = PreparedGenotypeEvaluation(_seal=_SEAL)
    if prepared_policies is not None:
        if type(prepared_policies) is not PreparedGenotypePolicies:
            raise ValueError('invalid prepared genotype policies')
        prepared_policies.require(policies)
    values = tuple(alleles.values()) if isinstance(alleles, Mapping) else alleles
    with _validation_session(isolated=True):
        _typed(context, HypothesisContext)
        _typed(values, tuple[AlleleEvidence, ...])
        _typed(genotypes, SuppliedGenotypeProperties | None)
        if prepared_gene_evidence is None:
            _typed(gene_evidence, tuple[DiseaseGeneEvidence, ...])
        else:
            if type(prepared_gene_evidence) is not PreparedGeneEvidence:
                raise ValueError('invalid prepared case gene evidence')
            prepared_gene_evidence.require(gene_evidence)
        if prepared_policies is None:
            for policy in policies.values():
                _validate_policy(policy)
        _typed(pairs, tuple[tuple[str, str], ...])
    _ids(tuple(a.allele_id for a in values))
    _ids(tuple(e.evidence_id for e in gene_evidence))
    if len(set(pairs)) != len(pairs) or any(a >= b for a, b in pairs):
        raise ValueError('invalid prepared pair universe')
    owner.context, owner.genotypes, owner.gene_evidence = context, genotypes, gene_evidence
    owner.alleles = MappingProxyType({a.allele_id: a for a in values})
    owner.policies = MappingProxyType(dict(policies))
    owner.pairs = pairs
    owner.allele_bytes = {a.allele_id: _canonical(a) for a in values}
    owner.genotypes_sha256 = _hash(genotypes)
    owner.gene_evidence_sha256 = _hash(gene_evidence) if prepared_gene_evidence is None else prepared_gene_evidence.sha256
    if prepared_policies is None:
        owner.policy_hashes = {name: _hash(value) for name, value in policies.items()}
        owner.loci = tuple(sorted(policies.get('loci', ()), key=_canonical))
        if 'loci' in policies:
            owner.policy_hashes['loci'] = _hash(owner.loci)
    else:
        owner.policy_hashes = dict(prepared_policies.hashes)
        owner.loci = prepared_policies.loci
    owner.properties = {} if genotypes is None else {p.allele_id: p for p in genotypes.allele_properties}
    owner.source_properties = {} if genotypes is None else {p.source_allele_id: p for p in genotypes.source_properties}
    owner.phases = {} if genotypes is None else {tuple(sorted(p.allele_ids)): p for p in genotypes.pair_phases}
    if genotypes is not None:
        _ids(tuple(p.allele_id for p in genotypes.allele_properties))
        _ids(tuple(p.source_allele_id for p in genotypes.source_properties))
        if (genotypes.policy_sha256 != POLICY_SHA256
                or len(owner.properties) != len(genotypes.allele_properties)
                or len(owner.source_properties) != len(genotypes.source_properties)
                or len(owner.phases) != len(genotypes.pair_phases)):
            raise ValueError('invalid prepared genotype coverage')
        for value in (genotypes.inventory_sha256, genotypes.input_sha256, genotypes.reference_manifest_sha256):
            _sha(value)
    if prepared_gene_evidence is None:
        owner.registry = {EvidenceReference(e.evidence_id, key, e.public_manifest_sha256): record
                          for e in gene_evidence for key, record in e.store.records.items()}
        owner.gene_stores = frozenset(e.evidence_id for e in gene_evidence)
        assertions = {}
        for public in gene_evidence:
            for assertion in public.store.assertions:
                assertions.setdefault((public.evidence_id, assertion.record_id), []).append(assertion)
        owner.assertions_by_record = MappingProxyType({key: tuple(rows) for key, rows in assertions.items()})
    else:
        owner.registry = prepared_gene_evidence.registry
        owner.gene_stores = prepared_gene_evidence.gene_stores
        owner.assertions_by_record = prepared_gene_evidence.assertions_by_record
        owner.gene_evidence_owner = prepared_gene_evidence
    owner.contexts, owner.context_source_membership, owner.calls = {}, {}, {}
    owner.validated_sources = set()
    owner.requirements, owner.requirement_hashes, owner.clinical, owner.model_members = {}, {}, {}, {}
    owner._clinical_reached = None
    owner.context_key(context)
    return owner
