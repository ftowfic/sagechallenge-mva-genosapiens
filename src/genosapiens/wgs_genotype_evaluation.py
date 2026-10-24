"""Bounded exact-owner windows supplying the existing scientific kernels.

The genotype source is explicitly a D28 fact root, not a pretend full
SuppliedGenotypeProperties object or hash. The scientific owners still perform
all model, requirement and clinical decisions. An allele window is authenticated
against the sealed enriched context and discarded, with its caches, on exit.
"""
from contextlib import contextmanager
from types import MappingProxyType

from genosapiens.evidence_model import AlleleEvidence, _typed
from genosapiens.evidence_serialization import _canonical
from genosapiens.genotype_evaluation import (
    PreparedGeneEvidence, PreparedGenotypeEvaluation, PreparedGenotypePolicies, prepare_genotype_evaluation, _SEAL,
)
from genosapiens.wgs_context_evidence import EnrichedContextEvidence, VerifiedContextEvidence
from genosapiens.wgs_evidence_store import EvidenceStoreError
from genosapiens.wgs_genotype_facts import ContextGenotypeFacts
from genosapiens.wgs_validation_scope import member_validation_scope

MAX_WINDOW_BYTES = 128 * 1024**2
MAX_INVARIANT_BYTES = 128 * 1024**2


class _FactLookup:
    def __init__(self, lookup):
        self.lookup = lookup

    def __getitem__(self, key):
        return self.lookup(key)


class PreparedBoundedGenotypeEvaluation(PreparedGenotypeEvaluation):
    def genotype_binding_fields(self):
        self.fact_owner._require_open()
        return {"genotype_evidence":{
            "schema":"genosapiens.wgs-genotype-source-member-binding","version":1,
            "facts":self.fact_owner.published,
            "scope":"complete-source-and-member-facts;individual-pair-phase-bound-separately"}}

    @contextmanager
    def window(self, keys, *, phase=None):
        with member_validation_scope(), self._window(keys, phase=phase) as owner:
            yield owner

    @contextmanager
    def _window(self, keys, *, phase=None):
        if (self._window_active or type(keys) is not tuple or not 1 <= len(keys) <= 2
                or keys != tuple(sorted(set(keys)))):
            raise EvidenceStoreError("bounded genotype member window is invalid")
        self.fact_owner._require_open()
        self.member_owner._require_open()
        self._window_active = True
        try:
            values, content, size = {},{},0
            for key in keys:
                allele = self.member_owner.member(key).allele
                _typed(allele,AlleleEvidence)
                encoded = _canonical(allele)
                size += len(encoded)
                if size > MAX_WINDOW_BYTES:
                    raise EvidenceStoreError("bounded genotype window exceeds capacity",capacity=True)
                values[key],content[key] = allele,encoded
            self.alleles, self.allele_bytes = MappingProxyType(values),content
            if phase is not None:
                if len(keys) != 2 or phase != self.fact_owner.phase(*keys):
                    raise EvidenceStoreError("bounded genotype phase disagrees with source facts")
                self.phases[keys] = phase
            yield self
        finally:
            self.alleles,self.allele_bytes = MappingProxyType({}),{}
            for cache in (self.calls,self.phases,self.requirements,self.requirement_hashes,self.clinical,self.model_members):
                cache.clear()
            self._window_active = False

    def require(self, **kwargs):
        if not self._window_active:
            raise EvidenceStoreError("bounded genotype window is not active")
        self.fact_owner._require_open()
        self.member_owner._require_open()
        return super().require(**kwargs)


def prepare_bounded_genotype_evaluation(*, members, facts, context, gene_evidence, policies,
                                       prepared_gene_evidence=None, prepared_policies=None):
    if (type(members) not in (EnrichedContextEvidence,VerifiedContextEvidence) or type(facts) is not ContextGenotypeFacts
            or members.writer is not facts.writer or members.context_id != facts.context_id
            or facts.published is None or members.seed.context != context):
        raise EvidenceStoreError("bounded genotype owners disagree")
    members._require_open()
    facts._require_open()
    # A prepared case resource lends its immutable indexes and logical hash.
    # Its full serialized size is not a context-local allocation. Inline gene
    # evidence still belongs to this operation's invariant budget below.
    size = 0
    if prepared_gene_evidence is not None:
        if type(prepared_gene_evidence) is not PreparedGeneEvidence:
            raise ValueError('invalid prepared case gene evidence')
        prepared_gene_evidence.require(gene_evidence)
    if prepared_policies is not None:
        if type(prepared_policies) is not PreparedGenotypePolicies:
            raise ValueError('invalid prepared genotype policies')
        prepared_policies.require(policies)
        size += sum(prepared_policies.sizes.values())
    for value in (context,*(() if prepared_gene_evidence is not None else (gene_evidence,)),
                  *(() if prepared_policies is not None else policies.values())):
        size += len(_canonical(value))
        if size > MAX_INVARIANT_BYTES:
            raise EvidenceStoreError("bounded genotype invariant capacity exceeded",capacity=True)
    # Reuse the ordinary immutable context/source/policy validation and indexes.
    # No genotype subset is ever presented as a full supplied genotype owner.
    base = prepare_genotype_evaluation(context=context,alleles=(),genotypes=None,
        gene_evidence=gene_evidence,policies=policies,pairs=(),prepared_gene_evidence=prepared_gene_evidence,
        prepared_policies=prepared_policies)
    owner = PreparedBoundedGenotypeEvaluation(_seal=_SEAL)
    owner.__dict__.update(base.__dict__)
    owner.member_owner,owner.fact_owner = members,facts
    owner.genotypes = facts.header
    del owner.genotypes_sha256
    owner.properties,owner.source_properties = _FactLookup(facts.allele),_FactLookup(facts.source)
    owner._window_active = False
    return owner
