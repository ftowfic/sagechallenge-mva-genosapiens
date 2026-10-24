"""Pure diagnosis-blind selection over authenticated source facts.

This owner defines scope, never scientific support, scoring or eligibility.
Acquisition adapters authenticate source facts before invoking this evaluator.
Phase 1 implements focused; exhaustive has no executable gate until M1.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass, field
import hashlib
import json
from pathlib import Path
from typing import Literal

from genosapiens.annotation_admission import INTRAGENIC_TERMS, NON_OVERLAP_TERMS, _gene
from genosapiens.wgs_ledger import canonical
from genosapiens.wgs_shards import read_bounded_metadata
from genosapiens import variant_retention as retention

POLICY_PATH = Path(__file__).resolve().parents[2] / 'config/track1/analysis-profiles.json'
POLICY_SHA256 = 'ea12ce48321da6add3bc3d929447415be0701c892b035b25d34f44a225b4b291'
MOLECULAR_CONSEQUENCES = retention.CONSEQUENCES
ROUTES = ('molecular', 'clinical', 'native', 'regulatory')
_MOLECULAR_SET = frozenset(MOLECULAR_CONSEQUENCES)
_KNOWN_CONSEQUENCES = INTRAGENIC_TERMS | NON_OVERLAP_TERMS
_NEGATIVE_VALIDITIES = frozenset(('', 'not yet assessed', 'disputed', 'disputed evidence',
    'refuted', 'refuted evidence', 'no known disease relationship'))
Disposition = Literal['selected', 'unresolved', 'outside_profile', 'not_applicable_upstream']


class AnalysisProfileError(ValueError):
    """Static source/policy failure; never expose case values in diagnostics."""


def _require(condition):
    if not condition:
        raise AnalysisProfileError('analysis profile facts or policy are invalid')


def _digest(value):
    _require(type(value) is str and len(value) == 64
             and all(c in '0123456789abcdef' for c in value))
    return value


def _text(value):
    _require(type(value) is str and bool(value) and value.strip() == value
             and not any(c in value for c in '\r\n\t'))
    return value


def _stable_gene(value):
    try:
        gene = _gene(value)
        _require(gene is not None)
        return gene.split('.')[0]
    except (ValueError, TypeError):
        raise AnalysisProfileError('analysis profile target identity is invalid') from None


@dataclass(frozen=True, slots=True)
class AnalysisProfilePolicy:
    content: bytes
    _content_owner: bytes = field(init=False, repr=False, compare=False)
    _content_sha256: str = field(init=False, repr=False, compare=False)

    def __post_init__(self):
        _require(type(self.content) is bytes and hashlib.sha256(self.content).hexdigest() == POLICY_SHA256)
        object.__setattr__(self, '_content_owner', self.content)
        object.__setattr__(self, '_content_sha256', POLICY_SHA256)

    @property
    def sha256(self):
        _require(self.content == self._content_owner)
        return self._content_sha256

    @property
    def document(self):
        return json.loads(self.content)


def load_analysis_profile_policy(path=POLICY_PATH, *, expected_sha256=POLICY_SHA256):
    try:
        content = read_bounded_metadata(path)
        _require(expected_sha256 == POLICY_SHA256)
        policy = AnalysisProfilePolicy(content)
        _require(tuple(policy.document['molecular_consequences']) == MOLECULAR_CONSEQUENCES)
        return policy
    except (OSError, KeyError, TypeError, ValueError):
        raise AnalysisProfileError('analysis profile policy is missing, stale or invalid') from None


def resolve_analysis_profile(policy, *, exhaustive=False):
    _require(type(policy) is AnalysisProfilePolicy and type(exhaustive) is bool)
    if exhaustive:
        raise AnalysisProfileError('profile_not_qualified')
    return 'focused'


@dataclass(frozen=True, slots=True)
class AnalysisProfileSelection:
    """One authenticated consumer scope; contains no release permission flag."""
    profile: str
    policy: AnalysisProfilePolicy

    def __post_init__(self):
        _require(type(self.policy) is AnalysisProfilePolicy)
        _require(self.profile in {'focused', 'exhaustive'})
        resolved = resolve_analysis_profile(self.policy, exhaustive=self.profile == 'exhaustive')
        _require(self.profile == resolved)


def authenticate_analysis_profile(profile='focused', *, policy=None, expected_policy_sha256=POLICY_SHA256):
    # Refuse unavailable/unknown modes before reading even the bounded policy.
    if profile == 'exhaustive':
        raise AnalysisProfileError('profile_not_qualified')
    _require(profile == 'focused' and expected_policy_sha256 == POLICY_SHA256)
    policy = load_analysis_profile_policy() if policy is None else policy
    _require(type(policy) is AnalysisProfilePolicy and policy.sha256 == expected_policy_sha256)
    return AnalysisProfileSelection(profile, policy)


@dataclass(frozen=True, slots=True)
class ProfileCatalogRecord:
    gene_id: str | None
    original_gene_id: str | None
    source_dataset: str
    source_disease_id: str | None
    assertion_id: str | None
    original_row_json: str
    catalog_eligible: bool
    reason: str


@dataclass(frozen=True, slots=True)
class ProfileCatalog:
    policy_sha256: str
    source_sha256: str
    records: tuple[ProfileCatalogRecord, ...]
    genes: tuple[str, ...]
    links: tuple[tuple[str, str], ...]
    coverage_issues: tuple[tuple[str, int], ...]
    source_binding_json: str
    _content_sha256: str = field(init=False, repr=False, compare=False)
    _gene_membership: frozenset[str] = field(init=False, repr=False, compare=False)

    def __post_init__(self):
        _digest(self.policy_sha256)
        _digest(self.source_sha256)
        _require(type(self.records) is tuple and all(type(r) is ProfileCatalogRecord for r in self.records)
                 and self.genes == tuple(sorted(set(self.genes)))
                 and self.links == tuple(sorted(set(self.links))))
        object.__setattr__(self, '_content_sha256', hashlib.sha256(self.as_bytes()).hexdigest())
        object.__setattr__(self, '_gene_membership', frozenset(self.genes))

    def as_bytes(self):
        return canonical(dict(policy_sha256=self.policy_sha256, source_sha256=self.source_sha256,
            records=tuple(asdict(r) for r in self.records), genes=self.genes, links=self.links,
            coverage_issues=self.coverage_issues, source_binding_json=self.source_binding_json))

    @property
    def sha256(self):
        return self._content_sha256

    def contains(self, gene_id):
        return _stable_gene(gene_id) in self._gene_membership


def _catalog_predicate(row, catalog):
    if not row.get('ensembl_gene_id') or not row.get('source_disease_id'):
        return False, 'missing_gene_or_native_disease'
    if (row.get('disease_specificity') not in catalog['disease_specificity']
            or row.get('association_scope') not in catalog['association_scope']
            or row.get('is_ontology_propagated') is not False):
        return False, 'outside_direct_assertion_scope'
    rules = catalog['source_predicates'].get(row.get('source_dataset'))
    if rules is None:
        return False, 'unknown_catalog_source'
    validity = row.get('validity') or row.get('confidence') or ''
    _require(type(validity) is str)
    validity = validity.lower()
    if validity not in rules['validity']:
        return False, ('excluded_source_validity' if validity in _NEGATIVE_VALIDITIES
                       else 'unknown_catalog_vocabulary')
    if 'association_type_prefix' in rules:
        association = row.get('association_type') or ''
        _require(type(association) is str)
        if not association.startswith(rules['association_type_prefix']):
            return False, 'outside_causal_source_role'
    if ('consequence_requirement' in rules
            and row.get('consequence_requirement') not in rules['consequence_requirement']):
        return False, 'outside_source_consequence_role'
    return True, 'direct_catalog_assertion'


def build_profile_catalog(rows, *, policy, source_sha256, source_binding=None):
    """Deterministic source-preserving catalog; I/O/authentication is caller-owned.

    Preserve all supplied records, including negative conflicts and model
    alternatives. No assertion is rewritten to increase scientific validity.
    """
    _require(type(policy) is AnalysisProfilePolicy)
    _digest(source_sha256)
    if source_binding is None:
        source_binding = dict(schema='genosapiens.profile-catalog-synthetic-source',
                              source_sha256=source_sha256)
    source_binding_json = canonical(source_binding).decode('ascii')
    records, genes, links, coverage = [], set(), set(), Counter()
    catalog_rules = policy.document['catalog']
    for row in rows:
        _require(type(row) is dict)
        original_gene = row.get('ensembl_gene_id')
        gene = None if not original_gene else _stable_gene(original_gene)
        source = _text(row.get('source_dataset'))
        disease = row.get('source_disease_id')
        if disease:
            _text(disease)
        eligible, reason = _catalog_predicate(row, catalog_rules)
        if reason.startswith('unknown_catalog_'):
            coverage[f'{source}:{reason}'] += 1
        record = ProfileCatalogRecord(gene, original_gene, source, disease,
            row.get('source_record_id'), canonical(row).decode('ascii'), eligible, reason)
        records.append(record)
        if eligible:
            genes.add(gene)
            links.add((gene, disease))
    return ProfileCatalog(policy.sha256, source_sha256,
        tuple(sorted(records, key=lambda r: r.original_row_json)), tuple(sorted(genes)),
        tuple(sorted(links)), tuple(sorted(coverage.items())), source_binding_json)


def build_profile_catalog_from_gene_evidence(store, *, policy):
    """Reuse the existing authenticated disease-gene resource owner, unchanged.

    The public context owns acquisition/authentication through the pinned v6
    reader. This adapter retains its complete bundle/contract/source identity
    and source-native rows, without replaying or changing scientific validity.
    """
    from genosapiens.disease_gene_evidence import DiseaseGeneEvidenceStore, _plain
    _require(type(store) is DiseaseGeneEvidenceStore)
    binding = _plain(store.source_binding)
    _require(binding.get('classification') == 'public'
             and binding.get('schema') == 'genosapiens.rare-disease-evidence-bundle')
    source = binding['outputs']['gene_assertion']
    _digest(source['sha256'])
    rows = tuple(_plain(record.source_fields) for record in store.records.values()
                 if record.table == 'gene_assertion')
    _require(len(rows) == source['row_count'])
    return build_profile_catalog(rows, policy=policy, source_sha256=source['sha256'],
        source_binding=dict(bundle=binding, scientific_source_evidence_sha256=store.evidence_sha256))


@dataclass(frozen=True, slots=True)
class ProfileRouteEvidence:
    target_id: str | None
    source_ref: str
    source_sha256: str
    transcript_id: str | None = None
    consequences: tuple[str, ...] = ()
    interpretation: Literal['supported', 'unresolved'] = 'supported'
    # Actual assertion state is preserved, never used as a clinical hit cutoff.
    assertion_state: str | None = None
    target_namespace: str | None = 'Ensembl'
    source_target_id: str | None = None
    source_target_kind: Literal['ensembl', 'missing', 'unmapped'] | None = None
    source_dataset: str | None = None
    exact_allele: retention.ExactAllele | None = None
    clinical_classifications: tuple[str, ...] = ()
    clinical_conditions: tuple[str, ...] = ()
    clinical_conflicts: tuple[str, ...] = ()
    review_status: str | None = None
    sift: str | None = None
    polyphen: str | None = None
    structural_native_identity: str | None = None
    functional: retention.FunctionalFact | None = None

    def __post_init__(self):
        _require((self.source_target_id is None) == (self.source_target_kind is None))
        if self.source_target_id is not None:
            from genosapiens.admission_links import classify_encode_target
            _text(self.source_target_id)
            _require(classify_encode_target(self.source_target_id) == self.source_target_kind)
            if self.source_target_kind == 'ensembl':
                _require(self.target_id is not None and self.target_namespace == 'Ensembl'
                         and _stable_gene(self.target_id) == _stable_gene(self.source_target_id))
            else:
                _require(self.target_id is None and self.target_namespace is None
                         and self.interpretation == 'unresolved')
        if self.target_id is not None:
            _text(self.target_namespace)
            if self.target_namespace == 'Ensembl':
                _stable_gene(self.target_id)
            else:
                _text(self.target_id)
        _text(self.source_ref)
        _digest(self.source_sha256)
        if self.transcript_id is not None:
            _text(self.transcript_id)
        _require(type(self.consequences) is tuple
                 and all(type(v) is str and bool(v) for v in self.consequences)
                 and self.interpretation in {'supported', 'unresolved'})
        _require(self.exact_allele is None or type(self.exact_allele) is retention.ExactAllele)
        _require(self.functional is None or type(self.functional) is retention.FunctionalFact)
        for values in (self.clinical_classifications, self.clinical_conditions, self.clinical_conflicts):
            _require(type(values) is tuple and all(type(v) is str and v for v in values))


@dataclass(frozen=True, slots=True)
class ProfileRouteFacts:
    route: Literal['molecular', 'clinical', 'native', 'regulatory']
    status: Literal['matched', 'no_match', 'unresolved', 'not_applicable']
    source_sha256: str
    evidence: tuple[ProfileRouteEvidence, ...] = ()
    reasons: tuple[str, ...] = ()
    source_refs: tuple[str, ...] = ()

    def __post_init__(self):
        _require(self.route in ROUTES and self.status in
                 {'matched', 'no_match', 'unresolved', 'not_applicable'})
        _digest(self.source_sha256)
        _require(type(self.evidence) is tuple and all(type(e) is ProfileRouteEvidence
                 and e.source_sha256 == self.source_sha256 for e in self.evidence))
        _require(self.status not in {'no_match', 'not_applicable'} or not self.evidence)
        _require(self.status != 'matched' or bool(self.evidence))
        _require(all(type(values) is tuple and all(type(value) is str and bool(value)
                     for value in values) for values in (self.reasons, self.source_refs)))


@dataclass(frozen=True, slots=True)
class ProfileCandidateFacts:
    identity_id: str
    original_input_sha256: str
    upstream_disposition: Literal['assess', 'exclude_common', 'defer_refcall']
    routes: tuple[ProfileRouteFacts, ...]
    technical_reasons: tuple[str, ...] = ()
    exact_allele: retention.ExactAllele | None = None
    source_occurrences: tuple[str, ...] = ()
    public_coverage: tuple[retention.LookupCoverage, ...] = ()
    sample_ids: tuple[str, ...] = ()

    def __post_init__(self):
        _text(self.identity_id)
        _digest(self.original_input_sha256)
        _require(self.upstream_disposition in {'assess', 'exclude_common', 'defer_refcall'}
                 and type(self.routes) is tuple and all(type(r) is ProfileRouteFacts for r in self.routes)
                 and tuple(r.route for r in self.routes) == ROUTES
                 and type(self.technical_reasons) is tuple
                 and all(type(v) is str and bool(v) for v in self.technical_reasons))
        _require(self.exact_allele is None or type(self.exact_allele) is retention.ExactAllele)
        _require(type(self.source_occurrences) is tuple
                 and all(type(v) is str and v for v in self.source_occurrences))
        _require(type(self.public_coverage) is tuple
                 and all(type(v) is retention.LookupCoverage for v in self.public_coverage))
        _require(type(self.sample_ids) is tuple and all(type(v) is str and v for v in self.sample_ids)
                 and len(set(self.sample_ids)) == len(self.sample_ids))


@dataclass(frozen=True, slots=True)
class ProfileRouteDecision:
    route: str
    target_id: str | None
    outcome: str
    reasons: tuple[str, ...]
    source_sha256: str
    source_refs: tuple[str, ...]
    transcript_id: str | None = None
    target_namespace: str | None = None


@dataclass(frozen=True, slots=True)
class ProfileDecision:
    identity_id: str
    original_input_sha256: str
    profile: str
    policy_sha256: str
    catalog_sha256: str
    upstream_disposition: str
    disposition: Disposition
    routes: tuple[ProfileRouteDecision, ...]
    technical_reasons: tuple[str, ...]
    retention_disposition: str
    retention_routes: tuple[str, ...]
    retention_policy_sha256: str

    @property
    def advancing(self):
        return self.retention_disposition == 'retained'

    @property
    def target_ids(self):
        return tuple(sorted({r.target_id for r in self.routes if r.target_id is not None
                             and r.outcome in {'selected', 'unresolved'}}))

    @property
    def target_keys(self):
        return tuple(sorted({(r.target_namespace, r.target_id) for r in self.routes
            if r.target_id is not None and r.outcome in {'selected', 'unresolved'}}))


def _profile_evidence_sources(facts):
    """Own source identifiers for one synchronous policy evaluation only."""
    def source_for(route, item):
        dataset = {'molecular':'Ensembl VEP', 'clinical':'clinical-query',
                   'regulatory':'regulatory-query', 'native':'supplied-native'}[route.route]
        identity = item.source_ref + '#' + hashlib.sha256(canonical(asdict(item))).hexdigest()
        return retention.SourceFact(item.source_dataset or dataset, route.source_sha256, identity)
    return {(id(route), id(item)): source_for(route, item)
            for route in facts.routes for item in route.evidence}


def shared_retention_facts(facts):
    """Expose the same immutable selection facts to both track adapters."""
    _require(type(facts) is ProfileCandidateFacts)
    return _shared_retention_facts(facts, _profile_evidence_sources(facts))


def _shared_retention_facts(facts, sources):
    transcripts, clinical, functional, structural, coverage = [], [], [], [], {}
    for route in facts.routes:
        dataset = {'molecular':'Ensembl VEP', 'clinical':'clinical-query',
                   'regulatory':'regulatory-query', 'native':'supplied-native'}[route.route]
        coverage[dataset, route.source_sha256] = retention.LookupCoverage(
            dataset, route.source_sha256, True, 'complete')
        for item in route.evidence:
            source = sources[id(route), id(item)]
            coverage[source.dataset, source.release_sha256] = retention.LookupCoverage(
                source.dataset, source.release_sha256, True, 'complete')
            target = (None if item.target_id is None or item.target_namespace != 'Ensembl'
                      else _stable_gene(item.target_id))
            if route.route == 'molecular':
                transcripts.append(retention.TranscriptFact(source, target, item.transcript_id,
                    item.consequences, 'relevant' if item.transcript_id and target else 'unresolved',
                    sift=item.sift, polyphen=item.polyphen))
            elif route.route == 'clinical':
                clinical.append(retention.ClinicalFact(source, item.exact_allele,
                    'exact' if item.exact_allele is not None else 'unmapped',
                    item.clinical_classifications, item.clinical_conditions, target,
                    item.review_status, item.clinical_conflicts))
            elif route.route == 'native' and item.structural_native_identity is not None:
                structural.append(retention.StructuralFact(source, item.structural_native_identity, True,
                    ('native_interpretation_unresolved',) if item.interpretation == 'unresolved' else ()))
            if item.functional is not None:
                functional.append(item.functional)
                original = item.functional.source
                coverage[original.dataset, original.release_sha256] = retention.LookupCoverage(
                    original.dataset, original.release_sha256, True, 'complete')
    for source in facts.public_coverage:
        coverage[source.dataset, source.release_sha256] = source
    return retention.RetentionFacts(facts.identity_id, facts.exact_allele, facts.upstream_disposition,
        facts.source_occurrences, tuple(transcripts), tuple(clinical), tuple(functional), tuple(structural),
        tuple(coverage.values()))


def evaluate_profile_candidate(facts, *, policy, catalog, profile='focused'):
    _require(type(facts) is ProfileCandidateFacts and type(policy) is AnalysisProfilePolicy
             and type(catalog) is ProfileCatalog and catalog.policy_sha256 == policy.sha256)
    _require(profile in {'focused', 'exhaustive'})
    if profile == 'exhaustive':
        raise AnalysisProfileError('profile_not_qualified')
    if catalog.coverage_issues:
        raise AnalysisProfileError('profile_catalog_coverage_unresolved')
    sources = _profile_evidence_sources(facts)
    shared = _shared_retention_facts(facts, sources)
    selected = retention.select_variant(shared, established_genes=catalog._gene_membership,
        catalog_sha256=catalog.sha256)
    selected_refs = {source for _, source in selected.evidence}
    decisions = []
    for route in facts.routes:
        if facts.upstream_disposition != 'assess':
            decisions.append(ProfileRouteDecision(route.route, None, 'not_applicable_upstream',
                (facts.upstream_disposition,), route.source_sha256, ()))
            continue
        if route.status in {'no_match', 'not_applicable'}:
            decisions.append(ProfileRouteDecision(route.route, None, 'outside_profile',
                ('authenticated_no_match' if route.status == 'no_match' else 'not_applicable', *route.reasons),
                route.source_sha256, route.source_refs))
        elif route.status == 'unresolved':
            decisions.append(ProfileRouteDecision(route.route, None, 'outside_profile',
                ('deferred_no_retention_evidence',), route.source_sha256,
                (*route.source_refs, *(e.source_ref for e in route.evidence))))
        for evidence in route.evidence:
            if route.route in {'molecular', 'regulatory'}:
                _require(evidence.target_id is None or evidence.target_namespace == 'Ensembl')
            target = (None if evidence.target_id is None else
                      _stable_gene(evidence.target_id) if evidence.target_namespace == 'Ensembl'
                      else evidence.target_id)
            qualifies = sources[id(route), id(evidence)] in selected_refs or (evidence.functional is not None
                and evidence.functional.source in selected_refs)
            # Retention is an allele decision. Existing supported associations
            # remain usable after another route retains the allele. They never
            # add a retention route or restore an otherwise deferred allele.
            supported_association = selected.retained and (
                route.route == 'regulatory' and evidence.interpretation == 'supported'
                and target in catalog._gene_membership
                or route.route == 'clinical' and evidence.exact_allele == facts.exact_allele
                and evidence.exact_allele is not None
                or route.route == 'native')
            outcome = 'selected' if qualifies or supported_association else 'outside_profile'
            reasons = (selected.routes if qualifies else ('supported_association_of_retained_allele',)
                       if supported_association else ('deferred_no_retention_evidence',))
            decisions.append(ProfileRouteDecision(route.route, target, outcome, reasons,
                route.source_sha256, (evidence.source_ref,), evidence.transcript_id,
                None if target is None else evidence.target_namespace))
    disposition = ('not_applicable_upstream' if facts.upstream_disposition != 'assess' else
                   'selected' if selected.retained else 'outside_profile')
    return ProfileDecision(facts.identity_id, facts.original_input_sha256, profile, policy.sha256,
        catalog.sha256, facts.upstream_disposition, disposition, tuple(decisions), facts.technical_reasons,
        selected.disposition, selected.routes, selected.policy_sha256)
