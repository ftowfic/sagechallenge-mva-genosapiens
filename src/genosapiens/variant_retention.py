"""Shared post-annotation retention over source-authenticated allele facts.

Adapters own normalization, QC, source authentication, and complete lookups.
This module performs no I/O and establishes neither pathogenicity nor mechanism.
Keep every source fact in the decision, including nonqualifying and conflicting facts.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import re
from typing import Literal


CONSEQUENCES = (
    'transcript_ablation', 'splice_acceptor_variant', 'splice_donor_variant',
    'stop_gained', 'frameshift_variant', 'stop_lost', 'start_lost',
    'transcript_amplification', 'inframe_insertion', 'inframe_deletion',
    'missense_variant', 'protein_altering_variant', 'splice_region_variant',
    'splice_donor_5th_base_variant', 'splice_donor_region_variant',
    'splice_polypyrimidine_tract_variant',
)
ROUTES = ('known_disease_allele', 'predicted_disease_gene_effect',
          'structural_variant', 'relevant_coding_splice_consequence')
POLICY_VERSION = 'shared-variant-retention/1.0.0'


def policy_document():
    """Return a new policy document without mutable shared configuration."""
    return {
        'schema': 'genosapiens.variant-retention-policy', 'version': POLICY_VERSION,
        'routes': ROUTES, 'combination': 'union_all_satisfied_routes',
        'consequences': CONSEQUENCES,
        'clinical_sources': ('ClinVar VCF', 'Open Targets evidence_eva'),
        'clinical_classifications': ('pathogenic', 'likely pathogenic',
                                     'pathogenic/likely pathogenic'),
        'clinical_review': 'preserve_all_review_levels_and_conflicts',
        'clinical_identity': 'exact_normalized_reference_contig_position_ref_alt',
        'clinical_condition': 'source_disease_required_gene_optional',
        'functional_sources': ('gencc_submissions', 'evidence_clingen'),
        'functional_assertion': 'exact_damaging_disease_relevant_nonnegated_with_citation_and_passage',
        'transcripts': 'all_applicable_source_records_with_explicit_relevance',
        'predictions': 'missense_sift_deleterious_or_absent_sift_polyphen_probably_damaging',
        'prediction_conflicts': 'preserve_without_inventing_functional_direction',
        'missing_evidence': 'deferred_no_retention_evidence_after_complete_required_lookup',
        'required_lookup_failure': 'stop_selection',
        'optional_unavailable': 'preserve_coverage_gap',
        'upstream_exclusion': 'never_override',
        'structural': 'scoped_native_event_without_sequence_proxy',
    }


def _canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode('ascii')


POLICY_SHA256 = hashlib.sha256(_canonical(policy_document())).hexdigest()


class RetentionError(ValueError):
    """Static errors never disclose case alleles or private paths."""


def _require(condition):
    if not condition:
        raise RetentionError('variant retention facts or dependencies are invalid')


def _text(value):
    _require(type(value) is str and bool(value) and value.strip() == value
             and not any(c in value for c in '\r\n\t'))


def _texts(value):
    _require(type(value) is tuple)
    for item in value:
        _text(item)


def _digest(value):
    _require(type(value) is str and re.fullmatch('[0-9a-f]{64}', value) is not None)


@dataclass(frozen=True, slots=True)
class ExactAllele:
    """An adapter-normalized sequence key on an authenticated reference."""
    reference_identity: str
    reference_sha256: str
    contig: str
    position: int
    reference: str
    alternate: str

    def __post_init__(self):
        _text(self.reference_identity)
        _digest(self.reference_sha256)
        _text(self.contig)
        _require(type(self.position) is int and self.position > 0)
        for sequence in (self.reference, self.alternate):
            _require(type(sequence) is str and re.fullmatch('[ACGT]+', sequence) is not None)
        _require(self.reference != self.alternate)


@dataclass(frozen=True, slots=True)
class SourceFact:
    dataset: str
    release_sha256: str
    record_id: str

    def __post_init__(self):
        _text(self.dataset)
        _digest(self.release_sha256)
        _text(self.record_id)


@dataclass(frozen=True, slots=True)
class LookupCoverage:
    dataset: str
    release_sha256: str
    required: bool
    status: Literal['complete', 'unavailable', 'failed', 'incomplete', 'corrupt']

    def __post_init__(self):
        _text(self.dataset)
        _digest(self.release_sha256)
        _require(type(self.required) is bool and self.status in
                 {'complete', 'unavailable', 'failed', 'incomplete', 'corrupt'})


@dataclass(frozen=True, slots=True)
class TranscriptFact:
    source: SourceFact
    gene_id: str | None
    transcript_id: str | None
    consequences: tuple[str, ...]
    relevance: Literal['relevant', 'irrelevant', 'unresolved']
    sift: str | None = None
    polyphen: str | None = None
    impact: str | None = None
    conflicts: tuple[str, ...] = ()

    def __post_init__(self):
        _require(type(self.source) is SourceFact)
        for value in (self.gene_id, self.transcript_id):
            if value is not None:
                _text(value)
        _texts(self.consequences)
        _texts(self.conflicts)
        _require(self.relevance in {'relevant', 'irrelevant', 'unresolved'})
        for value in (self.sift, self.polyphen, self.impact):
            _require(value is None or type(value) is str)


@dataclass(frozen=True, slots=True)
class ClinicalFact:
    source: SourceFact
    allele: ExactAllele | None
    mapping: Literal['exact', 'ambiguous', 'unmapped']
    classifications: tuple[str, ...]
    disease_ids: tuple[str, ...]
    gene_id: str | None = None
    review_status: str | None = None
    conflicts: tuple[str, ...] = ()

    def __post_init__(self):
        _require(type(self.source) is SourceFact and self.mapping in {'exact', 'ambiguous', 'unmapped'})
        _require(self.allele is None or type(self.allele) is ExactAllele)
        _require(self.mapping != 'exact' or self.allele is not None)
        for values in (self.classifications, self.disease_ids, self.conflicts):
            _texts(values)
        for value in (self.gene_id, self.review_status):
            if value is not None:
                _text(value)


@dataclass(frozen=True, slots=True)
class FunctionalFact:
    source: SourceFact
    allele: ExactAllele | None
    mapping: Literal['exact', 'ambiguous', 'unmapped']
    gene_id: str | None
    disease_ids: tuple[str, ...]
    effect: Literal['damaging', 'not_damaging', 'uncertain']
    negated: bool
    citation: str | None
    passage: str | None
    extraction_policy: str
    extraction_limits: tuple[str, ...] = ()
    conflicts: tuple[str, ...] = ()

    def __post_init__(self):
        _require(type(self.source) is SourceFact and self.mapping in {'exact', 'ambiguous', 'unmapped'})
        _require(self.allele is None or type(self.allele) is ExactAllele)
        _require(self.mapping != 'exact' or self.allele is not None)
        _require(self.effect in {'damaging', 'not_damaging', 'uncertain'} and type(self.negated) is bool)
        for value in (self.gene_id, self.citation):
            if value is not None:
                _text(value)
        _require(self.passage is None or type(self.passage) is str and bool(self.passage.strip()))
        _text(self.extraction_policy)
        for values in (self.disease_ids, self.extraction_limits, self.conflicts):
            _texts(values)


@dataclass(frozen=True, slots=True)
class StructuralFact:
    source: SourceFact
    native_identity: str
    in_scope: bool
    interpretation_limits: tuple[str, ...] = ()

    def __post_init__(self):
        _require(type(self.source) is SourceFact and type(self.in_scope) is bool)
        _text(self.native_identity)
        _texts(self.interpretation_limits)


@dataclass(frozen=True, slots=True)
class RetentionFacts:
    identity_id: str
    allele: ExactAllele | None
    upstream_disposition: Literal['assess', 'exclude_common', 'defer_refcall']
    occurrence_ids: tuple[str, ...]
    transcripts: tuple[TranscriptFact, ...]
    clinical: tuple[ClinicalFact, ...]
    functional: tuple[FunctionalFact, ...]
    structural: tuple[StructuralFact, ...]
    coverage: tuple[LookupCoverage, ...]

    def __post_init__(self):
        _text(self.identity_id)
        _require(self.allele is None or type(self.allele) is ExactAllele)
        _require(self.upstream_disposition in {'assess', 'exclude_common', 'defer_refcall'})
        _texts(self.occurrence_ids)
        _require(bool(self.occurrence_ids) and len(set(self.occurrence_ids)) == len(self.occurrence_ids))
        for values, kind in ((self.transcripts, TranscriptFact), (self.clinical, ClinicalFact),
                             (self.functional, FunctionalFact), (self.structural, StructuralFact),
                             (self.coverage, LookupCoverage)):
            _require(type(values) is tuple and all(type(value) is kind for value in values))
        coverage = {(r.dataset, r.release_sha256): r for r in self.coverage}
        _require(len(coverage) == len(self.coverage))
        for fact in (*self.transcripts, *self.clinical, *self.functional, *self.structural):
            key = fact.source.dataset, fact.source.release_sha256
            _require(key in coverage and coverage[key].status == 'complete')


@dataclass(frozen=True, slots=True)
class RetentionDecision:
    facts: RetentionFacts
    policy_version: str
    policy_sha256: str
    catalog_sha256: str
    disposition: Literal['retained', 'deferred_no_retention_evidence', 'not_applicable_upstream']
    routes: tuple[str, ...]
    evidence: tuple[tuple[str, SourceFact], ...]
    optional_coverage_gaps: tuple[LookupCoverage, ...]

    @property
    def retained(self):
        return self.disposition == 'retained'

    def as_bytes(self):
        return _canonical(asdict(self))


def prediction_label(value):
    """Reuse the existing SIFT/PolyPhen label rule without numeric cutoffs."""
    if value in (None, '', '.'):
        return None
    match = re.fullmatch(r'([a-z_]+)(?:\(([0-9]+(?:\.[0-9]+)?(?:e-[0-9]+)?)\))?', value)
    if not match or match[2] is not None and not 0 <= float(match[2]) <= 1:
        return 'malformed'
    return match[1]


def predicted_damage(transcript):
    """A damage prediction supplies retention evidence, never functional direction."""
    _require(type(transcript) is TranscriptFact)
    if 'missense_variant' not in transcript.consequences:
        return False
    sift, polyphen = prediction_label(transcript.sift), prediction_label(transcript.polyphen)
    return sift == 'deleterious' or sift is None and polyphen == 'probably_damaging'


def _clinical_pathogenic(fact):
    if fact.source.dataset == 'ClinVar VCF':
        # CLNSIG uses pipes for distinct aggregate labels. Do not split commas
        # into invented condition/assertion pairs or infer a positive conflict.
        values = (part.lower().replace('_', ' ') for value in fact.classifications
                  for part in value.split('|'))
    elif fact.source.dataset == 'Open Targets evidence_eva':
        values = (value.lower() for value in fact.classifications)
    else:
        return False
    return bool(set(values) & {'pathogenic', 'likely pathogenic', 'pathogenic/likely pathogenic'})


def select_variant(facts, *, established_genes, catalog_sha256):
    """Evaluate the four OR routes after the adapter closes required evidence.

    The catalog is public and diagnosis blind for Track 1. Track 2 restricts
    retrieval before calling this same predicate. Neither track passes ranking,
    phenotype, mechanism, or drug evidence to this function.
    """
    _require(type(facts) is RetentionFacts and type(established_genes) is frozenset)
    _digest(catalog_sha256)
    evidence = []
    gaps = tuple(row for row in facts.coverage if not row.required and row.status == 'unavailable')
    if facts.upstream_disposition != 'assess':
        disposition = 'not_applicable_upstream'
    else:
        if any(row.status in {'failed', 'incomplete', 'corrupt'}
               or row.required and row.status != 'complete' for row in facts.coverage):
            raise RetentionError('required variant selection evidence did not complete')
        if not facts.coverage:
            raise RetentionError('variant selection coverage is missing')
        for fact in facts.clinical:
            if (fact.mapping == 'exact' and fact.allele == facts.allele and facts.allele is not None
                    and fact.disease_ids and _clinical_pathogenic(fact)):
                evidence.append((ROUTES[0], fact.source))
        for fact in facts.functional:
            if (fact.source.dataset in {'gencc_submissions', 'evidence_clingen'}
                    and fact.mapping == 'exact' and fact.allele == facts.allele and facts.allele is not None
                    and fact.gene_id and fact.disease_ids and fact.effect == 'damaging'
                    and not fact.negated and fact.citation and fact.passage):
                evidence.append((ROUTES[0], fact.source))
        for fact in facts.transcripts:
            if (fact.gene_id not in established_genes or not fact.transcript_id
                    or fact.relevance != 'relevant'):
                continue
            if predicted_damage(fact):
                evidence.append((ROUTES[1], fact.source))
            if set(fact.consequences) & set(CONSEQUENCES):
                evidence.append((ROUTES[3], fact.source))
        for fact in facts.structural:
            if fact.in_scope and fact.native_identity == facts.identity_id:
                evidence.append((ROUTES[2], fact.source))
        disposition = 'retained' if evidence else 'deferred_no_retention_evidence'
    routes = tuple(route for route in ROUTES if any(item[0] == route for item in evidence))
    return RetentionDecision(facts, POLICY_VERSION, POLICY_SHA256, catalog_sha256,
                             disposition, routes, tuple(evidence), gaps)
