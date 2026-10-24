"""Compact molecular exceptions derived from authenticated profile facts.

These flags select work. They do not assert a diagnosis or replace scientific
evaluation. The caller must close the complete profile stream before selection.
"""
from dataclasses import asdict, dataclass
import hashlib

from genosapiens import variant_retention as retention
from genosapiens.evidence_serialization import _owned_canonical_bytes
from genosapiens.track1_analysis_profiles import ProfileCandidateFacts, shared_retention_facts


@dataclass(frozen=True, slots=True)
class ProtectionFact:
    kind: str
    disease_ids: tuple[str, ...]
    gene_id: str | None
    source_dataset: str
    source_sha256: str
    source_record_id: str
    uncertain: bool = False


@dataclass(frozen=True, slots=True)
class MolecularProtection:
    identity_id: str
    original_input_sha256: str
    facts_sha256: str
    retention_policy_sha256: str
    facts: tuple[ProtectionFact, ...]
    uncertain: bool


def extract_protection(facts):
    """Project source facts without queries, genotype hypotheses, or scoring."""
    if type(facts) is not ProfileCandidateFacts:
        raise ValueError('scientific protection requires typed profile facts')
    shared = shared_retention_facts(facts)
    # The existing kernel checks source coverage and supplies qualifying source
    # identities. Empty established genes prevent gene-only protection.
    selected = retention.select_variant(shared, established_genes=frozenset(),
        catalog_sha256=retention.POLICY_SHA256)
    known = {source for route, source in selected.evidence if route == 'known_disease_allele'}
    native = {source for route, source in selected.evidence if route == 'structural_variant'}
    projected = []
    uncertain = bool(selected.optional_coverage_gaps)

    def add(kind, value, conditions, gene, ambiguity):
        source = value.source
        projected.append(ProtectionFact(kind, tuple(sorted(set(conditions))), gene,
            source.dataset, source.release_sha256, source.record_id, ambiguity))

    if facts.upstream_disposition == 'assess':
        for value in shared.clinical:
            # A source can carry several assertions. Recheck the assertion's
            # exact mapping and classification rather than treating its presence
            # in the source set as authority for every assertion in that source.
            applies = (value.source in known and value.mapping == 'exact'
                and shared.allele is not None and value.allele == shared.allele
                and bool(value.disease_ids) and retention._clinical_pathogenic(value))
            ambiguity = bool(value.conflicts or value.mapping != 'exact')
            if applies:
                add('clinical', value, value.disease_ids, value.gene_id, ambiguity)
            uncertain |= ambiguity
        for value in shared.functional:
            applies = (value.source in known and value.source.dataset in
                {'gencc_submissions', 'evidence_clingen'} and value.mapping == 'exact'
                and shared.allele is not None and value.allele == shared.allele
                and value.gene_id and value.disease_ids and value.effect == 'damaging'
                and not value.negated and value.citation and value.passage)
            ambiguity = bool(value.conflicts or value.extraction_limits
                or value.mapping != 'exact' or value.effect == 'uncertain')
            if applies:
                add('functional', value, value.disease_ids, value.gene_id, ambiguity)
            uncertain |= ambiguity
        for value in shared.structural:
            if value.source in native and value.in_scope and value.native_identity == shared.identity_id:
                add('native', value, (), None, bool(value.interpretation_limits))
        # Supported native routes include classes that have no sequence proxy.
        # The roster must additionally require an applicable disease link.
        for route in facts.routes:
            if route.route == 'native' and route.status == 'matched':
                for value in route.evidence:
                    projected.append(ProtectionFact('native', (), value.target_id,
                        value.source_dataset or 'supplied-native', value.source_sha256,
                        value.source_ref, value.interpretation != 'supported'))
            if route.route in {'clinical', 'native'} and route.status == 'unresolved':
                uncertain = True
    else:
        uncertain = False
    unique = tuple(sorted(set(projected), key=lambda item: (
        item.kind, item.disease_ids, item.gene_id or '', item.source_dataset,
        item.source_sha256, item.source_record_id, item.uncertain)))
    return MolecularProtection(facts.identity_id, facts.original_input_sha256,
        hashlib.sha256(_owned_canonical_bytes(facts)).hexdigest(), retention.POLICY_SHA256,
        unique, uncertain)


def iter_validated_profile_facts(reader, *, metrics=None):
    """Yield decoded facts only within a fully authenticated policy replay."""
    from genosapiens.wgs_profile_projection import ProfileProjectionReader, _Closure, _iter_sealed_profile_rows
    from genosapiens.track1_analysis_profiles import evaluate_profile_candidate
    from genosapiens.wgs_ledger import canonical
    if type(reader) is not ProfileProjectionReader:
        raise ValueError('scientific protection requires an authenticated profile reader')
    reader.require_stable()
    closure = _Closure(None, reader.identity)
    for row, facts in _iter_sealed_profile_rows(reader.root, reader.manifest, reader.identity,
                                              reader.signatures, _with_facts=True):
        decision = evaluate_profile_candidate(facts, policy=reader.selection.policy,
            catalog=reader.catalog, profile=reader.selection.profile)
        if canonical(asdict(decision)) != canonical(row.decision):
            raise ValueError('scientific protection profile decision differs')
        closure.add(row)
        if metrics is not None:
            metrics['profile_rows'] = metrics.get('profile_rows', 0) + 1
        if row.decision['retention_policy_sha256'] != retention.POLICY_SHA256:
            raise ValueError('scientific protection retention policy differs')
        if (facts.identity_id, facts.original_input_sha256) != (row.original_ref[2], row.original_ref[4]):
            raise ValueError('scientific protection input identity differs')
        yield row, facts
    if closure.result() != reader.manifest['counts']:
        raise ValueError('scientific protection profile closure differs')
    reader.require_stable()


def iter_profile_protection(reader, *, include=None, metrics=None):
    """Reuse decoded facts while authenticating and replaying every source row."""
    from contextlib import closing
    with closing(iter_validated_profile_facts(reader, metrics=metrics)) as rows:
        for row, facts in rows:
            if include is None or include(row.original_ref):
                yield row.original_ref, extract_protection(facts)
