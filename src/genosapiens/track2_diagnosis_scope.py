"""Frozen, local diagnosis scope; interval inclusion never establishes causality."""
from collections import defaultdict
from contextlib import ExitStack
from dataclasses import asdict
import json
import math
from pathlib import Path
import re
import sqlite3
import time

from genosapiens import admission_links as links
from genosapiens.disease_gene_evidence import (
    DiseaseGeneEvidenceStore, _plain, assertion_native_disease_id,
)
from genosapiens.track1_analysis_profiles import build_profile_catalog, load_analysis_profile_policy
from genosapiens.track2_published_genotype import canonical, sha256


SCHEMA = 'genosapiens.track2-diagnosis-scope/v1'
POLICY = {
    'version': 'diagnosis-scope/v2', 'diagnosis': 'pinned_exact_identity_with_unresolved_legacy_conflicts',
    'subtypes': 'explicit_active_orpha_and_pinned_mondo_subclasses',
    'genes': 'existing_direct_catalog_predicate_and_positive_source_assertion',
    'regulatory': 'all_explicit_stable_target_links_without_score_or_tissue_filter',
    'coordinates': 'zero_based_half_open', 'hpo': 'context_only',
    'replacement': 'unique_orpha_21471_moved_to',
}
GENE = re.compile(r'ENSG[0-9]{11}(?:\.[0-9]+)?')
ACTIVE = {'Active', 'Active: Historic'}


class DiagnosisScopeError(ValueError):
    """Static messages; case diagnosis, selected genes and paths are not logged."""


def require(ok, message):
    if not ok:
        raise DiagnosisScopeError(message)


def resolve_diagnosis(store, diagnosis):
    """Resolve a native concept, not an assertion's mapped display view."""
    require(type(store) is DiseaseGeneEvidenceStore, 'disease evidence owner is required')
    require(isinstance(diagnosis, str) and re.fullmatch(r'(?:ORPHA|OMIM|MONDO):[0-9]+', diagnosis),
            'unsupported diagnosis identifier')
    catalog = defaultdict(list)
    for record in store.records.values():
        if record.table in {'disease', 'track2_concept'}:
            catalog[record.source_fields['disease_id']].append(record)
    exact = defaultdict(set)
    for relation in store.relations:
        if relation.kind == 'equivalent':
            exact[relation.subject_id].add(relation.object_id)
            exact[relation.object_id].add(relation.subject_id)

    def closure(start):
        found, pending = set(), list(start)
        while pending:
            item = pending.pop()
            if item not in found:
                found.add(item)
                pending.extend(exact[item] - found)
        return found

    equivalent = closure([diagnosis])
    replacements = set()
    # The source adapter deliberately does not turn replacements into scientific
    # hierarchy. Here only the explicit native identity replacement is followed.
    visited = set()
    while True:
        active = {key for key in equivalent if len(catalog[key]) == 1
                  and catalog[key][0].source_fields['status'] in ACTIVE}
        require(not any(len(catalog[key]) > 1 for key in equivalent), 'ambiguous diagnosis metadata')
        if active:
            orpha = {key for key in active if key.startswith('ORPHA:')}
            mondo = {key for key in active if key.startswith('MONDO:')}
            require(len(orpha) <= 1 and len(mondo) <= 1 and active == orpha | mondo,
                    'ambiguous diagnosis mapping')
            resolved = next(iter(orpha or mondo))
            break
        require(not equivalent & visited, 'cyclic diagnosis replacement')
        visited.update(equivalent)
        moves = [r for r in store.records.values() if r.table == 'disease_relation'
                 and r.source_fields['root_disease_id'] in equivalent
                 and r.source_fields['relation_type_id'] == '21471'
                 and r.source_fields['relation_type'] == 'Moved to'
                 and r.source_fields['unresolved_endpoint'] is False]
        destinations = {r.source_fields['target_disease_id'] for r in moves}
        require(len(destinations) == 1 and any(catalog[k] for k in equivalent),
                'diagnosis is unresolved or has no unique supported replacement')
        replacements.update(r.record_id for r in moves)
        equivalent = closure(destinations)
    require(catalog[resolved][0].source_fields['classification_level'] in {'Disorder', 'Subtype of disorder'},
            'broad diagnosis cannot define a molecular scope')
    concepts = set(equivalent)
    hierarchy_records = set()
    while True:
        children = [r for r in store.relations if r.kind == 'narrower'
                    and 'active_subtype_aggregation' in r.diagnostics and r.object_id in concepts]
        expanded = closure({r.subject_id for r in children}) | concepts
        hierarchy_records.update(i for r in children for i in r.record_ids)
        if expanded == concepts:
            break
        concepts = expanded
    evidence_ids = replacements | hierarchy_records | {
        r.record_id for key in concepts | visited for r in catalog[key]}
    relevant = [r for r in store.relations if r.subject_id in concepts | visited
                or r.object_id in concepts | visited]
    evidence_ids.update(i for r in relevant for i in r.record_ids)
    evidence = {'concept_ids': sorted(concepts), 'resolution': 'unique_active_concept',
        'records': [_plain(store.records[key]) for key in sorted(evidence_ids)],
        'relations': [_plain(asdict(r)) for r in relevant]}
    return resolved, concepts, evidence


def project_regulatory_targets(index, genes, *, maximum_intervals=100000, deadline_seconds=120):
    """One compact selected-gene projection, never read/decode provenance chunks.

    The existing index is keyed by element. Its membership scan occurs once per
    newly frozen scope; downstream selection reuses the resulting tiny intervals.
    Overflow/interruption raises, never publishes a truncated successful scope.
    """
    genes = tuple(sorted(set(genes)))
    require(genes and len(genes) <= 1000 and all(GENE.fullmatch(g) for g in genes),
            'regulatory projection target scope is invalid')
    require(type(maximum_intervals) is int and 0 < maximum_intervals <= 1000000
            and type(deadline_seconds) in (int, float) and math.isfinite(deadline_seconds)
            and 0 < deadline_seconds <= 3600, 'regulatory projection budget is invalid')
    index.require_stable()
    deadline = time.monotonic() + deadline_seconds
    db = index.connection
    db.set_progress_handler(lambda: int(time.monotonic() > deadline), 10000)
    result = []
    try:
        cursor = db.execute('''SELECT m.gene_id,e.chromosome,e.start,e.end,e.ordinal,
            e.current_id,e.legacy_id,m.source_rows FROM memberships m
            JOIN elements e ON e.legacy_id=m.legacy_id WHERE m.target_kind='ensembl'
            AND (''' + ' OR '.join('(m.gene_id=? OR m.gene_id GLOB ?)' for _ in genes) + ')',
            [value for g in genes for value in (g, g+'.*')])
        for gene, chrom, start, end, ordinal, current, element, count in cursor:
            require(len(result) < maximum_intervals and time.monotonic() <= deadline,
                    'regulatory scope projection incomplete: budget exceeded')
            result.append({'gene_id':gene.split('.')[0], 'source_target_id':gene,
                'chromosome':chrom, 'start':start, 'end':end,
                'coordinate_system':'zero_based_half_open', 'assembly':'GRCh38', 'kind':'regulatory',
                'element_ordinal':ordinal, 'element_id':element, 'anchor_id':current,
                'source_row_count':count, 'source_manifest_sha256':index.manifest_sha256})
    except sqlite3.Error:
        raise DiagnosisScopeError('regulatory scope projection incomplete') from None
    finally:
        db.set_progress_handler(lambda: (db._read_cache.release(), 0)[1], 100000)
        index.require_stable()
    return sorted(result, key=canonical)


def gene_intervals(pin, genes, reference):
    """Read selected stable IDs from the already published Open Targets target set."""
    import duckdb
    with ExitStack() as stack:
        files = links._Files(stack, reuse_published_outputs=True)
        release, contract, manifest, artifacts, _ = links._open_targets_metadata(files, pin)
        targets = [a for a in artifacts if a.dataset == 'target']
        require(bool(targets), 'required gene-coordinate resource is absent')
        outputs = {o.path:o for o in manifest.outputs}
        paths = []
        for artifact in targets:
            path = release / artifact.relative_path
            output = outputs[f'{contract.output_prefix}/{artifact.relative_path}']
            files.open(path, output.sha256, output.size_bytes)
            paths.append(str(path))
        _, fasta, proof = links._reference(files, reference)
        db = stack.enter_context(duckdb.connect(':memory:'))
        db.execute("SET threads=1")
        db.execute("SET memory_limit='256MB'")
        rows = db.execute('''SELECT id,genomicLocation.chromosome,genomicLocation.start,
            genomicLocation.end FROM read_parquet(?) WHERE id IN ('''
            + ','.join('?' for _ in genes) + ')', [paths, *genes]).fetchall()
        found = {}
        for gene, chrom, start, end in rows:
            # Same source-native one-based inclusive interpretation as the
            # existing common-allele interval owner; convert once to BED geometry.
            chrom = 'chrM' if chrom == 'MT' else 'chr' + str(chrom)
            require(gene not in found and chrom in fasta.index and type(start) is int
                    and type(end) is int and 0 < start <= end <= fasta.index[chrom].length,
                    'gene-coordinate mapping is ambiguous or invalid')
            found[gene] = {'gene_id':gene, 'chromosome':chrom, 'start':start-1, 'end':end,
                'coordinate_system':'zero_based_half_open', 'assembly':'GRCh38', 'kind':'gene',
                'source_manifest_sha256':pin.manifest_sha256, 'source_record_id':gene,
                'source_start':start, 'source_end':end}
        require(set(found) == set(genes), 'essential selected-gene coordinates are incomplete')
        return [found[g] for g in sorted(found)], proof


def build_diagnosis_scope(diagnosis, *, store, gene_source, regulatory_index, reference, retention_index=None):
    """Freeze scoped source facts before VCF interpretation; no phenotype input."""
    resolved, concepts, mapping = resolve_diagnosis(store, diagnosis)
    policy = load_analysis_profile_policy()
    candidates = [a for a in store.assertions if assertion_native_disease_id(store, a) in concepts]
    require(all(a.locus_id is None for a in candidates), 'native locus scope is unsupported')
    rows = [_plain(store.records[a.record_id].source_fields) for a in candidates]
    catalog = build_profile_catalog(rows, policy=policy, source_sha256=store.evidence_sha256,
                                    source_binding=_plain(store.source_binding))
    by_row = {r.original_row_json:r for r in catalog.records}
    assertions, genes = [], set()
    for assertion, row in zip(candidates, rows):
        record = by_row[canonical(row).decode('ascii')]
        included = record.catalog_eligible and assertion.polarity == 'positive'
        if included:
            genes.add(record.gene_id)
        assertions.append({'record_id':assertion.record_id, 'gene_id':record.gene_id,
            'included':included, 'scope_reason':record.reason if not record.catalog_eligible
                else 'established_direct_assertion' if included else 'source_validity_not_positive',
            'assertion':_plain(asdict(assertion)), 'source':row})
    require(bool(genes), 'diagnosis has no established directly supported gene scope')
    genes = sorted(genes)
    try:
        intervals, proof = gene_intervals(gene_source, genes, reference)
        require(regulatory_index.manifest['reference'] == proof, 'regulatory and gene reference bindings differ')
        intervals += project_regulatory_targets(regulatory_index, genes)
    except DiagnosisScopeError:
        raise
    except Exception:
        raise DiagnosisScopeError('diagnosis scope resource input is missing, stale or invalid') from None
    mapping['regulatory_source_coverage'] = regulatory_index.manifest['counts']
    extra_bindings = {}
    if retention_index is not None:
        from genosapiens.track2_retention import known_scope
        from genosapiens.track2_reference import scope_reference_mapping
        reference_mapping = scope_reference_mapping(reference)
        mapping['known_alleles'] = known_scope(retention_index, concepts=concepts, genes=genes,
            reference_identity=reference_mapping['scope_identity_id'],
            reference_sha256=reference.manifest_sha256, catalog_sha256=catalog.sha256)
        extra_bindings['retention_index'] = retention_index.verified.manifest_sha256
    return {'schema':SCHEMA, 'original_diagnosis':diagnosis, 'resolved_diagnosis':resolved,
        'reference_sha256':reference.manifest_sha256,
        'resource_bindings':{'disease_gene':store.evidence_sha256,
            'disease_gene_sources':sha256(canonical(_plain(store.source_binding))),
            'gene_coordinates':gene_source.manifest_sha256,
            'regulatory':regulatory_index.manifest_sha256,
            'reference_contract':reference.contract_sha256, **extra_bindings},
        'policy_sha256':sha256(canonical({'scope':POLICY,'catalog':policy.sha256})),
        'genes':genes, 'intervals':sorted(intervals, key=canonical),
        'assertions':sorted(assertions, key=lambda r:r['record_id']), 'mapping_evidence':mapping,
        'limitations':['Scope membership is not pathogenicity, genotype causality or molecular direction.',
            'Only established direct disease genes and explicit regulatory links are included; no pathway expansion.',
            'Missing/unmapped regulatory targets remain outside scope; no symbol or nearest-gene inference.',
            'Regulatory membership preserves the upstream source scope; it does not establish function or tissue relevance.',
            'SNV/indel interpretation only; SV, CNV, repeat and symbolic partners remain unassessed.']}


def freeze_diagnosis_scope(path, **inputs):
    """Publish only a complete scope; a second write cannot replace prior scope."""
    result = build_diagnosis_scope(**inputs)
    raw = canonical(result)
    with Path(path).open('xb') as target:
        target.write(raw)
    return result, sha256(raw)
