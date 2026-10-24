"""Pinned MONDO identities and subclass evidence for Track 2 only."""
from dataclasses import asdict, replace
from graphlib import TopologicalSorter
import json
from pathlib import Path
from types import MappingProxyType
import re

from genosapiens.track2_published_genotype import canonical, sha256
from genosapiens.track2_operator_inputs import require

PATH = 'config/track2/concept-supplement.json'
POLICY = 'track2-pinned-mondo-scope/1.0.0'
PIN_SHA256 = 'e2899d46433d20f99e6df1c1bca0b6430f4d15aa59357970d965ef8130d04942'


def exact_identifier(url):
    patterns = ((r'https?://omim.org/entry/(\d+)', 'OMIM:'),
                (r'https?://omim.org/phenotypicSeries/PS(\d+)', 'OMIMPS:'),
                (r'https?://www.orpha.net/ORDO/Orphanet_(\d+)', 'ORPHA:'),
                (r'http://purl.obolibrary.org/obo/MONDO_(\d+)', 'MONDO:'))
    for pattern, prefix in patterns:
        match = re.fullmatch(pattern, url)
        if match:
            return prefix + match[1]
    return None


def checked_supplement(data):
    require(data['schema'] == 'genosapiens.track2-concept-supplement/v1'
            and data['policy'] == POLICY, 'Unsupported Track 2 concept supplement.')
    def response(source):
        value = source['response']
        require(sha256(canonical(value)) == source['response_sha256'], 'Pinned concept response changed.')
        require(source['url'].startswith('https://www.ebi.ac.uk/ols4/api/ontologies/mondo'),
                'Unsupported concept evidence source.')
        return value
    for name in ('ontology_before', 'ontology_after'):
        require(response(data[name])['config']['versionIri'] == data['release'],
                'Concept source versions disagree.')
    concepts, mappings, parents = {}, {}, {}
    for record in data['records']:
        term = response(record['term'])
        cid = term['obo_id']
        require(re.fullmatch(r'MONDO:\d+', cid) and cid not in concepts
                and not term['is_obsolete'] and term['is_defining_ontology'],
                'Concept identity is ambiguous or inactive.')
        concepts[cid] = record
        aliases = {cid} | {identifier for url in term.get('annotation', {}).get('has exact match', [])
                          if (identifier := exact_identifier(url))}
        for alias in aliases:
            require(alias not in mappings, 'Ambiguous pinned exact concept mapping.')
            mappings[alias] = cid
        related = response(record['parents'])
        rows = related.get('_embedded', {}).get('terms', [])
        require(related['page']['totalElements'] == len(rows), 'Concept parent coverage is incomplete.')
        parents[cid] = {row['obo_id'] for row in rows}
    declared = set(data['coverage']['concept_ids'])
    require(declared == set(concepts), 'Concept coverage differs from its declaration.')
    descendants = response(data['descendants'])
    rows = descendants.get('_embedded', {}).get('terms', [])
    require(descendants['page']['totalElements'] == len(rows)
            and {r['obo_id'] for r in rows} | set(data['coverage']['roots']) == declared,
            'Concept descendant coverage is incomplete.')
    graph = {key: value & declared for key, value in parents.items()}
    tuple(TopologicalSorter(graph).static_order())
    reachable = set(data['coverage']['roots'])
    while expanded := {key for key, value in graph.items() if value & reachable} - reachable:
        reachable.update(expanded)
    require(reachable == declared, 'A concept lacks an explicit subtype path.')
    return concepts, mappings, graph


def installed_supplement():
    raw = (Path(__file__).resolve().parents[2]/PATH).read_bytes()
    require(sha256(raw) == PIN_SHA256, 'The pinned concept supplement changed.')
    data = json.loads(raw)
    checked_supplement(data)
    return data


def validate_concept_record(fields):
    from genosapiens.disease_gene_evidence import _plain
    fields = _plain(fields)
    require(set(fields) == {'disease_id', 'status', 'classification_level', 'release', 'record'}
            and re.fullmatch(r'MONDO:\d+', fields['disease_id'])
            and fields['status'] == 'Active'
            and fields['classification_level'] in {'Disorder', 'Subtype of disorder'}
            and fields['release'].startswith('http://purl.obolibrary.org/obo/mondo/releases/'),
            'Invalid retained Track 2 concept identity.')
    for source in fields['record'].values():
        require(sha256(canonical(source['response'])) == source['response_sha256'],
                'Retained concept response changed.')
    term = fields['record']['term']['response']
    require(term['obo_id'] == fields['disease_id'] and not term['is_obsolete'],
            'Retained concept identity differs from its source.')


def supplement_store(store, data=None):
    """Keep legacy conflicts as unresolved records, without changing assertions."""
    from genosapiens.disease_gene_evidence import SourceRecord, DiseaseEvidenceRelation, _freeze, _plain
    data = installed_supplement() if data is None else data
    concepts, mappings, graph = checked_supplement(data)
    identity = sha256(canonical(data))
    previous = store.source_binding.get('track2_concept_supplement')
    if previous is not None:
        require(previous['sha256'] == identity, 'A different concept supplement requires reconstruction.')
        return store
    records = dict(store.records)
    relations = []
    # Reject indirect legacy paths that would merge distinct pinned concepts.
    components = {alias: {alias} for alias in mappings}
    def component(key):
        return components.setdefault(key, {key})
    for alias, cid in mappings.items():
        merged = component(alias) | component(cid)
        for key in merged:
            components[key] = merged
    for relation in store.relations:
        if relation.kind == 'equivalent':
            merged = component(relation.subject_id) | component(relation.object_id)
            owners = {mappings[k] for k in merged if k in mappings}
            if len(owners) > 1:
                relations.append(replace(relation, kind='unknown', diagnostics=(*relation.diagnostics,
                    'track2_legacy_mapping_conflicts_with_pinned_mondo')))
                continue
            for key in merged:
                components[key] = merged
        relations.append(relation)
    for cid, record in concepts.items():
        key = 'track2-mondo:'+cid
        records[key] = SourceRecord(key, 'track2_concept', _freeze(dict(disease_id=cid,
            status='Active', classification_level='Disorder' if cid in data['coverage']['roots'] else 'Subtype of disorder',
            release=data['release'], record=record)))
        for alias, owner in mappings.items():
            if owner == cid and alias != cid:
                relations.append(DiseaseEvidenceRelation(cid, alias, 'equivalent', (key,), ('pinned_mondo_exact_mapping',)))
        for parent in sorted(graph[cid]):
            relations.append(DiseaseEvidenceRelation(cid, parent, 'narrower', (key, 'track2-mondo:'+parent),
                ('active_subtype_aggregation', 'pinned_mondo_subclass')))
    binding = _freeze({**_plain(store.source_binding), 'track2_concept_supplement': dict(
        sha256=identity, policy=POLICY, release=data['release'], coverage=data['coverage'])})
    return replace(store, records=MappingProxyType(records), relations=tuple(relations), source_binding=binding,
        evidence_sha256=sha256(canonical(dict(source=store.evidence_sha256, supplement=identity,
            policy=POLICY, relations=[_plain(asdict(r)) for r in relations]))))
