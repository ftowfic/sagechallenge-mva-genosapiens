"""Frozen public resource inputs for the operator's existing scoped producers."""

from contextlib import contextmanager, ExitStack
from dataclasses import dataclass
import json
from pathlib import Path

from genosapiens.admission_links import PinnedPublicSource
from genosapiens.disease_gene_evidence import DiseaseGeneEvidenceStore, rebuild_disease_gene_evidence
from genosapiens.evidence_model import _typed
from genosapiens.evidence_serialization import _plain
from genosapiens.genotype_hypotheses import HypothesisPolicy
from genosapiens.regulatory_index import open_regulatory_index
from genosapiens.track2_codec import decode
from genosapiens.track2_operator_inputs import require, small_json, write_new
from genosapiens.track2_published_genotype import canonical, sha256
from genosapiens.track2_reference import scope_reference_mapping, vep_index_binding, read_vep_index


SCHEMA = 'genosapiens.track2-operator-public-resources/v5'


def export_installed(repository, reference_root, destination, *, retention_index=None, retention_disease_lookup=None):
    """Project the already installed public owners, without resource rebuilds.

    This is installation work. It never opens a VCF or queries case findings.
    Population policy and genotype-specific adjudication are not invented by
    exporting resource locations; absent reviewed inputs remain unassessed.
    """
    from genosapiens import track1_public_validation as pins
    from genosapiens.ontology_bundle import load_ontology_bundle_contract, _parse_obo
    from genosapiens.bounded_json import read_bounded_bytes
    from genosapiens.reference_concordance import load_reference_identity
    from genosapiens.vep_annotation import OpenTargetsConstraintReference, VepReferenceCompatibility
    from genosapiens.vep_bundle import load_vep_bundle_contract
    repository=Path(repository).resolve();references=Path(reference_root).resolve()
    identity=repository/pins.REFERENCE_IDENTITY_PATH
    contract=load_reference_identity(identity)
    manifest=references/'manifests/grch38'/f'{contract.resource_id}.{sha256(identity.read_bytes())[:12]}.json'
    reference=PinnedPublicSource(identity,sha256(identity.read_bytes()),references,manifest,
        pins.REFERENCE_MANIFEST_SHA256)
    genes_contract=repository/pins.OPEN_TARGETS_CONTRACT_PATH
    genes_manifest=references/'manifests/open-targets/open-targets-platform-26.06.0cae4540cfea.json'
    genes_release=references/'releases/open-targets/26.06'
    genes=PinnedPublicSource(genes_contract,sha256(genes_contract.read_bytes()),
        references,genes_manifest,pins.OPEN_TARGETS_MANIFEST_SHA256)
    store={'contract_path':str(repository/pins.RARE_DISEASE_CONTRACT_PATH),
        'reference_root':str(references),'expected_bundle_sha256':pins.RARE_DISEASE_BUNDLE_SHA256}
    print('Bind the installed genotype policy and reference metadata.',flush=True)
    policy=pins.load_public_policy_family(repository,
        model_loci=pins.load_public_model_loci(repository,references)).hypothesis_policy
    regulatory=small_json(repository/'config/track1/regulatory-index-pin.json')
    require(regulatory['reference_contract_sha256']==reference.contract_sha256
        and regulatory['reference_manifest_sha256']==reference.manifest_sha256,
        'The installed regulatory index uses a different reference.')
    print('Read the pinned HPO term dictionary. No phenotype scoring corpus is needed.',flush=True)
    ontology=load_ontology_bundle_contract(repository/pins.ONTOLOGY_CONTRACT_PATH)
    artifact=ontology.artifact('hpo_ontology');hpo_path=references/artifact.path
    hpo_raw=read_bounded_bytes(hpo_path,maximum_bytes=32*1024*1024)
    require(len(hpo_raw)==artifact.size_bytes and sha256(hpo_raw)==artifact.sha256,
        'The installed HPO dictionary changed.')
    hpo_manifest=references/'manifests/ontologies/hpo-ontology-2026-06-23.json'
    hpo_metadata=small_json(hpo_manifest)
    require(hpo_metadata['classification']=='public' and any(
        row['path']==artifact.path and row['sha256']==artifact.sha256
        and row['size_bytes']==artifact.size_bytes for row in hpo_metadata['outputs']),
        'The HPO manifest does not bind the selected dictionary.')
    hpo=_parse_obo(hpo_path,ontology.hpo_ontology)
    clinical,clinvar,_,_=pins._public_consumer_descriptors(repository,references)
    vep={'contract_path':str(repository/pins.VEP_CONTRACT_PATH),'reference_root':str(references),
        'manifest_sha256':pins.VEP_BUNDLE_SHA256,
        'upstream_image_digest':load_vep_bundle_contract(repository/pins.VEP_CONTRACT_PATH).tool_image_digest,
        'constraint':_plain(OpenTargetsConstraintReference(genes_contract,genes_release,genes_manifest,genes.manifest_sha256)),
        'compatibility':_plain(VepReferenceCompatibility(identity,reference.contract_sha256,references,
            manifest,reference.manifest_sha256)),
        'fai':vep_index_binding(repository/pins.VEP_CONTRACT_PATH,references)}
    from genosapiens.rare_disease import RareDiseaseContract, _bundle_manifest_path
    from genosapiens.disease_gene_evidence import PUBLIC_TABLES
    rare_path=repository/pins.RARE_DISEASE_CONTRACT_PATH
    rare=RareDiseaseContract.from_mapping(small_json(rare_path))
    mounts=[manifest,genes_manifest,genes_release,clinvar.manifest_path,
        references/regulatory['index_path'],
        references/'manifests/ensembl/vep',references/'releases/ensembl/vep/116_GRCh38',
        _bundle_manifest_path(references,rare,sha256(rare_path.read_bytes())),
        *(references/item.path for item in contract.expected_outputs),
        *(references/rare.projection.output_prefix/(table+'.parquet') for table in PUBLIC_TABLES)]
    mounts.extend(references/row['path'] for row in small_json(clinvar.manifest_path)['outputs'])
    from genosapiens.installed_public_index import load_installed_index_pin
    index_pin = retention_index or load_installed_index_pin(repository)
    require(index_pin is not None, 'Prepare the shared public selection index before exporting Track 2 resources.')
    index_manifest = small_json(references/index_pin['manifest_path'])
    require(sha256((references/index_pin['manifest_path']).read_bytes()) == index_pin['manifest_sha256'],
        'The installed selection index manifest changed.')
    mounts.extend([references/index_pin['manifest_path'],
                   *(references/row['path'] for row in index_manifest['outputs'].values())])
    if retention_disease_lookup is not None:
        from genosapiens.public_disease_lookup import load_manifest
        _, lookup = load_manifest(references, retention_disease_lookup, index_pin['manifest_sha256'])
        mounts.extend([references/retention_disease_lookup['manifest_path'], lookup])
    from genosapiens.track2_functional_sources import installed_catalog
    functional_catalog = installed_catalog(repository, references)
    if functional_catalog is not None:
        functional_manifest = Path(functional_catalog['manifest_path'])
        mounts.extend([functional_manifest, functional_manifest.parent/'catalog.sqlite'])
        print('Bind the frozen functional source for gene-level local lookup.', flush=True)
    return publish_resources(destination,reference=reference,gene_source=genes,store=store,policy=policy,
        regulatory_root=references/regulatory['index_path'],regulatory_sha256=regulatory['manifest_sha256'],
        hpo={'ontology':_plain(hpo),
            'manifest_sha256':sha256(hpo_manifest.read_bytes())},vep=vep,clinvar=clinvar,open_targets=clinical,
        public_mounts=mounts, retention_index=index_pin, retention_index_root=references,
        retention_disease_lookup=retention_disease_lookup, functional_catalog=functional_catalog)


def validate_resource_roots(data):
    """Check consumer path semantics without opening or rehashing public outputs."""
    for name in ('reference', 'gene_source', 'clinvar', 'open_targets'):
        source = data.get(name)
        if source is None:
            continue
        root, manifest = Path(source['root']), Path(source['manifest_path'])
        require(root.is_absolute() and manifest.is_absolute()
                and '..' not in root.parts and '..' not in manifest.parts
                and manifest != root and manifest.is_relative_to(root),
                'A public source root must contain its manifest and release. Re-export the resource bindings.')


def publish_resources(destination, *, reference, gene_source, store, policy,
                      regulatory_root, regulatory_sha256,
                      hpo=None, vep=None, clinvar=None, open_targets=None, population=None, public_mounts=None,
                      retention_index=None, retention_index_root=None, retention_disease_lookup=None,
                      functional_catalog=None):
    """Public preparation API; no case data, source selection or ready flags."""
    require(type(reference) is PinnedPublicSource and type(gene_source) is PinnedPublicSource
            and (type(store) is DiseaseGeneEvidenceStore or type(store) is dict
                and set(store)=={'contract_path','reference_root','expected_bundle_sha256'})
            and type(policy) is HypothesisPolicy,
            'Resource publication requires the existing scientific source owners.')
    value = {'schema': SCHEMA, 'classification': 'public',
        'reference': _plain(reference), 'gene_source': _plain(gene_source),
        'store': _plain(store), 'policy': _plain(policy),
        'regulatory': {'root': str(regulatory_root), 'sha256': regulatory_sha256},
        'reference_mapping': scope_reference_mapping(reference), 'hpo': hpo, 'vep': vep,
        'clinvar': _plain(clinvar), 'open_targets': _plain(open_targets),
        'population': _plain(population),
        'retention_index': retention_index,
        'retention_index_root': str(retention_index_root or reference.root),
        'retention_disease_lookup': retention_disease_lookup,
        'functional_catalog': functional_catalog,
        'public_mounts':sorted(set(str(p) for p in (public_mounts if public_mounts is not None
            else (reference.root,gene_source.root,reference.contract_path,gene_source.contract_path,
                  reference.manifest_path,gene_source.manifest_path,regulatory_root))))}
    # Avoid overlapping bind mounts while retaining every selected resource.
    value['public_mounts']=[p for p in value['public_mounts'] if not any(
        p!=other and Path(p).is_relative_to(other) for other in value['public_mounts'])]
    def encode(item):
        if isinstance(item,Path):return str(item)
        if isinstance(item,bytes):return item.hex()
        if isinstance(item,(set,frozenset)):return sorted(item)
        raise TypeError('Unsupported public resource value')
    value=json.loads(json.dumps(value,default=encode))
    validate_resource_roots(value)
    write_new(destination, value)
    return value


@dataclass
class Resources:
    data: dict
    reference: PinnedPublicSource
    genes: PinnedPublicSource
    store: DiseaseGeneEvidenceStore
    policy: HypothesisPolicy
    index: object
    retention_index: object | None = None


@contextmanager
def open_resources(path, expected_sha256, *, diagnosis=None):
    raw = Path(path).read_bytes()
    require(sha256(raw) == expected_sha256, 'Installed public resource bindings changed.')
    data = small_json(path, maximum=128 * 1024 * 1024)
    require(set(data) == {'schema','classification','reference','gene_source','store','policy',
            'regulatory','reference_mapping','hpo','vep','clinvar','open_targets','population','public_mounts',
            'retention_index','retention_index_root','retention_disease_lookup','functional_catalog'}
            and data['schema'] == SCHEMA and data['classification'] == 'public',
            'Installed public resource record is invalid.')
    reference = decode(PinnedPublicSource, data['reference'])
    require(data['reference_mapping'] == scope_reference_mapping(reference),
            'Installed reference-role mappings changed. Re-export the resource bindings.')
    genes = decode(PinnedPublicSource, data['gene_source'])
    if set(data['store'])=={'contract_path','reference_root','expected_bundle_sha256'}:
        from genosapiens.disease_gene_evidence import read_pinned_disease_gene_evidence
        require(diagnosis is not None,'The installed disease sources require an explicit diagnosis.')
        store=read_pinned_disease_gene_evidence(**data['store'],diagnosis=diagnosis,track2_concepts=True)
    else:
        supplied = decode(DiseaseGeneEvidenceStore, data['store'])
        store = rebuild_disease_gene_evidence(supplied)
        require(store.evidence_sha256 == supplied.evidence_sha256,
                'Frozen disease evidence differs from its source assertions.')
        if diagnosis is not None:
            from genosapiens.track2_concept_supplement import supplement_store
            store = supplement_store(store)
    policy = _typed(data['policy'], HypothesisPolicy, decode=True)
    with ExitStack() as stack:
        index = stack.enter_context(open_regulatory_index(data['regulatory']['root'],
            expected_manifest_sha256=data['regulatory']['sha256']))
        retention = None
        if data['retention_index'] is not None:
            from genosapiens.installed_public_index import validate_index_pin
            from genosapiens.common_allele_public_indexes import (
                verify_common_allele_public_indexes, CommonAllelePublicIndexSession)
            from genosapiens.bounded_wgs_config import PublicQueryLimits
            pin = validate_index_pin(data['retention_index'])
            verified = verify_common_allele_public_indexes(reference_root=data['retention_index_root'],
                manifest_path=pin['manifest_path'], expected_manifest_sha256=pin['manifest_sha256'],
                verify_sources=False, reuse_qualified_reconciliation=True)
            retention = stack.enter_context(CommonAllelePublicIndexSession(verified,
                limits=PublicQueryLimits(1024, 10, 64 * 1024 * 1024, 300)))
            retention.retention_coverage()
            if data['retention_disease_lookup'] is not None:
                from genosapiens.public_disease_lookup import attach_disease_lookup
                stack.enter_context(attach_disease_lookup(retention, data['retention_index_root'], data['retention_disease_lookup']))
            require(retention.retention_coverage().get('disease_lookup_version') == 1
                or retention.disease_lookup is not None,
                'Prepare the public disease lookup before Track 2 construction. Resume with the current launcher.')
        yield Resources(data, reference, genes, store, policy, index, retention)


def interpret(selection, scope, resources, output, runtime, *, annotation_root=None):
    from genosapiens.track2_annotation_artifacts import FileAnnotationInput
    from genosapiens.track2_scoped_interpretation import ScopedSources, selected_alleles, interpret_selection
    from genosapiens.vep_annotation import (VepRequest, run_vep_annotation,
        OpenTargetsConstraintReference, VepReferenceCompatibility)
    annotation = None
    config = resources.data['vep']
    alleles = selected_alleles(selection)
    if config is not None and alleles:
        index_bytes = read_vep_index(config)
        requests = tuple(VepRequest(a) for a in alleles)
        from genosapiens.track2_annotation_checkpoint import dependency_binding, load_or_run
        destination = (Path(output) if annotation_root is None else Path(annotation_root)) / 'annotation'
        destination.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
        def produce(path):
            run_vep_annotation(requests, contract_path=config['contract_path'],
                reference_root=config['reference_root'], bundle_manifest_sha256=config['manifest_sha256'],
                runtime_image=runtime['image'], upstream_image_digest=config['upstream_image_digest'],
                output_directory=path,
                constraint_reference=decode(OpenTargetsConstraintReference, config['constraint']),
                reference_compatibility=decode(VepReferenceCompatibility, config['compatibility']))
        from genosapiens.track2_annotation_recovery import validate_recovery
        def recover(current):
            validate_recovery(current, requests=requests, index_bytes=index_bytes, image=runtime['image'])
        artifacts = load_or_run(destination, binding=dependency_binding(requests, config, runtime,
            Path(__file__).resolve().parents[2]), produce=produce, recover=recover)
        require(read_vep_index(config) == index_bytes, 'The annotation index changed during annotation.')
        artifacts['reference.fai'] = index_bytes
        annotation = FileAnnotationInput(artifacts, requests, sha256(artifacts['binding.json']))
    population = None
    if resources.data['population'] is not None:
        from genosapiens.genotype_hypotheses import HypothesisPopulationInputs
        record=dict(resources.data['population'])
        record['policy_bytes']=bytes.fromhex(record['policy_bytes'])
        record['public_curation_bytes']={key:bytes.fromhex(value)
            for key,value in record['public_curation_bytes'].items()}
        population = decode(HypothesisPopulationInputs,record)
    sources = ScopedSources(annotation, resources.index, resources.reference,
        *(None if resources.data[k] is None else decode(PinnedPublicSource, resources.data[k])
          for k in ('clinvar','open_targets')), population)
    require(resources.retention_index is not None,
        'Complete local clinical and curated selection evidence is required before interpretation.')
    if not selection.selection['selected_alleles']:
        return interpret_selection(selection, scope=scope, store=resources.store, sources=sources)
    require(sources.clinvar is not None and sources.open_targets is not None,
        'Install the ClinVar and Open Targets sources before interpretation.')
    from genosapiens.wgs_clinical_session import WgsClinicalResourceSession
    from genosapiens.bounded_wgs_config import load_bounded_wgs_execution_config
    from dataclasses import replace
    raw = (Path(__file__).resolve().parents[2]/'config/track1/bounded-wgs-execution-v5.json').read_bytes()
    config = load_bounded_wgs_execution_config(raw, sha256(raw))
    with WgsClinicalResourceSession(indexes=resources.retention_index.verified,
            open_targets=sources.open_targets, clinvar=sources.clinvar,
            reference=resources.reference, execution_config=config,
            scratch_path=Path(output).resolve()/'clinical-session') as session:
        from genosapiens.wgs_allele_query_cache import adapter_identity, clinical_source_bindings
        bound = selection.execution['binding']
        session.use_case_cache(path=session.scratch_path/'case-queries.sqlite', binding=canonical({
            'schema':'genosapiens.wgs-case-allele-query-cache','version':1,
            'consumer_executable_sha256':sha256(Path(__file__).read_bytes()),
            'case_input_sha256':bound['source_sha256'],
            'reference_role':bound['reference_role'],'reference_identity_id':bound['reference_identity_id'],
            'reference_manifest_sha256':resources.reference.manifest_sha256,
            'population_sources':{'availability':'not_queried_in_this_session'},
            'clinical_sources':clinical_source_bindings(session),
            'query_semantics_sha256':sha256(raw),'adapter_sha256':adapter_identity()}))
        return interpret_selection(selection, scope=scope, store=resources.store,
            sources=replace(sources, clinical_session=session))


def phenotype(path, resources, digest):
    from genosapiens.track2_phenotype import ingest_optional_phenotype
    if path is None:
        return ingest_optional_phenotype()
    from genosapiens.ontology_bundle import ParsedOntology, OntologyTerm
    record = resources.data['hpo']
    require(record is not None, 'Install the frozen HPO ontology before using supplied HPO terms.')
    raw=record['ontology']
    ontology=ParsedOntology({k:decode(OntologyTerm,v) for k,v in raw['terms'].items()},
        raw['active_aliases'],raw['parent_edges'],raw['obsolete_mappings'])
    return ingest_optional_phenotype(path, expected_input_sha256=digest,
        hpo_ontology=ontology,
        ontology_manifest_sha256=record['manifest_sha256'])
