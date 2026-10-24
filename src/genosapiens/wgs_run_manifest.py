"""Current run inputs and direct stage dependencies, without attempt histories."""
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import stat
import tempfile

from genosapiens.bounded_json import read_bounded_bytes
from genosapiens.wgs_ledger import canonical

SCHEMA = 'genosapiens.track1-current-run'
VERSION = 1
PREPARATION_CONTRACT = 'track1-early-scope/v1'
LEGACY_PREPARATION_MESSAGE = ('This run predates the early scope contract. Use a fresh output directory with this package. '
    'Preserve the original run and package for separate recovery.')
STAGES = ('preparation', 'annotation', 'selection', 'associations', 'science', 'reports')
PARENTS = dict(preparation=(), annotation=('preparation',), selection=('annotation',),
    associations=('selection',), science=('associations',), reports=('science',))
MAX_BYTES = 2 * 1024**2
MAX_OUTPUT_BYTES = 8 * 1024**2

class RunManifestError(ValueError):
    """Static errors omit private paths and values."""


def digest(content):
    return hashlib.sha256(content).hexdigest()


def _require(condition, message='current run manifest is invalid'):
    if not condition:
        raise RunManifestError(message)


def _sha(value):
    return type(value) is str and len(value) == 64 and all(c in '0123456789abcdef' for c in value)


def atomic_current(path, value):
    """Replace one current document after a complete durable write."""
    path = Path(path)
    content = canonical(value)
    _require(len(content) <= MAX_BYTES and path.parent.is_dir())
    _require(not any(p.is_symlink() for p in (path, *path.parents)))
    descriptor, name = tempfile.mkstemp(prefix='.current-', dir=path.parent)
    try:
        os.fchmod(descriptor, 0o600)
        with os.fdopen(descriptor, 'wb') as stream:
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(name, path)
        directory = os.open(path.parent, os.O_DIRECTORY | os.O_RDONLY)
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
    finally:
        if os.path.exists(name):
            os.unlink(name)


def load_current(path):
    path = Path(path)
    _require(not any(p.is_symlink() for p in (path, *path.parents)))
    metadata = path.stat()
    _require(stat.S_ISREG(metadata.st_mode) and metadata.st_uid == os.getuid()
        and stat.S_IMODE(metadata.st_mode) in (0o400, 0o600))
    content = read_bounded_bytes(path, maximum_bytes=MAX_BYTES)
    try:
        value = json.loads(content)
        _require(type(value) is dict and canonical(value) == content
            and value.get('schema') == SCHEMA and value.get('version') == VERSION,
            'current run format is unsupported; rebuild intake from the original inputs')
        _require(type(value.get('inputs')) is dict and type(value.get('stages')) is dict
            and not set(value['stages']) - set(STAGES))
        return value
    except (KeyError, TypeError, json.JSONDecodeError):
        raise RunManifestError('current run manifest is invalid') from None


def resolved_inputs(value, supplied, *, host=False):
    """Resolve omitted inputs and reject explicit conflicts before any writes."""
    saved = value['inputs']
    resolved = {}
    for name in ('vcf', 'hpo', 'sample', 'concordance'):
        key = 'host_' + name if host and name in ('vcf', 'hpo', 'concordance') else name
        previous = saved.get(key, saved.get(name))
        current = supplied.get(name)
        if current is not None:
            if name != 'sample':
                current = str(Path(current).resolve(strict=True))
            _require(current == previous, 'explicit input conflicts with the current run manifest')
        resolved[name] = previous
    _require(all(resolved.get(name) for name in ('vcf', 'hpo', 'sample')),
        'current run input binding is incomplete')
    return resolved


@dataclass(frozen=True, slots=True)
class RestartPlan:
    restart_stage: str | None
    reused_stages: tuple[str, ...]
    invalid_stages: tuple[str, ...]
    reason: str


def dependency_binding(stage, direct, completed):
    _require(stage in STAGES and type(direct) is dict and direct
        and all(type(k) is str and _sha(v) for k, v in direct.items()))
    parents = {parent: completed[parent]['output_sha256'] for parent in PARENTS[stage]}
    _require(all(_sha(value) for value in parents.values()))
    return dict(stage=stage, direct=dict(sorted(direct.items())), parents=parents)


def restart_plan(value, contracts, *, output_root):
    """Inspect only bounded completion metadata and each named output manifest."""
    _require(set(contracts) == set(STAGES))
    root = Path(output_root).resolve(strict=True)
    completed = value['stages']
    reused, invalid = [], []
    reason = 'all stages have matching completed dependencies'
    for stage in STAGES:
        previous = completed.get(stage)
        if any(parent in invalid for parent in PARENTS[stage]):
            invalid.append(stage)
            continue
        if previous is None:
            invalid.append(stage)
            if len(invalid) == 1:
                reason = 'stage has no complete output'
            continue
        _require(set(previous) == {'dependencies', 'outputs', 'output_sha256', 'status'}
            and previous['status'] == 'complete' and _sha(previous['output_sha256']))
        expected = dependency_binding(stage, contracts[stage], completed)
        if previous['dependencies'] != expected:
            invalid.append(stage)
            if len(invalid) == 1:
                reason = 'direct code, policy, source, or parent dependency changed'
            continue
        outputs = previous['outputs']
        _require(type(outputs) is list and outputs)
        for output in outputs:
            _require(set(output) == {'path', 'sha256', 'size_bytes'} and _sha(output['sha256']))
            relative = Path(output['path'])
            _require(not relative.is_absolute() and '..' not in relative.parts and str(relative) == output['path'])
            path = root / relative
            _require(not any(p.is_symlink() for p in (path, *path.parents)))
            if not path.is_file():
                invalid.append(stage)
                if len(invalid) == 1:
                    reason = 'completed stage output is missing'
                break
            content = read_bounded_bytes(path, maximum_bytes=MAX_OUTPUT_BYTES,
                expected_size_bytes=output['size_bytes'])
            _require(digest(content) == output['sha256'], 'completed stage output integrity check failed')
        else:
            _require(digest(canonical(outputs)) == previous['output_sha256'])
            reused.append(stage)
    return RestartPlan(invalid[0] if invalid else None, tuple(reused), tuple(invalid), reason)


def completed_stage(stage, direct, completed, *, output_root, paths):
    root = Path(output_root).resolve(strict=True)
    outputs = []
    for path in paths:
        path = Path(path)
        _require(not any(p.is_symlink() for p in (path, *path.parents)) and path.is_relative_to(root))
        content = read_bounded_bytes(path, maximum_bytes=MAX_OUTPUT_BYTES)
        outputs.append(dict(path=str(path.relative_to(root)), sha256=digest(content), size_bytes=len(content)))
    _require(outputs and len({v['path'] for v in outputs}) == len(outputs))
    return dict(status='complete', dependencies=dependency_binding(stage, direct, completed),
        outputs=outputs, output_sha256=digest(canonical(outputs)))

# Each stage binds its implementation and its direct public inputs. Package
# identity remains provenance and never substitutes for these dependencies.
CODE_MODULES = {
    'preparation': ('wgs_preparation', 'wgs_filter', 'wgs_prepared_case', 'vcf_qc',
        'vcf_contracts', 'vcf_context', 'vcf_header', 'allele_identity',
        'common_allele_filter', 'population_clinical_session', 'wgs_admission', 'track1_early_scope'),
    'annotation': ('wgs_annotation', 'vep_annotation', 'vep_bundle', 'wgs_annotation_coverage',
        'wgs_native_annotation', 'vep_native_annotation'),
    'selection': ('variant_retention', 'retention_sources', 'curated_alleles',
        'wgs_profile_projection', 'track1_analysis_profiles', 'wgs_clinical_session',
        'common_allele_public_indexes', 'wgs_analysis_eligibility', 'wgs_admission_projection'),
    'associations': ('wgs_association_projection', 'track1_association_scope', 'candidate_admission',
        'annotation_admission', 'regulatory_index', 'admission_links', 'wgs_context_pipeline'),
    'science': ('wgs_contexts', 'wgs_context_admission', 'genotype_hypotheses',
        'genotype_properties', 'wgs_pair_scores', 'wgs_context_scores', 'wgs_hypothesis_evidence', 'model_requirements', 'wgs_assessment_checkpoint', 'wgs_scientific_products', 'wgs_shards',
        'wgs_scientific_artifact', 'wgs_scientific_checkpoint', 'track1_epcr'),
    'reports': ('wgs_science_reporting', 'wgs_reporting', 'wgs_primary_publication', 'wgs_handoff', 'track1_methods',
        'track1_preflight', 'track1_publication'),
}

# These direct helpers own scientific behavior used inside each stage. Later
# consumers bind parent outputs, so report helpers cannot invalidate annotation.
STAGE_HELPERS = {
    'preparation': ('variant_normalization', 'variant_inventory', 'population_annotation',
        'installed_public_index',
        'phenotype_ingestion', 'reference_concordance', 'wgs_source_routing',
        'wgs_public_filter', 'common_allele_exception_bundle', 'gnomad_faf95_adapter'),
    'annotation': ('annotation_admission', 'grch38_reference', 'variant_normalization',
        'wgs_annotation_prefix', 'wgs_read_cache', 'wgs_shards'),
    'selection': ('disease_gene_evidence', 'regulatory_index', 'admission_links',
        'wgs_allele_query_cache', 'wgs_admission_screening', 'wgs_public_filter',
        'wgs_admission_diagnostics', 'wgs_clinvar_reference_sites', 'tabix'),
    'associations': ('wgs_association_codec', 'wgs_admission_graph', 'disease_gene_evidence',
        'disease_profiles', 'hpo_ic', 'rare_disease', 'wgs_context_payloads'),
    'science': ('genotype_evaluation', 'genotype_model', 'evidence_fusion', 'phenotype_scoring',
        'partner_qualification', 'partner_relevance', 'population_annotation', 'population_policy',
        'population_filter_candidate', 'clinical_annotation', 'track1_sensitivity',
        'wgs_case_context', 'wgs_case_evidence', 'wgs_case_resources', 'wgs_completed_science',
        'wgs_context_evidence', 'wgs_context_payloads', 'wgs_context_runtime', 'wgs_evidence_replay',
        'wgs_evidence_store', 'wgs_gene_evidence', 'wgs_genotype_evaluation', 'wgs_packed_outcomes',
        'wgs_pair_evidence', 'wgs_pair_outcomes', 'wgs_phenotype_evidence', 'wgs_population_evidence',
        'wgs_population_work_filter', 'wgs_reference_evidence', 'wgs_score_decisions',
        'wgs_score_sources', 'wgs_population_input', 'wgs_scientific',
        'track1_science_scope', 'track1_science_protection', 'track1_science_roster',
        'wgs_context_domain', 'wgs_science_aggregation', 'wgs_science_checkpoints',
        'wgs_science_coordinator', 'wgs_science_execution', 'wgs_science_failures',
        'wgs_science_workers', 'wgs_science_parallel', 'wgs_science_publication',
        'wgs_science_options', 'wgs_science_runtime', 'wgs_science_scope_evidence',
        'wgs_science_selection', 'wgs_science_resources', 'wgs_science_hpo', 'wgs_resource_rows', 'wgs_ic_evidence', 'wgs_science_graph', 'wgs_science_transport', 'wgs_validation_scope', 'wgs_case_order', 'wgs_case_identities',
        'wgs_case_science_facts', 'wgs_selection_inputs', 'wgs_source_closure'),
    'reports': ('case_manifest', 'challenge_evaluator', 'track1_artifacts', 'track1_handoff',
        'track1_submission', 'submission_rendering', 'wgs_logical_evidence', 'wgs_provenance'),
}

# Mixed orchestration modules are bound at function boundaries. This prevents a
# report renderer edit from invalidating completed upstream scientific producers.
STAGE_FUNCTIONS = {
    'preparation': {'wgs_public_runner': ('_preparation', '_load_preparation_resources', '_indexes'),
        'common_allele_public_indexes': ('verify_common_allele_public_indexes', '_bound_query',
            'query_many', 'query', 'assert_stable')},
    'annotation': {'wgs_public_analysis': ('_pinned_annotation_coverage', '_resource_session_sha256',
        '_acquire_admission_projection')},
    'selection': {'wgs_public_analysis': ('_admission_inputs', '_intake_binding',
        '_admission_projection_jobs', '_admission_projection_binding', '_admission_stage_binding',
        '_publish_admission_projection', '_publish_admission_projection_owned',
        '_screening_session', '_build_with_population_input')},
    'associations': {'wgs_public_analysis': ('_disease_subset',)},
    'science': {'wgs_public_analysis': ('_science_disease_subset', '_build', '_build_with_session', '_build_with_query_cache', '_build_scoped_with_population_input', 'checkpoint_identity', '__post_init__')},
    'reports': {'wgs_public_analysis': ('_publish_scientific_products', '_selected_technical_assessments',
        'reuse_current_scientific', 'assemble_primary_publication'),
        'wgs_analysis': ('_checkpoint_products', '_load_checkpointed_products', 'analyze_bounded_case'),
        'wgs_public_runner': ('_analysis',)},
}


def function_identity(read, module, names):
    import ast
    tree = ast.parse(read('python/genosapiens/' + module + '.py'))
    found = {node.name:node for node in ast.walk(tree) if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and node.name in names}
    _require(set(found) == set(names), 'required stage entry is missing')
    constants = [node for node in tree.body if isinstance(node, (ast.Assign, ast.AnnAssign, ast.Import, ast.ImportFrom))]
    return digest(canonical(dict(functions={name:ast.dump(found[name],include_attributes=False) for name in sorted(found)},
        module_bindings=[ast.dump(node,include_attributes=False) for node in constants])))


def require_scoped_upstream_equivalence(read):
    """The new entry cannot change upstream decisions under the legacy identity."""
    import ast
    tree = ast.parse(read('python/genosapiens/wgs_public_analysis.py'))
    names = {'_build_with_population_input', '_build_scoped_with_population_input'}
    methods = {n.name: n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name in names}
    if '_build_scoped_with_population_input' not in methods:
        return
    _require(set(methods) == names, 'legacy upstream entry is missing')
    original = methods['_build_with_population_input'].body
    scoped = methods['_build_scoped_with_population_input'].body
    instrumentation = ast.parse('import time\nscience_setup_started = time.monotonic()').body
    _require([ast.dump(n) for n in scoped[:2]] == [ast.dump(n) for n in instrumentation],
        'scoped upstream instrumentation changed')
    def prefix(nodes):
        result = []
        for node in nodes:
            result.append(ast.dump(node))
            if (isinstance(node, ast.Expr) and isinstance(node.value, ast.Call)
                    and isinstance(node.value.func, ast.Attribute) and node.value.func.attr == 'stage_complete'
                    and node.value.args and isinstance(node.value.args[0], ast.Constant)
                    and node.value.args[0].value == 'associations'):
                return result
        raise RunManifestError('scoped upstream completion boundary is missing')
    _require(prefix(original) == prefix(scoped[2:]),
        'scoped upstream decisions differ and require a new upstream dependency')


def installation_contracts(contract, versions, files, read):
    """Use authenticated installation metadata without reading resource payloads."""
    shared = ('evidence_model', 'evidence_serialization', 'wgs_ledger', 'bounded_wgs_config',
        'bounded_json', 'wgs_query_limits', 'wgs_runtime', 'wgs_failures')
    def code(stage):
        if stage == 'selection':
            require_scoped_upstream_equivalence(read)
        rows = {}
        for name in (*shared, *CODE_MODULES[stage], *STAGE_HELPERS[stage]):
            path = 'python/genosapiens/' + name + '.py'
            _require(path in files, 'required stage implementation is missing')
            rows[path] = files[path]['sha256']
        for module, names in STAGE_FUNCTIONS[stage].items():
            rows[module] = function_identity(read, module, names)
        return digest(canonical(rows))
    selectors = {
        'preparation': dict(reference=versions.reference_manifest_sha256,
            early_scope_policy=contract.scope_policy_sha256, early_scope_manifest=contract.scope_manifest_sha256,
            policy=contract.common_allele_policy_sha256, population=contract.population_manifest_sha256,
            clinical=contract.clinical_manifest_sha256, exceptions=contract.exception_manifest_sha256,
            precision=contract.precision_manifest_sha256, ontology=contract.ontology_manifest_sha256),
        'annotation': dict(reference=versions.reference_manifest_sha256),
        'selection': dict(policy=contract.analysis_profile_policy_sha256,
            evidence=contract.common_allele_public_index_manifest_sha256),
        'associations': dict(evidence=contract.common_allele_public_index_manifest_sha256),
        'science': dict(policy=versions.policy_bindings.sha256,
            public_resources=digest(canonical(versions.resource_manifest_sha256s)),
            work=contract.focused_work_contract_sha256),
        'reports': dict(schema=contract.canonical_handoff_schema_sha256),
    }
    config_scopes = {
        'preparation': ('config/reference-identities/', 'config/track1/vcf-evidence-',
            'config/public-resources/track1-common-allele-', 'config/public-resources/gnomad-'),
        'science': ('config/track1/public-ranking-', 'config/track1/public-population-',
            'config/track1/epcr-', 'config/public-resources/grch38-model-loci'),
        'annotation': ('config/public-resources/ensembl-vep-',),
        'selection': ('config/public-resources/rare-disease', 'config/track1/regulatory-index-pin.json'),
        'associations': ('config/track1/regulatory-index-pin.json',),
    }
    for stage in STAGES:
        selectors[stage].update(code=code(stage), limits=contract.execution_limits_sha256,
            runtime=contract.runtime_image_digest.removeprefix('sha256:'))
        for name, row in files.items():
            if any(name.startswith(prefix) for prefix in config_scopes.get(stage, ())):
                selectors[stage][name] = row['sha256']
    return selectors

class CurrentRun:
    """Update atomic stage completion records while the output lock is held."""
    def __init__(self, path, contracts):
        self.path = Path(path)
        self.root = self.path.parent.parent
        self.value = load_current(path)
        self.contracts = contracts
        self.plan = restart_plan(self.value, contracts, output_root=self.root)

    def save(self):
        self.value['contracts'] = self.contracts
        atomic_current(self.path, self.value)

    def invalidate(self):
        for stage in self.plan.invalid_stages:
            self.value['stages'].pop(stage, None)
        self.save()

    def complete(self, stage, *paths):
        self.value['stages'][stage] = completed_stage(stage, self.contracts[stage],
            self.value['stages'], output_root=self.root, paths=paths)
        self.save()

    def reusable(self, stage):
        return stage in self.plan.reused_stages

    def publication_identity(self):
        # Preserve every accepted result when any scientific dependency changes.
        return digest(canonical(self.contracts))

    def code_identity(self, stage):
        return digest(canonical(self.contracts[stage]))


def require_resume_archive_complete(root):
    """Do not continue after an interrupted checkpoint archive."""
    pending = Path(root) / '.operator/resume-archive-pending.json'
    _require(not pending.exists() and not pending.is_symlink(),
        'resume archive is incomplete; preserve the run and its recovery archive before continuing')


def preserve_previous_invocation(root):
    """Preserve the previous launcher evidence while holding its invocation lock."""
    root = Path(root)
    require_resume_archive_complete(root)
    metadata = root / '.operator'
    previous = metadata / 'current'
    store = metadata / 'recovery'
    _require(not any(p.is_symlink() for p in (previous, store, *previous.parents)))
    if not previous.exists():
        return
    _require(previous.is_dir())
    store.mkdir(mode=0o700, exist_ok=True)
    value = store.stat()
    _require(value.st_uid == os.getuid() and stat.S_IMODE(value.st_mode) == 0o700)
    archive = Path(tempfile.mkdtemp(prefix='invocation-', dir=store))
    record = dict(schema='genosapiens.invocation-evidence-archive', version=1,
        status='incomplete', source='.operator/current', destination='current')
    atomic_current(archive / 'archive.json', record)
    pending = metadata / 'resume-archive-pending.json'
    atomic_current(pending, dict(archive=str(archive.relative_to(root))))
    os.rename(previous, archive / 'current')
    for path in (archive, store, metadata):
        descriptor = os.open(path, os.O_DIRECTORY | os.O_RDONLY | os.O_NOFOLLOW)
        try:
            os.fsync(descriptor)
        finally:
            os.close(descriptor)
    record['status'] = 'complete'
    atomic_current(archive / 'archive.json', record)
    pending.unlink()
    descriptor = os.open(metadata, os.O_DIRECTORY | os.O_RDONLY | os.O_NOFOLLOW)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _archive_invalid_work(run, targets):
    """Move invalidated work without copying payloads or changing scientific identities.

    The invocation lock must be held. A pending marker blocks another resume if
    any move or manifest update fails. The archive is retained evidence, not a
    standalone runnable copy of the original output.
    """
    root = run.root
    require_resume_archive_complete(root)
    targets = list(dict.fromkeys(targets))
    for target in targets:
        _require(target.is_relative_to(root) and not any(p.is_symlink() for p in (target, *target.parents)),
            'disposable stage output has an unsafe path')
    targets = [target for target in targets if target.exists()]
    metadata = root / '.operator'
    store = metadata / 'recovery'
    _require(metadata.is_dir() and not any(p.is_symlink() for p in (store, *store.parents)))
    # Read only small control records. Never read checkpoint payloads to archive them.
    snapshots = []
    for path, limit in ((run.path, MAX_BYTES),
            (root / 'analysis/controlled-analysis-result.json', 16384),
            (root / 'analysis/failure-analysis.json', 16384),
            (metadata / 'last-error.json', MAX_BYTES)):
        _require(not any(p.is_symlink() for p in (path, *path.parents)))
        if path.exists():
            snapshots.append((path, read_bounded_bytes(path, maximum_bytes=limit)))
    store.mkdir(mode=0o700, exist_ok=True)
    metadata_mode = store.stat()
    _require(stat.S_ISDIR(metadata_mode.st_mode) and metadata_mode.st_uid == os.getuid()
        and stat.S_IMODE(metadata_mode.st_mode) == 0o700,
        'resume recovery directory must be private and owned by the run owner')
    archive = Path(tempfile.mkdtemp(prefix='resume-', dir=store))
    record = dict(schema='genosapiens.resume-work-archive', version=1,
        status='incomplete', restart_stage=run.plan.restart_stage,
        invalid_stages=list(run.plan.invalid_stages),
        moved=[dict(source=str(path.relative_to(root)), archive=f'work-{index:06d}')
               for index, path in enumerate(targets)],
        snapshots=[dict(source=str(path.relative_to(root)), archive=f'metadata-{index:06d}.json',
                        sha256=digest(content)) for index, (path, content) in enumerate(snapshots)])
    for item, (_, content) in zip(record['snapshots'], snapshots, strict=True):
        with (archive / item['archive']).open('xb') as stream:
            os.fchmod(stream.fileno(), 0o600)
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
    atomic_current(archive / 'archive.json', record)
    def sync(path):
        descriptor = os.open(path, os.O_DIRECTORY | os.O_RDONLY | os.O_NOFOLLOW)
        try:
            os.fsync(descriptor)
        finally:
            os.close(descriptor)
    sync(store)
    pending = metadata / 'resume-archive-pending.json'
    atomic_current(pending, dict(archive=str(archive.relative_to(root))))
    for path, item in zip(targets, record['moved'], strict=True):
        os.rename(path, archive / item['archive'])
        sync(path.parent)
        sync(archive)
    run.invalidate()
    record['status'] = 'complete'
    atomic_current(archive / 'archive.json', record)
    pending.unlink()
    sync(metadata)


def discard_invalid_work(run):
    """Archive invalidated products while keeping reusable upstream work in place."""
    require_resume_archive_complete(run.root)
    root = run.root
    invalid = set(run.plan.invalid_stages)
    work = root / 'analysis/analysis-work'
    _require(not any(p.is_symlink() for p in (work / 'state', work, *work.parents)),
        'disposable stage output has an unsafe path')
    targets = []
    preparation_changed = (run.value.get('contracts', {}).get('preparation') != run.contracts['preparation'])
    annotation_changed = (preparation_changed or run.value.get('contracts', {}).get('annotation') != run.contracts['annotation'])
    if 'preparation' in invalid and preparation_changed:
        targets.append(root / 'preparation')
        targets.append(root / '.operator/preparation-versions.json')
    if 'annotation' in invalid and annotation_changed:
        targets.append(work)
    else:
        if 'selection' in invalid:
            targets += [work / name for name in ('analysis/profile', 'eligibility')]
            state = work / 'state'
            if state.is_dir():
                targets += [p for p in state.iterdir() if not p.name.startswith(('annotation.sqlite', 'retained.sqlite'))]
        if 'associations' in invalid:
            targets.append(work / 'associations')
        if 'science' in invalid:
            targets += [work / name for name in ('scientific', 'order', 'primary', 'selection', 'product-checkpoint')]
            # Assessments are science outputs stored beside reusable discovery.
            # Archive the child only when its parent is not already invalidated.
            if 'associations' not in invalid:
                targets.append(work / 'associations/assessments')
            state = work / 'state'
            if state.is_dir():
                targets += [p for p in state.iterdir() if p.name.startswith(('products', 'context', 'assessment', 'bounded-products'))]
        if 'reports' in invalid:
            targets += [work / 'reports', work / 'product-checkpoint', work / 'state/bounded-products.sqlite']
    if invalid:
        targets.append(root / '.operator/final-check.json')
    if invalid:
        _archive_invalid_work(run, targets)
    else:
        run.invalidate()
