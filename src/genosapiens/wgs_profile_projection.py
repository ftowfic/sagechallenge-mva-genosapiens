"""Ordered profile decisions over an authenticated immutable admission parent.

Original producer and consumer identities remain separate. This owner persists
scope evidence, not scientific admission or a claim that excluded inputs are benign.
"""
from __future__ import annotations

from collections import Counter
from contextlib import ExitStack
from dataclasses import asdict, dataclass, replace
import hashlib
import json
import os
from pathlib import Path
import shutil
import sqlite3
import tempfile
import uuid

from genosapiens.track1_analysis_profiles import (
    AnalysisProfileSelection, ProfileCatalog, ProfileCandidateFacts,
    ProfileRouteEvidence, ProfileRouteFacts, ROUTES, evaluate_profile_candidate,
)
from genosapiens.wgs_admission_projection import AuthenticatedAdmissionOwner
from genosapiens.wgs_ledger import canonical, MAX_ROW_BYTES, MAX_SHARD_BYTES, MAX_SHARD_ROWS
from genosapiens.wgs_shards import ScientificShard, _ShardWriter, encode_row, read_bounded_metadata
from genosapiens.wgs_read_cache import ReadFileCache
from genosapiens.variant_retention import POLICY_SHA256 as RETENTION_POLICY_SHA256, POLICY_VERSION

SCHEMA = 'genosapiens.wgs-profile-projection'
DISPOSITIONS = ('selected', 'unresolved', 'outside_profile', 'not_applicable_upstream')
UPSTREAM = ('assess', 'exclude_common', 'defer_refcall')
RETENTION_DISPOSITIONS = ('retained', 'deferred_no_retention_evidence', 'not_applicable_upstream')


class ProfileProjectionError(ValueError):
    """Static incomplete, foreign, corrupt or publication failures."""


def _require(value):
    if not value:
        raise ProfileProjectionError('profile projection binding or closure is invalid')


def _sha(content):
    return hashlib.sha256(content).hexdigest()


def profile_adapter_identity():
    """Bind selection to its direct adapters without a package-wide identity."""
    root = Path(__file__).parent
    names = ('wgs_profile_projection.py', 'track1_analysis_profiles.py', 'variant_retention.py',
             'retention_sources.py', 'curated_alleles.py', 'wgs_clinical_session.py',
             'common_allele_public_indexes.py')
    return _sha(canonical({name: _sha((root / name).read_bytes()) for name in names}))


def _digest(value):
    _require(type(value) is str and len(value) == 64 and all(c in '0123456789abcdef' for c in value))
    return value


def _unique(pairs):
    result = {}
    for key, value in pairs:
        _require(key not in result)
        result[key] = value
    return result


def _json(content):
    value = json.loads(content, object_pairs_hook=_unique)
    _require(canonical(value) == content)
    return value


def _signature(path):
    path = Path(path)
    _require(not any(p.is_symlink() for p in (path, *path.parents)) and path.is_file())
    value = path.stat()
    return value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns


def _durable(path):
    descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


@dataclass(frozen=True, slots=True)
class PublicProfilePilotCohort:
    """Complete declared finite domain; never a full-parent qualification claim."""
    identity_bytes: bytes

    @property
    def identity(self):
        return _json(self.identity_bytes)

    @property
    def ordinals(self):
        return tuple(self.identity['original_ordinals'])


@dataclass(frozen=True, slots=True)
class CompletedProfileDiscoveryInputs:
    """Exact public restart inputs, separately bound by the invocation authority."""
    identity_bytes: bytes
    plan_sha256: str

    def __post_init__(self):
        _digest(self.plan_sha256)
        _require(type(self.identity_bytes) is bytes and len(self.identity_bytes) <= 32768)
        value = self.identity
        optional = {'profile_producer_executable_stack_sha256'} & set(value)
        _require(set(value) == optional | ({'restart_stage'} if self.profile_only else set()) | {'producer_candidate_sha256', 'producer_executable_stack_sha256',
            'consumer_executable_stack_sha256', 'profile_manifest', 'association_manifest',
            'profile_adapter_sha256', 'discovery_binding_sha256', 'acceptance_amendment'})
        for key, item in value.items():
            if key == 'restart_stage':
                _require(item == 'association_discovery')
            elif key.endswith('_sha256'):
                _digest(item)
            else:
                _require(type(item) is dict and set(item) == {'path', 'sha256', 'size_bytes'})
                _digest(item['sha256'])
                path = Path(item['path'])
                _require(path.is_absolute() and '..' not in path.parts
                    and type(item['size_bytes']) is int and item['size_bytes'] > 0)
        profile = Path(value['profile_manifest']['path'])
        association = Path(value['association_manifest']['path'])
        _require(profile.name == association.name == 'manifest.json'
            and profile.parent.name == 'profile' and profile.parent.parent.name == 'analysis'
            and association.parent.name == 'associations'
            and profile.parent.parent.parent == association.parent.parent)

    @property
    def profile_producer_executable_sha256(self):
        # A later discovery producer can reuse an older byte-exact profile.
        value = self.identity
        return value.get('profile_producer_executable_stack_sha256',
                         value['producer_executable_stack_sha256'])

    @property
    def profile_only(self):
        return self.identity.get('restart_stage') == 'association_discovery'

    @property
    def identity(self):
        return _json(self.identity_bytes)

    @property
    def sha256(self):
        return _sha(canonical(dict(inputs=self.identity, plan_sha256=self.plan_sha256)))

    def require_authority(self, authority, *, consumer_executable_sha256):
        value = self.identity
        _require(authority.is_finite_profile_pilot
            and authority.binding('public_profile_pilot_context_inputs') == self.sha256
            and authority.plan_sha256 == self.plan_sha256
            and authority.executable_stack_sha256 == consumer_executable_sha256
            and value['consumer_executable_stack_sha256'] == consumer_executable_sha256
            and authority.authorization_receipt_sha256 == value['acceptance_amendment']['sha256'])
        root = Path('/mnt/reference/scratch/public-build') / (
            'track1-bounded-wgs-prefilter-replacement-' + value['producer_candidate_sha256'])
        _require(Path(value['profile_manifest']['path']).parent.parent.parent
            == root/'pilot-01/analysis/analysis-work')
        for key in (('profile_manifest',) if self.profile_only else ('profile_manifest', 'association_manifest')):
            path = Path(value[key]['path'])
            _signature(path)
            _require(bool(os.statvfs(path).f_flag & os.ST_RDONLY))
        from genosapiens.track1_roster import validate_context_restart_acceptance
        validate_context_restart_acceptance(self)


@dataclass(frozen=True, slots=True)
class PublicProfilePilotInputs:
    """Public candidate-harness inputs; ordinary/operator callers expose none."""
    admission_owner: AuthenticatedAdmissionOwner
    selection_bytes: bytes
    dense_bytes: bytes
    reconciled_directory: Path
    completed_context_inputs: CompletedProfileDiscoveryInputs | None = None
    completed_assessment_inputs: object | None = None
    completed_scientific_inputs: object | None = None

    @property
    def binding_sha256(self):
        return _sha(canonical(dict(selection_sha256=_sha(self.selection_bytes),
            dense_selection_sha256=_sha(self.dense_bytes),
            admission_completion_sha256=self.admission_owner.completion.sha256)))


def load_public_profile_pilot_ordinals(repository, *, selection_bytes, dense_bytes):
    """Validate the pinned cohort union before opening the admission owner."""
    from genosapiens.track1_roster import load_analysis_profile_acceptance
    declaration = load_analysis_profile_acceptance(repository)
    pilot = declaration['pilot']
    _require(_sha(selection_bytes) == pilot['selection']['sha256']
        and _sha(dense_bytes) == pilot['dense_selection']['sha256'])
    random = tuple(json.loads(selection_bytes)['random_ordinals'])
    dense = tuple(row['ordinal'] for row in json.loads(dense_bytes)['samples'])
    _require(len(random) == len(set(random)) == pilot['strata'][0]['sample_size']
        and len(dense) == len(set(dense)) == pilot['strata'][1]['sample_size']
        and all(type(ordinal) is int and ordinal > 0 for ordinal in (*random, *dense)))
    # The original random and dense cohorts share 50 inputs. Each original
    # ordinal participates once in the scientific union, preserving both lists.
    ordinals = tuple(sorted(set(random) | set(dense)))
    _require(len(ordinals) == 36286)
    return ordinals


def authenticate_public_profile_pilot_cohort(repository, admission_owner, *, selection_bytes, dense_bytes):
    """Bind the original frozen pilot to an authentic complete public parent.

    This has no CLI flag and cannot redefine the cohort from ranking results.
    The whole admission owner remains unchanged; only its point-read domain is
    finite. Original ordinals, occurrence IDs and producer bindings survive.
    """
    from genosapiens.track1_roster import PROFILE_ACCEPTANCE_SHA256
    _require(type(admission_owner) is AuthenticatedAdmissionOwner
        and Path(admission_owner.path).is_relative_to('/mnt/reference/scratch/public-build'))
    ordinals = load_public_profile_pilot_ordinals(repository,
        selection_bytes=selection_bytes, dense_bytes=dense_bytes)
    admission_owner._assert_stable()
    identity = dict(schema='genosapiens.track1-public-profile-pilot-cohort', version=2,
        purpose='finite-scientific-output-pilot', classification='public',
        full_parent_qualification=False, acceptance_sha256=PROFILE_ACCEPTANCE_SHA256,
        selection_sha256=_sha(selection_bytes), dense_selection_sha256=_sha(dense_bytes),
        original_ordinals=ordinals, ordered_ordinals_sha256=_sha(canonical(ordinals)),
        rows=len(ordinals),
        binding_sha256=admission_owner.binding_sha256,
        completion_sha256=admission_owner.completion.sha256)
    return PublicProfilePilotCohort(canonical(identity))


def _cohort(identity):
    cohort = identity.get('finite_cohort')
    if cohort is None:
        return None
    _require(type(cohort) is dict and set(cohort) == {'schema','version','purpose','classification',
        'full_parent_qualification','acceptance_sha256','selection_sha256','dense_selection_sha256',
        'original_ordinals','ordered_ordinals_sha256','rows','binding_sha256','completion_sha256'}
        and cohort['schema'] == 'genosapiens.track1-public-profile-pilot-cohort'
        and cohort['version'] == 2 and cohort['purpose'] == 'finite-scientific-output-pilot'
        and cohort['classification'] == 'public' and cohort['full_parent_qualification'] is False
        and cohort['binding_sha256'] == identity['parent']['binding_sha256']
        and cohort['completion_sha256'] == identity['parent']['completion_sha256'])
    values = cohort['original_ordinals']
    _require(type(values) is list and 0 < len(values) <= 50000
        and len(values) == cohort['rows'] and values == sorted(set(values))
        and all(type(v) is int and 1 <= v <= identity['parent']['rows'] for v in values)
        and _sha(canonical(values)) == cohort['ordered_ordinals_sha256'])
    for key in ('acceptance_sha256','selection_sha256','dense_selection_sha256','ordered_ordinals_sha256'):
        _digest(cohort[key])
    return cohort


def _compact_domain(owner, identity):
    cohort = _cohort(identity)
    return owner.iter_screening_records(original_ordinals=None if cohort is None
        else tuple(cohort['original_ordinals']))


class ReconciledPublicProfilePilotFacts:
    """Bounded reuse of a complete reconciled original-cohort experiment.

    Every saved fact is checked against the live original compact owner and
    current policy; complete ledger authentication finishes before discovery.
    This cannot supply a partial prefix or authorize a full-source profile.
    """
    def __init__(self, directory, *, repository, admission_owner, pilot_cohort,
            selection, catalog, route_sources, regulatory_index_sha256):
        import ast
        from genosapiens.track1_roster import load_analysis_profile_acceptance, PROFILE_ACCEPTANCE_SHA256
        _require(type(pilot_cohort) is PublicProfilePilotCohort
            and type(admission_owner) is AuthenticatedAdmissionOwner)
        self.root = Path(directory)
        self.owner, self.cohort = admission_owner, pilot_cohort
        self.ordinals = pilot_cohort.ordinals
        self.selection, self.catalog = selection, catalog
        self.route_sources = dict(route_sources)
        self.index_sha256 = regulatory_index_sha256
        declaration = load_analysis_profile_acceptance(repository)
        _require(pilot_cohort.identity['acceptance_sha256'] == PROFILE_ACCEPTANCE_SHA256)
        expected = {Path(row['path']).name: row for row in declaration['prospective_amendment']['selection_receipts']}
        self.signatures = {}
        values = {}
        for name in ('plan.json','result.json','exit.json','reconciliation.json'):
            path = self.root/name
            signature = _signature(path)
            content = read_bounded_metadata(path)
            _require((_sha(content), len(content)) == (expected[name]['sha256'], expected[name]['size_bytes']))
            self.signatures[path] = signature
            values[name] = json.loads(content)
        result, terminal, reconciliation = (values[n] for n in ('result.json','exit.json','reconciliation.json'))
        plan = values['plan.json']
        descriptor_path = Path(plan['admission_descriptor'])
        _require(descriptor_path.is_relative_to('/mnt/reference/scratch/public-build'))
        descriptor_signature = _signature(descriptor_path)
        descriptor_bytes = read_bounded_metadata(descriptor_path)
        _require(_sha(descriptor_bytes) == plan['admission_descriptor_sha256'])
        descriptor = json.loads(descriptor_bytes)
        _require(descriptor['binding_sha256'] == admission_owner.binding_sha256
            and descriptor['database']['path'] == str(admission_owner.path))
        self.signatures[descriptor_path] = descriptor_signature
        _require(result['status'] == 'complete' and result['unique_inputs'] == len(self.ordinals)
            and terminal['returncode'] == 0 and terminal['timed_out'] is False
            and terminal['docker_state']['OOMKilled'] is False and reconciliation['status'] == 'passed'
            and reconciliation['unique_inputs'] == result['unique_inputs']
            and reconciliation['decisions_sha256'] == result['decisions_sha256']
            and result['policy_sha256'] == selection.policy.sha256
            and result['catalog_sha256'] == catalog.sha256
            and result['index_manifest_sha256'] == regulatory_index_sha256
            and plan['selection_sha256'] == pilot_cohort.identity['selection_sha256']
            and plan['dense_sample_inventory_sha256'] == pilot_cohort.identity['dense_selection_sha256'])
        # Exact unchanged source owners plus explicit unchanged functions in the
        # extended projection owner. This is not general compatibility dispatch.
        for relative in ('src/genosapiens/track1_analysis_profiles.py', 'src/genosapiens/admission_links.py',
                'src/genosapiens/wgs_clinical_session.py',
                'src/genosapiens/annotation_admission.py', 'src/genosapiens/wgs_admission_screening.py'):
            _require(_sha(read_bounded_metadata(Path(repository)/relative)) == plan['source_sha256s'][relative])
        relative = 'src/genosapiens/wgs_profile_projection.py'
        producer_bytes = read_bounded_metadata(self.root/'frozen-source'/relative)
        _require(_sha(producer_bytes) == plan['source_sha256s'][relative])
        producer, consumer = (ast.parse(content) for content in
            (producer_bytes, read_bounded_metadata(Path(repository)/relative)))
        for name in ('_facts','derive_profile_candidate_facts','profile_route_sources'):
            before = next(node for node in producer.body if isinstance(node, ast.FunctionDef) and node.name == name)
            after = next(node for node in consumer.body if isinstance(node, ast.FunctionDef) and node.name == name)
            _require(ast.dump(before, include_attributes=False) == ast.dump(after, include_attributes=False))
        # The revised owners reuse completion instead of replaying payloads.
        # Saved scientific facts remain compatible only with the exact codecs,
        # query semantics and policy definitions of their frozen producer.
        for relative, changed_definitions in (
                ('src/genosapiens/wgs_admission_projection.py', {'AuthenticatedAdmissionOwner'}),
                ('src/genosapiens/regulatory_index.py', {'open_regulatory_index'}),
                ('src/genosapiens/wgs_analysis_eligibility.py', {'scientific_assessment_domain'})):
            producer_bytes = read_bounded_metadata(self.root/'frozen-source'/relative)
            _require(_sha(producer_bytes) == plan['source_sha256s'][relative])
            trees = [ast.parse(content) for content in
                (producer_bytes, read_bounded_metadata(Path(repository)/relative))]
            for tree in trees:
                _require(changed_definitions <= {getattr(node, 'name', None) for node in tree.body})
                tree.body = [node for node in tree.body
                    if not isinstance(node, (ast.FunctionDef, ast.ClassDef))
                    or node.name not in changed_definitions]
            _require(ast.dump(trees[0], include_attributes=False) == ast.dump(trees[1], include_attributes=False))
        self.expected_digest = result['decisions_sha256']
        self.path = self.root/'decisions.jsonl'
        self.signatures[self.path] = _signature(self.path)
        self.stream = self.path.open('rb')
        self.cache = ReadFileCache(self.stream.fileno())
        self.digest = hashlib.sha256()
        self.count = 0
        self.finished = False

    def _next(self):
        raw = self.stream.readline(MAX_ROW_BYTES+1)
        _require(raw and len(raw) <= MAX_ROW_BYTES and raw.endswith(b'\n')
            and self.count < len(self.ordinals))
        self.digest.update(raw)
        row = _json(raw)
        ordinal = self.ordinals[self.count]
        _require(row['ordinal'] == ordinal)
        compact = self.owner._compact(ordinal)
        _require(row['identity_id'] == compact.ref[2]
            and row['occurrences'] == len(compact.facts.source_allele_ids)
            and row['upstream_disposition'] == compact.decision.disposition)
        self.count += 1
        self.cache.release()
        if compact.decision.disposition != 'assess':
            _require(row['disposition'] == 'not_applicable_upstream'
                and 'facts' not in row and 'decision' not in row)
            return compact, None
        facts = _facts(row['facts'])
        _require((facts.identity_id, facts.original_input_sha256, facts.upstream_disposition) ==
            (compact.ref[2], compact.ref[4], compact.decision.disposition)
            and all(r.source_sha256 == self.route_sources[r.route] for r in facts.routes))
        decision = evaluate_profile_candidate(facts, policy=self.selection.policy,
            catalog=self.catalog, profile=self.selection.profile)
        _require(canonical(asdict(decision)) == canonical(row['decision'])
            and row['disposition'] == decision.disposition)
        return compact, facts

    def __call__(self, compact):
        while True:
            original, facts = self._next()
            _require(original.ref[1] <= compact.ref[1])
            if original.ref == compact.ref:
                _require(facts is not None)
                return facts
            _require(facts is None)

    def finish(self):
        while self.count < len(self.ordinals):
            _, facts = self._next()
            _require(facts is None)
        _require(not self.stream.read(1) and self.digest.hexdigest() == self.expected_digest
            and all(_signature(path) == signature for path, signature in self.signatures.items()))
        self.owner._assert_stable()
        self.finished = True

    def close(self):
        try:
            self.cache.release(force=True)
        finally:
            self.stream.close()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()


def _binding(owner, selection, catalog, regulatory_index_sha256, consumer_executable_sha256, route_sources,
        *, pilot_cohort=None, completed_context_inputs=None):
    _require(type(owner) is AuthenticatedAdmissionOwner
             and type(selection) is AnalysisProfileSelection and type(catalog) is ProfileCatalog
             and catalog.policy_sha256 == selection.policy.sha256
             and not catalog.coverage_issues and type(route_sources) is dict
             and set(route_sources) == set(ROUTES))
    owner._assert_stable()
    identity = dict(schema=SCHEMA, version=3, classification='controlled',llm_policy='no-llm',
        profile=selection.profile,
        policy_sha256=selection.policy.sha256, catalog_sha256=catalog.sha256,
        regulatory_index_sha256=_digest(regulatory_index_sha256),
        consumer_executable_sha256=_digest(consumer_executable_sha256),
        adapter_sha256=profile_adapter_identity(),
        route_sources={key: _digest(route_sources[key]) for key in ROUTES},
        parent=dict(binding_sha256=owner.binding_sha256, completion_sha256=owner.completion.sha256,
            producer_executable_sha256=owner.parents['consumer_executable_sha256'],
            source_manifest_sha256=owner.parents['source_manifest_sha256'],
            reference_sha256=owner.eligibility_manifest.reference_manifest_sha256,
            input_rows_sha256=owner.input_rows_sha256, rows=owner.row_count,
            occurrences=owner.occurrence_count,
            input_dispositions=dict(owner.eligibility_manifest.input_dispositions),
            occurrence_dispositions=dict(owner.eligibility_manifest.occurrence_dispositions)))
    if pilot_cohort is not None:
        _require(type(pilot_cohort) is PublicProfilePilotCohort)
        identity['finite_cohort'] = pilot_cohort.identity
        _cohort(identity)
    if completed_context_inputs is not None:
        from genosapiens.track1_roster import validate_context_restart_acceptance
        _require(type(completed_context_inputs) is CompletedProfileDiscoveryInputs
            and pilot_cohort is not None)
        original = completed_context_inputs.identity
        compatibility = validate_context_restart_acceptance(completed_context_inputs)
        _require(original['consumer_executable_stack_sha256'] == consumer_executable_sha256)
        identity['consumer_executable_sha256'] = completed_context_inputs.profile_producer_executable_sha256
        identity['adapter_sha256'] = original['profile_adapter_sha256']
        identity['finite_cohort']['acceptance_sha256'] = compatibility['producer_profile_acceptance_sha256']
    return identity


def load_installed_profile_resources(repository, references, *, gene_store, reference):
    """Load required pinned resources; runtime never builds or guesses an index.

    The pin is published only after installed-source qualification. Missing or
    mismatching resources stop analysis before admission or case discovery.
    Preparation retains its separate, independently reusable resource owner.
    """
    from genosapiens.admission_links import PinnedPublicSource
    from genosapiens.regulatory_index import open_regulatory_index
    from genosapiens.track1_analysis_profiles import authenticate_analysis_profile,build_profile_catalog_from_gene_evidence
    index = None
    try:
        repository,references=Path(repository),Path(references)
        pin_path=repository/'config/track1/regulatory-index-pin.json'
        signature=_signature(pin_path)
        # This configuration pin is published as ordinary JSON by the host
        # driver. Canonical bytes remain required for profile/index artifacts;
        # here enforce unique keys and the exact typed fields below instead.
        pin=json.loads(read_bounded_metadata(pin_path),object_pairs_hook=_unique)
        _require(set(pin)=={'schema','version','index_path','manifest_sha256','source_contract_sha256',
            'source_manifest_sha256s','semantic_acceptance_sha256','reference_contract_sha256','reference_manifest_sha256'}
            and pin['schema']=='genosapiens.track1-regulatory-index-pin'
            and type(pin['version']) is int and pin['version']==1
            and pin['reference_contract_sha256']==reference.contract_sha256
            and pin['reference_manifest_sha256']==reference.manifest_sha256)
        relative=Path(pin['index_path'])
        _require(not relative.is_absolute() and len(relative.parts)==3
            and relative.parts[:2]==('derived','track1-regulatory-index'))
        _digest(relative.parts[-1])
        for key in ('manifest_sha256','source_contract_sha256','semantic_acceptance_sha256',
                    'reference_contract_sha256','reference_manifest_sha256'):
            _digest(pin[key])
        source_ids=('encode-registry-v4-ccre-grch38','encode-re2g-registry-v4-grch38')
        _require(type(pin['source_manifest_sha256s']) is dict
            and set(pin['source_manifest_sha256s'])==set(source_ids))
        contract=repository/'config/public-resources/variant-class-resources.json'
        _require(_sha(read_bounded_metadata(contract))==pin['source_contract_sha256'])
        sources=[]
        for name in source_ids:
            digest=_digest(pin['source_manifest_sha256s'][name])
            path=references/'manifests/variant-class-resources'/(name+'.json')
            _require(_sha(read_bounded_metadata(path))==digest)
            sources.append(PinnedPublicSource(contract,pin['source_contract_sha256'],references,path,digest))
        root=references/relative
        _require(not any(path.is_symlink() for path in (root,*root.parents)))
        selection=authenticate_analysis_profile()
        catalog=build_profile_catalog_from_gene_evidence(gene_store,policy=selection.policy)
        _require(not catalog.coverage_issues)
        index=open_regulatory_index(root,expected_manifest_sha256=pin['manifest_sha256'])
        _require(index.manifest['semantic_acceptance_sha256']==pin['semantic_acceptance_sha256']
            and index.manifest['reference']['reference_contract_sha256']==reference.contract_sha256
            and index.manifest['reference']['reference_manifest_sha256']==reference.manifest_sha256
            and sorted((row['resource_id'],row['contract_sha256'],row['manifest_sha256'])
                for row in index.manifest['sources'])==sorted((name,pin['source_contract_sha256'],
                pin['source_manifest_sha256s'][name]) for name in source_ids)
            and _signature(pin_path)==signature)
        return catalog,index,tuple(sources)
    except (Exception,KeyboardInterrupt):
        if index is not None:
            index.close()
        raise ProfileProjectionError('required installed profile resources are missing, stale or corrupt') from None


@dataclass(frozen=True, slots=True)
class ProfileDecisionRow:
    original_ref: tuple
    source_allele_ids: tuple[str, ...]
    screening_sha256: str
    upstream_disposition: str
    facts: dict
    decision: dict
    deferred_inventory: dict | None = None


def _deferred_inventory(row):
    """Compact, source-addressable inventory in the existing ordered shards."""
    if row.decision['retention_disposition'] != 'deferred_no_retention_evidence':
        return None
    terms = sorted({term for route in row.facts['routes'] for evidence in route['evidence']
                    for term in evidence['consequences']})
    return dict(schema='genosapiens.deferred-allele-inventory', version=1,
        input_ordinal=row.original_ref[1], identity_id=row.original_ref[2], allele_id=row.original_ref[3],
        sample_ids=row.facts['sample_ids'], source_occurrences=row.source_allele_ids,
        original_input_ref=row.original_ref, exact_allele=row.facts['exact_allele'],
        annotation_summary=dict(consequences=terms, transcript_records=sum(
            bool(e['transcript_id']) for route in row.facts['routes'] for e in route['evidence'])),
        reason='deferred_no_retention_evidence', policy_version=POLICY_VERSION,
        policy_sha256=row.decision['retention_policy_sha256'],
        profile_policy_sha256=row.decision['policy_sha256'], catalog_sha256=row.decision['catalog_sha256'],
        route_sources={r['route']:r['source_sha256'] for r in row.facts['routes']},
        public_coverage=row.facts['public_coverage'])


def profile_association_applicability(row, associations):
    """Project the sealed gate's target routes, without scientific reinterpretation."""
    from genosapiens import candidate_admission as scientific
    _require(type(row) is ProfileDecisionRow and row.decision['disposition'] in {'selected','unresolved'}
             and type(associations) is tuple
             and all(type(value) is scientific._Association for value in associations))
    subjects=(row.original_ref[3],) if row.original_ref[3] is not None else row.source_allele_ids
    eligible=tuple(route for route in row.decision['routes'] if route['outcome'] in {'selected','unresolved'})
    native=any(route['route']=='native' for route in eligible)
    result=[]
    for value in associations:
        _require(value.subject in subjects)
        if value.state!='match':
            result.append(True)
            continue
        if value.adapter=='annotation':
            if native:
                result.append(True)
                continue
            route='molecular'
        elif value.adapter=='encode':route='regulatory'
        elif value.adapter in {'open_targets','clinvar','curated'}:route='clinical'
        elif value.adapter=='repeat':route='native'
        else:raise ProfileProjectionError('profile association adapter is invalid')
        target=(value.target.split('.')[0] if value.namespace=='Ensembl' and value.target else value.target)
        result.append(any(candidate['route']==route and (
            candidate['target_id'] is not None
            and (candidate['target_namespace'],candidate['target_id'])==(value.namespace,target)
            or candidate['target_id'] is None and route=='native'
            or candidate['target_id'] is None and route=='clinical'
        ) for candidate in eligible))
    return tuple(result)


def _facts(document):
    from genosapiens.variant_retention import ExactAllele, FunctionalFact, SourceFact, LookupCoverage
    def evidence(e):
        e = dict(e)
        for key in ('consequences', 'clinical_classifications', 'clinical_conditions', 'clinical_conflicts'):
            e[key] = tuple(e[key])
        if e['exact_allele'] is not None:
            e['exact_allele'] = ExactAllele(**e['exact_allele'])
        if e['functional'] is not None:
            f = dict(e['functional'])
            f['source'] = SourceFact(**f['source'])
            if f['allele'] is not None:
                f['allele'] = ExactAllele(**f['allele'])
            for key in ('disease_ids', 'extraction_limits', 'conflicts'):
                f[key] = tuple(f[key])
            e['functional'] = FunctionalFact(**f)
        return ProfileRouteEvidence(**e)
    value = dict(document)
    value['routes'] = tuple(ProfileRouteFacts(**{**r,
        'reasons': tuple(r['reasons']), 'source_refs': tuple(r['source_refs']), 'evidence': tuple(
        evidence(e) for e in r['evidence'])}) for r in value['routes'])
    value['technical_reasons'] = tuple(value['technical_reasons'])
    value['source_occurrences'] = tuple(value['source_occurrences'])
    value['public_coverage'] = tuple(LookupCoverage(**v) for v in value['public_coverage'])
    value['sample_ids'] = tuple(value['sample_ids'])
    if value['exact_allele'] is not None:
        value['exact_allele'] = ExactAllele(**value['exact_allele'])
    return ProfileCandidateFacts(**value)


@dataclass(frozen=True, slots=True, repr=False)
class _ValidatedNativeAdmission:
    row: object
    intake: object


class FreshAdmissionProfileIntake:
    """Borrow live source owners; persist bounded gate facts before rich admission.

    The disposable SQLite file is a staging owner, never a published profile or
    recovery authority. The ordinary profile writer independently reconciles
    its facts against the sealed original admission and eligibility ledgers.
    """

    def __init__(self, path, *, selection, catalog, clinical_session, regulatory_index, reference,
                 annotation_manifest_sha256, native_source_sha256, execution_config):
        from genosapiens.wgs_admission_projection import _settings
        _require(type(selection) is AnalysisProfileSelection and type(catalog) is ProfileCatalog
                 and catalog.policy_sha256 == selection.policy.sha256 and not catalog.coverage_issues)
        self.selection,self.catalog=selection,catalog
        self.clinical,self.regulatory,self.reference=clinical_session,regulatory_index,reference
        self.route_sources=profile_route_sources(annotation_manifest_sha256=annotation_manifest_sha256,
            native_source_sha256=native_source_sha256,clinical_session=clinical_session,
            regulatory_index=regulatory_index)
        self._route_owner = _ProfileRouteContext(clinical_session, regulatory_index, self.route_sources)
        self.identity=dict(profile=selection.profile,policy_sha256=selection.policy.sha256,
            catalog_sha256=catalog.sha256,regulatory_index_sha256=regulatory_index.manifest_sha256,
            route_sources=self.route_sources,adapter_sha256=profile_adapter_identity())
        self._identity_snapshot = _checked_metadata(self.identity)
        self.path=Path(path)
        _require(not self.path.exists() and not any(p.is_symlink() for p in (self.path,*self.path.parents)))
        self.descriptor=os.open(self.path,os.O_CREAT|os.O_EXCL|os.O_NOFOLLOW|os.O_RDWR,0o600)
        try:
            self.db=sqlite3.connect(self.path)
            _settings(self.db,execution_config,readonly=False)
            self.db.execute('CREATE TABLE facts(ordinal INTEGER PRIMARY KEY, identity TEXT UNIQUE, input_sha TEXT, content BLOB, digest TEXT)')
            self.db.execute('CREATE TABLE binding(content BLOB)')
        except BaseException:
            if hasattr(self,'db'):self.db.close()
            os.close(self.descriptor);self.path.unlink();raise
        self.binding=self.screening=None
        self.rows=self.occurrences=self.read_rows=0
        self._last_read_ordinal=0
        self.upstream=Counter()
        self.finished=self.closed=False
        self.metrics=dict(profile_inputs=0,sequence_proof_calls=0,native_proof_calls=0,
            rich_inputs=0,outside_inputs=0,upstream_inputs=0,facts_bytes=0)
        self.signature=_signature(self.path)

    def _check(self):
        _require(not self.closed and _signature(self.path)==self.signature)
        if self.binding is not None:
            metadata=self.db.execute('SELECT content FROM binding').fetchall()
            _require(self.binding == self._record_context.binding and metadata == [(self._binding_content,)])
        _require(_metadata_unchanged(self.identity, self._identity_snapshot))
        self._route_owner.require(self.route_sources, self.clinical, self.regulatory,
            self.reference, self.route_sources['native'])

    def bind(self, binding, screening):
        from genosapiens.wgs_admission_projection import _ProfileRecordContext, _row_parents
        _require(not self.closed and self.binding is None and not self.finished)
        self._check()
        context = _ProfileRecordContext(binding, screening)
        parents = context.parents
        _require(self.route_sources['molecular']==parents['annotation_stage_sha256']
                 and self.route_sources['native']==_row_parents(parents)[2])
        self.binding,self.screening=binding,screening
        self.parents=parents
        self._record_context = context
        self._binding_content = canonical(dict(self.identity, admission_binding_sha256=context.binding_sha256))
        self.db.execute('INSERT INTO binding VALUES(?)', (self._binding_content,))
        self.db.commit()
        self.signature=_signature(self.path)

    def prepare_native(self, pending, session):
        from genosapiens.wgs_admission_projection import PendingNativeAdmission, EncodedAdmissionRow
        _require(not self.finished and not self.closed and self.binding is not None
                 and type(pending) is PendingNativeAdmission)
        self._check()
        enriched,intake=session.enrich_intake(pending.decode(self.binding))
        return _ValidatedNativeAdmission(EncodedAdmissionRow.encode(enriched,(),
            binding=self.binding,screening=self.screening),intake)

    def finalize_window(self, values, batch):
        from genosapiens.wgs_admission_diagnostics import boundary
        from genosapiens.wgs_admission_projection import (EncodedAdmissionRow, ValidatedSequenceAdmission,
            _finish_profile_sequence_admission, _owned_admission_proofs_bytes, _pack_admission_blob)
        from genosapiens.wgs_admission import sequence_annotation_intake_provenance
        from genosapiens.wgs_analysis_eligibility import EligibilityBatch, evaluate_profile_intake
        with boundary('profile_window_validate'):
            _require(not self.finished and not self.closed and self.binding is not None
                     and type(values) is tuple and 0<len(values)<=64 and type(batch) is EligibilityBatch
                     and len(values)==len(batch.decisions))
            self._check()
        from genosapiens.wgs_admission_projection import _ProfileRecord
        rows,decisions=[],[]
        with ExitStack() as scope:
            decoded, alleles = [], []
            with boundary('profile_window_decode'):
                for value, upstream in zip(values, batch.decisions, strict=True):
                    intake = value.row if type(value) is _ValidatedNativeAdmission else value
                    record = scope.enter_context(_ProfileRecord(intake, self._record_context))
                    decoded.append(record)
                    if upstream.disposition == 'assess' and record.value.allele is not None:
                        alleles.append(record.value.allele)
            with boundary('profile_clinical_query', blocking=True, reason='lookup_failed'):
                clinical_window = scope.enter_context(self.clinical._retention_window(tuple(alleles)))
            for value,upstream,record in zip(values,batch.decisions,decoded,strict=True):
                enriched, verified = record.value, record.verified
                with boundary('profile_row_convert'):
                    arguments=dict(upstream_decision=upstream,selection=self.selection,catalog=self.catalog,
                        route_sources=self.route_sources,clinical_session=self.clinical,
                        regulatory_index=self.regulatory,reference=self.reference,
                        _clinical_window=clinical_window, _record_owner=record, _route_owner=self._route_owner)
                    if type(value) is ValidatedSequenceAdmission:
                        provenance = sequence_annotation_intake_provenance(verified,
                            annotation_stage_manifest_sha256=self.parents['annotation_stage_sha256'],
                            bundle_manifest_sha256=self.parents['bundle_manifest_sha256'])
                        row,upstream,facts,decision=_finish_profile_sequence_admission(enriched,verified,provenance,
                            binding=self.binding,screening=self.screening,original_row=value.row,**arguments)
                        self.metrics['sequence_proof_calls']+=int(decision.disposition in ('selected','unresolved'))
                    else:
                        native=type(value) is _ValidatedNativeAdmission
                        row=value.row if native else value
                        _require(type(row) is EncodedAdmissionRow and row.proof_count==0)
                        facts,decision=evaluate_profile_intake(enriched,molecular_sources=(),**arguments)
                        if native and decision.disposition in ('selected','unresolved'):
                            proofs=value.intake.proofs()
                            content=_owned_admission_proofs_bytes(proofs)
                            row=replace(row,proof_blob=_pack_admission_blob(content),proof_size=len(content),
                                proof_sha256=_sha(content),proof_count=len(proofs))
                            self.metrics['native_proof_calls']+=1
                        upstream=replace(upstream,proof_sha256=row.proof_sha256)
                    _require(upstream.retained_ordinal==self.rows+1 and row.input_sha256==upstream.input_sha256)
                    content=canonical(asdict(facts));_require(len(content)<=MAX_ROW_BYTES)
                with boundary('profile_row_persist'):
                    self.db.execute('INSERT INTO facts VALUES(?,?,?,?,?)',
                        (upstream.retained_ordinal,upstream.identity_id,upstream.input_sha256,content,_sha(content)))
                self.rows+=1;self.occurrences+=len(upstream.source_allele_ids)
                self.upstream[upstream.disposition]+=1
                self.metrics['profile_inputs']+=1;self.metrics['facts_bytes']+=len(content)
                self.metrics['rich_inputs']+=int(decision.disposition in ('selected','unresolved'))
                self.metrics['outside_inputs']+=int(decision.disposition=='outside_profile')
                self.metrics['upstream_inputs']+=int(decision.disposition=='not_applicable_upstream')
                rows.append(row);decisions.append(upstream)
            with boundary('profile_row_persist', blocking=True):
                self.db.commit()
                self.signature=_signature(self.path)
            return tuple(rows),EligibilityBatch(tuple(decisions))

    def finish_intake(self, *, expected_inputs, expected_occurrences, input_dispositions):
        _require(not self.closed and self.binding is not None and not self.finished
                 and _signature(self.path)==self.signature
                 and (expected_inputs,expected_occurrences)==(self.parents['identities'],self.parents['occurrences'])
                 and (self.rows,self.occurrences)==(expected_inputs,expected_occurrences)
                 and all(self.upstream.get(k,0)==input_dispositions[k] for k in UPSTREAM))
        self._check()
        self.regulatory.require_stable()
        self.db.commit();os.fsync(self.descriptor)
        self.signature=_signature(self.path)
        self.finished=True

    def __call__(self, compact):
        value, = self.batch((compact,))
        return value

    def batch(self, compacts):
        """Read at most 64 exact ordinals under one live staging owner.

        The cursor retains one bounded row. Boundary checks cover external file
        changes. Each consumed row still authenticates its digest and identity.
        Staging format and SQLite durability settings remain unchanged.
        """
        _require(not self.closed and self.finished and type(compacts) is tuple
                 and 0 < len(compacts) <= 64
                 and all(c.decision.disposition == 'assess' for c in compacts))
        ordinals = tuple(c.ref[1] for c in compacts)
        _require(ordinals == tuple(sorted(set(ordinals)))
                 and ordinals[0] > self._last_read_ordinal)
        self._check()
        marks = ','.join('?' for _ in ordinals)
        cursor = self.db.execute('SELECT ordinal,identity,input_sha,content,digest FROM facts '
            f'WHERE ordinal IN ({marks}) ORDER BY ordinal', ordinals)
        try:
            for compact in compacts:
                found = cursor.fetchone()
                _require(found is not None and found[:3] == (compact.ref[1],compact.ref[2],compact.ref[4])
                         and len(found[3]) <= MAX_ROW_BYTES and _sha(found[3]) == found[4])
                facts = _facts(_json(found[3]))
                _require((facts.identity_id,facts.original_input_sha256,facts.upstream_disposition) ==
                         (compact.ref[2],compact.ref[4],compact.decision.disposition))
                self.read_rows += 1
                self._last_read_ordinal = compact.ref[1]
                yield facts
            _require(cursor.fetchone() is None)
        finally:
            cursor.close()
        self._check()

    def finish(self):
        _require(not self.closed and self.finished and self.read_rows==self.upstream.get('assess',0)
                 and _signature(self.path)==self.signature)
        self._check()
        self.regulatory.require_stable()

    def close(self):
        if not self.closed:
            self.closed=True
            try:self.db.close()
            finally:
                os.close(self.descriptor)
                self.path.unlink()

    def __enter__(self):
        _require(not self.closed);return self

    def __exit__(self, kind, value, traceback):
        from genosapiens.wgs_admission_diagnostics import close_preserving
        close_preserving(self.close, 'writer_close', primary=value)


class AdmissionProfileFactsProvider:
    """Borrow one live admission owner and the original annotation stream.

    Non-transcript CSQ records are not in LedgerAdmissionInput. Read those from
    the existing authenticated annotation owner, never from broad proof blobs.
    The same derive_profile_candidate_facts function also serves fresh intake.
    """
    def __init__(self, *, admission_owner, annotation_source, source_manifest,
                 working_manifest, clinical_session, regulatory_index, reference,
                 prepared_index=None, prepared_state_root=None):
        from genosapiens.wgs_annotation import BoundedAnnotationSource, iter_annotation_stage_encoded_records
        _require(type(admission_owner) is AuthenticatedAdmissionOwner
                 and type(annotation_source) is BoundedAnnotationSource
                 and annotation_source.result.manifest.sha256 == admission_owner.parents['annotation_stage_sha256']
                 and source_manifest.sha256 == admission_owner.parents['source_manifest_sha256']
                 and working_manifest.sha256 == admission_owner.parents['working_manifest_sha256'])
        admission_owner._assert_stable()
        self.owner = admission_owner
        self.clinical, self.regulatory, self.reference = clinical_session, regulatory_index, reference
        from genosapiens.wgs_admission_projection import _row_parents
        self.route_sources = profile_route_sources(
            annotation_manifest_sha256=annotation_source.result.manifest.sha256,
            native_source_sha256=_row_parents(admission_owner.parents)[2],
            clinical_session=clinical_session, regulatory_index=regulatory_index)
        self.records = iter_annotation_stage_encoded_records(annotation_source.result,
            plan_root=annotation_source.plan_root, plan=annotation_source.plan,
            source_manifest=source_manifest, working_manifest=working_manifest,
            output_root=annotation_source.output_root,
            resource_session_sha256=admission_owner.parents['resource_session_sha256'],
            bundle_manifest_sha256=admission_owner.parents['bundle_manifest_sha256'],
            execution_config=admission_owner.config, prepared_index=prepared_index,
            prepared_state_root=prepared_state_root)
        self.pending = None
        self.previous = 0
        self.previous_input = 0
        self.closed = False
        self.finished = False
        self.metrics = dict(annotation_rows=0, decoded_annotation_rows=0, profile_inputs=0)

    def __enter__(self):
        _require(not self.closed and not self.finished)
        return self

    def _next(self):
        from genosapiens.wgs_annotation import _load_request
        record = next(self.records, None)
        if record is None:
            self.pending = None
            return
        request = _load_request(record.encoded.canonical_request)
        _require(request.retained_ordinal > self.previous
                 and request.parent_working_set_manifest_sha256 == self.owner.parents['working_manifest_sha256'])
        self.previous = request.retained_ordinal
        self.pending = record, request
        self.metrics['annotation_rows'] += 1

    def __call__(self, compact):
        _require(not self.closed and not self.finished and compact.decision.disposition == 'assess')
        value, = self.owner.get_inputs((compact.ref,), purpose='profile')
        return self._derive(compact, value)

    def batch(self, compacts):
        """Read and query one bounded window before deriving ordered facts."""
        _require(not self.closed and not self.finished and 0 < len(compacts) <= 64
                 and all(row.decision.disposition == 'assess' for row in compacts))
        values = self.owner.get_inputs(tuple(row.ref for row in compacts), purpose='profile')
        _require(len(values) == len(compacts))
        alleles = tuple(value.allele for value in values if value.allele is not None)
        if alleles:
            self.clinical.retention_batch(alleles)
        return tuple(self._derive(compact, value) for compact, value in zip(compacts, values))

    def _derive(self, compact, value):
        _require(value.retained_ordinal > self.previous_input)
        self.previous_input = value.retained_ordinal
        self.metrics['profile_inputs'] += 1
        if self.pending is None:
            self._next()
        while self.pending is not None and self.pending[1].retained_ordinal < value.retained_ordinal:
            self._next()
        sources = ()
        if value.transport == 'sequence_vep':
            _require(self.pending is not None and
                     (self.pending[1].retained_ordinal, self.pending[1].identity_id,
                      self.pending[1].source_allele_ids) ==
                     (value.retained_ordinal, value.identity_id, value.source_allele_ids))
            record = self.pending[0].decode()
            self.metrics['decoded_annotation_rows'] += 1
            from genosapiens.wgs_annotation import UnassessedVepRecord
            if type(record) is UnassessedVepRecord:
                _require(value.annotation_status == 'unassessed'
                    and value.unresolved_reasons == (record.reason,))
            else:
                _require(value.annotation_status == 'annotated')
                key = record.request.request_id, 1
                transcripts = record.parsed.transcripts[key]
                actual = value.allele.transcript_consequences or ()
                _require(not transcripts or actual[-len(transcripts):] == transcripts)
                sources = (*(t.source for t in actual), *record.parsed.non_transcripts[key])
            self._next()
        return derive_profile_candidate_facts(value, original_input_sha256=compact.ref[4],
            upstream_disposition=compact.decision.disposition, molecular_sources=sources,
            route_sources=self.route_sources, clinical_session=self.clinical,
            regulatory_index=self.regulatory, reference=self.reference)

    def finish(self):
        """Reconcile the complete original stream before a gate may publish."""
        _require(not self.closed)
        if not self.finished:
            # Includes upstream rows which needed neither profile evidence nor
            # rich proof decoding. Original reader verifies shard/index tails.
            while self.pending is not None:
                self._next()
            for _ in self.records:
                self.metrics['annotation_rows'] += 1
            _require(self.metrics['annotation_rows'] == self.owner.request_count)
            self.finished = True
        self.owner._assert_stable()
        self.regulatory.require_stable()

    def __exit__(self, kind, value, traceback):
        try:
            if kind is None:
                self.finish()
        finally:
            self.records.close()
            self.closed = True


def profile_route_sources(*, annotation_manifest_sha256, native_source_sha256,
                          clinical_session, regulatory_index):
    """Bind gate facts to the same live scientific source owners."""
    from genosapiens.wgs_clinical_session import WgsClinicalResourceSession
    from genosapiens.regulatory_index import RegulatoryIndex
    _require(type(clinical_session) is WgsClinicalResourceSession
             and clinical_session._stack is not None and type(regulatory_index) is RegulatoryIndex
             and regulatory_index.connection._read_cache is not None)
    regulatory_index.require_stable()
    return dict(molecular=_digest(annotation_manifest_sha256), native=_digest(native_source_sha256),
        clinical=_sha(canonical(dict(index=clinical_session.indexes.manifest_sha256,
            closure=clinical_session.indexes.closure_sha256,
            open_targets=clinical_session.open_targets.manifest_sha256,
            clinvar=clinical_session.clinvar.manifest_sha256,
            reference=clinical_session.reference.manifest_sha256))),
        regulatory=regulatory_index.manifest_sha256)


class _ProfileRouteContext:
    """Bind immutable route identities once without caching external file checks."""
    def __init__(self, clinical, regulatory, sources):
        from genosapiens.wgs_clinical_session import _retention_binding
        self.clinical, self.regulatory = clinical, regulatory
        self.binding = _retention_binding(clinical)
        self.sources = tuple(sorted(sources.items()))

    def require(self, sources, clinical, regulatory, reference, native_source):
        from genosapiens.wgs_clinical_session import _retention_binding
        _require(clinical is self.clinical and regulatory is self.regulatory
            and _retention_binding(clinical) == self.binding
            and reference == clinical.reference and tuple(sorted(sources.items())) == self.sources
            and sources['native'] == native_source and sources['regulatory'] == regulatory.manifest_sha256)
        clinical._active()
        clinical._index.assert_stable()
        regulatory.require_stable()


def derive_profile_candidate_facts(value, *, original_input_sha256, upstream_disposition,
        molecular_sources, route_sources, clinical_session, regulatory_index, reference, _clinical_window=None,
        _route_owner=None):
    """Shared intake/reuse facts; no broad annotation-association proof replay.

    molecular_sources is the complete CSQ source tuple from the existing
    authenticated annotation reader/intake, including non-transcript records.
    No caller may substitute a most-severe consequence or a gene-only summary.
    """
    from genosapiens import admission_links as links, annotation_admission as annotation
    from genosapiens import vep_annotation as vep
    from genosapiens.allele_identity import MosaicAllele, MitochondrialAllele, StructuralAllele
    from genosapiens.variant_retention import ExactAllele
    from genosapiens.retention_sources import RetentionEvidence
    from genosapiens.evidence_model import SourcedRecord
    from genosapiens.evidence_serialization import _plain
    from genosapiens.regulatory_index import indexed_regulatory_report
    from genosapiens.wgs_admission import LedgerAdmissionInput
    _require(type(value) is LedgerAdmissionInput and type(route_sources) is dict
             and set(route_sources) == set(ROUTES) and upstream_disposition == 'assess'
             and type(molecular_sources) is tuple
             and all(type(source) is SourcedRecord for source in molecular_sources))
    _digest(original_input_sha256)
    _require(reference == clinical_session.reference)
    if _route_owner is None:
        expected = profile_route_sources(annotation_manifest_sha256=route_sources['molecular'],
            native_source_sha256=value.native_source_contract_sha256,
            clinical_session=clinical_session, regulatory_index=regulatory_index)
        _require(expected == route_sources)
    else:
        _require(type(_route_owner) is _ProfileRouteContext)
        _route_owner.require(route_sources, clinical_session, regulatory_index, reference,
            value.native_source_contract_sha256)
    original_locator = canonical(dict(input_sha256=original_input_sha256,
        ordinal=value.retained_ordinal, identity_id=value.identity_id)).decode().strip()
    def locator(source):
        return canonical(dict(dataset=source.source_dataset, release=source.source_release,
            record_id=source.source_record_id, source_sha256=source.source_sha256,
            manifest_sha256=source.public_manifest_sha256)).decode().strip()
    def evidence(route, target, source_ref, **kwargs):
        return ProfileRouteEvidence(target, source_ref, route_sources[route], **kwargs)
    routes = {}
    allele = value.allele
    native = None if allele is None else links._native(allele).identity
    if _clinical_window is None:
        lookup = clinical_session.retention_batch(() if allele is None else (allele,))
        retention_evidence = (RetentionEvidence(None, (), (), (), 'native_identity_unresolved')
                              if allele is None else lookup[0])
    else:
        from genosapiens.wgs_clinical_session import _RetentionWindow
        _require(type(_clinical_window) is _RetentionWindow)
        retention_evidence = (RetentionEvidence(None, (), (), (), 'native_identity_unresolved')
            if allele is None else _clinical_window.take(clinical_session, allele))
    _require(type(retention_evidence) is RetentionEvidence)
    exact = retention_evidence.allele
    if value.transport == 'sequence_vep' and value.annotation_status == 'unassessed':
        _require(not molecular_sources and bool(value.unresolved_reasons))
        routes['molecular'] = ProfileRouteFacts('molecular', 'unresolved', route_sources['molecular'],
            reasons=value.unresolved_reasons, source_refs=(original_locator,))
    elif value.transport == 'sequence_vep':
        _require(value.annotation_status == 'annotated' and bool(molecular_sources))
        molecular = []
        for source in molecular_sources:
            csq = source.source_fields['csq']
            terms = tuple(vep._value(csq['Consequence']).split('&'))
            transcript = (vep._value(csq['Feature']) or None) if csq['Feature_type'] == 'Transcript' else None
            target = annotation._gene(vep._value(csq['Gene']))
            unknown = bool(set(terms) - (annotation.INTRAGENIC_TERMS | annotation.NON_OVERLAP_TERMS))
            if transcript or unknown:
                molecular.append(evidence('molecular', target, locator(source),
                    target_namespace='Ensembl' if target else None,
                    transcript_id=transcript, consequences=terms,
                    sift=csq.get('SIFT'), polyphen=csq.get('PolyPhen'),
                    interpretation='unresolved' if unknown or not target or not transcript else 'supported'))
        routes['molecular'] = ProfileRouteFacts('molecular', 'matched' if molecular else 'no_match',
            route_sources['molecular'], tuple(molecular),
            ('complete_original_CSQ_sources',), tuple(locator(source) for source in molecular_sources))
    else:
        _require(not molecular_sources)
        routes['molecular'] = ProfileRouteFacts('molecular','not_applicable',route_sources['molecular'])
    native_applicable = value.transport != 'sequence_vep' or type(native) in (MosaicAllele, MitochondrialAllele)
    routes['native'] = ProfileRouteFacts('native', 'matched' if native_applicable else 'not_applicable',
        route_sources['native'], (evidence('native', None, original_locator, target_namespace=None,
            structural_native_identity=value.identity_id if type(native) is StructuralAllele else None,
            interpretation='unresolved' if allele is None or value.transport=='unresolved' else 'supported'),)
            if native_applicable else ())
    if allele is None:
        routes['clinical'] = ProfileRouteFacts('clinical','unresolved',route_sources['clinical'],
            reasons=('native_identity_unresolved_no_sequence_proxy',), source_refs=(original_locator,))
        routes['regulatory'] = ProfileRouteFacts('regulatory','unresolved',route_sources['regulatory'],
            reasons=('native_geometry_unresolved_no_guessed_point',), source_refs=(original_locator,))
    else:
        clinical_evidence = []
        for fact in retention_evidence.clinical:
            clinical_evidence.append(evidence('clinical', fact.gene_id, fact.source.record_id,
                target_namespace='Ensembl' if fact.gene_id else None, source_dataset=fact.source.dataset,
                exact_allele=fact.allele, clinical_classifications=fact.classifications,
                clinical_conditions=fact.disease_ids, clinical_conflicts=fact.conflicts,
                review_status=fact.review_status))
        for fact in retention_evidence.functional:
            clinical_evidence.append(evidence('clinical', fact.gene_id, fact.source.record_id,
                target_namespace='Ensembl' if fact.gene_id else None, source_dataset=fact.source.dataset,
                exact_allele=fact.allele, functional=fact))
        routes['clinical'] = ProfileRouteFacts('clinical',
            'unresolved' if retention_evidence.mapping_reason else 'matched' if clinical_evidence else 'no_match',
            route_sources['clinical'], tuple(clinical_evidence),
            (retention_evidence.mapping_reason,) if retention_evidence.mapping_reason else
            ('complete_indexed_clinical_and_curated_lookup',))
        regulatory, no_links, unknown = [], [], False
        # Exhaust the actual adapter even if a molecular or clinical route passes.
        for proof in indexed_regulatory_report((allele,), index=regulatory_index, reference=reference,metadata_only=True):
            records = [source for source in proof.source_refs if source.source_dataset == links.ENCODE_RE2G]
            metadata=next((source.source_fields for source in proof.source_refs
                if source.source_dataset=='genosapiens.focused-association-scope'
                and source.source_fields.get('kind')=='regulatory-membership-metadata'),None)
            if records or metadata is not None:
                element = next(source for source in proof.source_refs if source.source_dataset == links.ENCODE_CCRE)
                raw = records[0].source_fields['record'] if records else dict(
                    source_target_id=metadata['raw_target'],target_kind=metadata['target_kind'])
                if metadata is not None:
                    _require(metadata['index_manifest_sha256']==regulatory_index.manifest_sha256
                        and metadata['element_ordinal']==element.source_fields['row_ordinal'])
                regulatory.append(evidence('regulatory', proof.target_id,
                    canonical(dict(index_manifest_sha256=regulatory_index.manifest_sha256,
                        element_ordinal=element.source_fields['row_ordinal'],
                        raw_target=raw['source_target_id'], first_source_ordinal=(records[0].source_fields['row_ordinal']
                            if records else metadata['first_source_ordinal']),
                        source_rows=len(records) if records else metadata['source_rows'])).decode().strip(),
                    target_namespace=proof.target_namespace,
                    source_target_id=raw['source_target_id'], source_target_kind=raw['target_kind'],
                    interpretation='supported' if proof.target_id else 'unresolved'))
            elif 'overlapping-element-without-matched-target-link' in proof.reasons:
                no_links.extend(locator(source) for source in proof.source_refs)
            elif proof.outcome == 'unresolved':
                unknown = True
        routes['regulatory'] = ProfileRouteFacts('regulatory',
            'unresolved' if unknown else 'matched' if regulatory else 'no_match',
            route_sources['regulatory'], tuple(regulatory),
            ('overlapping_element_genuinely_without_source_link',) if no_links else (), tuple(no_links))
    technical = {flag for obligation in value.technical_obligations for flag in obligation.technical_flags}
    if any(obligation.read_validation_required for obligation in value.technical_obligations):
        technical.add('mandatory_selected_read_validation_unperformed')
    return ProfileCandidateFacts(value.identity_id, original_input_sha256, upstream_disposition,
        tuple(routes[route] for route in ROUTES), tuple(sorted(technical)), exact, value.source_allele_ids,
        retention_evidence.coverage, tuple(sorted({o.sample_name for o in value.technical_obligations})))


class _Closure:
    """Bounded ordered digests plus disk-backed uniqueness; no source-sized sets."""
    def __init__(self, path, identity):
        self.identity = identity
        self.cohort = _cohort(identity)
        self.db = None if path is None else sqlite3.connect(path)
        if self.db is not None:
            self._open_uniqueness(path)
        self.rows = self.occurrences = self.advancing_rows = self.advancing_occurrences = 0
        self._initialize_counts()

    def _open_uniqueness(self, path):
        self.db.executescript('''PRAGMA mmap_size=0; PRAGMA temp_store=FILE;
            PRAGMA cache_size=-8192; PRAGMA journal_mode=DELETE;
            CREATE TABLE identities(id TEXT PRIMARY KEY);
            CREATE TABLE occurrences(id TEXT PRIMARY KEY);''')
        self.descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
        self.cache = ReadFileCache(self.descriptor)
        self.db.set_progress_handler(lambda: (self.cache.release(), 0)[1], 10000)

    def _initialize_counts(self):
        self.counts = Counter({key: 0 for key in DISPOSITIONS})
        self.retention_counts = Counter({key: 0 for key in RETENTION_DISPOSITIONS})
        self.retention_occurrences = Counter({key: 0 for key in RETENTION_DISPOSITIONS})
        self.upstream = Counter({key: 0 for key in UPSTREAM})
        self.upstream_occurrences = Counter({key: 0 for key in UPSTREAM})
        self.reasons = Counter()
        self.refs = hashlib.sha256()
        self.selected = hashlib.sha256()
        self.advancing_ordinals = hashlib.sha256()
        self.occurrence_digest = hashlib.sha256()
        self.targets = Counter()
        self.regulatory_class_inputs = Counter({kind: 0 for kind in ('ensembl','missing','unmapped')})
        self.regulatory_linked_unknown_inputs = 0

    def add(self, row):
        ref = row.original_ref
        _require(type(ref) is tuple and len(ref) == 6
                 and ref[0] == self.identity['parent']['binding_sha256']
                 and type(ref[1]) is int and ref[1] == (self.rows+1 if self.cohort is None
                     else self.cohort['original_ordinals'][self.rows] if self.rows < self.cohort['rows'] else None)
                 and type(ref[2]) is str and bool(ref[2])
                 and (ref[3] is None or type(ref[3]) is str and bool(ref[3]))
                 and type(ref[5]) is int and 0 < ref[5] <= MAX_ROW_BYTES)
        _digest(ref[4]); _digest(row.screening_sha256)
        _require(type(row.source_allele_ids) is tuple and bool(row.source_allele_ids)
                 and all(type(v) is str and bool(v) for v in row.source_allele_ids)
                 and row.upstream_disposition in UPSTREAM)
        if self.db is not None:
            self.db.execute('INSERT INTO identities VALUES (?)', (ref[2],))
            self.db.executemany('INSERT INTO occurrences VALUES (?)', ((v,) for v in row.source_allele_ids))
        decision = row.decision
        _require(decision['identity_id'] == ref[2] and decision['original_input_sha256'] == ref[4]
                 and decision['upstream_disposition'] == row.upstream_disposition
                 and decision['profile'] == self.identity['profile']
                 and decision['policy_sha256'] == self.identity['policy_sha256']
                 and decision['catalog_sha256'] == self.identity['catalog_sha256']
                 and decision['disposition'] in DISPOSITIONS)
        _require((row.upstream_disposition != 'assess') ==
                 (decision['disposition'] == 'not_applicable_upstream'))
        retained = decision['retention_disposition']
        _require(retained in RETENTION_DISPOSITIONS
                 and decision['retention_policy_sha256'] == RETENTION_POLICY_SHA256
                 and (retained == 'retained') == (decision['disposition'] == 'selected')
                 and (retained == 'not_applicable_upstream') == (row.upstream_disposition != 'assess')
                 and bool(decision['retention_routes']) == (retained == 'retained')
                 and tuple(row.facts['source_occurrences']) == row.source_allele_ids
                 and canonical(row.deferred_inventory) == canonical(_deferred_inventory(row)))
        count = len(row.source_allele_ids)
        self.rows += 1; self.occurrences += count
        self.refs.update(canonical(ref))
        self.occurrence_digest.update(canonical((ref[1], row.source_allele_ids)))
        self.counts[decision['disposition']] += 1
        self.retention_counts[retained] += 1
        self.retention_occurrences[retained] += count
        self.upstream[row.upstream_disposition] += 1
        self.upstream_occurrences[row.upstream_disposition] += count
        self.reasons.update(reason for route in decision['routes'] for reason in route['reasons'])
        regulatory_classes = {e['source_target_kind'] for route in row.facts['routes']
            if route['route']=='regulatory' for e in route['evidence']
            if e['source_target_kind'] is not None}
        self.regulatory_class_inputs.update(regulatory_classes)
        self.regulatory_linked_unknown_inputs += int(bool(regulatory_classes & {'missing','unmapped'}))
        if decision['disposition'] in ('selected', 'unresolved'):
            self.selected.update(canonical(ref))
            self.advancing_ordinals.update((str(ref[1])+'\n').encode('ascii'))
            self.advancing_rows += 1; self.advancing_occurrences += count
            self.targets.update(set((r['target_namespace'], r['target_id']) for r in decision['routes']
                if r['target_id'] is not None and r['outcome'] in ('selected', 'unresolved')))
        if self.db is not None and self.rows % 1024 == 0:
            self.db.commit(); self.cache.release()

    def result(self):
        parent = self.identity['parent'] if self.cohort is None else self.cohort
        _require(self.rows == parent['rows'])
        if self.cohort is None:
            _require(self.occurrences == parent['occurrences']
                and dict(self.upstream) == parent['input_dispositions']
                and dict(self.upstream_occurrences) == parent['occurrence_dispositions'])
        return dict(rows=self.rows, occurrences=self.occurrences, advancing_rows=self.advancing_rows,
            advancing_occurrences=self.advancing_occurrences, dispositions=dict(self.counts),
            retention_dispositions=dict(self.retention_counts), retention_occurrences=dict(self.retention_occurrences),
            upstream_dispositions=dict(self.upstream), upstream_occurrences=dict(self.upstream_occurrences),
            reasons=dict(sorted(self.reasons.items())),
            regulatory_target_class_inputs=dict(self.regulatory_class_inputs),
            regulatory_linked_unknown_inputs=self.regulatory_linked_unknown_inputs,
            target_counts=[dict(namespace=k[0], target=k[1], rows=v) for k,v in sorted(self.targets.items())],
            ordered_input_refs_sha256=self.refs.hexdigest(), advancing_refs_sha256=self.selected.hexdigest(),
            advancing_ordinals_sha256=self.advancing_ordinals.hexdigest(),
            ordered_occurrences_sha256=self.occurrence_digest.hexdigest())

    def close(self):
        if self.db is None:
            return
        try:
            self.db.close()
        finally:
            self.cache.release(force=True); os.close(self.descriptor)


def _profile_fact_domain(admission_owner, identity, facts_provider):
    """Bound original-input reads and public queries to at most 64 identities."""
    from itertools import islice
    domain = iter(_compact_domain(admission_owner, identity))
    while window := tuple(islice(domain, 64)):
        assess = tuple(row for row in window if row.decision.disposition == 'assess')
        if type(facts_provider) in (AdmissionProfileFactsProvider, FreshAdmissionProfileIntake) and assess:
            facts = iter(facts_provider.batch(assess))
        else:
            facts = (facts_provider(row) for row in assess)
        for compact in window:
            yield compact, next(facts) if compact.decision.disposition == 'assess' else None
        # Exhaust the provider so its final ownership check precedes publication.
        _require(next(facts, None) is None)


def build_profile_projection(destination, *, admission_owner, selection, catalog,
        regulatory_index_sha256, consumer_executable_sha256, route_sources, facts_provider,
        maximum_shard_rows=MAX_SHARD_ROWS, maximum_shard_bytes=MAX_SHARD_BYTES, pilot_cohort=None):
    """Publish complete scope closure; facts_provider never receives old proof blobs."""
    scratch = writer = closure = None
    lock = None
    published = None
    try:
        destination = Path(destination)
        _require(not destination.exists() and not any(p.is_symlink() for p in (destination, *destination.parents)))
        _require(callable(facts_provider) and type(maximum_shard_rows) is int
                 and 0 < maximum_shard_rows <= MAX_SHARD_ROWS and type(maximum_shard_bytes) is int
                 and 0 < maximum_shard_bytes <= MAX_SHARD_BYTES)
        identity = _binding(admission_owner, selection, catalog, regulatory_index_sha256,
                            consumer_executable_sha256, route_sources, pilot_cohort=pilot_cohort)
        if type(facts_provider) is AdmissionProfileFactsProvider:
            _require(facts_provider.owner is admission_owner and facts_provider.route_sources == route_sources
                     and facts_provider.regulatory.manifest_sha256 == regulatory_index_sha256)
        if type(facts_provider) is FreshAdmissionProfileIntake:
            _require(facts_provider.finished and facts_provider.binding==admission_owner.binding
                     and facts_provider.selection==selection and facts_provider.catalog==catalog
                     and facts_provider.route_sources==route_sources
                     and facts_provider.regulatory.manifest_sha256==regulatory_index_sha256)
        if type(facts_provider) is ReconciledPublicProfilePilotFacts:
            _require(pilot_cohort is facts_provider.cohort and admission_owner is facts_provider.owner
                and selection == facts_provider.selection and catalog == facts_provider.catalog
                and route_sources == facts_provider.route_sources
                and regulatory_index_sha256 == facts_provider.index_sha256)
        lock = destination.parent / ('.profile-lock-'+destination.name)
        lock.mkdir(mode=0o700)
        scratch = destination.parent / ('.profile-building-'+uuid.uuid4().hex)
        scratch.mkdir(mode=0o700)
        closure = _Closure(scratch/'closure.sqlite', identity)
        writer = _ShardWriter(scratch, 'profile_decisions', max_rows=maximum_shard_rows,
                              max_bytes=maximum_shard_bytes, filename_prefix='decisions')
        for compact, facts in _profile_fact_domain(admission_owner, identity, facts_provider):
            upstream = compact.decision.disposition
            if upstream != 'assess':
                facts = ProfileCandidateFacts(compact.ref[2], compact.ref[4], upstream,
                    tuple(ProfileRouteFacts(r, 'not_applicable', route_sources[r]) for r in ROUTES),
                    source_occurrences=compact.facts.source_allele_ids)
            _require(type(facts) is ProfileCandidateFacts
                     and (facts.identity_id, facts.original_input_sha256, facts.upstream_disposition) ==
                         (compact.ref[2], compact.ref[4], upstream)
                     and facts.source_occurrences == compact.facts.source_allele_ids
                     and all(route.source_sha256 == route_sources[route.route] for route in facts.routes))
            decision = evaluate_profile_candidate(facts, policy=selection.policy, catalog=catalog,
                                                  profile=selection.profile)
            row = ProfileDecisionRow(compact.ref, compact.facts.source_allele_ids,
                _sha(compact.as_bytes()), upstream, asdict(facts), asdict(decision))
            row = replace(row, deferred_inventory=_deferred_inventory(row))
            encoded = encode_row(row)
            _require(len(encoded.content) <= MAX_ROW_BYTES)
            closure.add(row)
            writer.write(row, key=str(compact.ref[1]), canonical_bytes=encoded)
        if type(facts_provider) in (AdmissionProfileFactsProvider,FreshAdmissionProfileIntake,ReconciledPublicProfilePilotFacts):
            facts_provider.finish()
        descriptors = writer.close()
        counts = closure.result()
        closure.close(); closure = None
        (scratch/'closure.sqlite').unlink()
        manifest = dict(**identity, status='complete', counts=counts,
                        shards=[asdict(d) for d in descriptors], stream_sha256=writer.stream.hexdigest())
        content = canonical(manifest)
        with (scratch/'manifest.json').open('xb') as target:
            target.write(content); target.flush(); os.fsync(target.fileno())
        # The producer closed counts, decisions, uniqueness and source refs
        # during writing. Check completed metadata before atomic publication.
        with ProfileProjectionReader(scratch, expected_manifest_sha256=_sha(content),
                selection=selection, catalog=catalog, expected_binding=identity,
                admission_owner=admission_owner, verification_parent=destination.parent):
            pass
        admission_owner._assert_stable()
        _durable(scratch)
        _require(not destination.exists())
        os.rename(scratch, destination); scratch = None; published = destination
        _durable(destination.parent)
        lock.rmdir(); lock = None
        published = None
        return _sha(content)
    except ProfileProjectionError:
        raise
    except (Exception, KeyboardInterrupt):
        raise ProfileProjectionError('profile projection failed before complete publication') from None
    finally:
        if writer is not None and writer._handle is not None:
            writer._handle.close()
        if closure is not None:
            closure.close()
        if scratch is not None:
            shutil.rmtree(scratch)
        if published is not None:
            shutil.rmtree(published)
        if lock is not None and lock.exists():
            lock.rmdir()


def _iter_sealed_profile_rows(root, manifest, identity, signatures, *, _with_facts=False):
    """One bounded decoder for gate replay and retained publication evidence.

    This authenticates sealed facts and decisions. Catalog predicate replay is
    additionally performed by ProfileProjectionReader while the catalog owner
    is live; publication never constructs a catalog from selected targets.
    """
    from genosapiens.track1_analysis_profiles import ProfileDecision, ProfileRouteDecision
    stream = hashlib.sha256()
    for value in manifest['shards']:
        descriptor = ScientificShard(**value)
        path = root/descriptor.relative_path
        _require(_signature(path) == signatures[path])
        digest = hashlib.sha256(); count = size = 0; first = last = None
        with path.open('rb') as source:
            cache = ReadFileCache(source.fileno())
            try:
                while True:
                    raw = source.readline(MAX_ROW_BYTES+1)
                    if not raw:
                        break
                    _require(len(raw) <= MAX_ROW_BYTES and raw.endswith(b'\n'))
                    digest.update(raw); stream.update(raw); size += len(raw); count += 1
                    _require(count <= descriptor.row_count and size <= descriptor.size_bytes)
                    document = _json(raw)
                    row = ProfileDecisionRow(**{**document,
                        'original_ref': tuple(document['original_ref']),
                        'source_allele_ids': tuple(document['source_allele_ids'])})
                    facts = _facts(row.facts)
                    _require(all(r.source_sha256 == identity['route_sources'][r.route] for r in facts.routes)
                        and (facts.identity_id, facts.original_input_sha256, facts.upstream_disposition) ==
                            (row.original_ref[2], row.original_ref[4], row.upstream_disposition))
                    decision = ProfileDecision(**{**row.decision,
                        'technical_reasons': tuple(row.decision['technical_reasons']),
                        'retention_routes': tuple(row.decision['retention_routes']),
                        'routes': tuple(ProfileRouteDecision(**{**route,
                            'reasons': tuple(route['reasons']), 'source_refs': tuple(route['source_refs'])})
                            for route in row.decision['routes'])})
                    _require(canonical(asdict(decision)) == canonical(row.decision)
                        and decision.technical_reasons == facts.technical_reasons
                        and set(route.route for route in decision.routes) == set(ROUTES)
                        and all(route.route in ROUTES and route.outcome in DISPOSITIONS
                            and route.source_sha256 == identity['route_sources'][route.route]
                            and type(route.reasons) is tuple and bool(route.reasons)
                            and all(type(reason) is str and bool(reason) for reason in route.reasons)
                            and all(type(ref) is str and bool(ref) for ref in route.source_refs)
                            for route in decision.routes))
                    key = str(row.original_ref[1]); first = key if first is None else first; last = key
                    cache.release()
                    yield (row, facts) if _with_facts else row
                _require((size, count, digest.hexdigest(), first, last) ==
                    (descriptor.size_bytes, descriptor.row_count, descriptor.sha256,
                     descriptor.first_key, descriptor.last_key))
            finally:
                cache.release(force=True)
        _require(_signature(path) == signatures[path])
    _require(stream.hexdigest() == manifest['stream_sha256'])


def _validate_completed_counts(manifest, identity):
    cohort = _cohort(identity)
    parent = identity['parent']
    counts = manifest['counts']
    _require(type(counts) is dict and type(counts['rows']) is int and counts['rows'] == (parent['rows'] if cohort is None else cohort['rows'])
        and type(counts['occurrences']) is int and counts['occurrences'] >= counts['rows']
        and sum(s['row_count'] for s in manifest['shards']) == counts['rows'])
    for key, keys, total in (('dispositions', DISPOSITIONS, 'rows'),
            ('retention_dispositions', RETENTION_DISPOSITIONS, 'rows'),
            ('retention_occurrences', RETENTION_DISPOSITIONS, 'occurrences'),
            ('upstream_dispositions', UPSTREAM, 'rows'), ('upstream_occurrences', UPSTREAM, 'occurrences')):
        values = counts[key]
        _require(type(values) is dict and set(values) == set(keys)
            and all(type(v) is int and v >= 0 for v in values.values()) and sum(values.values()) == counts[total])
    _require(type(counts['advancing_rows']) is int and counts['advancing_rows'] == counts['dispositions']['selected'] + counts['dispositions']['unresolved']
        and type(counts['advancing_occurrences']) is int
        and counts['advancing_rows'] <= counts['advancing_occurrences'] <= counts['occurrences'])
    _require(counts['dispositions']['unresolved'] == 0
        and counts['retention_dispositions']['retained'] == counts['advancing_rows']
        and counts['retention_occurrences']['retained'] == counts['advancing_occurrences']
        and counts['retention_dispositions']['deferred_no_retention_evidence'] == counts['dispositions']['outside_profile'])
    if cohort is None:
        _require(counts['occurrences'] == parent['occurrences']
            and counts['upstream_dispositions'] == parent['input_dispositions']
            and counts['upstream_occurrences'] == parent['occurrence_dispositions'])
    for key in ('ordered_input_refs_sha256','advancing_refs_sha256',
                'advancing_ordinals_sha256','ordered_occurrences_sha256'):
        _digest(counts[key])
    _digest(manifest['stream_sha256'])


def _profile_component(manifest, manifest_sha256):
    return dict(schema='genosapiens.wgs-analysis-profile-component', version=1,
        manifest_path='analysis/profile/manifest.json', manifest_sha256=_digest(manifest_sha256),
        identity={key: value for key, value in manifest.items()
                  if key not in {'status', 'counts', 'shards', 'stream_sha256'}},
        counts=manifest['counts'], stream_sha256=manifest['stream_sha256'])


def analysis_profile_component(reader):
    """Reference a fully replayed gate without retaining a disposable index."""
    _require(type(reader) is ProfileProjectionReader and reader.root.name == 'profile')
    reader.require_stable()
    # Return independent bounded metadata, not aliases to a mutable reader.
    return _json(canonical(_profile_component(reader.manifest, reader.manifest_sha256)))


def profile_component_parents(value):
    _require(type(value) is dict and value.get('schema') == 'genosapiens.wgs-analysis-profile-component'
             and type(value.get('version')) is int and value['version'] == 1)
    identity = value['identity']
    return dict(analysis_profile_manifest=_digest(value['manifest_sha256']),
        analysis_profile_policy=_digest(identity['policy_sha256']),
        analysis_profile_catalog=_digest(identity['catalog_sha256']),
        analysis_profile_regulatory_index=_digest(identity['regulatory_index_sha256']),
        analysis_profile_adapter=_digest(identity['adapter_sha256']))


def verify_analysis_profile_component(scope, value, *, full=True, scratch_directory=None, profile_root=None):
    """Authenticate retained sealed scope evidence and complete ordered closure.

    Scientific readers authenticate the exact selected domain certified by the
    original gate. They do not relabel upstream inputs or reopen their database.
    """
    closure = None
    try:
        root = Path(scope)/'analysis/profile' if profile_root is None else Path(profile_root)
        _require(root.is_absolute() and root.name == 'profile'
            and type(value) is dict and set(value) == {'schema','version','manifest_path',
                'manifest_sha256','identity','counts','stream_sha256'}
            and value['schema'] == 'genosapiens.wgs-analysis-profile-component'
            and type(value['version']) is int and value['version'] == 1
            and value['manifest_path'] == 'analysis/profile/manifest.json')
        manifest_signature = _signature(root/'manifest.json')
        content = read_bounded_metadata(root/'manifest.json')
        _require(_sha(content) == _digest(value['manifest_sha256']))
        manifest = _json(content)
        identity = value['identity']
        _require(type(identity) is dict and set(identity) - {'finite_cohort'} == {'schema','version','classification','llm_policy','profile',
            'policy_sha256','catalog_sha256','regulatory_index_sha256','consumer_executable_sha256',
            'adapter_sha256','route_sources','parent'}
            and set(manifest) == set(identity) | {'status','counts','shards','stream_sha256'}
            and manifest['status'] == 'complete' and identity['schema'] == SCHEMA
            and type(identity['version']) is int and identity['version'] == 3
            and identity['classification'] == 'controlled' and identity['llm_policy'] == 'no-llm'
            and identity['profile'] in {'focused','exhaustive'}
            and canonical(value) == canonical(_profile_component(manifest, _sha(content))))
        for key in ('policy_sha256','catalog_sha256','regulatory_index_sha256',
                    'consumer_executable_sha256','adapter_sha256'):
            _digest(identity[key])
        _require(type(identity['route_sources']) is dict and set(identity['route_sources']) == set(ROUTES))
        for digest in identity['route_sources'].values():
            _digest(digest)
        parent = identity['parent']
        _require(type(parent) is dict and set(parent) == {'binding_sha256','completion_sha256',
            'producer_executable_sha256','source_manifest_sha256','reference_sha256',
            'input_rows_sha256','rows','occurrences','input_dispositions','occurrence_dispositions'})
        for key in ('binding_sha256','completion_sha256','producer_executable_sha256',
                    'source_manifest_sha256','reference_sha256','input_rows_sha256'):
            _digest(parent[key])
        _require(type(parent['rows']) is int and parent['rows'] >= 0
            and type(parent['occurrences']) is int and parent['occurrences'] >= parent['rows'])
        for key, total in (('input_dispositions','rows'),('occurrence_dispositions','occurrences')):
            _require(type(parent[key]) is dict and set(parent[key]) == set(UPSTREAM)
                and all(type(v) is int and v >= 0 for v in parent[key].values())
                and sum(parent[key].values()) == parent[total])
        _validate_completed_counts(manifest, identity)
        _digest(manifest['stream_sha256'])
        _require(type(manifest['shards']) is list)
        shards = tuple(ScientificShard(**d) for d in manifest['shards'])
        _require(all(d.kind == 'profile_decisions' and d.relative_path.startswith('decisions-') for d in shards)
            and len({d.relative_path for d in shards}) == len(shards))
        signatures = {root/'manifest.json': manifest_signature}
        signatures.update({root/d.relative_path: _signature(root/d.relative_path) for d in shards})
        _require(_signature(root/'manifest.json') == manifest_signature)
        if full:
            with tempfile.TemporaryDirectory(prefix='.profile-verification-', dir=scratch_directory) as temporary:
                closure = _Closure(Path(temporary)/'closure.sqlite', identity)
                try:
                    for row in _iter_sealed_profile_rows(root, manifest, identity, signatures):
                        closure.add(row)
                    _require(closure.result() == manifest['counts'])
                finally:
                    closure.close(); closure = None
        _require(all(_signature(path) == signature for path, signature in signatures.items()))
        return manifest
    except (Exception, KeyboardInterrupt):
        raise ProfileProjectionError('retained profile component is incomplete, foreign or corrupt') from None


def profile_retained_files(scope, value, *, profile_root=None):
    root = Path(scope)/'analysis/profile' if profile_root is None else Path(profile_root)
    manifest = verify_analysis_profile_component(scope, value, full=False, profile_root=root)
    path = root/'manifest.json'
    yield path, 'analysis/profile/manifest.json', value['manifest_sha256'], path.stat().st_size
    for value in manifest['shards']:
        shard = ScientificShard(**value)
        yield root/shard.relative_path, 'analysis/profile/'+shard.relative_path, shard.sha256, shard.size_bytes


def scientific_profile_component(artifact, directory, *, profile_root=None):
    """Reconstruct the exact retained gate from explicit scientific parents."""
    parents = {key: value for key, value in artifact.parent_bindings if key.startswith('analysis_profile_')}
    if not parents:
        return None
    root = Path(directory).parent/'analysis/profile' if profile_root is None else Path(profile_root)
    content = read_bounded_metadata(root/'manifest.json')
    value = _profile_component(_json(content), _sha(content))
    _require(parents == profile_component_parents(value))
    verify_analysis_profile_component(root.parent, value, full=False, profile_root=root)
    return value


def verify_profile_assessment_binding(projections, value):
    """Selected scope and actual indexed discovery share one closed domain."""
    _require(type(projections) is tuple and len(projections) == 2)
    associations, assessments = projections
    identity, counts = value['identity'], value['counts']
    _require(associations.profile_manifest_sha256 == value['manifest_sha256']
        and associations.regulatory_index_manifest_sha256 == identity['regulatory_index_sha256']
        and associations.admission_completion_sha256 == identity['parent']['completion_sha256']
        and associations.original_inputs == identity['parent']['rows']
        and associations.row_count == assessments.assessed_inputs == counts['advancing_rows']
        and associations.subject_occurrences == counts['advancing_occurrences']
        and associations.assessed_ordinals_sha256 == counts['advancing_ordinals_sha256'])


def _checked_metadata(value):
    """Immutable type-preserving comparison tree, built only at authentication."""
    if type(value) is dict:
        return (dict, tuple((key, _checked_metadata(child)) for key, child in value.items()))
    if type(value) is list:
        return (list, tuple(_checked_metadata(child) for child in value))
    return (type(value), value)


def _metadata_unchanged(value, checked):
    kind, original = checked
    if type(value) is not kind:
        return False
    if kind is dict:
        return len(value) == len(original) and all(
            key in value and _metadata_unchanged(value[key], child) for key, child in original)
    if kind is list:
        return len(value) == len(original) and all(
            _metadata_unchanged(child, before) for child, before in zip(value, original))
    return value == original


class ProfileProjectionReader:
    """Bounded completed metadata reuse, with validation during requested reads."""
    def __init__(self, root, *, expected_manifest_sha256, selection, catalog,
            expected_binding, admission_owner=None, verification_parent=None, completed_context_inputs=None):
        try:
            adapter = profile_adapter_identity()
            if completed_context_inputs is not None:
                from genosapiens.track1_roster import validate_context_restart_acceptance
                _require(type(completed_context_inputs) is CompletedProfileDiscoveryInputs)
                validate_context_restart_acceptance(completed_context_inputs)
                original = completed_context_inputs.identity
                _require(original['profile_manifest']['sha256'] == expected_manifest_sha256
                    and expected_binding['consumer_executable_sha256'] == completed_context_inputs.profile_producer_executable_sha256)
                adapter = original['profile_adapter_sha256']
            self.root = Path(root)
            manifest_signature = _signature(self.root/'manifest.json')
            content = read_bounded_metadata(self.root/'manifest.json')
            _require(_sha(content) == expected_manifest_sha256)
            self._manifest_content = content
            self.manifest = _json(content)
            _require(set(self.manifest) == set(expected_binding) | {'status', 'counts', 'shards', 'stream_sha256'})
            self.identity = {key: self.manifest[key] for key in expected_binding}
            self._identity_content = canonical(self.identity)
            _require(self.identity == expected_binding and self.identity['schema'] == SCHEMA
                     and self.identity['version'] == 3 and self.manifest['status'] == 'complete'
                     and self.identity['classification'] == 'controlled' and self.identity['llm_policy'] == 'no-llm'
                     and self.identity['adapter_sha256'] == adapter
                     and self.identity['profile'] == selection.profile
                     and self.identity['policy_sha256'] == selection.policy.sha256
                     and self.identity['catalog_sha256'] == catalog.sha256)
            self.selection, self.catalog = selection, catalog
            self.shards = tuple(ScientificShard(**d) for d in self.manifest['shards'])
            _require(all(d.kind == 'profile_decisions' and d.relative_path.startswith('decisions-')
                         for d in self.shards)
                     and len({d.relative_path for d in self.shards}) == len(self.shards))
            self.manifest_sha256 = expected_manifest_sha256
            _require(_signature(self.root/'manifest.json') == manifest_signature)
            self.signatures = {self.root/'manifest.json': manifest_signature}
            self.signatures.update({self.root/d.relative_path: _signature(self.root/d.relative_path) for d in self.shards})
            _require(all((self.root/d.relative_path).stat().st_size == d.size_bytes for d in self.shards))
            if admission_owner is not None:
                admission_owner._assert_stable()
                _require(type(admission_owner) is AuthenticatedAdmissionOwner
                    and self.identity['parent']['binding_sha256'] == admission_owner.binding_sha256
                    and self.identity['parent']['completion_sha256'] == admission_owner.completion.sha256)
            _validate_completed_counts(self.manifest, self.identity)
            self._checked_manifest = _checked_metadata(self.manifest)
            self._checked_identity = _checked_metadata(self.identity)
            self._manifest_owner, self._identity_owner = self.manifest, self.identity
            self.require_stable()
        except (Exception, KeyboardInterrupt):
            raise ProfileProjectionError('profile projection is incomplete, foreign or corrupt') from None

    def require_stable(self):
        _require(self.manifest is self._manifest_owner and self.identity is self._identity_owner
                 and _metadata_unchanged(self.manifest, self._checked_manifest)
                 and _metadata_unchanged(self.identity, self._checked_identity))
        _require(all(_signature(path) == signature for path,signature in self.signatures.items()))

    def iter_rows(self):
        self.require_stable()
        # Uniqueness is certified by the producer. Recompute ordered closure
        # while the caller actually consumes rows, without another disk index.
        closure = _Closure(None, self.identity)
        for row, facts in _iter_sealed_profile_rows(self.root, self.manifest, self.identity,
                                                  self.signatures, _with_facts=True):
            decision = evaluate_profile_candidate(facts, policy=self.selection.policy,
                catalog=self.catalog, profile=self.selection.profile)
            _require(canonical(asdict(decision)) == canonical(row.decision))
            closure.add(row)
            yield row
        _require(closure.result() == self.manifest['counts'])
        self.require_stable()

    def __enter__(self):
        self.require_stable(); return self

    def iter_deferred_inventory(self):
        """Read the compact inventory without opening original input objects."""
        for row in self.iter_rows():
            if row.deferred_inventory is not None:
                yield _json(canonical(row.deferred_inventory))

    def __exit__(self, *_):
        self.require_stable()


def _stage_completed_files(source, destination, *, manifest_sha256, shards, progress=None):
    """Copy one manifest-listed closure with bounded hashing and manifest-last commit."""
    from genosapiens.track1_publication import _publish_immutable_report_shard
    source, destination = Path(source), Path(destination)
    _require(source != destination and not destination.exists()
        and not any(p.is_symlink() for p in (destination, *destination.parents)))
    descriptors = tuple(shards)
    paths = tuple(d.relative_path for d in descriptors)
    _require(len(set(paths)) == len(paths) and all(
        Path(p).name == p and p not in {'manifest.json', '.', '..'} for p in paths))
    manifest = source/'manifest.json'
    signatures = {source/d.relative_path: _signature(source/d.relative_path) for d in descriptors}
    signatures[manifest] = _signature(manifest)
    _require(all(signatures[source/d.relative_path][2] == d.size_bytes for d in descriptors))
    total = sum(d.size_bytes for d in descriptors) + signatures[manifest][2]
    destination.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix='.'+destination.name+'-staging-', dir=destination.parent))
    published = False
    copied_files = copied_bytes = 0
    def report():
        from genosapiens.wgs_runtime import analysis_progress
        analysis_progress(progress, 'completed_input_reuse', 'staging',
            counts={'staging_files_completed':copied_files, 'staging_files_total':len(descriptors)+1,
                'staging_bytes_completed':copied_bytes, 'staging_bytes_total':total},
            unit='bytes', completed=copied_bytes, total=total)
    report()
    try:
        for descriptor in descriptors:
            path = source/descriptor.relative_path
            _require(_signature(path) == signatures[path])
            _publish_immutable_report_shard(path, temporary/descriptor.relative_path,
                expected_sha256=descriptor.sha256, expected_size=descriptor.size_bytes, root=temporary)
            _require(_signature(path) == signatures[path])
            copied_files += 1
            copied_bytes += descriptor.size_bytes
            report()
        _durable(temporary)
        _require(all(_signature(p) == signature for p,signature in signatures.items()))
        _publish_immutable_report_shard(manifest, temporary/'manifest.json',
            expected_sha256=manifest_sha256, expected_size=signatures[manifest][2], root=temporary)
        _require(_signature(manifest) == signatures[manifest])
        _durable(temporary)
        _require(not destination.exists())
        os.rename(temporary, destination)
        published = True
        _durable(destination.parent)
        copied_files += 1
        copied_bytes += signatures[manifest][2]
        report()
        return total
    except BaseException as error:
        from genosapiens.wgs_failures import run_cleanup
        run_cleanup((('context_store', lambda: shutil.rmtree(destination if published else temporary)),),
            primary=error)
        raise


def stage_completed_profile_projection(destination, *, inputs, selection, catalog,
        expected_binding, admission_owner, progress=None):
    """Preserve producer metadata and stage only its authenticated dependency list."""
    _require(type(inputs) is CompletedProfileDiscoveryInputs)
    descriptor = inputs.identity['profile_manifest']
    source = Path(descriptor['path']).parent
    _require(_signature(source/'manifest.json')[2] == descriptor['size_bytes'])
    arguments = dict(expected_manifest_sha256=descriptor['sha256'], selection=selection,
        catalog=catalog, expected_binding=expected_binding, admission_owner=admission_owner,
        completed_context_inputs=inputs)
    reader = ProfileProjectionReader(source, **arguments)
    _stage_completed_files(source, destination, manifest_sha256=descriptor['sha256'], shards=reader.shards, progress=progress)
    reader.require_stable()
    return ProfileProjectionReader(destination, **arguments)


def iter_selected_admission_rows(reader, admission_owner):
    """Read advancing original inputs only; never request old admission proofs."""
    _require(type(reader) is ProfileProjectionReader and type(admission_owner) is AuthenticatedAdmissionOwner
             and reader.identity['parent']['binding_sha256'] == admission_owner.binding_sha256
             and reader.identity['parent']['completion_sha256'] == admission_owner.completion.sha256)
    for row in reader.iter_rows():
        if row.decision['disposition'] in ('selected', 'unresolved'):
            value, = admission_owner.get_inputs((row.original_ref,), purpose='discovery')
            yield row, value
    admission_owner._assert_stable()
