"""One live scientific resource producer with independent context derivatives."""
from contextlib import ExitStack
from copy import copy
from pathlib import Path
import tempfile
import time

from genosapiens.wgs_case_resources import ScientificCaseResources
from genosapiens.wgs_evidence_store import EvidenceWriter, EvidenceStoreError
from genosapiens.evidence_serialization import _owned_canonical_bytes as canonical

SCIENCE_INDEX_CACHE_KIB = 512 * 1024

def _containers(value):
    """Copy mutable metadata containers; preserve privately owned frozen values."""
    if type(value) is dict:
        return {key: _containers(child) for key, child in value.items()}
    if type(value) is list:
        return [_containers(child) for child in value]
    if type(value) is tuple:
        return tuple(_containers(child) for child in value)
    return value


class ScientificResourceTemplate:
    """Never accepts persisted validation flags or a foreign resource writer.

    The ordinary factory validates and publishes all resources once. Contexts
    receive private writer owners and policy tables. Only immutable prepared
    scientific inputs and the case-scoped phenotype score cache are shared.
    """
    def __init__(self, factory, *, scratch_directory):
        self.stack = ExitStack()
        self.closed = False
        self.started = time.monotonic()
        self.pid = __import__('os').getpid()
        self.directory = Path(tempfile.mkdtemp(prefix='context-science-resources-', dir=scratch_directory))
        try:
            self.writer = self.stack.enter_context(EvidenceWriter(self.directory,
                scratch_directory=scratch_directory, _skip_verification=True))
            self.writer.connection.execute(f'PRAGMA cache_size=-{SCIENCE_INDEX_CACHE_KIB}')
            self.resources, self.opened = factory(self.writer, self.stack)
            self.phase_seconds = dict(getattr(factory, 'phase_seconds', {}))
            if type(self.resources) is not ScientificCaseResources:
                raise EvidenceStoreError('scientific resource producer has a foreign type')
            self.resources.require(self.writer)
            from genosapiens.wgs_science_sources import ScientificSourceTemplate
            source_started = time.monotonic()
            self.source_records = self.stack.enter_context(ScientificSourceTemplate(self.writer, self.resources))
            self.phase_seconds['source_registry_preparation'] = time.monotonic()-source_started
            self.population_tables = tuple(self.resources.population.tables)
            self.tables = (*self.population_tables, *self.source_records.tables)
            self.binding = canonical(self.resources.header)
            self.parents = canonical(self.resources.parent_bindings)
            self.owners = self._owners()
            self.writer.connection.commit()
            self.evidence_prefix = self.writer.produced_prefix(extra_tables=self.tables)
            self.index_path = Path(self.writer.connection.execute('PRAGMA database_list').fetchone()[2])
            self.elapsed_seconds = time.monotonic() - self.started
            self.phase_seconds['other_initialization'] = max(0, self.elapsed_seconds-sum(self.phase_seconds.values()))
        except BaseException:
            self.close()
            raise

    def _owners(self):
        r = self.resources
        p = r.phenotypes
        profile = p.profile
        return (r, r.genes, r.genes.prepared, r.phenotypes, r.population,
            p.view, p.observations, profile, p.ic_store, p.ic_resource, p.prepared,
            p.prepared.profile_store, p.prepared.resource_binding,
            *(() if profile is None else (profile.evidence_id, profile.public_manifest_sha256, profile.store)),
            r.genotype_policies, r.hypothesis_policy, r.fusion_policy, r.prior_policy, r.diagnosis, r.reference_binding)

    def bind(self, writer, stack):
        if self.closed or type(writer) is not EvidenceWriter or type(stack) is not ExitStack:
            raise EvidenceStoreError('scientific resource template is closed or foreign')
        owners = self._owners()
        if (len(owners) != len(self.owners)
                or any(a is not b for a, b in zip(owners, self.owners, strict=True))
                or canonical(self.resources.header) != self.binding
                or canonical(self.resources.parent_bindings) != self.parents
                or tuple(self.resources.population.tables) != self.population_tables):
            raise EvidenceStoreError('scientific resource template dependencies changed')
        self.resources.require(self.writer)
        self.source_records.require()
        writer.connection.execute(f'PRAGMA cache_size=-{SCIENCE_INDEX_CACHE_KIB}')
        writer.copy_produced_prefix(self.writer, extra_tables=self.tables)
        resources = copy(self.resources)
        resources.writer, resources.connection = writer, writer.connection
        resources.header = _containers(self.resources.header)
        resources.parent_bindings = dict(self.resources.parent_bindings)
        resources.assessment_domain = _containers(self.resources.assessment_domain)
        resources.analysis_eligibility = _containers(self.resources.analysis_eligibility)
        resources.analysis_profile = _containers(self.resources.analysis_profile)
        resources.genes = copy(self.resources.genes)
        resources.genes.writer, resources.genes.connection = writer, writer.connection
        resources.phenotypes = copy(self.resources.phenotypes)
        resources.phenotypes.writer, resources.phenotypes.connection = writer, writer.connection
        resources.phenotypes.ic_resource = copy(self.resources.phenotypes.ic_resource)
        resources.phenotypes.ic_resource.writer = writer
        resources.phenotypes.ic_resource.connection = writer.connection
        resources.population = copy(self.resources.population)
        resources.population.connection = writer.connection
        resources.population.tables = list(self.population_tables)
        resources.population.closed = False
        stack.callback(resources.population.close)
        resources.require(writer)
        resources._source_prefix = self.source_records.bind(writer, resources)
        return resources, self.opened

    def attach_worker(self, stack, open_reference):
        """Replace inherited handles before using the frozen resource graph.

        Fork preserves immutable Python values through copy-on-write. The shared
        SQLite file is opened read-only by each child. Parent cleanup callbacks
        never execute in children, and children never finalize the shared writer.
        """
        import os
        import sqlite3
        from genosapiens.wgs_science_sources import _version
        if os.getpid() == self.pid or self.closed:
            raise EvidenceStoreError('resource worker attachment requires a new process')
        old = self.writer.connection
        old.close()
        connection = sqlite3.connect(self.index_path.as_uri() + '?mode=ro', uri=True)
        connection.execute('PRAGMA query_only=ON')
        connection.execute('PRAGMA cache_size=-65536')
        stack.callback(connection.close)
        self.writer.connection = connection
        for owner in (self.resources, self.resources.genes, self.resources.phenotypes,
                self.resources.phenotypes.ic_resource, self.resources.population):
            owner.connection = connection
        self.evidence_prefix.connection = connection
        self.evidence_prefix.version = self.evidence_prefix._version()
        self.source_records.version = _version(connection)
        self.opened = stack.enter_context(open_reference())
        self.pid = os.getpid()
        self.evidence_prefix.require()
        self.source_records.require()

    def close(self):
        if not self.closed:
            self.closed = True
            self.stack.close()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()


def prepare_live_gene_evidence(evidence):
    """Build scientific indexes without revalidating producer-owned records."""
    import hashlib
    from types import MappingProxyType
    from genosapiens.genotype_evaluation import PreparedGeneEvidence, _GENE_SEAL
    from genosapiens.disease_gene_evidence import PreparedDiseaseGeneStore
    from genosapiens.evidence_model import EvidenceReference
    from genosapiens.candidate_admission import _relation_edges
    owner = PreparedGeneEvidence(_seal=_GENE_SEAL)
    owner.evidence = tuple(evidence)
    digest = hashlib.sha256(b'[')
    size = 1
    for ordinal, public in enumerate(owner.evidence):
        if ordinal:
            digest.update(b',')
            size += 1
        content = canonical(public)[:-1]
        digest.update(content)
        size += len(content)
    digest.update(b']\n')
    owner.sha256, owner.canonical_size = digest.hexdigest(), size + 2
    owner.registry = MappingProxyType({EvidenceReference(e.evidence_id, key, e.public_manifest_sha256): record
        for e in owner.evidence for key, record in e.store.records.items()})
    owner.gene_stores = frozenset(e.evidence_id for e in owner.evidence)
    assertions = {}
    for public in owner.evidence:
        for assertion in public.store.assertions:
            assertions.setdefault((public.evidence_id, assertion.record_id), []).append(assertion)
    owner.assertions_by_record = MappingProxyType({key: tuple(rows) for key, rows in assertions.items()})
    owner.stores = MappingProxyType({public.evidence_id: PreparedDiseaseGeneStore(public.store)
        for public in owner.evidence})
    owner.relation_edges = MappingProxyType({key: tuple(rows)
        for key, rows in _relation_edges(owner.evidence).items()})
    owner._ready = True
    return owner
