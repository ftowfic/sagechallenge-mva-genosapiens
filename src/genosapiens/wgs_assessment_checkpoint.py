"""Explicit complete context snapshot for a declared public pilot continuation.

This is a stage-start boundary, never a checkpoint of a partial evidence tree.
The producer sealer performs closure once. Consumers authenticate the immutable
descriptor and copy/hash its database into private scratch before using it.
"""
from contextlib import closing
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import sqlite3
import time

from genosapiens.wgs_ledger import canonical
from genosapiens.wgs_shards import read_bounded_metadata

SCHEMA = 'genosapiens.complete-assessment-input'
TABLES = frozenset(('batches', 'inputs', 'assessment_subjects', 'contexts',
    'seed_payloads', 'baseline_seeds', 'seeds', 'members', 'decisions',
    'public_payloads', 'public_sources', 'context_payloads'))
BINDINGS = frozenset(('producer_candidate_sha256', 'producer_executable_stack_sha256',
    'producer_plan_sha256', 'execution_config_sha256', 'admission_completion_sha256',
    'admission_binding_sha256', 'eligibility_manifest_sha256', 'profile_manifest_sha256',
    'association_manifest_sha256', 'policy_sha256', 'parent_source_sha256',
    'parent_working_sha256', 'native_contract_sha256'))


class AssessmentCheckpointError(ValueError):
    """Static recovery errors; no input identifiers in diagnostics."""


def require(value):
    if not value:
        raise AssessmentCheckpointError('complete assessment checkpoint is invalid')


def digest(value):
    return type(value) is str and len(value) == 64 and all(c in '0123456789abcdef' for c in value)


def signature(path):
    path = Path(path)
    require(path.is_absolute() and '..' not in path.parts
        and not any(p.is_symlink() for p in (path, *path.parents)) and path.is_file())
    s = path.stat()
    return s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns


def file_sha(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def durable(path):
    with Path(path).open('rb') as stream:
        os.fsync(stream.fileno())


def readonly(path):
    signature(path)
    connection = sqlite3.connect(Path(path).as_uri() + '?mode=ro', uri=True)
    connection.enable_load_extension(False)
    connection.execute('PRAGMA trusted_schema=OFF')
    connection.execute('PRAGMA cache_size=-32768')
    connection.execute('PRAGMA temp_store=FILE')
    connection.execute('PRAGMA mmap_size=0')
    connection.execute('PRAGMA query_only=ON')
    return connection


def seal_snapshot(source, destination, *, bindings, expected_inputs,
                  expected_ordinals_sha256, expected_baselines, deadline_seconds=600):
    """Back up committed SQLite/WAL state without modifying the live producer.

    Caller binds the declared complete producer boundary and its expectations;
    progress alone is not that authority. Failure never writes a manifest.
    """
    require(set(bindings) == BINDINGS and all(digest(v) for v in bindings.values()))
    require(type(expected_inputs) is int and expected_inputs >= 0
        and type(expected_baselines) is int and expected_baselines >= 0
        and digest(expected_ordinals_sha256) and 0 < deadline_seconds <= 1800)
    source, destination = Path(source), Path(destination)
    before = signature(source)
    destination.mkdir(mode=0o700)
    target = destination/'contexts.sqlite'
    deadline = time.monotonic() + deadline_seconds
    def progress(*_):
        require(time.monotonic() < deadline)
    with closing(readonly(source)) as old, closing(sqlite3.connect(target)) as new:
        os.chmod(target, 0o600)
        new.execute('PRAGMA cache_size=-32768')
        new.execute('PRAGMA temp_store=FILE')
        old.backup(new, pages=1024, progress=progress)
        new.execute('PRAGMA journal_mode=DELETE')
        new.set_progress_handler(lambda: int(time.monotonic() >= deadline), 10000)
        tables = {r[0] for r in new.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        require(tables in (TABLES, TABLES | {'context_public_baseline'}))
        if 'context_public_baseline' in tables:
            # Science may have built this disposable projection from the
            # complete public_sources/public_payloads tables. It is not a
            # checkpoint input or authority. Remove only in the private copy;
            # the normal closed-context owner rebuilds it on restart.
            new.execute('DROP TABLE context_public_baseline')
        require(new.execute('PRAGMA integrity_check').fetchall() == [('ok',)])
        require(new.execute('PRAGMA foreign_key_check').fetchone() is None)
        require(new.execute('SELECT count(*) FROM inputs').fetchone()[0] == expected_inputs)
        ordinals, count = hashlib.sha256(), 0
        for ordinal, binding in new.execute('SELECT s.retained_ordinal,i.projection_binding_sha256 '
                'FROM assessment_subjects s JOIN inputs i USING(retained_ordinal) ORDER BY s.retained_ordinal'):
            require(binding == bindings['admission_binding_sha256'])
            ordinals.update((str(ordinal)+'\n').encode()); count += 1
        require(count == expected_inputs and ordinals.hexdigest() == expected_ordinals_sha256)
        batches = 0
        expected = tuple(bindings[k] for k in ('parent_source_sha256', 'parent_working_sha256',
            'native_contract_sha256', 'policy_sha256'))
        for ordinal, *parents in new.execute('SELECT * FROM batches ORDER BY batch_ordinal'):
            batches += 1
            require(ordinal == batches and tuple(parents) == expected)
        require(bool(batches) == bool(expected_inputs))
        require(new.execute('SELECT count(*) FROM baseline_seeds').fetchone()[0] == expected_baselines)
        require(new.execute('SELECT 1 FROM contexts c WHERE NOT EXISTS '
            '(SELECT 1 FROM seeds s WHERE s.context_id=c.context_id) LIMIT 1').fetchone() is None)
        counts = dict(inputs=count, baselines=expected_baselines, batches=batches)
        for table in ('contexts', 'decisions', 'public_sources'):
            counts[table] = new.execute('SELECT count(*) FROM '+table).fetchone()[0]
        schema_sha = hashlib.sha256(canonical(list(new.execute(
            "SELECT type,name,tbl_name,sql FROM sqlite_master ORDER BY type,name")))).hexdigest()
        new.commit()
    require(signature(source) == before)
    durable(target)
    manifest = dict(schema=SCHEMA, version=1, stage='assessment_publication',
        bindings=dict(bindings), counts=counts, assessed_ordinals_sha256=expected_ordinals_sha256,
        schema_sha256=schema_sha, database=dict(name='contexts.sqlite', size_bytes=target.stat().st_size,
            sha256=file_sha(target)))
    content = canonical(manifest)
    with (destination/'manifest.json').open('xb') as stream:
        stream.write(content); stream.flush(); os.fsync(stream.fileno())
    descriptor = os.open(destination, os.O_RDONLY | os.O_DIRECTORY)
    try: os.fsync(descriptor)
    finally: os.close(descriptor)
    return manifest


@dataclass(frozen=True, slots=True)
class CompletedAssessmentInputs:
    """Bound by a prospective public plan and exact compatibility amendment."""
    identity_bytes: bytes
    plan_sha256: str

    def __post_init__(self):
        require(type(self.identity_bytes) is bytes and len(self.identity_bytes) <= 32768
            and digest(self.plan_sha256))
        value = self.identity
        require(set(value) == {'manifest', 'consumer_executable_stack_sha256', 'acceptance_amendment'})
        require(digest(value['consumer_executable_stack_sha256']))
        for key in ('manifest', 'acceptance_amendment'):
            row = value[key]
            require(set(row) == {'path', 'sha256', 'size_bytes'} and digest(row['sha256'])
                and type(row['size_bytes']) is int and 0 < row['size_bytes'] <= 8 * 1024**2
                and Path(row['path']).is_absolute() and '..' not in Path(row['path']).parts)

    @property
    def identity(self):
        value = json.loads(self.identity_bytes)
        require(canonical(value) == self.identity_bytes)
        return value

    @property
    def sha256(self):
        return hashlib.sha256(canonical(dict(inputs=self.identity, plan_sha256=self.plan_sha256))).hexdigest()

    def metadata(self):
        row = self.identity['manifest']; path = Path(row['path']); before = signature(path)
        content = read_bounded_metadata(path)
        require(len(content) == row['size_bytes'] and hashlib.sha256(content).hexdigest() == row['sha256']
            and signature(path) == before)
        value = json.loads(content)
        require(canonical(value) == content and set(value) == {'schema', 'version', 'stage', 'bindings',
            'counts', 'assessed_ordinals_sha256', 'schema_sha256', 'database'}
            and value['schema'] == SCHEMA and value['version'] == 1 and value['stage'] == 'assessment_publication'
            and set(value['bindings']) == BINDINGS and all(digest(v) for v in value['bindings'].values())
            and set(value['counts']) == {'inputs', 'baselines', 'batches', 'contexts', 'decisions', 'public_sources'}
            and all(type(v) is int and v >= 0 for v in value['counts'].values())
            and digest(value['assessed_ordinals_sha256']) and digest(value['schema_sha256'])
            and set(value['database']) == {'name', 'size_bytes', 'sha256'}
            and value['database']['name'] == 'contexts.sqlite' and digest(value['database']['sha256'])
            and type(value['database']['size_bytes']) is int and value['database']['size_bytes'] > 0)
        return value

    def require_authority(self, authority, *, consumer_executable_sha256):
        identity = self.identity
        require(authority.is_finite_profile_pilot
            and authority.binding('public_profile_pilot_assessment_inputs') == self.sha256
            and authority.plan_sha256 == self.plan_sha256
            and authority.executable_stack_sha256 == consumer_executable_sha256
            and identity['consumer_executable_stack_sha256'] == consumer_executable_sha256
            and authority.authorization_receipt_sha256 == identity['acceptance_amendment']['sha256'])
        entry = identity['acceptance_amendment']; content = read_bounded_metadata(entry['path'])
        require(len(content) == entry['size_bytes'] and hashlib.sha256(content).hexdigest() == entry['sha256'])
        compatibility = json.loads(content)['assessment_restart']
        require(compatibility == dict(inputs={k:v for k,v in identity.items() if k != 'acceptance_amendment'},
            producer_bindings=self.metadata()['bindings'], checkpoint_schema=SCHEMA, checkpoint_version=1,
            equivalence='same-scientific-policy-and-complete-domain'))
        for path in (Path(identity['manifest']['path']), Path(identity['manifest']['path']).parent/'contexts.sqlite'):
            signature(path)
            require(path.is_relative_to('/mnt/reference/scratch/public-build')
                and bool(os.statvfs(path).f_flag & os.ST_RDONLY))

    def stage(self, target, *, expected_bindings):
        """Hash during the necessary copy; original snapshot stays immutable."""
        value = self.metadata()
        require(all(value['bindings'].get(k) == v for k,v in expected_bindings.items()))
        source = Path(self.identity['manifest']['path']).parent/'contexts.sqlite'
        before = signature(source); target = Path(target)
        require(before[2] == value['database']['size_bytes']
            and target.is_absolute() and not any(p.is_symlink() for p in (target, *target.parents)))
        result = hashlib.sha256(); total = 0
        with source.open('rb') as stream, target.open('xb') as out:
            os.chmod(target, 0o600)
            for block in iter(lambda: stream.read(1024 * 1024), b''):
                result.update(block); total += len(block); out.write(block)
            out.flush(); os.fsync(out.fileno())
        require(signature(source) == before and total == before[2]
            and result.hexdigest() == value['database']['sha256'])
        return value


def iter_recovered_contexts(inputs, *, admission_reader, profile_reader, assessment_owner,
        scratch_path, execution_config, execution_config_sha256, population_input, models,
        native_contract_sha256, progress=None, science_selector=None):
    """Restart the stage with a private copy, then use the ordinary merge owner."""
    from genosapiens.wgs_contexts import _iter_stored_contexts, _ContextTransactions
    from genosapiens.wgs_context_payloads import Payloads
    from genosapiens.population_policy import ContextKey
    from genosapiens.wgs_runtime import analysis_progress
    require(type(inputs) is CompletedAssessmentInputs)
    expected = dict(admission_completion_sha256=admission_reader.completion.sha256,
        admission_binding_sha256=admission_reader.binding_sha256,
        eligibility_manifest_sha256=admission_reader.eligibility_manifest.sha256,
        profile_manifest_sha256=profile_reader.manifest_sha256,
        association_manifest_sha256=assessment_owner.associations.sha256,
        execution_config_sha256=execution_config_sha256,
        policy_sha256=assessment_owner.policy_sha256, native_contract_sha256=native_contract_sha256,
        parent_source_sha256=admission_reader.parents['source_manifest_sha256'],
        parent_working_sha256=admission_reader.parents['working_manifest_sha256'])
    path = Path(scratch_path)
    require(not path.exists() and not path.is_symlink())
    connection = None
    try:
        analysis_progress(progress, 'completed_input_reuse', 'assessment_snapshot')
        value = inputs.stage(path, expected_bindings=expected)
        require(value['counts']['inputs'] == assessment_owner.associations.subject_count
            and value['assessed_ordinals_sha256'] == assessment_owner.associations.assessed_ordinals_sha256)
        connection = sqlite3.connect(path)
        connection.enable_load_extension(False)
        connection.execute('PRAGMA foreign_keys=ON')
        connection.execute('PRAGMA trusted_schema=OFF')
        connection.execute('PRAGMA cache_size=-65536')
        connection.execute('PRAGMA temp_store=FILE')
        connection.execute('PRAGMA mmap_size=0')
        require(hashlib.sha256(canonical(list(connection.execute(
            'SELECT type,name,tbl_name,sql FROM sqlite_master ORDER BY type,name')))).hexdigest()
            == value['schema_sha256'])
        connection.execute('PRAGMA journal_mode=WAL')
        payloads = Payloads(connection)
        payloads.transactions = _ContextTransactions(connection, path)
        # The producer graph contributes these exact grouping fields for all
        # reached keys, including keys that will not yield a nonempty context.
        for key, content in connection.execute('SELECT context_id,grouping_bytes FROM contexts ORDER BY context_id'):
            require(type(content) is bytes and len(content) <= execution_config.contexts.max_context_bytes)
            row = json.loads(content)
            require(row['context_id'] == key)
            if row['native_disease_id'] is not None and row['inheritance'] is not None:
                for model in models:
                    population_input.add(ContextKey(row['native_disease_id'], row['gene_id'],
                        row['locus_id'], row['inheritance'], model))
        analysis_progress(progress, 'context_assembly', 'assessment_publication')
        assessment_owner.publish(connection, payloads, admission_reader)
        analysis_progress(progress, 'context_assembly', 'context_merge',
            counts={'contexts_obtained':value['counts']['contexts'], 'contexts_merged':0, 'contexts_emitted':0})
        yield from _iter_stored_contexts(connection, payloads, admission_reader, execution_config,
            scratch_path=path, progress=progress, total_contexts=value['counts']['contexts'], science_selector=science_selector)
    finally:
        if connection is not None:
            connection.close()
        for owned in (path, path.with_name(path.name+'-wal'), path.with_name(path.name+'-shm')):
            owned.unlink(missing_ok=True)
