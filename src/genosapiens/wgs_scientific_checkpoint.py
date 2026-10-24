"""Explicit finite-pilot inputs at the sealed scientific publication boundary.

The snapshot is an input inventory, never a validated scientific output. Its
consumer must bind completed production and run the required publication checks.
"""
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import time

from genosapiens.wgs_ledger import canonical, publish_bytes
from genosapiens.wgs_shards import read_bounded_metadata
from genosapiens.wgs_scientific_artifact import ScientificArtifact, scientific_artifact_from_bytes, scientific_artifact_complete

SCHEMA = 'genosapiens.complete-scientific-input'
COMPLETED_CONSUMPTION_POLICY = 'completed-science-publication-v1'
_REUSE_ISSUER = object()
MAX_METADATA = 8 * 1024**2
PRODUCER_KEYS = frozenset(('candidate_sha256', 'plan_sha256', 'executable_stack_sha256',
    'runtime_image_sha256', 'authority_sha256', 'release_sha256'))


class ScientificCheckpointError(ValueError):
    pass


@dataclass(frozen=True, slots=True, init=False)
class CompletedScientificReuse:
    """Operation-local producer authority, never a persisted validation flag."""
    scientific_sha256: str
    producer_binding: str
    consumer_binding: str

    def __init__(self, scientific_sha256, producer_binding, consumer_binding, *, _issuer=None):
        require(_issuer is _REUSE_ISSUER and all(digest(v) for v in
            (scientific_sha256, producer_binding, consumer_binding)))
        object.__setattr__(self, 'scientific_sha256', scientific_sha256)
        object.__setattr__(self, 'producer_binding', producer_binding)
        object.__setattr__(self, 'consumer_binding', consumer_binding)

    def require_artifact(self, artifact):
        require(scientific_artifact_complete(artifact) and artifact.sha256 == self.scientific_sha256)

    def cover_producer_roots(self, reader, artifact):
        """Reuse exact completed DAG coverage; do not claim semantic replay."""
        from genosapiens.wgs_evidence_store import EvidenceReader
        self.require_artifact(artifact)
        require(type(reader) is EvidenceReader and reader.manifest == artifact.evidence
            and reader.verified and reader._root_coverage and reader._read_scope_depth > 0)
        require(reader.connection.execute('SELECT count(*) FROM roots').fetchone()[0] == reader.manifest.root_count)
        reader.connection.execute('INSERT OR IGNORE INTO normalized_completed_root_coverage SELECT name FROM roots')


def _completed_producer_reuse(artifact, *, consumer_executable_sha256):
    """Called by scientific production only after its final descriptor is sealed."""
    require(scientific_artifact_complete(artifact))
    return CompletedScientificReuse(artifact.sha256, artifact.sha256, consumer_executable_sha256, _issuer=_REUSE_ISSUER)


def _completed_product_reuse(artifact, binding):
    """Issued only after existing committed product-checkpoint authentication."""
    from genosapiens.wgs_ledger import StageRunBinding
    require(scientific_artifact_complete(artifact) and type(binding) is StageRunBinding and binding.stage == 'bounded_products')
    return CompletedScientificReuse(artifact.sha256, binding.sha256,
        binding.executable_stack_sha256, _issuer=_REUSE_ISSUER)


def require(value):
    if not value:
        raise ScientificCheckpointError('complete scientific input is invalid')


def digest(value):
    return type(value) is str and len(value) == 64 and all(c in '0123456789abcdef' for c in value)


def signature(path):
    path = Path(path)
    require(path.is_absolute() and '..' not in path.parts
        and not any(p.is_symlink() for p in (path, *path.parents)) and path.is_file())
    value = path.stat()
    return value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns


def retained_files(scientific, directory, *, eligibility_root=None, profile_root=None):
    """Enumerate sealed leaves, including selected-index nodes and dependencies.

    Context version 2 binds its selected reconstruction root into the ordinary
    evidence store. Keeping every manifest shard therefore retains the index
    and its referenced member/outcome arrays without an extra store traversal.
    """
    from genosapiens.wgs_analysis_eligibility import scientific_eligibility_component, eligibility_retained_files
    from genosapiens.wgs_profile_projection import scientific_profile_component, profile_retained_files
    from genosapiens.wgs_association_projection import assessment_retained_files
    require(scientific_artifact_complete(scientific))
    directory = Path(directory)
    yield (directory / 'scientific-artifact.json', 'scientific/scientific-artifact.json',
        scientific.sha256, len(scientific.as_bytes()))
    for shard in (*scientific.evidence.node_shards, *scientific.evidence.root_shards):
        yield directory / shard.relative_path, 'scientific/' + shard.relative_path, shard.sha256, shard.size_bytes
    eligibility = scientific_eligibility_component(scientific, directory, eligibility_root=eligibility_root)
    profile = scientific_profile_component(scientific, directory, profile_root=profile_root)
    require(eligibility is not None and profile is not None)
    yield from eligibility_retained_files(directory.parent, eligibility, eligibility_root=eligibility_root)
    yield from profile_retained_files(directory.parent, profile, profile_root=profile_root)
    yield from assessment_retained_files(directory.parent, dict(scientific.parent_bindings))


def _copy(source, target, expected_sha256, size, deadline):
    before = signature(source)
    require(before[2] == size and time.monotonic() < deadline)
    target.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
    result, count = hashlib.sha256(), 0
    fd = os.open(source, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    with os.fdopen(fd, 'rb') as reader:
        value = os.fstat(reader.fileno())
        require((value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns) == before)
        fd = os.open(target, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600)
        with os.fdopen(fd, 'wb') as writer:
            while block := reader.read(min(1024**2, size - count + 1)):
                count += len(block)
                require(count <= size and time.monotonic() < deadline)
                result.update(block)
                writer.write(block)
            require(count == size and result.hexdigest() == expected_sha256 and signature(source) == before)
            writer.flush()
            os.fsync(writer.fileno())


def seal_snapshot(destination, *, scientific, directory, producer, eligibility_root=None,
                  profile_root=None, deadline_seconds=1800):
    """Copy exact sealed leaves; never hardlink a shard held by an active reader."""
    require(type(producer) is dict and set(producer) == PRODUCER_KEYS
        and all(digest(value) for value in producer.values()) and 0 < deadline_seconds <= 7200)
    destination = Path(destination)
    require(destination.is_absolute() and '..' not in destination.parts
        and not any(p.is_symlink() for p in (destination, *destination.parents)))
    destination.mkdir(mode=0o700)
    deadline = time.monotonic() + deadline_seconds
    try:
        rows = []
        for source, relative, sha256, size in retained_files(scientific, directory,
                eligibility_root=eligibility_root, profile_root=profile_root):
            _copy(source, destination / relative, sha256, size, deadline)
            rows.append(dict(path=relative, sha256=sha256, size_bytes=size))
        rows.sort(key=lambda row: row['path'])
        require(len({row['path'] for row in rows}) == len(rows))
        content = canonical(dict(schema=SCHEMA, version=1, stage='scientific_publication',
            producer=producer, scientific_sha256=scientific.sha256, files=rows))
        require(len(content) <= MAX_METADATA)
        publish_bytes(destination, 'manifest.json', content)
        return hashlib.sha256(content).hexdigest()
    except BaseException as error:
        from genosapiens.wgs_failures import run_cleanup
        run_cleanup((('context_store', lambda: shutil.rmtree(destination)),), primary=error)
        raise


@dataclass(frozen=True, slots=True)
class CompletedScientificInputs:
    identity_bytes: bytes
    plan_sha256: str

    def __post_init__(self):
        require(type(self.identity_bytes) is bytes and len(self.identity_bytes) <= 32768 and digest(self.plan_sha256))
        value = self.identity
        require(set(value) == {'manifest', 'consumer_executable_stack_sha256', 'acceptance_amendment'}
            and digest(value['consumer_executable_stack_sha256']))
        for key in ('manifest', 'acceptance_amendment'):
            row = value[key]
            require(type(row) is dict and set(row) == {'path', 'sha256', 'size_bytes'}
                and digest(row['sha256']) and type(row['size_bytes']) is int and 0 < row['size_bytes'] <= MAX_METADATA
                and type(row['path']) is str and Path(row['path']).is_absolute() and '..' not in Path(row['path']).parts)
        require(Path(value['manifest']['path']).name == 'manifest.json')

    @property
    def identity(self):
        value = json.loads(self.identity_bytes)
        require(canonical(value) == self.identity_bytes)
        return value

    @property
    def sha256(self):
        return hashlib.sha256(canonical(dict(inputs=self.identity, plan_sha256=self.plan_sha256))).hexdigest()

    def _bound(self, row):
        before = signature(row['path'])
        content = read_bounded_metadata(row['path'])
        require(len(content) == row['size_bytes'] and hashlib.sha256(content).hexdigest() == row['sha256']
            and signature(row['path']) == before)
        return content

    def metadata(self):
        content = self._bound(self.identity['manifest'])
        value = json.loads(content)
        require(canonical(value) == content and set(value) == {'schema', 'version', 'stage', 'producer', 'scientific_sha256', 'files'}
            and value['schema'] == SCHEMA and type(value['version']) is int and value['version'] == 1
            and value['stage'] == 'scientific_publication' and digest(value['scientific_sha256'])
            and type(value['producer']) is dict and set(value['producer']) == PRODUCER_KEYS
            and all(digest(v) for v in value['producer'].values())
            and type(value['files']) is list and value['files'])
        paths = []
        for row in value['files']:
            require(type(row) is dict and set(row) == {'path', 'sha256', 'size_bytes'}
                and digest(row['sha256']) and type(row['size_bytes']) is int and 0 < row['size_bytes'] < 2**63
                and type(row['path']) is str)
            path = PurePosixPath(row['path'])
            require(str(path) == row['path'] and not path.is_absolute() and '..' not in path.parts
                and '\\' not in row['path'] and len(path.parts) >= 2
                and path.parts[0] in {'scientific', 'eligibility', 'analysis', 'associations'})
            paths.append(row['path'])
        require(paths == sorted(set(paths)) and 'scientific/scientific-artifact.json' in paths)
        return value

    def require_authority(self, authority, *, consumer_executable_sha256):
        value = self.identity
        require(authority.is_finite_profile_pilot
            and authority.binding('public_profile_pilot_scientific_inputs') == self.sha256
            and authority.plan_sha256 == self.plan_sha256
            and authority.executable_stack_sha256 == consumer_executable_sha256
            and value['consumer_executable_stack_sha256'] == consumer_executable_sha256
            and authority.authorization_receipt_sha256 == value['acceptance_amendment']['sha256'])
        amendment = json.loads(self._bound(value['acceptance_amendment']))
        metadata = self.metadata()
        expected = dict(
            inputs={key: item for key, item in value.items() if key != 'acceptance_amendment'},
            producer=metadata['producer'], checkpoint_schema=SCHEMA, checkpoint_version=1,
            equivalence='same-scientific-policy-and-complete-domain')
        restart = amendment['scientific_restart']
        if 'consumption_policy' in restart:
            expected['consumption_policy'] = COMPLETED_CONSUMPTION_POLICY
        require(restart == expected)
        path = Path(value['manifest']['path'])
        require(path.is_relative_to('/mnt/reference/scratch/public-build') and bool(os.statvfs(path).f_flag & os.ST_RDONLY))
        for row in metadata['files']:
            leaf = path.parent / row['path']
            require(signature(leaf)[2] == row['size_bytes'] and bool(os.statvfs(leaf).f_flag & os.ST_RDONLY))

    def completed_reuse(self, authority, *, consumer_executable_sha256, artifact):
        """Issue only under the prospectively bound completed-consumer contract."""
        self.require_authority(authority, consumer_executable_sha256=consumer_executable_sha256)
        self.require_completed_publication()
        metadata = self.metadata()
        require(scientific_artifact_complete(artifact) and metadata['scientific_sha256'] == artifact.sha256)
        return CompletedScientificReuse(artifact.sha256, self.identity['manifest']['sha256'],
            consumer_executable_sha256, _issuer=_REUSE_ISSUER)

    def require_completed_publication(self):
        amendment = json.loads(self._bound(self.identity['acceptance_amendment']))
        require(amendment['scientific_restart'].get('consumption_policy') == COMPLETED_CONSUMPTION_POLICY)

    def stage(self, destination, *, deadline_seconds=1800):
        """Create fresh private inputs; the normal reader owns scientific replay."""
        require(0 < deadline_seconds <= 7200)
        metadata = self.metadata()
        source = Path(self.identity['manifest']['path']).parent
        destination = Path(destination)
        require(destination.is_absolute() and '..' not in destination.parts
            and not any(p.is_symlink() for p in (destination, *destination.parents)))
        destination.mkdir(mode=0o700)
        deadline = time.monotonic() + deadline_seconds
        try:
            for row in metadata['files']:
                _copy(source / row['path'], destination / row['path'], row['sha256'], row['size_bytes'], deadline)
            content = read_bounded_metadata(destination / 'scientific/scientific-artifact.json')
            scientific = scientific_artifact_from_bytes(content, expected_sha256=metadata['scientific_sha256'])
            expected = sorted((dict(path=relative, sha256=sha256, size_bytes=size)
                for _, relative, sha256, size in retained_files(scientific, destination / 'scientific')),
                key=lambda row: row['path'])
            require(expected == metadata['files'])
            return scientific
        except BaseException as error:
            from genosapiens.wgs_failures import run_cleanup
            run_cleanup((('context_store', lambda: shutil.rmtree(destination)),), primary=error)
            raise

    def stage_work_root(self, work_root, *, deadline_seconds=1800):
        """Keep the ordinary product-checkpoint layout and existing telemetry."""
        root = Path(work_root)
        require(root.is_absolute() and root.is_dir() and '..' not in root.parts
            and not any(path.is_symlink() for path in (root, *root.parents)))
        names = ('scientific', 'eligibility', 'analysis', 'associations')
        require(all(not (root / name).exists() and not (root / name).is_symlink() for name in names))
        staging = root / '.scientific-input-staging'
        scientific = self.stage(staging, deadline_seconds=deadline_seconds)
        moved = []
        try:
            require({path.name for path in staging.iterdir()} == set(names))
            for name in names:
                os.rename(staging / name, root / name)
                moved.append(root / name)
            staging.rmdir()
            return scientific
        except BaseException as error:
            from genosapiens.wgs_failures import run_cleanup
            paths = (*moved, *((staging,) if staging.exists() else ()))
            run_cleanup(tuple(('context_store', lambda path=path: shutil.rmtree(path)) for path in paths), primary=error)
            raise
