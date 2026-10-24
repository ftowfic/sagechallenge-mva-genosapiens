"""Explicit, byte-preserving migration of stopped operator checkpoints.

The inventory authorizes exact identity substitutions. It never rewrites a
producer descriptor, a scientific binding, or a stage database. Ordinary new
commits continue to use their own publication identities.
"""

from contextlib import ExitStack, closing, contextmanager
import fcntl
import hashlib
import json
import os
from pathlib import Path
import sqlite3
import stat
import time
import uuid

from genosapiens.wgs_ledger import canonical


RECORD = '.operator/output-migration.json'
INCOMPLETE = '.track1-migration-incomplete'
SCHEMA = 'genosapiens.track1-output-migration'
CHUNK = 1024 * 1024


class MigrationError(RuntimeError):
    """Static operator error that does not expose case contents."""


def _require(ok, message='migration checkpoint or file binding is invalid'):
    if not ok:
        raise MigrationError(message)


def _signature(value):
    return (value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns)


def _directory(path):
    value = path.lstat()
    _require(stat.S_ISDIR(value.st_mode) and value.st_uid == os.getuid()
             and stat.S_IMODE(value.st_mode) == 0o700)
    return value


def _file(path):
    value = path.lstat()
    _require(stat.S_ISREG(value.st_mode) and value.st_uid == os.getuid()
             and value.st_nlink == 1 and stat.S_IMODE(value.st_mode) in (0o400, 0o500, 0o600))
    return value


def _identity(path, kind):
    value = _file(path)
    if kind == 'file':
        return _signature(value)
    _require(kind == 'publication' and stat.S_IMODE(value.st_mode) == 0o600)
    parent = _directory(path.parent)
    return (parent.st_dev, parent.st_ino, *_signature(value))


def _small(path, limit=65536):
    before = _file(path)
    _require(before.st_size <= limit)
    fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    with os.fdopen(fd, 'rb') as stream:
        _require(_signature(os.fstat(stream.fileno())) == _signature(before))
        content = stream.read(limit + 1)
        _require(len(content) == before.st_size
                 and _signature(os.fstat(stream.fileno())) == _signature(before))
    _require(_signature(_file(path)) == _signature(before))
    return content


def _json(content):
    def unique(pairs):
        value = dict(pairs)
        _require(len(value) == len(pairs))
        return value
    value = json.loads(content, object_pairs_hook=unique)
    _require(type(value) is dict and canonical(value) == content)
    return value


def _sha(content):
    return hashlib.sha256(content).hexdigest()


def _sync(path):
    fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def _write(path, content, mode=0o600):
    fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, mode)
    with os.fdopen(fd, 'wb') as stream:
        stream.write(content)
        stream.flush()
        os.fsync(stream.fileno())


def _walk(root):
    """Stream directory entries without retaining a file list in memory."""
    with os.scandir(root) as entries:
        for entry in entries:
            path = Path(entry.path)
            value = path.lstat()
            if stat.S_ISDIR(value.st_mode):
                _directory(path)
                yield path, True
                yield from _walk(path)
            else:
                _file(path)
                yield path, False


def _database(path):
    _file(path)
    for suffix in ('-journal', '-wal', '-shm'):
        sidecar = Path(str(path) + suffix)
        _require(not sidecar.exists() and not sidecar.is_symlink(),
                 'database recovery is required at the original location before migration')
    db = sqlite3.connect(path.as_uri() + '?mode=ro&immutable=1', uri=True)
    db.enable_load_extension(False)
    db.execute('PRAGMA query_only=ON')
    db.execute('PRAGMA trusted_schema=OFF')
    db.execute('PRAGMA cache_size=-4096')
    db.execute('PRAGMA mmap_size=0')
    return db


@contextmanager
def _lock(path):
    before = _file(path)
    fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    try:
        _require(_signature(os.fstat(fd)) == _signature(before))
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise MigrationError('the run is active; stop its container before migration') from None
        _require(_signature(_file(path)) == _signature(before))
        yield
    finally:
        os.close(fd)


class _MigrationMap:
    def __init__(self, root):
        self.root = root
        self.record_path = root / RECORD
        self.record_signature = _signature(_file(self.record_path))
        value = _json(_small(self.record_path))
        _require(set(value) == {'schema', 'version', 'status', 'case_sha256',
            'package_sha256', 'root_identity', 'inventory', 'inventory_identity',
            'inventory_sha256', 'files', 'bytes'})
        _require(value['schema'] == SCHEMA and value['version'] == 1 and value['status'] == 'complete')
        relative = Path(value['inventory'])
        _require(len(relative.parts) == 4 and relative.parts[:2] == ('.operator', 'migrations')
                 and relative.name == 'inventory.sqlite'
                 and relative.parts[2].startswith('migration-')
                 and relative.parts[2] not in ('.', '..'))
        self.inventory = root / relative
        self.value = value
        self.stable()
        with closing(_database(self.inventory)) as db:
            _require(db.execute('PRAGMA application_id').fetchone() == (1196648786,))
            _require(db.execute('PRAGMA user_version').fetchone() == (1,))
            _require(db.execute('SELECT case_sha256 FROM binding').fetchall() == [(value['case_sha256'],)])

    def stable(self):
        root_stat = _directory(self.root)
        _require(not (self.root / INCOMPLETE).exists())
        _require([root_stat.st_dev, root_stat.st_ino] == self.value['root_identity']
                 and _signature(_file(self.record_path)) == self.record_signature
                 and list(_signature(_file(self.inventory))) == self.value['inventory_identity']
                 and _sha(_small(self.root / '.operator/case.json')) == self.value['case_sha256'])
        for parent in self.inventory.parents:
            _directory(parent)
            if parent == self.root:
                break

    def matches(self, path, recorded, kind, sha256, size_bytes):
        self.stable()
        relative = path.relative_to(self.root).as_posix()
        with closing(_database(self.inventory)) as db:
            row = db.execute('SELECT f.destination_identity,f.sha256,f.size_bytes,d.destination_identity '
                'FROM substitutions s JOIN files f ON f.path=s.path '
                'JOIN directories d ON d.path=s.parent '
                'WHERE s.path=? AND s.kind=? AND s.recorded=?',
                (relative, kind, canonical(list(recorded)))).fetchone()
        if row is None:
            return False
        identity, digest, size, directory = row
        actual = _identity(path, kind)
        expected = tuple(json.loads(identity))
        if kind == 'publication':
            expected = (*json.loads(directory)[:2], *expected)
        self.stable()
        return (actual == expected and (sha256 is None or sha256 == digest)
                and (size_bytes is None or size_bytes == size))


_MAPS = {}


def _map(root):
    # Bound cached controls to immutable file identities. Never cache failed reads.
    existing = _MAPS.get(root)
    if existing is not None:
        existing.stable()
        return existing
    value = _MigrationMap(root)
    if len(_MAPS) >= 8:
        _MAPS.clear()
    _MAPS[root] = value
    return value


def require_migration_ready(root):
    root = Path(root)
    _require(not (root / INCOMPLETE).exists() and not (root / INCOMPLETE).is_symlink(),
             'output migration is incomplete; retain it and use the original run')
    path = root / RECORD
    if path.exists() or path.is_symlink():
        _map(root).stable()


def matches_migrated_identity(path, recorded, *, kind='file', sha256=None, size_bytes=None):
    """Accept only an exact original identity or a checked destination substitute."""
    path = Path(path)
    try:
        _require(path.is_absolute() and path.resolve() == path)
        _require(type(recorded) in (list, tuple) and len(recorded) == (5 if kind == 'file' else 7))
        if _identity(path, kind) == tuple(recorded):
            return True
        for parent in path.parents:
            if (parent / '.operator').is_dir():
                require_migration_ready(parent)
                if (parent / RECORD).exists():
                    return _map(parent).matches(path, recorded, kind, sha256, size_bytes)
                return False
        return False
    except (OSError, ValueError, sqlite3.Error, MigrationError, KeyError, TypeError):
        return False


def _check_stage_head(stage_root, db, binding, latest, event):
    """Validate a current head or one pending publication without writing it.

StageCommitState commits SQLite before publishing the event and head. Its
ordinary destination recovery publishes this exact pending event after copying.
"""
    from genosapiens.wgs_ledger import JournalEntry, load_model

    def pointer(entry):
        return dict(schema='genosapiens.wgs-stage-journal-head', version=1,
                    stage=binding.stage, ordinal=entry.ordinal, sha256=entry.sha256)

    def exists(path):
        return path.exists() or path.is_symlink()

    head_path = stage_root / f'{binding.stage}-head.json'
    head = _json(_small(head_path, 4096)) if exists(head_path) else None
    if head is not None:
        _require(type(head.get('version')) is int and type(head.get('ordinal')) is int,
                 'stage journal head has invalid fields')
    next_path = stage_root / f'{binding.stage}-{event.ordinal + 1:06d}.journal.json'
    _require(not exists(next_path), 'stage journal is ahead of its database checkpoint')
    pending = head != pointer(event)
    if pending:
        if event.ordinal == 1:
            _require(head is None, 'first stage event has a conflicting journal head')
        else:
            previous = db.execute('SELECT ordinal,event_bytes,event_sha256 FROM events WHERE ordinal=?',
                                  (event.ordinal - 1,)).fetchone()
            _require(previous is not None, 'stage journal predecessor is missing')
            predecessor = load_model(previous[1], JournalEntry)
            _require(predecessor.ordinal == event.ordinal - 1 == previous[0]
                     and predecessor.sha256 == previous[2] == event.previous_entry_sha256
                     and predecessor.stage == binding.stage
                     and predecessor.stage_parent_sha256 == binding.stage_parent_sha256
                     and predecessor.parents == binding.parents
                     and head == pointer(predecessor),
                     'stage journal head is not current or one recoverable event behind')
            previous_path = stage_root / f'{binding.stage}-{predecessor.ordinal:06d}.journal.json'
            _require(_small(previous_path) == previous[1],
                     'stage journal predecessor differs from its database checkpoint')
    latest_path = stage_root / f'{binding.stage}-{event.ordinal:06d}.journal.json'
    if exists(latest_path):
        _require(_small(latest_path) == latest[1],
                 'stage journal event differs from its database checkpoint')
    else:
        _require(pending, 'current stage journal head has no matching journal event')


def _stage_substitutions(source, stage_root, state_path, inventory):
    """Check persisted commitments without opening a stage for writing."""
    from genosapiens.wgs_ledger import StageRunBinding, JournalEntry, load_model
    with closing(_database(state_path)) as db:
        binding = load_model(db.execute('SELECT canonical_bytes FROM run_binding WHERE singleton=1').fetchone()[0], StageRunBinding)
        latest = db.execute('SELECT ordinal,event_bytes,event_sha256 FROM events ORDER BY ordinal DESC LIMIT 1').fetchone()
        _require(latest is not None)
        event = load_model(latest[1], JournalEntry)
        _require(event.ordinal == latest[0] and event.sha256 == latest[2]
            and event.stage == binding.stage and event.stage_parent_sha256 == binding.stage_parent_sha256
            and event.parents == binding.parents)
        _check_stage_head(stage_root, db, binding, latest, event)
        for name, commit, content in db.execute('SELECT name,commit_ordinal,metadata_bytes FROM committed_outputs'):
            data = _json(content)
            _require(set(data) == {'name', 'sha256', 'size_bytes', 'publication'}
                     and data['name'] == name and Path(name).name == name and name not in ('.', '..'))
            saved = db.execute('SELECT outputs_bytes FROM commits WHERE ordinal=?', (commit,)).fetchone()
            _require(saved is not None and name in json.loads(saved[0]))
            _substitution(source, stage_root / name, data['publication'], 'publication', inventory,
                          data['sha256'], data['size_bytes'])


def _substitution(source, path, recorded, kind, inventory, digest=None, size=None):
    _require(path.is_relative_to(source))
    _require(matches_migrated_identity(path, recorded, kind=kind, sha256=digest, size_bytes=size))
    inventory.execute('INSERT INTO substitutions VALUES (?,?,?,?,?,?)',
        (path.relative_to(source).as_posix(), path.parent.relative_to(source).as_posix(),
         kind, canonical(list(recorded)), digest, size))


def _operator_preflight(source, inventory, package_sha, package_root):
    """Authenticate original controls without reading VCFs or public resources."""
    from genosapiens.wgs_operator_invocation import _package, CASE_SCHEMA
    from genosapiens.wgs_operator_maintenance import compatible_package_binding
    from genosapiens.wgs_preparation import load_completed_source_input, resolve_completed_source
    from genosapiens.wgs_ledger import load_model, PreparedCaseManifest
    case_content = _small(source / '.operator/case.json')
    case = _json(case_content)
    contract, release, _, _, _ = _package(package_root, package_sha)
    _require(case['schema'] == CASE_SCHEMA and case['version'] == 2
             and case['runtime_uid'] == os.getuid() and case['runtime_gid'] == os.getgid()
             and case['output_root'] == str(source)
             and compatible_package_binding(package_root, case['package_bindings_sha256'], package_sha))
    _require(_sha(_small(Path(case['hpo_input_path']))) == case['hpo_input_sha256']
             and _sha(_small(Path(case['concordance_path']))) == case['concordance_sha256'])
    receipt = _json(_small(source / 'preparation/controlled-preparation-result.json'))
    _require(receipt['status'] == 'passed' and receipt['base_release_sha256'] == release.sha256,
             'migration requires completed preparation')
    manifest = load_model(_small(source / 'preparation/prepared/prepared-case-manifest.json'), PreparedCaseManifest)
    path = Path(manifest.completed_source_input_path)
    _require(path.is_relative_to(source))
    completed = load_completed_source_input(_small(path))
    _require(completed.sha256 == manifest.completed_source_input.sha256
             and completed.producer.source.input_sha256 == case['input_sha256']
             and completed.producer.runtime_image_digest == contract.runtime_image_digest)
    for path in (source / 'preparation/state/source.sqlite', source / 'preparation/state/source-audit.sqlite',
                 source / 'preparation/state/filter.sqlite', source / 'analysis/analysis-work/state/annotation.sqlite',
                 source / 'analysis/analysis-work/state/retained.sqlite'):
        if path.exists():
            with closing(_database(path)):
                pass
    resolve_completed_source(completed, expected_source=completed.producer.source,
                             audit_requested=completed.audit_requested)
    for state in (completed.source_state, completed.source_audit_state):
        if state is not None:
            _substitution(source, Path(completed.source_state_root) / state.relative_path,
                          state.file_identity, 'file', inventory)
    _require_supported_checkpoint(source)
    work = source / 'analysis/analysis-work'
    for root, state in (
        (source / 'preparation/source', source / 'preparation/state/source.sqlite'),
        (source / 'preparation/source', source / 'preparation/state/source-audit.sqlite'),
        (source / 'preparation/prefilter', source / 'preparation/state/filter.sqlite'),
        (work / 'annotation', work / 'state/annotation.sqlite'),
    ):
        if state.exists():
            _stage_substitutions(source, root, state, inventory)
    return _sha(case_content)


def _require_supported_checkpoint(source):
    work = source / 'analysis/analysis-work'
    for name in ('state/admission.sqlite', 'state/profile-intake.sqlite'):
        _require(not (work / name).exists(), 'migration currently supports runs before admission starts')
    for name in ('scientific', 'profile', 'assessment', 'admission'):
        path = work / name
        _require(not path.exists() or path.is_dir() and not any(path.iterdir()),
                 'migration currently supports preparation and annotation checkpoints only')


def _copy_file(source, target, expected, *, progress=None):
    before = _file(source)
    _require(_signature(before) == tuple(expected))
    digest = hashlib.sha256()
    fd = os.open(source, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    with os.fdopen(fd, 'rb') as reader:
        _require(_signature(os.fstat(reader.fileno())) == tuple(expected))
        fd = os.open(target, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600)
        with os.fdopen(fd, 'wb') as writer:
            size = 0
            while block := reader.read(CHUNK):
                size += len(block)
                _require(size <= before.st_size)
                digest.update(block)
                writer.write(block)
                if progress:
                    progress(size)
            _require(size == before.st_size and _signature(os.fstat(reader.fileno())) == tuple(expected))
            writer.flush()
            os.fchmod(writer.fileno(), stat.S_IMODE(before.st_mode))
            os.fsync(writer.fileno())
    os.utime(target, ns=(before.st_atime_ns, before.st_mtime_ns), follow_symlinks=False)
    _sync(target)
    _require(_signature(_file(source)) == tuple(expected))
    return digest.hexdigest()


def _hash_file(path, *, progress=None):
    before = _file(path)
    digest = hashlib.sha256()
    fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    with os.fdopen(fd, 'rb') as stream:
        _require(_signature(os.fstat(stream.fileno())) == _signature(before))
        size = 0
        while block := stream.read(CHUNK):
            digest.update(block)
            size += len(block)
            if progress:
                progress(size)
        _require(_signature(os.fstat(stream.fileno())) == _signature(before))
    _require(_signature(_file(path)) == _signature(before))
    return digest.hexdigest()


def _publish_directory(source, destination):
    """Linux no-replace rename prevents clobbering a concurrent destination."""
    import ctypes
    library = ctypes.CDLL(None, use_errno=True)
    rename = getattr(library, 'renameat2', None)
    _require(rename is not None, 'this platform lacks atomic migration publication')
    rename.argtypes = (ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_uint)
    rename.restype = ctypes.c_int
    if rename(-100, os.fsencode(source), -100, os.fsencode(destination), 1) != 0:
        raise OSError(ctypes.get_errno(), 'migration destination publication failed')


def migrate_output(source, destination, *, package_sha, package_root, progress=None):
    """Copy one stopped output. Failed staging remains explicitly incomplete."""
    source, destination = Path(source), Path(destination)
    _require(source.is_absolute() and destination.is_absolute()
             and source.resolve() == source and destination.resolve() == destination
             and not destination.is_relative_to(source) and not source.is_relative_to(destination))
    _directory(source)
    _require(destination.parent.is_dir() and not destination.exists() and not destination.is_symlink(),
             'migration requires a new destination in an existing directory')
    require_migration_ready(source)
    identity = 'migration-' + uuid.uuid4().hex
    staging = destination.parent / ('.' + destination.name + '-' + identity + '-incomplete')
    with ExitStack() as locks:
        locks.enter_context(_lock(source / '.operator/invocation.lock'))
        # Also exclude orphaned stage writers after a supervisor dies.
        for path, directory in _walk(source):
            if not directory and path.suffix == '.lock' and path != source / '.operator/invocation.lock':
                locks.enter_context(_lock(path))
        staging.mkdir(mode=0o700)
        _write(staging / INCOMPLETE, canonical(dict(schema=SCHEMA, status='incomplete')))
        control = staging / '.operator/migrations' / identity
        for path in (staging / '.operator', staging / '.operator/migrations', control):
            path.mkdir(mode=0o700)
        inventory_path = control / 'inventory.sqlite'
        _write(inventory_path, b'')
        with closing(sqlite3.connect(inventory_path)) as inventory:
            inventory.executescript('''
                PRAGMA application_id=1196648786;
                PRAGMA user_version=1;
                PRAGMA journal_mode=DELETE;
                PRAGMA synchronous=FULL;
                PRAGMA cache_size=-4096;
                CREATE TABLE binding(case_sha256 TEXT NOT NULL);
                CREATE TABLE files(path TEXT PRIMARY KEY, sha256 TEXT,
                    size_bytes INTEGER NOT NULL, source_identity BLOB NOT NULL, destination_identity BLOB);
                CREATE TABLE directories(path TEXT PRIMARY KEY, source_identity BLOB NOT NULL, destination_identity BLOB);
                CREATE TABLE substitutions(path TEXT, parent TEXT, kind TEXT, recorded BLOB, sha256 TEXT, size_bytes INTEGER,
                    PRIMARY KEY(path,kind,recorded));
            ''')
            root_before = _signature(_directory(source))
            inventory.execute('INSERT INTO directories VALUES (?,?,?)',
                ('.', canonical(root_before), canonical(_signature(_directory(staging)))))
            # Capture identities before checkpoint authentication. In particular,
            # a source database must not change between attestation and copying.
            for path, directory in _walk(source):
                relative = path.relative_to(source).as_posix()
                before = _signature(path.lstat())
                if directory:
                    inventory.execute('INSERT INTO directories VALUES (?,?,NULL)',
                        (relative, canonical(before)))
                else:
                    inventory.execute('INSERT INTO files VALUES (?,NULL,?,?,NULL)',
                        (relative, before[2], canonical(before)))
            inventory.commit()
            case_sha = _operator_preflight(source, inventory, package_sha, package_root)
            inventory.execute('INSERT INTO binding VALUES (?)', (case_sha,))
            count = total = 0
            last_progress = 0
            def report(phase, amount=0, force=False):
                nonlocal last_progress
                if progress and (force or time.monotonic() - last_progress >= 5):
                    progress(phase, count, total + amount)
                    last_progress = time.monotonic()
            for path, directory in _walk(source):
                relative = path.relative_to(source).as_posix()
                target = staging / relative
                table = 'directories' if directory else 'files'
                saved = inventory.execute(f'SELECT source_identity FROM {table} WHERE path=?', (relative,)).fetchone()
                _require(saved is not None)
                before = tuple(json.loads(saved[0]))
                _require(_signature(path.lstat()) == before)
                if directory:
                    target.mkdir(mode=0o700, exist_ok=True)
                    inventory.execute('UPDATE directories SET destination_identity=? WHERE path=?',
                        (canonical(_signature(_directory(target))), relative))
                    continue
                digest = _copy_file(path, target, before, progress=lambda size: report('copying', size))
                _require(_hash_file(target, progress=lambda size: report('checking', size)) == digest,
                         'migration destination content differs from the original')
                inventory.execute('UPDATE files SET sha256=?,destination_identity=? WHERE path=?',
                    (digest, canonical(_signature(_file(target))), relative))
                for expected_hash, expected_size in inventory.execute(
                        'SELECT sha256,size_bytes FROM substitutions WHERE path=?', (relative,)):
                    _require((expected_hash is None or expected_hash == digest)
                             and (expected_size is None or expected_size == before[2]))
                count += 1
                total += before[2]
                if count % 128 == 0:
                    inventory.commit()
                report('checked')
            # Compare the full inventory again to detect additions, removals and
            # changes while copying. This reads metadata, not scientific bytes.
            observed = observed_directories = 0
            for path, directory in _walk(source):
                relative = path.relative_to(source).as_posix()
                table = 'directories' if directory else 'files'
                row = inventory.execute(f'SELECT source_identity FROM {table} WHERE path=?', (relative,)).fetchone()
                _require(row is not None and tuple(json.loads(row[0])) == _signature(path.lstat()))
                if not directory:
                    observed += 1
                else:
                    observed_directories += 1
            _require(observed == count and _signature(_directory(source)) == root_before)
            _require(inventory.execute('SELECT count(*) FROM files').fetchone() == (observed,)
                     and inventory.execute('SELECT count(*) FROM directories').fetchone() == (observed_directories + 1,))
            _require(inventory.execute('SELECT count(*) FROM substitutions s LEFT JOIN files f ON s.path=f.path WHERE f.path IS NULL').fetchone() == (0,))
            for relative, identity in inventory.execute('SELECT path,destination_identity FROM files'):
                _require(identity is not None and list(_signature(_file(staging / relative))) == json.loads(identity))
            inventory.commit()
        inventory_path.chmod(0o400)
        _sync(inventory_path)
        inventory_sha = _hash_file(inventory_path)
        old = staging / RECORD
        if old.exists():
            os.rename(old, control / 'previous-migration.json')
        root_stat = _directory(staging)
        record = dict(schema=SCHEMA, version=1, status='complete', case_sha256=case_sha,
            package_sha256=package_sha, root_identity=[root_stat.st_dev, root_stat.st_ino],
            inventory=inventory_path.relative_to(staging).as_posix(),
            inventory_identity=list(_signature(_file(inventory_path))), inventory_sha256=inventory_sha,
            files=count, bytes=total)
        _write(staging / RECORD, canonical(record))
        # Flush copied directory entries before the single destination rename.
        for path, directory in _walk(staging):
            if directory:
                _sync(path)
        (staging / INCOMPLETE).unlink()
        _sync(staging)
        try:
            _publish_directory(staging, destination)
        except BaseException:
            if staging.exists():
                _write(staging / INCOMPLETE, canonical(dict(schema=SCHEMA, status='incomplete')))
                _sync(staging)
            raise
        _sync(destination.parent)
    if progress:
        progress('complete', count, total)
    return record


def main(argv=None):
    import argparse
    import signal
    import sys
    parser = argparse.ArgumentParser(description='Copy a stopped Track 1 output and preserve supported checkpoints.')
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--destination', type=Path, required=True)
    args = parser.parse_args(argv)
    os.umask(0o077)
    for sig in (signal.SIGINT, signal.SIGTERM):
        signal.signal(sig, signal.default_int_handler)
    try:
        root = Path(__file__).resolve().parents[2]
        expected = os.environ['GENOSAPIENS_PACKAGE_SHA256']
        from genosapiens.wgs_operator_invocation import _package
        contract, _, _, _, _ = _package(root, expected)
        _require(os.environ.get('GENOSAPIENS_RUNTIME_IMAGE', '').split('@')[-1] == contract.runtime_image_digest)
        import socket
        _require(os.getuid() > 0 and not any(name != 'lo' for _, name in socket.if_nameindex())
                 and os.statvfs('/').f_flag & os.ST_RDONLY, 'use the offline migration launcher')
        migrate_output(args.source, args.destination, package_sha=expected, package_root=root,
            progress=lambda phase, files, size: print(f'Migration {phase}: {files} files, {size} bytes.', flush=True))
        print('Migration complete. Resume with the destination output and unchanged inputs.', flush=True)
        return 0
    except KeyboardInterrupt:
        print('Migration interrupted. The original remains available. The incomplete copy cannot resume.', file=sys.stderr)
    except MigrationError as error:
        print(str(error), file=sys.stderr)
    except (OSError, ValueError, RuntimeError, KeyError, TypeError, sqlite3.Error):
        print('Migration failed. Check the original checkpoint, permissions, and destination storage.', file=sys.stderr)
    return 1


if __name__ == '__main__':
    raise SystemExit(main())
