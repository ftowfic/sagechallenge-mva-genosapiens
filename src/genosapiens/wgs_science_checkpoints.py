"""Independent, durable context artifacts under the existing case output lock.

Only a finished evidence manifest and its context root can be committed. Pending
directories are failure evidence, never completion authority. Scope selection is
deliberately absent from the computation identity.
"""
from dataclasses import asdict, dataclass
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import tempfile

from genosapiens.bounded_json import read_bounded_bytes
from genosapiens.wgs_evidence_store import EvidenceManifest, EvidenceReader, manifest_from_bytes
from genosapiens.wgs_ledger import canonical

SCHEMA = 'genosapiens.scientific-context-checkpoint'
MAX_BYTES = 8 * 1024**2


def _sha(value):
    return type(value) is str and len(value) == 64 and all(c in '0123456789abcdef' for c in value)


def _require(condition):
    if not condition:
        raise ValueError('scientific context checkpoint is invalid or incompatible')


def _path(path):
    path = Path(path).absolute()
    _require(not any(item.is_symlink() for item in (path, *path.parents)))
    return path


def _sync_directory(path):
    descriptor = os.open(path, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _write(path, content):
    _require(type(content) is bytes and 0 < len(content) <= MAX_BYTES)
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600)
    with os.fdopen(descriptor, 'wb') as stream:
        stream.write(content)
        stream.flush()
        os.fsync(stream.fileno())


@dataclass(frozen=True, slots=True)
class ContextIdentity:
    context_id: str
    original_ordinal: int
    membership_sha256: str
    computation_sha256: str

    def __post_init__(self):
        _require(type(self.context_id) is str and 0 < len(self.context_id) <= 256
            and type(self.original_ordinal) is int and 0 < self.original_ordinal < 2**32
            and _sha(self.membership_sha256) and _sha(self.computation_sha256))

    @property
    def sha256(self):
        return hashlib.sha256(canonical(asdict(self))).hexdigest()


@dataclass(frozen=True, slots=True)
class ContextCheckpoint:
    identity: ContextIdentity
    directory: Path
    evidence: EvidenceManifest
    context_name: str
    context_sha256: str
    completion_sha256: str
    case_name: str | None = None
    case_sha256: str | None = None


class ContextCheckpointStore:
    def __init__(self, root):
        self._produced_prefix = None
        self._resource_bundle = None
        self._fresh_output = self._fresh_checkpoint = None
        self.omitted_fresh_replays = 0
        self.root = _path(root)
        self.root.mkdir(mode=0o700, parents=True, exist_ok=True)
        _require(self.root.is_dir())

    def begin(self, identity):
        _require(type(identity) is ContextIdentity)
        _path(self.root)
        return Path(tempfile.mkdtemp(prefix='.pending-' + identity.sha256[:16] + '-', dir=self.root))

    def _destination(self, identity):
        _require(type(identity) is ContextIdentity)
        return _path(self.root / ('context-' + identity.sha256))

    def _read(self, directory, identity):
        directory = _path(directory)
        path = _path(directory / 'complete.json')
        content = read_bounded_bytes(path, maximum_bytes=MAX_BYTES)
        value = json.loads(content)
        extra = {'case_name', 'case_sha256'} if type(value) is dict and value.get('version') == 2 else set()
        _require(type(value) is dict and set(value) == {'schema', 'version', 'identity',
            'evidence_sha256', 'context_name', 'context_sha256', *extra}
            and canonical(value) == content and value['schema'] == SCHEMA
            and type(value['version']) is int and value['version'] in (1, 2)
            and value['identity'] == asdict(identity) and _sha(value['evidence_sha256'])
            and _sha(value['context_sha256']) and type(value['context_name']) is str
            and 0 < len(value['context_name']) <= 512)
        if extra:
            _require(type(value['case_name']) is str and value['case_name'] ==
                'checkpoint-case-' + str(identity.original_ordinal) and _sha(value['case_sha256']))
        raw = read_bounded_bytes(_path(directory / 'evidence.json'), maximum_bytes=MAX_BYTES)
        evidence = manifest_from_bytes(raw, expected_sha256=value['evidence_sha256'])
        return ContextCheckpoint(identity, directory, evidence, value['context_name'],
            value['context_sha256'], hashlib.sha256(content).hexdigest(), value.get('case_name'), value.get('case_sha256'))

    @contextmanager
    def using_produced_prefix(self, prefix):
        from genosapiens.wgs_evidence_store import _ProducedEvidencePrefix
        _require(type(prefix) is _ProducedEvidencePrefix and self._produced_prefix is None)
        prefix.require()
        self._produced_prefix = prefix
        try:
            yield self
        finally:
            self._produced_prefix = None
            self._resource_bundle = None

    @contextmanager
    def evidence_writer(self, pending, *, scratch_directory):
        from genosapiens.wgs_evidence_store import EvidenceWriter
        from genosapiens.wgs_evidence_shared import CheckpointResourcePrefix
        pending = _path(pending)
        _require(pending.parent == self.root and pending.name.startswith('.pending-'))
        if self._produced_prefix is not None and self._resource_bundle is None:
            self._resource_bundle = CheckpointResourcePrefix(self.root, self._produced_prefix)
        self._fresh_output = None
        with EvidenceWriter(pending, scratch_directory=scratch_directory,
                _checkpoint_prefix=self._resource_bundle, _skip_verification=True) as writer:
            yield writer
        if writer.manifest is not None and not writer.failed:
            self._fresh_output = (pending, writer.manifest)

    def accept_fresh(self, identity, checkpoint):
        """Consume live publication ownership, never a persisted trust flag."""
        if self._fresh_checkpoint is None or checkpoint is not self._fresh_checkpoint:
            return None
        self._fresh_checkpoint = None
        durable = self._read(self._destination(identity), identity)
        _require(durable == checkpoint)
        self.omitted_fresh_replays += 1
        return durable

    def open_evidence(self, checkpoint, *, scratch_directory):
        return EvidenceReader(checkpoint.evidence, checkpoint.directory, scratch_directory=scratch_directory,
            _produced_prefix=self._produced_prefix)

    @staticmethod
    def authenticate(checkpoint, *, scratch_directory, _produced_prefix=None):
        """Check source bytes and context identity with the ordinary evidence reader."""
        with EvidenceReader(checkpoint.evidence, checkpoint.directory, scratch_directory=scratch_directory,
                _produced_prefix=_produced_prefix) as reader, reader.read_scope():
            # One complete authentication operation brackets all root and
            # membership reads. Each node still checks fresh bytes. Full-byte
            # boundary checks include shards unused by these header reads.
            reference = reader.reference(checkpoint.context_name, kind='scientific-case-context-v1')
            _require(reference.sha256 == checkpoint.context_sha256)
            header = reader.reconstruct(reference)
            _require(header['schema'] == 'genosapiens.wgs-scientific-case-context'
                and header['context_id'] == checkpoint.identity.context_id
                and header['ordinal'] == checkpoint.identity.original_ordinal)
            if checkpoint.case_name is not None:
                from genosapiens.wgs_case_resources import bound_reference
                from genosapiens.wgs_completed_science import child
                case_ref = reader.reference(checkpoint.case_name, kind='complete-scientific-case-v1')
                _require(case_ref.sha256 == checkpoint.case_sha256)
                case = reader.reconstruct(case_ref)
                _require(case['version'] == 3 and case['context_count'] == 1 and case['resources'] == header['resources'])
                from itertools import islice
                domain = tuple(islice(reader.iter_array(bound_reference(reader, case['context_domain'], 'scientific-context-domain-v1')), 2))
                _require(domain == (checkpoint.identity.original_ordinal,))
                primary = reader.reconstruct(bound_reference(reader, case['primary'], 'scientific-primary-evidence-v1'))
                rows = iter(reader.iter_array(bound_reference(reader, primary['contexts'], 'scientific-case-contexts-v1')))
                row = next(rows, None)
                _require(row is not None and next(rows, None) is None
                    and row['context'] == dict(name=checkpoint.context_name, kind='scientific-case-context-v1',
                        sha256=checkpoint.context_sha256))
                _verify_membership(reader, header, checkpoint.identity.membership_sha256)

    def load(self, identity, *, scratch_directory):
        destination = self._destination(identity)
        if not destination.exists():
            return None
        # An existing final directory with a missing or damaged marker is an
        # error. Never replace it with a new computation or treat it as pending.
        result = self._read(destination, identity)
        self.authenticate(result, scratch_directory=scratch_directory, _produced_prefix=self._produced_prefix)
        return result

    def commit(self, pending, identity, evidence, *, context_name, context_sha256, scratch_directory,
               case_name=None, case_sha256=None):
        pending = _path(pending)
        _require(pending.parent == self.root and pending.name.startswith('.pending-')
            and pending.is_dir() and type(identity) is ContextIdentity and type(evidence) is EvidenceManifest
            and type(context_name) is str and _sha(context_sha256))
        _require((case_name is None and case_sha256 is None) or
            (case_name == 'checkpoint-case-' + str(identity.original_ordinal) and _sha(case_sha256)))
        raw = evidence.as_bytes()
        value = dict(schema=SCHEMA, version=1 if case_name is None else 2, identity=asdict(identity),
            evidence_sha256=hashlib.sha256(raw).hexdigest(), context_name=context_name,
            context_sha256=context_sha256)
        if case_name is not None:
            value.update(case_name=case_name, case_sha256=case_sha256)
        _write(pending / 'evidence.json', raw)
        # Fresh output belongs to this operation's finished writer. The user
        # authorized omitting its immediate independent replay. Imported or
        # restored output still receives ordinary evidence authentication.
        checkpoint = ContextCheckpoint(identity, pending, evidence, context_name, context_sha256, '', case_name, case_sha256)
        fresh = self._fresh_output == (pending, evidence)
        self._fresh_output = None
        if fresh:
            self.omitted_fresh_replays += 1
        else:
            self.authenticate(checkpoint, scratch_directory=scratch_directory, _produced_prefix=self._produced_prefix)
        _write(pending / 'complete.json', canonical(value))
        _sync_directory(pending)
        destination = self._destination(identity)
        _require(not destination.exists())
        os.rename(pending, destination)
        _sync_directory(self.root)
        result = self._read(destination, identity)
        if fresh:
            self._fresh_checkpoint = result
        return result


def _verify_membership(reader, context, expected):
    """Bind the complete pre-filter membership, including excluded members."""
    from genosapiens.evidence_model import _typed
    from genosapiens.wgs_admission import LedgerAdmissionInput
    from genosapiens.wgs_case_resources import bound_reference
    from genosapiens.wgs_completed_science import child
    from genosapiens.wgs_context_evidence import MAX_MEMBER_BYTES
    members = reader.reconstruct(bound_reference(reader, context['members'], 'context-evidence-v1'))
    reader.connection.execute('CREATE TEMP TABLE checkpoint_members(id TEXT PRIMARY KEY,digest TEXT) WITHOUT ROWID')
    count = 0
    try:
        if members.get('work_filter') is not None:
            work = reader.reconstruct(bound_reference(reader, members['work_filter'], 'context-population-work-filter-v1'))
            source = reader.iter_array(bound_reference(reader, work['members'], 'population-work-member-decisions-v1'),
                max_item_bytes=MAX_MEMBER_BYTES)
            pairs = ((row['allele_id'], row['original_input_sha256']) for row in source)
        else:
            def originals():
                for row in reader.iter_array(child(reader, members['members'], 'context-enriched-members-v1'),
                        max_item_bytes=MAX_MEMBER_BYTES + 5 * 1024**2):
                    value = _typed(row['admitted'], LedgerAdmissionInput, decode=True)
                    _require(value.allele is not None)
                    yield value.allele.allele_id, hashlib.sha256(value.as_bytes()).hexdigest()
            pairs = originals()
        for key, digest in pairs:
            count += 1
            _require(count <= 100_000 and type(key) is str and 0 < len(key) <= 1024 and _sha(digest))
            reader.connection.execute('INSERT INTO checkpoint_members VALUES (?,?)', (key, digest))
        digest = hashlib.sha256()
        for row in reader.connection.execute('SELECT id,digest FROM checkpoint_members ORDER BY id'):
            digest.update(canonical(row))
        _require(digest.hexdigest() == expected)
    finally:
        reader.connection.execute('DROP TABLE checkpoint_members')
