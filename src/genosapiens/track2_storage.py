"""Minimal diagnosis-package staging and existing controlled retention ownership.

These are storage controls, not participant authorization or a release gate.
The controlled launcher must authorize access before mounting any input.
"""

from contextlib import contextmanager
from dataclasses import dataclass
import hashlib
import os
from pathlib import Path
import stat

from genosapiens.governance import (
    ArtifactClassification, ArtifactRecord, CaseInventory, GovernanceError,
    load_inventory, validate_relative_classification, write_inventory,
)


class Track2StorageError(ValueError):
    """Static error safe for an operational log."""


def _storage_error(message, operation, error):
    """Expose an operation and numeric errno, never paths or source values."""
    seen = set()
    current = error
    detail = f'operation={operation}'
    while current is not None and id(current) not in seen and len(seen) < 8:
        seen.add(id(current))
        if isinstance(current, OSError) and type(current.errno) is int:
            detail += f', errno={current.errno}'
            break
        current = current.__cause__ or current.__context__
    return Track2StorageError(f'{message} ({detail})')


@dataclass(frozen=True)
class ControlledStore:
    reference_root: Path
    inventory_path: Path
    case_id: str
    output_root: Path | None = None
    input_paths: tuple[Path, ...] = ()
    path_bindings: tuple[tuple[Path, Path], ...] = ()

    def __post_init__(self):
        # The host validates this location. Later workers use it only as an
        # inventory namespace and need not mount the public reference tree.
        root = Path(self.reference_root).absolute()
        object.__setattr__(self, 'reference_root', root)
        if self.output_root is not None:
            object.__setattr__(self, 'output_root', Path(self.output_root).absolute())
        object.__setattr__(self, 'input_paths', tuple(Path(p).absolute() for p in self.input_paths))
        object.__setattr__(self, 'path_bindings', tuple((Path(guest).absolute(), Path(host).absolute())
            for guest, host in self.path_bindings))
        object.__setattr__(self, 'inventory_path', self.path(self.inventory_path))
        self.inventory()

    def path(self, value):
        """Validate the location supplied by this operation, without root rules.

        Docker mounts determine worker access. The inventory records the actual
        location and controlled classification independently of directory names.
        """
        path = Path(value).absolute()
        try:
            if path.resolve() != path:
                raise ValueError
            if any(parent.is_symlink() for parent in (path, *path.parents)):
                raise ValueError
        except (ValueError, OSError):
            raise Track2StorageError('controlled artifact path is inaccessible or contains a symbolic link') from None
        return path

    def host_path(self, value):
        """Translate only selected mounts; never resolve a host path inside Docker."""
        path = self.path(value)
        for guest, host in self.path_bindings:
            if path == guest or (guest == self.output_root and path.is_relative_to(guest)):
                return host / path.relative_to(guest)
        return path

    def relative_path(self, value):
        path = self.host_path(value)
        try:
            relative = path.relative_to(self.reference_root).as_posix()
            validate_relative_classification(relative, ArtifactClassification.CONTROLLED)
            return relative
        except (ValueError, GovernanceError):
            return 'controlled-work/operator-paths/' + hashlib.sha256(str(path).encode()).hexdigest()

    def record_for(self, value, *, kind, digest=None, size=None):
        """Retain the actual host location, including paths outside reference storage."""
        from genosapiens.governance import _file_identity
        path = self.path(value)
        host = self.host_path(path)
        relative = self.relative_path(path)
        if digest is None:
            digest, size = _file_identity(path)
        return ArtifactRecord(hashlib.sha256(f'controlled\0{relative}\0{digest}'.encode()).hexdigest()[:24],
            relative, ArtifactClassification.CONTROLLED, kind, digest, size, 'delete',
            operator_path=str(host) if self.reference_root/relative != host else None)

    def inventory(self):
        try:
            inventory = (load_inventory(self.inventory_path) if self.inventory_path.exists()
                         else CaseInventory(self.case_id, ()))
            if inventory.case_id != self.case_id:
                raise ValueError
            return inventory
        except (ValueError, OSError, GovernanceError) as error:
            raise _storage_error('controlled inventory identity is invalid', 'read_inventory', error) from None

    def record(self, paths, *, kind):
        """Inventory only this operation's small generated outputs, including partials."""
        operation = 'read_inventory'
        try:
            inventory = self.inventory()
            records = {record.relative_path: record for record in inventory.artifacts}
            for value in paths:
                operation = 'read_artifact'
                path = self.path(value)
                if path == self.inventory_path or not stat.S_ISREG(path.lstat().st_mode):
                    raise Track2StorageError('controlled output inventory could not be retained (reason=invalid_artifact_type)')
                record = self.record_for(path, kind=kind)
                previous = records.get(record.relative_path)
                if previous is not None and previous != record:
                    raise Track2StorageError('controlled output inventory could not be retained (reason=record_conflict)')
                records[record.relative_path] = record
            operation = 'write_inventory'
            write_inventory(CaseInventory(self.case_id, tuple(records[key] for key in sorted(records))),
                            self.inventory_path)
        except Track2StorageError:
            raise
        except (ValueError, OSError, GovernanceError) as error:
            raise _storage_error('controlled output inventory could not be retained', operation, error) from None

    def model_artifacts(self):
        """Use the existing NO_LLM guard for every inventoried case derivative."""
        from genosapiens.llm_guard import LLMArtifact, LLMArtifactKind, LLMPolicy
        return {record.artifact_id: LLMArtifact(record.artifact_id,
            record.operator_path or str(self.reference_root / record.relative_path), record.sha256,
            LLMArtifactKind.PATIENT_DERIVED_OUTPUT, LLMPolicy.NO_LLM)
            for record in self.inventory().artifacts
            if record.classification is ArtifactClassification.CONTROLLED and record.status == 'present'}


@contextmanager
def output_storage(output, *, classification, controlled_store=None, kind):
    """Deny untracked controlled writes; keep partial output private and inventoried."""
    target = Path(output)
    if classification not in {'public', 'controlled'}:
        raise Track2StorageError('output classification is invalid')
    if classification == 'public':
        if controlled_store is not None:
            raise Track2StorageError('public output does not use a controlled inventory')
        yield target
        return
    if type(controlled_store) is not ControlledStore:
        raise Track2StorageError('controlled output requires its retention inventory')
    target = controlled_store.path(target)
    if (target == controlled_store.inventory_path
            or controlled_store.inventory_path.is_relative_to(target)
            or target.exists() or target.is_symlink()):
        raise Track2StorageError('controlled output must use a new isolated location')
    previous_umask = os.umask(0o077)
    try:
        yield target
    finally:
        try:
            if target.exists():
                paths = [target] if target.is_file() else sorted(target.rglob('*'))
                files = []
                for path in paths:
                    controlled_store.path(path)
                    if path.is_symlink():
                        raise Track2StorageError('controlled output contains an invalid link')
                    if path.is_file():
                        path.chmod(0o600)
                        files.append(path)
                    elif path.is_dir():
                        path.chmod(0o700)
                    else:
                        raise Track2StorageError('controlled output contains a non-regular artifact')
                if target.is_dir():
                    target.chmod(0o700)
                controlled_store.record(files, kind=kind)
        finally:
            os.umask(previous_umask)


def stage_diagnosis_input(source, output, *, expected_input_sha256, expected_review_sha256,
                          allow_synthetic=False, controlled_store=None, draft=False):
    """Copy exactly the validated package closure, never a source directory mount.

    Extra source files (including raw VCF/HPO, logs and private notes) are omitted.
    The staging receipt is separate from the package mounted by synthesis.
    """
    from genosapiens.track2_diagnosis_package import load_diagnosis_package
    from genosapiens.track2_published_genotype import canonical, sha256, _source_file

    if not expected_input_sha256 or (not draft and not expected_review_sha256):
        raise Track2StorageError('immutable input and review selection are required')
    bundle = load_diagnosis_package(source, allow_synthetic=allow_synthetic, require_review=not draft,
        expected_input_sha256=expected_input_sha256, expected_review_sha256=expected_review_sha256)
    if bundle.review_sha256 != expected_review_sha256:
        raise Track2StorageError('selected diagnosis review changed before staging')
    hashes = {'package.json': bundle.input_sha256,
              'acceptance.json': bundle.automatic_acceptance_sha256,
              **({'review.json': bundle.review_sha256} if bundle.review_sha256 is not None else {}),
              **{record['path']: record['sha256'] for record in bundle.data['artifacts'].values()}}
    contents = {name: _source_file(Path(source), name).read_bytes() for name in hashes}
    if any(sha256(raw) != hashes[name] for name, raw in contents.items()):
        raise Track2StorageError('selected diagnosis input changed during staging')
    with output_storage(output, classification=bundle.data['classification'],
                        controlled_store=controlled_store, kind='track2-minimal-input') as target:
        target.mkdir(parents=True, exist_ok=False, mode=0o700)
        package = target / 'package'
        package.mkdir(mode=0o700)
        for name, raw in contents.items():
            with (package / name).open('xb') as stream:
                stream.write(raw)
        from genosapiens.track2_exclusions import copy_exclusions, DIRECTORY, ROWS, MANIFEST
        descriptor=bundle.components['selection']['coverage']['filter_exclusions']
        copy_exclusions(Path(source)/DIRECTORY,package/DIRECTORY,descriptor)
        hashes[DIRECTORY+'/'+ROWS]=descriptor['sha256']
        hashes[DIRECTORY+'/'+MANIFEST]=sha256(canonical(descriptor))
        from genosapiens.track2_deferred import copy_inventory, DIRECTORY as DEFERRED_DIRECTORY
        descriptor=bundle.components['selection']['coverage']['deferred_inventory']
        copy_inventory(Path(source)/DEFERRED_DIRECTORY,package/DEFERRED_DIRECTORY,bundle.components['selection'])
        hashes[DEFERRED_DIRECTORY+'/'+ROWS]=descriptor['sha256']
        hashes[DEFERRED_DIRECTORY+'/'+MANIFEST]=sha256(canonical(descriptor))
        load_diagnosis_package(package, allow_synthetic=allow_synthetic, require_review=not draft,
            expected_input_sha256=bundle.input_sha256, expected_review_sha256=bundle.review_sha256)
        receipt = {'schema': 'genosapiens.track2-input-staging/v1',
            'input_origin': 'diagnosis_guided_vcf', 'classification': bundle.data['classification'],
            'evidence_origin': bundle.data['evidence_origin'], 'files': hashes}
        (target / 'staging.json').write_bytes(canonical(receipt))
    return receipt
