"""Shared confined immutable I/O for current Track 1 publications.

Scientific intake and prepared evidence are owned by the bounded ledger contracts.
"""
from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from enum import StrEnum
import hashlib
import os
from pathlib import Path, PurePosixPath
import re
import stat
import uuid

from genosapiens.governance import ArtifactClassification, GovernanceError, validate_relative_classification

SHA256 = re.compile(r"^[0-9a-f]{64}$")
IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:-]*$")


class Track1ArtifactError(RuntimeError):
    """Raised with static messages when publication evidence is inconsistent."""


class AnalysisMode(StrEnum):
    PRIMARY = "primary"
    DIAGNOSIS_INFORMED_SENSITIVITY = "diagnosis-informed-sensitivity"
    DOMINANT_TERM_SENSITIVITY = "dominant-term-sensitivity"
    UNCONFIRMED_REFERENCE_INVESTIGATION = "unconfirmed-reference-investigation"
    VALIDATION_ONLY = "validation-only"


def _sha256(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def _require_sha256(value: str) -> None:
    if type(value) is not str or SHA256.fullmatch(value) is None:
        raise Track1ArtifactError("invalid publication SHA-256 binding")


@dataclass(frozen=True, slots=True)
class ArtifactBinding:
    role: str
    relative_path: str
    sha256: str
    size_bytes: int

    def __post_init__(self) -> None:
        path = PurePosixPath(self.relative_path)
        if (
            IDENTIFIER.fullmatch(self.role) is None
            or path.is_absolute()
            or path == PurePosixPath(".")
            or ".." in path.parts
            or not isinstance(self.size_bytes, int)
            or isinstance(self.size_bytes, bool)
            or self.size_bytes < 0
        ):
            raise Track1ArtifactError("invalid publication artifact binding")
        _require_sha256(self.sha256)


def _artifact(role: str, relative_path: str, content: bytes) -> ArtifactBinding:
    return ArtifactBinding(role, relative_path, _sha256(content), len(content))


def _controlled_path(
    reference_root: str | Path,
    path: str | Path,
) -> tuple[Path, Path]:
    try:
        root = Path(reference_root).resolve(strict=True)
        target = Path(path)
        if not target.is_absolute() or not root.is_dir():
            raise ValueError
        relative = target.relative_to(root)
        validate_relative_classification(
            relative.as_posix(), ArtifactClassification.CONTROLLED
        )
        current = root
        for part in relative.parts:
            current = current / part
            if current.is_symlink():
                raise ValueError
            if current.exists() and current != target and not current.is_dir():
                raise ValueError
        if not target.resolve(strict=False).is_relative_to(root):
            raise ValueError
        return root, target
    except (GovernanceError, OSError, ValueError) as error:
        raise Track1ArtifactError("invalid controlled publication path") from error


def _read_regular_file(path: Path, root: Path, required_mode: int) -> bytes:
    relative = path.relative_to(root)
    directory_flags = (
        os.O_RDONLY
        | getattr(os, "O_CLOEXEC", 0)
        | getattr(os, "O_DIRECTORY", 0)
        | getattr(os, "O_NOFOLLOW", 0)
    )
    file_flags = (
        os.O_RDONLY
        | getattr(os, "O_CLOEXEC", 0)
        | getattr(os, "O_NONBLOCK", 0)
        | getattr(os, "O_NOFOLLOW", 0)
    )
    directory_fd = os.open(root, directory_flags)
    try:
        for part in relative.parts[:-1]:
            child_fd = os.open(part, directory_flags, dir_fd=directory_fd)
            os.close(directory_fd)
            directory_fd = child_fd
        file_fd = os.open(relative.name, file_flags, dir_fd=directory_fd)
        try:
            metadata = os.fstat(file_fd)
            if (
                not stat.S_ISREG(metadata.st_mode)
                or stat.S_IMODE(metadata.st_mode) != required_mode
            ):
                raise ValueError
            with os.fdopen(file_fd, "rb", closefd=False) as source:
                content = source.read(metadata.st_size + 1)
            if len(content) != metadata.st_size:
                raise ValueError
            return content
        finally:
            os.close(file_fd)
    finally:
        os.close(directory_fd)


def _existing_content(path: Path, content: bytes, root: Path) -> bool:
    try:
        return _read_regular_file(path, root, 0o400) == content
    except (OSError, ValueError):
        return False


def _publish_immutable(path: Path, content: bytes, root: Path) -> bool:
    return _publish_immutable_validated(path, content, root, _controlled_path)


def _publish_immutable_validated(
    path: Path, content: bytes, root: Path,
    validate_path: Callable[[Path, Path], object],
) -> bool:
    """Share immutable I/O while each caller enforces its own path boundary."""
    validate_path(root, path)
    if path.exists() or path.is_symlink():
        if not _existing_content(path, content, root):
            raise Track1ArtifactError("published publication artifact is immutable")
        return True
    try:
        path.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
        validate_path(root, path)
        partial = path.parent / f".{path.name}.{uuid.uuid4().hex}.partial"
        try:
            with partial.open("xb") as output:
                os.chmod(partial, 0o400)
                output.write(content)
                output.flush()
                os.fsync(output.fileno())
            try:
                os.link(partial, path)
            except FileExistsError:
                if not _existing_content(path, content, root):
                    raise Track1ArtifactError("published publication artifact is immutable")
            directory_fd = os.open(path.parent, os.O_RDONLY)
            try:
                os.fsync(directory_fd)
            finally:
                os.close(directory_fd)
        finally:
            partial.unlink(missing_ok=True)
    except Track1ArtifactError:
        raise
    except OSError as error:
        raise Track1ArtifactError("cannot publish immutable publication artifact") from error
    if not _existing_content(path, content, root):
        raise Track1ArtifactError("published publication artifact failed verification")
    return False
