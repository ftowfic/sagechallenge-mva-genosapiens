"""Explicit terminal preservation before bounded preparation or analysis resume."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import re
import stat

from genosapiens.wgs_ledger import (
    Artifact, LedgerError, StageCommitState, StageRunBinding, canonical,
    publish_bytes, stage_lock,
)
from genosapiens.wgs_runtime import (
    MAX_RUNTIME_DOCUMENT_BYTES,
    WgsRuntimeError,
    load_runtime_phase_checkpoint,
    load_stage_runtime_evidence,
    load_stage_runtime_sample,
)


class WgsResumeError(RuntimeError):
    """Static refusal to resume an unauthenticated stage prefix."""


_SCHEMAS = {
    "preparation": "genosapiens.wgs-bounded-preparation-result",
    "analysis": "genosapiens.wgs-bounded-analysis-result",
}
_PARTIAL_ROOTS = {
    "preparation": (
        "source", "prefilter", "working-set", "phenotype", "prepared", "journals",
    ),
    "analysis": ("analysis-work", "controlled-work"),
}
_SHA = re.compile(r"[0-9a-f]{64}")


@dataclass(frozen=True, slots=True)
class ResumeArchive:
    stage: str
    attempt_number: int
    root: Path
    manifest: Artifact


def archive_operator_terminal_for_resume(workspace, *, stage, expected_release_sha256,
                                         attempt):
    """Preserve the small terminal files used by ``track1 run``.

    The invocation lock is held by the caller. Copy and seal both files before
    removing either current name; leave all journals, SQLite state, shards and
    previous invocation outcomes in place. No resource or partial-product scan.
    """
    from genosapiens.evidence_model import _typed
    from genosapiens.wgs_analysis import AnalysisRunReceipt
    from genosapiens.wgs_preparation import PreparationRunReceipt
    from genosapiens.wgs_failures import (
        MAX_FAILURE_EVIDENCE_BYTES, load_failure_evidence, terminal_failure_binding,
    )

    root, attempt = Path(workspace), Path(attempt)
    if (stage not in _SCHEMAS or _SHA.fullmatch(expected_release_sha256) is None
            or not root.is_absolute() or attempt.parent != root.parent / '.operator'
            or not attempt.name.startswith('attempt-')):
        raise WgsResumeError('operator resume archive binding is invalid')
    _private_directory(root)
    _private_directory(attempt)
    with stage_lock(root, 'operator-' + stage):
        path = root / f'controlled-{stage}-result.json'
        content = _private_file(path, 65536)
        kind = PreparationRunReceipt if stage == 'preparation' else AnalysisRunReceipt
        receipt = _typed(_json(content), kind, decode=True)
        if (receipt.as_bytes() != content or receipt.status not in {'failed', 'interrupted'}
                or receipt.base_release_sha256 != expected_release_sha256):
            raise WgsResumeError('only a matching failed or interrupted stage may resume')
        files = [(path, content)]
        failure_path = root / f'failure-{stage}.json'
        if failure_path.exists() or failure_path.is_symlink():
            failure_content = _private_file(failure_path, MAX_FAILURE_EVIDENCE_BYTES)
            code, heads = terminal_failure_binding(root, stage=stage, receipt_content=content)
            diagnostic = load_failure_evidence(failure_content, receipt_content=content,
                expected_failure_code=code, expected_journal_heads=heads)
            if diagnostic['stage'] != stage:
                raise WgsResumeError('operator resume diagnostic stage differs')
            files.append((failure_path, failure_content))
        history = attempt / ('previous-' + stage)
        history.mkdir(mode=0o700, exist_ok=True)
        _private_directory(history)
        def preserve(name, data):
            destination = history / name
            if destination.exists() or destination.is_symlink():
                if _private_file(destination, len(data)) != data:
                    raise WgsResumeError('operator resume archive already differs')
            else:
                publish_bytes(history, name, data)
        for source, data in files:
            preserve(source.name, data)
        preserve('archive.json', canonical(dict(
            schema='genosapiens.operator-resume-archive', version=1,
            classification='controlled', llm_policy='no-llm', stage=stage,
            base_release_sha256=expected_release_sha256, previous_status=receipt.status,
            files=[dict(name=source.name, sha256=hashlib.sha256(data).hexdigest(),
                        size_bytes=len(data)) for source, data in files])))
        # Remove the receipt last. A stop between removals leaves an authentic
        # receipt with an optional missing diagnostic, and a complete archive.
        for source, data in reversed(files):
            if _private_file(source, len(data)) != data:
                raise WgsResumeError('operator resume terminal changed while archiving')
            source.unlink()
        descriptor = os.open(root, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(descriptor)
        finally:
            os.close(descriptor)
        return history


def archive_stage_orphans_for_resume(
    stage_root: str | Path,
    state_path: str | Path,
    binding: StageRunBinding,
    *,
    lock,
    partition_id: str,
    names: tuple[str, ...],
) -> tuple[str, ...]:
    """Retire only declared uncommitted publications under the stage lock."""
    try:
        with StageCommitState(stage_root, state_path, binding, lock=lock) as state:
            start = state.open_attempt(partition_id)
            if start is None:
                raise WgsResumeError("no incomplete WGS stage attempt to preserve")
            moved = state.archive_orphans(names=names, attempt_id=start.attempt_id)
            state.fail(start, failure_code="process_lost", interrupted=True)
            return moved
    except LedgerError:
        raise WgsResumeError("WGS stage checkpoint is inconsistent") from None


def _private_directory(path: Path) -> None:
    metadata = path.lstat()
    if (
        path.is_symlink()
        or not stat.S_ISDIR(metadata.st_mode)
        or metadata.st_uid != os.getuid()
        or stat.S_IMODE(metadata.st_mode) != 0o700
    ):
        raise WgsResumeError("resume directory is unsafe")


def _private_file(path: Path, maximum_bytes: int | None = None) -> bytes:
    try:
        metadata = path.lstat()
        if (
            path.is_symlink()
            or not stat.S_ISREG(metadata.st_mode)
            or metadata.st_uid != os.getuid()
            or stat.S_IMODE(metadata.st_mode) not in {0o400, 0o500, 0o600}
            or metadata.st_nlink != 1
            or maximum_bytes is not None and metadata.st_size > maximum_bytes
        ):
            raise ValueError
        with path.open("rb") as stream:
            content = stream.read(
                (maximum_bytes + 1) if maximum_bytes is not None
                else metadata.st_size + 1
            )
        if maximum_bytes is not None and len(content) > maximum_bytes:
            raise ValueError
        return content
    except (OSError, ValueError):
        raise WgsResumeError("resume artifact is unsafe") from None


def _json(content: bytes) -> dict:
    def unique(pairs):
        value = dict(pairs)
        if len(value) != len(pairs):
            raise ValueError
        return value
    try:
        value = json.loads(content, object_pairs_hook=unique)
    except (UnicodeError, json.JSONDecodeError, ValueError):
        raise WgsResumeError("resume artifact is invalid") from None
    if type(value) is not dict:
        raise WgsResumeError("resume artifact is invalid")
    return value


def _identity(path: Path, content: bytes) -> dict:
    return {
        "name": path.name,
        "sha256": hashlib.sha256(content).hexdigest(),
        "size_bytes": len(content),
        "mode": stat.S_IMODE(path.lstat().st_mode),
    }


def _partial_metadata(path: Path) -> dict:
    """Inventory an owned partial file without rereading scientific bytes."""
    metadata = path.lstat()
    if (
        path.is_symlink()
        or not stat.S_ISREG(metadata.st_mode)
        or metadata.st_uid != os.getuid()
        or stat.S_IMODE(metadata.st_mode) not in {0o400, 0o500, 0o600}
        or metadata.st_nlink != 1
    ):
        raise WgsResumeError("partial resume product is unsafe")
    return {
        "size_bytes": metadata.st_size,
        "mode": stat.S_IMODE(metadata.st_mode),
        "device": metadata.st_dev,
        "inode": metadata.st_ino,
        "mtime_ns": metadata.st_mtime_ns,
        "ctime_ns": metadata.st_ctime_ns,
    }


def _partial_products(workspace: Path, stage: str) -> tuple[dict, ...]:
    records = []
    for name in _PARTIAL_ROOTS[stage]:
        root = workspace / name
        if not root.exists():
            continue
        _private_directory(root)
        for path in sorted(root.rglob("*")):
            relative = path.relative_to(workspace).as_posix()
            metadata = path.lstat()
            if path.is_symlink() or metadata.st_uid != os.getuid():
                raise WgsResumeError("partial resume product is unsafe")
            if path.is_dir():
                _private_directory(path)
                records.append({"path": relative, "kind": "directory"})
            elif stat.S_ISREG(metadata.st_mode) and metadata.st_nlink == 1:
                records.append({
                    "path": relative,
                    "kind": "file",
                    **_partial_metadata(path),
                })
            else:
                raise WgsResumeError("partial resume product is unsafe")
    return tuple(records)


def archive_terminal_for_resume(
    workspace: str | Path,
    *,
    stage: str,
    expected_release_sha256: str,
) -> ResumeArchive:
    """Authenticate and preserve one failed terminal without moving stage work."""

    root = Path(workspace)
    if (
        stage not in _SCHEMAS
        or not root.is_absolute()
        or _SHA.fullmatch(expected_release_sha256) is None
    ):
        raise WgsResumeError("resume request is invalid")
    _private_directory(root)
    with stage_lock(root, "operator-" + stage):
        receipt_path = root / f"controlled-{stage}-result.json"
        receipt_content = _private_file(receipt_path, 65536)
        receipt = _json(receipt_content)
        if (
            receipt.get("schema") != _SCHEMAS[stage]
            or receipt.get("version") != (3 if stage == "preparation" else 1)
            or receipt.get("classification") != "controlled"
            or receipt.get("llm_policy") != "no-llm"
            or receipt.get("base_release_sha256") != expected_release_sha256
            or receipt.get("status") not in {"failed", "interrupted"}
        ):
            raise WgsResumeError("only a failed or interrupted bounded stage may resume")
        runtime_path = root / f"runtime-{stage}.json"
        evidence_path = root / f"operator-{stage}-evidence.json"
        stdout_path = root / "operator-logs" / f"{stage}.stdout"
        stderr_path = root / "operator-logs" / f"{stage}.stderr"
        required = (
            (receipt_path, receipt_content),
            (runtime_path, _private_file(runtime_path, 65536)),
            (evidence_path, _private_file(evidence_path, 65536)),
            (stdout_path, _private_file(stdout_path, 64 * 1024**2)),
            (stderr_path, _private_file(stderr_path, 64 * 1024**2)),
        )
        try:
            runtime = load_stage_runtime_evidence(
                required[1][1], expected_stage=stage,
            )
        except WgsRuntimeError:
            raise WgsResumeError("resume runtime evidence is incompatible") from None
        evidence = _json(required[2][1])
        if (
            runtime.exit_code == 0
            or evidence.get("schema")
            != "genosapiens.wgs-bounded-operator-evidence"
            or evidence.get("stage") != stage
            or evidence.get("base_release_sha256") != expected_release_sha256
            or evidence.get("receipt_status") != receipt["status"]
        ):
            raise WgsResumeError("resume terminal evidence is incompatible")
        expected = {
            path.relative_to(root).as_posix(): hashlib.sha256(content).hexdigest()
            for path, content in (required[0], required[1], required[3], required[4])
        }
        artifacts = evidence.get("artifacts")
        if type(artifacts) is not list:
            raise WgsResumeError("resume terminal inventory is invalid")
        actual = {
            item.get("relative_path"): item.get("sha256")
            for item in artifacts if type(item) is dict
        }
        if (
            len(actual) != len(artifacts)
            or any(actual.get(path) != digest for path, digest in expected.items())
        ):
            raise WgsResumeError("resume terminal inventory changed")
        sample_path = root / f".runtime-sample-{stage}.json"
        checkpoint_path = root / f".runtime-phases-{stage}.json"
        optional = []
        from genosapiens.wgs_failures import (
            MAX_FAILURE_EVIDENCE_BYTES,
            load_failure_evidence,
            terminal_failure_binding,
        )

        failure_path = root / f"failure-{stage}.json"
        failure_content = _private_file(
            failure_path, MAX_FAILURE_EVIDENCE_BYTES,
        )
        try:
            expected_code, expected_heads = terminal_failure_binding(
                root, stage=stage, receipt_content=receipt_content,
            )
            load_failure_evidence(
                failure_content,
                receipt_content=receipt_content,
                expected_failure_code=expected_code,
                expected_journal_heads=expected_heads,
            )
        except ValueError:
            raise WgsResumeError("resume failure evidence is invalid") from None
        expected[failure_path.name] = hashlib.sha256(failure_content).hexdigest()
        if actual.get(failure_path.name) != expected[failure_path.name]:
            raise WgsResumeError("resume failure evidence changed")
        optional.append((failure_path, failure_content))
        if sample_path.exists() or sample_path.is_symlink():
            sample_content = _private_file(
                sample_path, MAX_RUNTIME_DOCUMENT_BYTES,
            )
            if actual.get(sample_path.name) != hashlib.sha256(sample_content).hexdigest():
                raise WgsResumeError("resume runtime sample changed")
            try:
                sample = load_stage_runtime_sample(
                    sample_content,
                    expected_stage=stage,
                    expected_attempt_id=runtime.attempt_id,
                    expected_execution_config_sha256=runtime.execution_config_sha256,
                    expected_exit_code=runtime.exit_code,
                )
            except WgsRuntimeError:
                raise WgsResumeError("resume runtime sample is invalid") from None
            if sample.metrics != runtime.metrics:
                raise WgsResumeError("resume runtime sample is incompatible")
            optional.append((sample_path, sample_content))
            expected[sample_path.name] = hashlib.sha256(sample_content).hexdigest()
        elif runtime.sample_sha256 is not None:
            raise WgsResumeError("resume runtime sample is missing")
        checkpoint_sha256 = runtime.metrics.phase_checkpoint_sha256
        if checkpoint_sha256 is not None:
            checkpoint_content = _private_file(
                checkpoint_path, MAX_RUNTIME_DOCUMENT_BYTES,
            )
            if (
                hashlib.sha256(checkpoint_content).hexdigest()
                != checkpoint_sha256
                or actual.get(checkpoint_path.name) != checkpoint_sha256
            ):
                raise WgsResumeError("resume runtime phase checkpoint changed")
            try:
                load_runtime_phase_checkpoint(
                    checkpoint_content,
                    expected_stage=stage,
                    expected_attempt_id=runtime.attempt_id,
                    expected_execution_config_sha256=runtime.execution_config_sha256,
                )
            except WgsRuntimeError:
                raise WgsResumeError("resume runtime phase checkpoint is invalid") from None
            optional.append((checkpoint_path, checkpoint_content))
            expected[checkpoint_path.name] = checkpoint_sha256
        elif checkpoint_path.exists() or checkpoint_path.is_symlink():
            raise WgsResumeError("resume runtime phase checkpoint is unexpected")
        if set(actual) != set(expected):
            raise WgsResumeError("resume terminal inventory is not exact")
        history_root = root / "resume-history"
        if not history_root.exists():
            history_root.mkdir(mode=0o700)
        _private_directory(history_root)
        names = sorted(
            path.name for path in history_root.iterdir()
            if path.name.startswith(stage + "-terminal-attempt-")
        )
        expected_names = [
            f"{stage}-terminal-attempt-{index:02d}"
            for index in range(1, len(names) + 1)
        ]
        if names != expected_names:
            raise WgsResumeError("resume history sequence is invalid")
        attempt_number = len(names) + 1
        history = history_root / (
            f"{stage}-terminal-attempt-{attempt_number:02d}"
        )
        history.mkdir(mode=0o700)
        partial = _partial_products(root, stage)
        terminal = []
        for path, content in (*required, *optional):
            destination = history / path.name
            os.rename(path, destination)
            terminal.append(_identity(destination, content))
        manifest_content = canonical({
            "schema": "genosapiens.wgs-bounded-resume-archive",
            "version": 1,
            "classification": "controlled",
            "llm_policy": "no-llm",
            "stage": stage,
            "attempt_number": attempt_number,
            "base_release_sha256": expected_release_sha256,
            "previous_status": receipt["status"],
            "terminal_artifacts": terminal,
            "preserved_partial_products": partial,
        })
        manifest = publish_bytes(history, "resume-archive.json", manifest_content)
        return ResumeArchive(stage, attempt_number, history, manifest)
