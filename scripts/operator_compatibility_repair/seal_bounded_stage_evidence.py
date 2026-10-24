#!/usr/bin/env python3
"""Seal bounded receipt, runtime evidence, and independent operator logs."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import stat
import sys

from genosapiens.wgs_runtime import WgsRuntimeError


SCHEMAS = {
    "preparation": "genosapiens.wgs-bounded-preparation-result",
    "analysis": "genosapiens.wgs-bounded-analysis-result",
}


def _canonical(value: object) -> bytes:
    return (
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
        + "\n"
    ).encode("ascii")


def _read_private(path: Path, *, maximum_bytes: int) -> bytes:
    metadata = path.lstat()
    if (
        path.is_symlink()
        or not stat.S_ISREG(metadata.st_mode)
        or metadata.st_uid != os.getuid()
        or stat.S_IMODE(metadata.st_mode) != 0o600
        or metadata.st_nlink != 1
        or metadata.st_size > maximum_bytes
    ):
        raise ValueError
    return path.read_bytes()


def _json(content: bytes) -> dict:
    def unique(pairs):
        value = dict(pairs)
        if len(value) != len(pairs):
            raise ValueError
        return value

    value = json.loads(content, object_pairs_hook=unique)
    if type(value) is not dict:
        raise ValueError
    return value


def _write_private(path: Path, content: bytes) -> None:
    descriptor = os.open(
        path,
        os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
        0o600,
    )
    try:
        os.fchmod(descriptor, 0o600)
        view = memoryview(content)
        while view:
            written = os.write(descriptor, view)
            if written <= 0:
                raise OSError
            view = view[written:]
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def seal_stage_evidence(
    workspace: Path,
    *,
    stage: str,
    expected_release_sha256: str,
) -> Path:
    if (
        stage not in SCHEMAS
        or not workspace.is_absolute()
        or workspace.is_symlink()
        or len(expected_release_sha256) != 64
        or any(value not in "0123456789abcdef" for value in expected_release_sha256)
    ):
        raise ValueError
    metadata = workspace.lstat()
    if (
        not stat.S_ISDIR(metadata.st_mode)
        or stat.S_IMODE(metadata.st_mode) != 0o700
        or metadata.st_uid != os.getuid()
    ):
        raise ValueError
    receipt_path = workspace / f"controlled-{stage}-result.json"
    receipt_content = _read_private(receipt_path, maximum_bytes=65536)
    receipt = _json(receipt_content)
    if (
        receipt.get("schema") != SCHEMAS[stage]
        or receipt.get("version") != (3 if stage == "preparation" else 1)
        or receipt.get("classification") != "controlled"
        or receipt.get("llm_policy") != "no-llm"
        or receipt.get("base_release_sha256") != expected_release_sha256
        or receipt.get("status") not in {"passed", "failed", "interrupted"}
    ):
        raise ValueError
    runtime_path = workspace / f"runtime-{stage}.json"
    runtime_content = _read_private(runtime_path, maximum_bytes=65536)
    from genosapiens.wgs_runtime import (
        MAX_RUNTIME_DOCUMENT_BYTES,
        load_runtime_phase_checkpoint,
        load_stage_runtime_evidence,
        load_stage_runtime_sample,
    )
    runtime = load_stage_runtime_evidence(
        runtime_content,
        expected_stage=stage,
    )
    if (
        (receipt["status"] == "passed") != (runtime.exit_code == 0)
    ):
        raise ValueError
    paths = [
        receipt_path,
        runtime_path,
        workspace / "operator-logs" / f"{stage}.stdout",
        workspace / "operator-logs" / f"{stage}.stderr",
    ]
    if receipt["status"] != "passed":
        from genosapiens.wgs_failures import (
            MAX_FAILURE_EVIDENCE_BYTES,
            load_failure_evidence,
            terminal_failure_binding,
        )

        failure_path = workspace / f"failure-{stage}.json"
        expected_code, expected_heads = terminal_failure_binding(
            workspace, stage=stage, receipt_content=receipt_content,
        )
        load_failure_evidence(
            _read_private(
                failure_path, maximum_bytes=MAX_FAILURE_EVIDENCE_BYTES,
            ),
            receipt_content=receipt_content,
            expected_failure_code=expected_code,
            expected_journal_heads=expected_heads,
        )
        paths.append(failure_path)
    sample_path = workspace / f".runtime-sample-{stage}.json"
    sample_sha256 = runtime.sample_sha256
    if sample_path.exists():
        sample_content = _read_private(
            sample_path, maximum_bytes=MAX_RUNTIME_DOCUMENT_BYTES,
        )
        if hashlib.sha256(sample_content).hexdigest() != sample_sha256:
            raise ValueError
        sample = load_stage_runtime_sample(
            sample_content,
            expected_stage=stage,
            expected_attempt_id=runtime.attempt_id,
            expected_execution_config_sha256=runtime.execution_config_sha256,
            expected_exit_code=runtime.exit_code,
        )
        if sample.metrics != runtime.metrics:
            raise ValueError
        paths.append(sample_path)
    elif sample_sha256 is not None:
        raise ValueError
    checkpoint_sha256 = runtime.metrics.phase_checkpoint_sha256
    if checkpoint_sha256 is not None:
        checkpoint_path = workspace / f".runtime-phases-{stage}.json"
        checkpoint_content = _read_private(
            checkpoint_path, maximum_bytes=MAX_RUNTIME_DOCUMENT_BYTES,
        )
        if hashlib.sha256(checkpoint_content).hexdigest() != checkpoint_sha256:
            raise ValueError
        load_runtime_phase_checkpoint(
            checkpoint_content,
            expected_stage=stage,
            expected_attempt_id=runtime.attempt_id,
            expected_execution_config_sha256=runtime.execution_config_sha256,
        )
        paths.append(checkpoint_path)
    records = []
    for path in paths:
        content = _read_private(path, maximum_bytes=64 * 1024**2)
        records.append({
            "relative_path": path.relative_to(workspace).as_posix(),
            "sha256": hashlib.sha256(content).hexdigest(),
            "size_bytes": len(content),
            "mode": 0o600,
            "classification": "controlled",
        })
    external_dependencies = []
    source_input_sha256 = receipt.get("completed_source_input_sha256") if stage == "preparation" else None
    if source_input_sha256 is not None:
        from genosapiens.wgs_preparation import load_preparation_source_binding, completed_source_readonly_paths
        binding = load_preparation_source_binding(workspace, expected_sha256=source_input_sha256)
        if receipt.get("source_origin") == "completed_input":
            external_dependencies = [{"path": path, "access": "read_only", "retention": "required"}
                for path in completed_source_readonly_paths(binding)]
    output = workspace / f"operator-{stage}-evidence.json"
    _write_private(output, _canonical({
        "schema": "genosapiens.wgs-bounded-operator-evidence",
        "version": 1,
        "classification": "controlled",
        "stage": stage,
        "base_release_sha256": expected_release_sha256,
        "receipt_status": receipt["status"],
        "artifacts": records,
        "completed_source_input_sha256": source_input_sha256,
        "external_dependencies": external_dependencies,
    }))
    return output


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=tuple(SCHEMAS))
    parser.add_argument("workspace", type=Path)
    parser.add_argument("release")
    arguments = parser.parse_args(argv)
    try:
        seal_stage_evidence(
            arguments.workspace,
            stage=arguments.stage,
            expected_release_sha256=arguments.release,
        )
        return 0
    except (OSError, UnicodeError, ValueError, TypeError, WgsRuntimeError):
        print("Bounded stage evidence sealing failed; existing evidence is preserved.", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
