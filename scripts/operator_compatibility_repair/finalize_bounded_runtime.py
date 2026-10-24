#!/usr/bin/env python3
"""Finalize bounded runtime evidence before removing a stopped container."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import stat
import sys

from genosapiens.wgs_runtime import WgsRuntimeError, finalize_bounded_runtime


_CONTAINER_ID = re.compile(r"[0-9a-f]{64}")


def _require_bounded_docker_limits(record: dict) -> None:
    """Authenticate the stopped container's effective v5 launch controls."""
    host = record.get("HostConfig")
    config = record.get("Config")
    if type(host) is not dict or type(config) is not dict:
        raise ValueError
    tmpfs = host.get("Tmpfs")
    options = None if type(tmpfs) is not dict else tmpfs.get("/tmp")
    flags = set(options.split(",")) if type(options) is str else set()
    sizes = {value for value in flags if value.startswith("size=")}
    security = host.get("SecurityOpt")
    if (
        host.get("Memory") != 24 * 1024**3
        or host.get("MemorySwap") != 24 * 1024**3
        or host.get("NanoCpus") != 12_000_000_000
        or host.get("PidsLimit") != 512
        or host.get("NetworkMode") != "none"
        or host.get("ReadonlyRootfs") is not True
        or host.get("Privileged") is not False
        or host.get("CapDrop") != ["ALL"]
        or host.get("CapAdd") not in (None, [])
        or type(security) is not list
        or "no-new-privileges" not in security
        or config.get("User") != "65532:65532"
        or not {"rw", "nosuid", "nodev"}.issubset(flags)
        or len(sizes) != 1
        or not sizes <= {"size=1g", "size=1024m", "size=1073741824"}
    ):
        raise ValueError


def load_terminal_docker_state(path: Path, *, expected_container_id: str):
    """Read one saved post-exit Docker inspect before container removal."""

    try:
        metadata = path.lstat()
        if (
            path.is_symlink()
            or not stat.S_ISREG(metadata.st_mode)
            or metadata.st_size <= 0
            or metadata.st_size > 256 * 1024
            or _CONTAINER_ID.fullmatch(expected_container_id) is None
        ):
            raise ValueError
        content = path.read_bytes()
        value = json.loads(content)
        if type(value) is not list or len(value) != 1 or type(value[0]) is not dict:
            raise ValueError
        record = value[0]
        state = record.get("State")
        _require_bounded_docker_limits(record)
        if (
            record.get("Id") != expected_container_id
            or type(state) is not dict
            or state.get("Status") != "exited"
            or type(state.get("ExitCode")) is not int
            or isinstance(state.get("ExitCode"), bool)
            or state["ExitCode"] < 0
            or state["ExitCode"] > 255
            or type(state.get("OOMKilled")) is not bool
            or not isinstance(state.get("FinishedAt"), str)
            or not state["FinishedAt"]
        ):
            raise ValueError
        return state["ExitCode"], state["OOMKilled"]
    except (OSError, UnicodeError, json.JSONDecodeError, ValueError, TypeError):
        raise WgsRuntimeError("terminal Docker state is unavailable or invalid") from None


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("preparation", "analysis"), required=True)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--attempt-id", required=True)
    parser.add_argument("--container-inspect", type=Path, required=True)
    parser.add_argument("--container-id", required=True)
    parser.add_argument("--execution-config-sha256", required=True)
    arguments = parser.parse_args(argv)
    try:
        exit_code, oom_killed = load_terminal_docker_state(
            arguments.container_inspect,
            expected_container_id=arguments.container_id,
        )
        finalize_bounded_runtime(
            stage=arguments.stage,
            workspace=arguments.workspace,
            attempt_id=arguments.attempt_id,
            exit_code=exit_code,
            oom_killed=oom_killed,
            execution_config_sha256=arguments.execution_config_sha256,
        )
        return 0
    except WgsRuntimeError as error:
        print(str(error), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
