#!/usr/bin/env python3
"""Observe one bounded stage inside its whole-container cgroup."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

from genosapiens.wgs_runtime import WgsRuntimeError, observe_bounded_stage


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("preparation", "analysis"), required=True)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--attempt-id", required=True)
    parser.add_argument("--execution-config-sha256", required=True)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    arguments = parser.parse_args(argv)
    command = tuple(arguments.command)
    if command and command[0] == "--":
        command = command[1:]
    try:
        exit_code, _ = observe_bounded_stage(
            stage=arguments.stage,
            workspace=arguments.workspace,
            attempt_id=arguments.attempt_id,
            execution_config_sha256=arguments.execution_config_sha256,
            command=command,
            require_bounded_limits=True,
        )
        return exit_code
    except WgsRuntimeError as error:
        print(str(error), file=sys.stderr)
        return 125


if __name__ == "__main__":
    raise SystemExit(main())
