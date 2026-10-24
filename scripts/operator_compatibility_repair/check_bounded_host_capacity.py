#!/usr/bin/env python3
"""Check public host capacity facts before a bounded controlled stage."""

from __future__ import annotations

import argparse
import hashlib
import os
from pathlib import Path
import stat
import sys

from genosapiens.bounded_wgs_config import load_bounded_wgs_execution_config
from genosapiens.wgs_runtime import WgsRuntimeError, require_fixed_host_capacity


def _read_bounded(path: Path, maximum_bytes: int) -> bytes:
    descriptor = -1
    try:
        descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
        metadata = os.fstat(descriptor)
        if not stat.S_ISREG(metadata.st_mode):
            raise OSError
        content = bytearray()
        while len(content) <= maximum_bytes:
            block = os.read(
                descriptor, min(8192, maximum_bytes + 1 - len(content)),
            )
            if not block:
                break
            content.extend(block)
        if not content or len(content) > maximum_bytes:
            raise OSError
        return bytes(content)
    except OSError:
        raise WgsRuntimeError("bounded host input is unavailable") from None
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _available_memory(path: Path) -> int:
    values = {}
    try:
        for line in _read_bounded(path, 1024 * 1024).decode("ascii").splitlines():
            fields = line.split()
            if len(fields) == 3 and fields[2] == "kB":
                values[fields[0].removesuffix(":")] = int(fields[1]) * 1024
        return values["MemAvailable"]
    except (OSError, UnicodeError, ValueError, KeyError):
        raise WgsRuntimeError("host memory evidence is unavailable") from None


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--config-sha256", required=True)
    parser.add_argument("--meminfo", type=Path, default=Path("/proc/meminfo"))
    # Accept historical driver arguments without probing or requiring disk capacity.
    parser.add_argument("--workspace", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--projected-scratch-bytes", type=int, help=argparse.SUPPRESS)
    arguments = parser.parse_args(argv)
    try:
        content = _read_bounded(arguments.config, 1024 * 1024)
        if hashlib.sha256(content).hexdigest() != arguments.config_sha256:
            raise WgsRuntimeError("bounded runtime configuration changed")
        config = load_bounded_wgs_execution_config(
            content, arguments.config_sha256,
        )
        require_fixed_host_capacity(
            config,
            host_available_memory_bytes=_available_memory(arguments.meminfo),
        )
        return 0
    except (OSError, ValueError, WgsRuntimeError) as error:
        print(str(error), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
