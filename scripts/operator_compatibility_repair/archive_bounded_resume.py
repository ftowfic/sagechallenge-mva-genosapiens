#!/usr/bin/env python3
"""Authenticate and preserve one bounded terminal before explicit resume."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

from genosapiens.wgs_resume import WgsResumeError, archive_terminal_for_resume


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("preparation", "analysis"))
    parser.add_argument("workspace", type=Path)
    parser.add_argument("release")
    arguments = parser.parse_args(argv)
    try:
        archive_terminal_for_resume(
            arguments.workspace,
            stage=arguments.stage,
            expected_release_sha256=arguments.release,
        )
        return 0
    except WgsResumeError as error:
        print(str(error), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
