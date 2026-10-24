"""Offline loading of the externally provisioned public scorer, not preflight."""

from __future__ import annotations

import hashlib
from pathlib import Path
import re
import sys
from types import MappingProxyType, ModuleType

from genosapiens.bounded_json import BoundedJSONError, load_bounded_json, read_bounded_bytes


class EvaluatorPinError(ValueError):
    """The public evaluator or its reviewed pin could not be verified."""


def load_pinned_evaluator(source_path: str | Path, *, pin_path: str | Path) -> ModuleType:
    """Execute verified bytes without downloads, installation, or a local scorer copy.

    The pin is trusted, reviewed project configuration, not caller-supplied data.
    Use only in an offline public-test sandbox. The returned upstream API does not
    enforce our privacy, reference, eligibility, or submission-readiness gates.
    """
    try:
        pin, pin_bytes = load_bounded_json(pin_path, maximum_bytes=65536)
        source = pin["source"]
        repository = "SageBio/rare-disease-real-kid-mva-hackathon-2026"
        revision = source["revision"]
        if (
            pin["schema"] != "genosapiens.challenge-evaluator-pin"
            or source["repository"] != repository
            or source["repository_type"] != "space"
            or source["path"] != "evaluation.py"
            or re.fullmatch(r"[0-9a-f]{40}", revision) is None
            or re.fullmatch(r"[0-9a-f]{64}", source["sha256"]) is None
            or re.fullmatch(r"[0-9a-f]{40}", source["git_blob_sha1"]) is None
            or type(source["size_bytes"]) is not int
            or not 0 < source["size_bytes"] <= 65536
            or source["url"] != (
                f"https://huggingface.co/spaces/{repository}/resolve/{revision}/evaluation.py"
            )
        ):
            raise ValueError
    except (BoundedJSONError, KeyError, TypeError, ValueError):
        raise EvaluatorPinError("public evaluator pin is invalid") from None

    try:
        content = read_bounded_bytes(source_path, maximum_bytes=source["size_bytes"])
    except BoundedJSONError:
        raise EvaluatorPinError("public evaluator source is unavailable or exceeds its pin") from None
    digest = hashlib.sha256(content).hexdigest()
    blob_digest = hashlib.sha1(
        f"blob {len(content)}\0".encode("ascii") + content, usedforsecurity=False
    ).hexdigest()
    if (len(content), digest, blob_digest) != (
        source["size_bytes"], source["sha256"], source["git_blob_sha1"]
    ):
        raise EvaluatorPinError("public evaluator bytes do not match the inspected pin")

    module = ModuleType(f"_genosapiens_public_evaluator_{digest}")
    module.__file__ = "pinned-public-evaluator/evaluation.py"
    # Dataclasses with postponed annotations consult sys.modules at definition time.
    previous = sys.modules.get(module.__name__)
    sys.modules[module.__name__] = module
    try:
        # Compile the hashed buffer, never reopen a potentially replaced source path.
        exec(compile(content, module.__file__, "exec"), module.__dict__)
    finally:
        if previous is None:
            sys.modules.pop(module.__name__, None)
        else:
            sys.modules[module.__name__] = previous
    module.__evaluator_provenance__ = MappingProxyType({
        "pin_sha256": hashlib.sha256(pin_bytes).hexdigest(),
        "revision": revision,
        "source_sha256": digest,
        "source_git_blob_sha1": blob_digest,
        "source_size_bytes": len(content),
        "source_url": source["url"],
    })
    return module
