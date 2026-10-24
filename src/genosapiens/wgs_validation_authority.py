"""Truth-free public worker authority before a controlled-use release exists.

The coordinator still owns scientific invocation budgets. This record binds one
already issued runner view to its finite policy family; it cannot authorize a
competitive publication, controlled input, rescue, or an accepted release.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import os
from pathlib import Path
import re
import stat
from typing import Literal

from genosapiens.evidence_model import _typed
from genosapiens.track1_policy_bindings import PolicyBindings
from genosapiens.track1_validation_plan import RankerCaseView, ranker_case_view_as_json
from genosapiens.wgs_preparation import CompletedSourceArtifact, COMPLETED_SOURCE_SCHEMA
from genosapiens.wgs_ledger import canonical, _json, MAX_MANIFEST_BYTES


SCHEMA = "genosapiens.wgs-public-validation-authority"
_SHA = re.compile(r"[0-9a-f]{64}")
_FORBIDDEN = frozenset({"data", "controlled", "controlled-work"})


class WgsValidationAuthorityError(ValueError):
    """Static rejection of a public validation binding."""


def _sha(value):
    if type(value) is not str or _SHA.fullmatch(value) is None:
        raise WgsValidationAuthorityError("invalid public validation hash")


def _public_path(value, *, output=False):
    path = Path(value)
    parts = path.parts
    fixture = any(parts[i:i + 3] == ("tests", "fixtures", "track1")
                  for i in range(len(parts) - 2))
    public_build = any(parts[i:i + 2] == ("scratch", "public-build")
                       for i in range(len(parts) - 1))
    if (not path.is_absolute() or path.as_posix() != value or ".." in parts
            or _FORBIDDEN.intersection(parts) or (not public_build if output else not fixture)
            or any(parent.is_symlink() for parent in (path, *path.parents))):
        raise WgsValidationAuthorityError("public validation path is outside its boundary")
    return path


@dataclass(frozen=True, slots=True)
class PublicValidationInput:
    kind: Literal["vcf", "hpo_text"]
    path: str
    sha256: str
    size_bytes: int

    def __post_init__(self):
        _sha(self.sha256)
        _public_path(self.path)
        if (self.kind not in ("vcf", "hpo_text") or type(self.size_bytes) is not int
                or self.size_bytes <= 0):
            raise WgsValidationAuthorityError("invalid public validation input")

    def verify(self):
        path = _public_path(self.path)
        descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
        with os.fdopen(descriptor, "rb") as stream:
            before = os.fstat(stream.fileno())
            if not stat.S_ISREG(before.st_mode) or before.st_size != self.size_bytes:
                raise WgsValidationAuthorityError("public validation input changed")
            digest = hashlib.sha256()
            for block in iter(lambda: stream.read(65536), b""):
                digest.update(block)
            after = os.fstat(stream.fileno())
        stable = lambda value: (value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns)
        if stable(before) != stable(after) or digest.hexdigest() != self.sha256:
            raise WgsValidationAuthorityError("public validation input changed")


@dataclass(frozen=True, slots=True)
class PublicValidationAuthority:
    case_id: str
    runner_view_sha256: str
    runner_view: RankerCaseView
    inputs: tuple[PublicValidationInput, ...]
    runner_resource_sha256s: tuple[str, ...]
    policy_family_sha256: str
    candidate_bindings: tuple[PolicyBindings, ...]
    reference_identity_id: str
    reference_manifest_sha256: str
    resource_manifest_sha256s: tuple[str, ...]
    executable_stack_sha256: str
    runtime_image_digest: str
    output_root: str
    preparation_run_options_sha256: str
    completed_source_input: CompletedSourceArtifact | None = None
    schema: Literal["genosapiens.wgs-public-validation-authority"] = SCHEMA
    version: Literal[3] = 3
    classification: Literal["public"] = "public"
    controlled_use_authorized: Literal[False] = False
    competitive_output_allowed: Literal[False] = False

    def __post_init__(self):
        for value in (self.runner_view_sha256, self.policy_family_sha256,
                self.reference_manifest_sha256, self.executable_stack_sha256,
                self.preparation_run_options_sha256,
                *self.runner_resource_sha256s, *self.resource_manifest_sha256s):
            _sha(value)
        if self.completed_source_input is not None:
            if (type(self.completed_source_input) is not CompletedSourceArtifact
                    or self.completed_source_input.schema != COMPLETED_SOURCE_SCHEMA
                    or sum(row.kind == "vcf" for row in self.inputs) != 1):
                raise WgsValidationAuthorityError("completed source authority is invalid")
            _public_path(self.completed_source_input.path, output=True)
        _public_path(self.output_root, output=True)
        if (self.schema != SCHEMA or type(self.version) is not int or self.version != 3
                or self.classification != "public" or self.controlled_use_authorized is not False
                or self.competitive_output_allowed is not False
                or type(self.case_id) is not str or not self.case_id
                or type(self.reference_identity_id) is not str or not self.reference_identity_id
                or type(self.runtime_image_digest) is not str
                or re.fullmatch(r"sha256:[0-9a-f]{64}", self.runtime_image_digest) is None
                or not self.inputs or sum(row.kind == "hpo_text" for row in self.inputs) != 1
                or not any(row.kind == "vcf" for row in self.inputs)
                or len({row.sha256 for row in self.inputs}) != len(self.inputs)
                or not self.runner_resource_sha256s or not self.resource_manifest_sha256s
                or self.resource_manifest_sha256s != tuple(sorted(set(self.resource_manifest_sha256s)))
                or not self.candidate_bindings
                or any(type(row) is not PolicyBindings for row in self.candidate_bindings)
                or tuple(row.sha256 for row in self.candidate_bindings)
                    != tuple(sorted({row.sha256 for row in self.candidate_bindings}))):
            raise WgsValidationAuthorityError("invalid public validation authority")
        # Preparation is shared only across candidates with identical preparation
        # and diagnostic policies. Ranking remains bound separately downstream.
        fields = ("common_allele_policy_sha256", "common_allele_exception_manifest_sha256",
                  "common_allele_public_index_manifest_sha256", "gnomad_projection_precision_sha256",
                  "bounded_wgs_execution_sha256", "diagnostic_policy_sha256")
        if any(len({getattr(row, field) for row in self.candidate_bindings}) != 1 for field in fields):
            raise WgsValidationAuthorityError("public candidates disagree on preparation authority")
        if (type(self.runner_view) is not RankerCaseView
                or hashlib.sha256(ranker_case_view_as_json(self.runner_view)).hexdigest() != self.runner_view_sha256
                or tuple((row.kind, row.sha256) for row in self.inputs)
                    != tuple((row.kind.value, row.sha256) for row in self.runner_view.inputs)
                or self.runner_resource_sha256s
                    != tuple(row.descriptor_sha256 for row in self.runner_view.resources)):
            raise WgsValidationAuthorityError("public validation runner view changed")
        if len(canonical(asdict(self))) > MAX_MANIFEST_BYTES:
            raise WgsValidationAuthorityError("public validation authority exceeds capacity")

    @property
    def policy_bindings(self):
        return self.candidate_bindings[0]

    def as_bytes(self):
        return canonical(asdict(self))

    @property
    def sha256(self):
        return hashlib.sha256(self.as_bytes()).hexdigest()

    def require_execution(self, *, original_input_path, hpo_input_sha256, output_path,
                          executable_stack_sha256, runtime_image_digest,
                          verify_input_bytes: bool = True,
                          preparation_run_options_sha256: str | None = None):
        original = _public_path(str(original_input_path))
        output = _public_path(str(output_path), output=True)
        row = next((row for row in self.inputs if row.kind == "vcf" and Path(row.path) == original), None)
        if (row is None or next(item.sha256 for item in self.inputs if item.kind == "hpo_text") != hpo_input_sha256
                or not output.is_relative_to(Path(self.output_root))
                or executable_stack_sha256 != self.executable_stack_sha256
                or runtime_image_digest != self.runtime_image_digest
                or type(verify_input_bytes) is not bool
                or preparation_run_options_sha256 is not None and
                preparation_run_options_sha256 != self.preparation_run_options_sha256):
            raise WgsValidationAuthorityError("public validation execution differs from authority")
        if verify_input_bytes:
            for item in self.inputs:
                item.verify()


def load_public_validation_authority(content, *, expected_sha256):
    _sha(expected_sha256)
    if type(content) is not bytes or hashlib.sha256(content).hexdigest() != expected_sha256:
        raise WgsValidationAuthorityError("public validation authority hash changed")
    try:
        result = _typed(_json(content), PublicValidationAuthority, decode=True)
        if result.as_bytes() != content:
            raise ValueError
        return result
    except (ValueError, TypeError, RuntimeError, KeyError):
        raise WgsValidationAuthorityError("invalid public validation authority bytes") from None
