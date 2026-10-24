"""Installed executable, reference, and policy versions for ordinary analysis."""

from dataclasses import asdict, dataclass
import hashlib
import json
import re

from genosapiens.track1_policy_bindings import PolicyBindings, decode_policy_bindings
from genosapiens.track1_release import BaseRelease, Track1ReleaseError, load_base_release


SCHEMA = "genosapiens.track1-installation-versions"


def canonical(value):
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode("ascii")


def digest(content):
    return hashlib.sha256(content).hexdigest()


@dataclass(frozen=True, slots=True)
class InstalledVersions:
    reference_identity_id: str
    reference_manifest_sha256: str
    reference_bindings_sha256: str
    assessment_scope_sha256: str
    policy_bindings: PolicyBindings
    evaluator_pin_sha256: str
    executable_stack_sha256: str
    runtime_image_digest: str
    resource_manifest_sha256s: tuple[str, ...]
    evidence_interface_sha256: str
    canonical_handoff_schema_sha256: str
    canonical_handoff_path_sha256: str
    path_selection_policy_sha256: str
    schema: str = SCHEMA
    version: int = 1
    def __post_init__(self):
        if (self.schema != SCHEMA or type(self.version) is not int or self.version != 1
                or type(self.policy_bindings) is not PolicyBindings
                or not self.reference_identity_id
                or re.fullmatch(r"sha256:[0-9a-f]{64}", self.runtime_image_digest) is None
                or self.resource_manifest_sha256s != tuple(sorted(set(self.resource_manifest_sha256s)))
                or not self.resource_manifest_sha256s):
            raise Track1ReleaseError("invalid installation versions")
        for name, value in asdict(self).items():
            if name.endswith("_sha256") and (type(value) is not str or re.fullmatch(r"[0-9a-f]{64}", value) is None):
                raise Track1ReleaseError("invalid installation identity")
        if any(type(value) is not str or re.fullmatch(r"[0-9a-f]{64}", value) is None
               for value in self.resource_manifest_sha256s):
            raise Track1ReleaseError("invalid installation resource identity")

    require_analysis_profile = BaseRelease.require_analysis_profile
    require_focused_work = BaseRelease.require_focused_work
    ranking_policy_sha256 = BaseRelease.ranking_policy_sha256
    epcr_policy_sha256 = BaseRelease.epcr_policy_sha256

    def as_bytes(self):
        return canonical(asdict(self))

    @property
    def sha256(self):
        return digest(self.as_bytes())


def load_run_versions(content, *, expected_sha256):
    """Read current versions or immutable historical accepted-release evidence."""
    try:
        data = json.loads(content)
        if type(data) is not dict or data.get("schema") != SCHEMA:
            return load_base_release(content, expected_sha256=expected_sha256)
        if (type(content) is not bytes or digest(content) != expected_sha256
                or set(data) != set(InstalledVersions.__dataclass_fields__)
                or canonical(data) != content):
            raise ValueError
        data["policy_bindings"] = decode_policy_bindings(data["policy_bindings"])
        data["resource_manifest_sha256s"] = tuple(data["resource_manifest_sha256s"])
        result = InstalledVersions(**data)
        if result.as_bytes() != content:
            raise ValueError
        return result
    except (TypeError, ValueError, KeyError):
        raise Track1ReleaseError("invalid run versions") from None
