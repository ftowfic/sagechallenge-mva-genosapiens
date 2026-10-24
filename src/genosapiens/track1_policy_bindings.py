"""The single current immutable Track 1 policy/resource/limit contract."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import re
from genosapiens.vcf_contracts import POLICY_SHA256 as VCF_POLICY_SHA256


SHA256 = re.compile(r"^[0-9a-f]{64}$")


class Track1PolicyBindingError(ValueError):
    pass


POLICY_BINDING_SCHEMA = "genosapiens.track1-policy-bindings"
POLICY_BINDING_VERSION = 3


@dataclass(frozen=True, slots=True)
class PolicyBindings:
    admission_policy_sha256: str
    hypothesis_policy_sha256: str
    fusion_policy_sha256: str
    diagnosis_prior_policy_sha256: str
    dominant_sensitivity_policy_sha256: str
    diagnostic_policy_sha256: str
    epcr_policy_sha256: str
    common_allele_policy_sha256: str
    common_allele_exception_manifest_sha256: str
    common_allele_public_index_manifest_sha256: str
    gnomad_projection_precision_sha256: str
    bounded_wgs_execution_sha256: str
    schema: str = POLICY_BINDING_SCHEMA
    version: int = POLICY_BINDING_VERSION
    vcf_evidence_policy_sha256: str = VCF_POLICY_SHA256

    def __post_init__(self) -> None:
        hashes = (
            self.admission_policy_sha256,
            self.hypothesis_policy_sha256,
            self.fusion_policy_sha256,
            self.diagnosis_prior_policy_sha256,
            self.dominant_sensitivity_policy_sha256,
            self.diagnostic_policy_sha256,
            self.epcr_policy_sha256,
            self.common_allele_policy_sha256,
            self.common_allele_exception_manifest_sha256,
            self.common_allele_public_index_manifest_sha256,
            self.gnomad_projection_precision_sha256,
            self.bounded_wgs_execution_sha256,
        )
        if (
            self.schema != POLICY_BINDING_SCHEMA
            or type(self.version) is not int
            or self.version != POLICY_BINDING_VERSION
            or self.vcf_evidence_policy_sha256 != VCF_POLICY_SHA256
            or any(type(value) is not str or SHA256.fullmatch(value) is None for value in hashes)
        ):
            raise Track1PolicyBindingError("invalid bounded Track 1 policy binding")

    @property
    def sha256(self) -> str:
        content = (
            json.dumps(asdict(self), sort_keys=True, separators=(",", ":")) + "\n"
        ).encode("ascii")
        return hashlib.sha256(content).hexdigest()

    @property
    def component_sha256s(self) -> tuple[str, ...]:
        """Every separately published current policy/resource/limit component."""

        return (
            self.common_allele_policy_sha256,
            self.common_allele_exception_manifest_sha256,
            self.common_allele_public_index_manifest_sha256,
            self.gnomad_projection_precision_sha256,
            self.bounded_wgs_execution_sha256,
            self.vcf_evidence_policy_sha256,
        )


POLICY_BINDING_FIELDS = frozenset(PolicyBindings.__dataclass_fields__)


def decode_policy_bindings(value: object) -> PolicyBindings:
    """Reject obsolete, partial and extended contracts without version dispatch."""

    if type(value) is not dict:
        raise Track1PolicyBindingError("invalid Track 1 policy binding version")
    fields = frozenset(value)
    try:
        if fields == POLICY_BINDING_FIELDS:
            return PolicyBindings(**value)
    except (TypeError, ValueError):
        pass
    raise Track1PolicyBindingError("invalid Track 1 policy binding version")
