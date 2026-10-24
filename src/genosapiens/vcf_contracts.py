"""Versioned source evidence and bounded VCF reading contracts."""

from dataclasses import asdict, dataclass
from enum import StrEnum
import hashlib
import json
import math
import re


SOURCE_PAYLOAD_VERSION = 2

DEFAULT_MAX_HEADER_BYTES = 8 * 1024 * 1024
DEFAULT_MAX_RECORD_BYTES = 4 * 1024 * 1024
POLICY = {
    "schema": "genosapiens.vcf-evidence-policy.v1",
    "phase": "separate-GT-PS-and-PGT-PID; source-bound-producer-contract-for-PID",
    "measurements": "independent-AD-and-FORMAT-DP; unexplained-AD-total-above-DP-unresolved",
    "balance": "ALT-AD-divided-by-FORMAT-DP; no-clipping; existing-thresholds",
    "quality": "track1-exact-record-PASS-and-existing-sample-class-checks",
    "source_states": "absent-key,omitted-value,dot,empty,malformed,present",
    "reads": "bounded-before-decode-split-fetch; cumulative-header; per-record-body-and-reference",
    "numbers": "original-tokens-and-bytes-authoritative; typed-binary64-is-a-conversion",
}
POLICY_BYTES = (json.dumps(POLICY, sort_keys=True, separators=(",", ":")) + "\n").encode("ascii")
POLICY_SHA256 = hashlib.sha256(POLICY_BYTES).hexdigest()


class VcfContractError(ValueError):
    """Reject an incompatible evidence contract without source data in logs."""


class FieldState(StrEnum):
    ABSENT = "absent-key"
    OMITTED = "omitted-value"
    DOT = "dot"
    EMPTY = "empty"
    MALFORMED = "malformed"
    PRESENT = "present"


class PhaseState(StrEnum):
    ABSENT = "absent"
    INCOMPLETE = "incomplete"
    CONFLICTING = "conflicting"
    UNRESOLVED = "unresolved"
    USABLE = "usable"


@dataclass(frozen=True, slots=True)
class FormatFieldEvidence:
    key: str
    raw_value: str | None
    state: FieldState
    consumed: bool

    def __post_init__(self):
        if (type(self.key) is not str or not self.key or ":" in self.key
                or type(self.consumed) is not bool
                or self.raw_value is not None and type(self.raw_value) is not str
                or self.state != field_state(self.raw_value)):
            raise VcfContractError("invalid FORMAT source state")


def field_state(value):
    if value is None:
        return FieldState.OMITTED
    if value == ".":
        return FieldState.DOT
    if value == "":
        return FieldState.EMPTY
    return FieldState.PRESENT


def expected_cardinality(number, alternate_count, ploidy=None):
    """Return None when the declaration or unknown ploidy leaves count open."""
    if type(alternate_count) is not int or alternate_count < 0:
        raise VcfContractError("invalid allele cardinality")
    if number == "A":
        return alternate_count
    if number == "R":
        return alternate_count + 1
    if number == "G":
        if ploidy is None:
            return None
        if type(ploidy) is not int or ploidy <= 0:
            raise VcfContractError("invalid genotype cardinality")
        return math.comb(alternate_count + ploidy, ploidy)
    if number == ".":
        return None
    if type(number) is str and number.isascii() and number.isdigit():
        return int(number)
    raise VcfContractError("invalid field cardinality")


@dataclass(frozen=True, slots=True)
class VcfReadLimits:
    max_header_bytes: int = DEFAULT_MAX_HEADER_BYTES
    max_record_bytes: int = DEFAULT_MAX_RECORD_BYTES
    max_reference_bytes: int = DEFAULT_MAX_RECORD_BYTES

    def __post_init__(self):
        if any(type(value) is not int or value <= 0 for value in asdict(self).values()):
            raise VcfContractError("invalid VCF read limits")
        if self.max_reference_bytes > self.max_record_bytes:
            raise VcfContractError("VCF reference limit exceeds record limit")


@dataclass(frozen=True, slots=True)
class PidPhaseContract:
    """Explicit producer semantics, never inferred from header declarations."""

    input_sha256: str
    producer: str
    version: str
    documentation: str
    policy_sha256: str = POLICY_SHA256
    namespace: str = "PGT/PID"
    block_scope: str = "input/sample/contig"
    gt_relationship: str = "complete-concordant-GT-multiset"

    def __post_init__(self):
        if (type(self.input_sha256) is not str
                or re.fullmatch(r"[0-9a-f]{64}", self.input_sha256) is None
                or any(type(value) is not str or not value.strip()
                       for value in (self.producer, self.version, self.documentation))
                or self.policy_sha256 != POLICY_SHA256
                or self.namespace != "PGT/PID"
                or self.block_scope != "input/sample/contig"
                or self.gt_relationship != "complete-concordant-GT-multiset"):
            raise VcfContractError("invalid PID phase contract")

    @property
    def sha256(self):
        return hashlib.sha256(json.dumps(asdict(self), sort_keys=True,
                                        separators=(",", ":")).encode("ascii")).hexdigest()


@dataclass(frozen=True, slots=True)
class ReferenceDeclaration:
    raw_value: str
    resolution: str

    def __post_init__(self):
        if (type(self.raw_value) is not str or not self.raw_value
                or any(c in self.raw_value for c in "\r\n\t")
                or self.resolution not in {"matched-operational-basename", "unresolved-metadata"}):
            raise VcfContractError("invalid reference declaration evidence")


def require_current_policy(value):
    if value != POLICY_SHA256:
        raise VcfContractError("incompatible VCF evidence policy; fresh affected output required")


def bounded_line(source, limit, *, header, error_type=VcfContractError):
    """Enforce the byte budget before any decoding or token allocation."""
    if type(limit) is not int or limit < 0:
        raise error_type("invalid VCF read capacity")
    content = source.readline(limit + 1)
    if len(content) > limit:
        label = "header" if header else "record"
        raise error_type(f"VCF {label} exceeds its bounded size (capacity failure)")
    if content and not content.endswith(b"\n"):
        raise error_type("VCF input is truncated")
    return content


def bounded_vcf_lines(source, limits, *, error_type=VcfContractError):
    header, used = True, 0
    while True:
        raw = bounded_line(source, limits.max_header_bytes - used if header
                           else limits.max_record_bytes, header=header, error_type=error_type)
        if not raw:
            return
        if header:
            used += len(raw)
            # A malformed first body line must also respect the record budget.
            if not raw.startswith(b"#"):
                if len(raw) > limits.max_record_bytes:
                    raise error_type("VCF record exceeds its bounded size (capacity failure)")
                header = False
            elif raw.startswith(b"#CHROM"):
                header = False
        yield raw
