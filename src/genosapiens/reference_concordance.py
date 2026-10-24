"""Reference identity contracts and the mandatory read-rescue safety gate."""

from __future__ import annotations

from collections.abc import Callable, Iterator, Mapping
from contextlib import ExitStack, contextmanager
from dataclasses import asdict, dataclass
import gzip
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
from typing import Any
from urllib.parse import unquote, urlparse
import uuid
import warnings
from genosapiens.operator_warnings import OperatorInputWarning
import zlib

from genosapiens.vcf_contracts import ReferenceDeclaration, VcfReadLimits, bounded_vcf_lines

from genosapiens.manifest import (
    ManifestError,
    PublicResourceManifest,
)


REFERENCE_IDENTITY_SCHEMA = "genosapiens.reference-identity"
CONCORDANCE_REPORT_SCHEMA = "genosapiens.reference-concordance"
PUBLIC_VALIDATION_CONCORDANCE_REPORT_SCHEMA = (
    "genosapiens.public-validation-reference-concordance"
)
CONTIG_ORDER_ADJUDICATION_SCHEMA = "genosapiens.contig-order-adjudication"
CONSTRUCTION_REPORT_SCHEMA = "genosapiens.reference-construction-report"
IDENTIFIER = re.compile(r"^[a-z0-9][a-z0-9._-]*$")
MD5 = re.compile(r"^[0-9a-f]{32}$")
SHA256 = re.compile(r"^[0-9a-f]{64}$")
REFERENCE_ALLELE = re.compile(r"^[ACGTRYSWKMBDHVN]+$")
PUBLIC_GIAB_REFERENCE_BASENAMES = {
    "GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta",
    "GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz",
    "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
    "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz",
}
ORDER_ADJUDICATION_RESOLUTIONS = {
    "vcf-header-reordered-without-record-change",
    "reconstructed-reference-order-corrected",
}
EXPECTED_SOURCE_ROLES = {
    "masked_component_fasta_gzip",
    "masked_component_fai",
    "masked_component_gzi",
    "unmasked_decoy_source_fasta_gzip",
    "unmasked_decoy_source_fai",
    "mask_provenance_bed",
}
EXPECTED_REFERENCE_ROLES = {
    "challenge_alignment",
    "challenge_tool_view",
    "public_annotation_benchmark",
}
EXPECTED_OUTPUT_ROLES = {
    "public_native_fasta",
    "public_native_fai",
    "public_native_dict",
    "public_native_m5",
    "challenge_tool_fasta",
    "challenge_tool_fai",
    "challenge_tool_dict",
    "challenge_tool_m5",
    "challenge_tool_bwa_0123",
    "challenge_tool_bwa_amb",
    "challenge_tool_bwa_ann",
    "challenge_tool_bwa_bwt",
    "challenge_tool_bwa_pac",
    "challenge_alignment_fasta",
    "challenge_alignment_fai",
    "challenge_alignment_dict",
    "challenge_alignment_m5",
    "challenge_alignment_bwa_0123",
    "challenge_alignment_bwa_amb",
    "challenge_alignment_bwa_ann",
    "challenge_alignment_bwa_bwt",
    "challenge_alignment_bwa_pac",
    "contig_aliases",
    "mask_delta",
    "construction_report",
}


class ReferenceConcordanceError(RuntimeError):
    """Raised when reference identity or rescue authorization is invalid."""


def _reject_unknown(data: Mapping[str, Any], allowed: set[str], context: str) -> None:
    unknown = sorted(set(data) - allowed)
    if unknown:
        raise ReferenceConcordanceError(
            f"unknown {context} fields: {', '.join(unknown)}"
        )


def _required_string(data: Mapping[str, Any], key: str, context: str) -> str:
    value = data.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ReferenceConcordanceError(f"{context}.{key} must be a non-empty string")
    return value


def _required_integer(data: Mapping[str, Any], key: str, context: str) -> int:
    value = data.get(key)
    if not isinstance(value, int) or isinstance(value, bool) or value < 0:
        raise ReferenceConcordanceError(
            f"{context}.{key} must be a non-negative integer"
        )
    return value


def _require_digest(value: str, pattern: re.Pattern[str], label: str) -> None:
    if not pattern.fullmatch(value):
        raise ReferenceConcordanceError(f"{label} is invalid")


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as source:
            while chunk := source.read(1024 * 1024):
                digest.update(chunk)
    except OSError as error:
        raise ReferenceConcordanceError("cannot hash required file") from error
    return digest.hexdigest()


@dataclass(frozen=True, slots=True)
class IdentitySourceArtifact:
    role: str
    url: str
    path: str
    purpose: str
    size_bytes: int
    md5: str
    sha256: str
    publisher_md5: str | None

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "IdentitySourceArtifact":
        context = "source_artifact"
        _reject_unknown(
            data,
            {
                "role",
                "url",
                "path",
                "purpose",
                "size_bytes",
                "md5",
                "sha256",
                "publisher_md5",
            },
            context,
        )
        publisher_md5 = data.get("publisher_md5")
        if publisher_md5 is not None and not isinstance(publisher_md5, str):
            raise ReferenceConcordanceError(
                "source_artifact.publisher_md5 must be a string or null"
            )
        artifact = cls(
            role=_required_string(data, "role", context),
            url=_required_string(data, "url", context),
            path=_required_string(data, "path", context),
            purpose=_required_string(data, "purpose", context),
            size_bytes=_required_integer(data, "size_bytes", context),
            md5=_required_string(data, "md5", context),
            sha256=_required_string(data, "sha256", context),
            publisher_md5=publisher_md5,
        )
        if not artifact.url.startswith("https://"):
            raise ReferenceConcordanceError("source artifact URL must use HTTPS")
        path = PurePosixPath(artifact.path)
        if path.is_absolute() or ".." in path.parts or path == PurePosixPath("."):
            raise ReferenceConcordanceError(
                "source artifact path must be safe and relative"
            )
        if artifact.purpose not in {"construction", "validation", "provenance"}:
            raise ReferenceConcordanceError("source artifact purpose is invalid")
        _require_digest(artifact.md5, MD5, "source artifact MD5")
        _require_digest(artifact.sha256, SHA256, "source artifact SHA-256")
        if artifact.publisher_md5 is not None:
            _require_digest(
                artifact.publisher_md5, MD5, "source artifact publisher MD5"
            )
            if artifact.publisher_md5 != artifact.md5:
                raise ReferenceConcordanceError(
                    "source artifact local and publisher MD5 values differ"
                )
        return artifact


@dataclass(frozen=True, slots=True)
class IdentityOutput:
    role: str
    path: str

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "IdentityOutput":
        context = "expected_output"
        _reject_unknown(data, {"role", "path"}, context)
        output = cls(
            role=_required_string(data, "role", context),
            path=_required_string(data, "path", context),
        )
        path = PurePosixPath(output.path)
        if path.is_absolute() or ".." in path.parts or path == PurePosixPath("."):
            raise ReferenceConcordanceError("expected output path must be safe and relative")
        return output


@dataclass(frozen=True, slots=True)
class ReferenceRole:
    role: str
    identity_id: str
    naming: str
    sequence_relationship: str
    fasta_output_role: str
    fai_output_role: str
    dictionary_output_role: str
    m5_output_role: str

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "ReferenceRole":
        context = "reference_role"
        expected = {
            "role",
            "identity_id",
            "naming",
            "sequence_relationship",
            "fasta_output_role",
            "fai_output_role",
            "dictionary_output_role",
            "m5_output_role",
        }
        _reject_unknown(data, expected, context)
        role = cls(**{key: _required_string(data, key, context) for key in expected})
        if role.role not in EXPECTED_REFERENCE_ROLES:
            raise ReferenceConcordanceError(f"unknown typed reference role: {role.role}")
        if not IDENTIFIER.fullmatch(role.identity_id):
            raise ReferenceConcordanceError("reference role identity_id is invalid")
        expected_values = {
            "challenge_alignment": (
                "unprefixed",
                "constructed-masked-v2-plus-hs38d1",
                "challenge_alignment",
            ),
            "challenge_tool_view": (
                "chr-prefixed",
                "sequence-identical-to-challenge-alignment",
                "challenge_tool",
            ),
            "public_annotation_benchmark": (
                "resource-native",
                "coordinate-compatible-not-sequence-identical",
                "public_native",
            ),
        }
        naming, relationship, output_prefix = expected_values[role.role]
        if role.naming != naming or role.sequence_relationship != relationship:
            raise ReferenceConcordanceError(
                f"typed reference role semantics are invalid: {role.role}"
            )
        actual_outputs = (
            role.fasta_output_role,
            role.fai_output_role,
            role.dictionary_output_role,
            role.m5_output_role,
        )
        expected_outputs = tuple(
            f"{output_prefix}_{suffix}" for suffix in ("fasta", "fai", "dict", "m5")
        )
        if actual_outputs != expected_outputs:
            raise ReferenceConcordanceError(
                f"typed reference role outputs are invalid: {role.role}"
            )
        return role


@dataclass(frozen=True, slots=True)
class ReferenceRoleBinding:
    role: str
    identity_id: str

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "ReferenceRoleBinding":
        context = "reference_binding"
        _reject_unknown(data, {"role", "identity_id"}, context)
        binding = cls(
            role=_required_string(data, "role", context),
            identity_id=_required_string(data, "identity_id", context),
        )
        if binding.role not in EXPECTED_REFERENCE_ROLES:
            raise ReferenceConcordanceError(
                f"unknown typed reference role: {binding.role}"
            )
        if not IDENTIFIER.fullmatch(binding.identity_id):
            raise ReferenceConcordanceError("reference binding identity_id is invalid")
        return binding

    def resolve(
        self,
        contract: "ReferenceIdentityContract",
        *,
        required_role: str,
    ) -> ReferenceRole:
        if self.role != required_role:
            raise ReferenceConcordanceError(
                f"reference binding must use role {required_role}"
            )
        role = contract.reference_role(self.role)
        if self.identity_id != role.identity_id:
            raise ReferenceConcordanceError(
                f"reference binding identity does not match role {self.role}"
            )
        return role


@dataclass(frozen=True, slots=True)
class MaskSentinel:
    contig: str
    position: int

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "MaskSentinel":
        _reject_unknown(data, {"contig", "position"}, "mask_sentinel")
        sentinel = cls(
            contig=_required_string(data, "contig", "mask_sentinel"),
            position=_required_integer(data, "position", "mask_sentinel"),
        )
        if sentinel.position <= 0 or not sentinel.contig.startswith("chr"):
            raise ReferenceConcordanceError("mask sentinel is invalid")
        return sentinel


@dataclass(frozen=True, slots=True)
class IdentityTransformation:
    schema: str
    masked_component_records: int
    masked_component_terminal: str
    decoy_suffix_records: int
    decoy_name_suffix: str
    source_prefix: str
    scope: str
    preserve_sequence_bytes: bool
    preserve_record_order: bool
    mask_differences: str
    mask_sentinels: tuple[MaskSentinel, ...]

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "IdentityTransformation":
        context = "transformation"
        allowed = {
            "schema",
            "masked_component_records",
            "masked_component_terminal",
            "decoy_suffix_records",
            "decoy_name_suffix",
            "source_prefix",
            "scope",
            "preserve_sequence_bytes",
            "preserve_record_order",
            "mask_differences",
            "mask_sentinels",
        }
        _reject_unknown(data, allowed, context)
        for key in ("preserve_sequence_bytes", "preserve_record_order"):
            if not isinstance(data.get(key), bool):
                raise ReferenceConcordanceError(f"transformation.{key} must be a boolean")
        sentinel_data = data.get("mask_sentinels")
        if not isinstance(sentinel_data, list) or not sentinel_data or not all(
            isinstance(item, Mapping) for item in sentinel_data
        ):
            raise ReferenceConcordanceError(
                "transformation.mask_sentinels must be a non-empty object array"
            )
        transformation = cls(
            schema=_required_string(data, "schema", context),
            masked_component_records=_required_integer(
                data, "masked_component_records", context
            ),
            masked_component_terminal=_required_string(
                data, "masked_component_terminal", context
            ),
            decoy_suffix_records=_required_integer(data, "decoy_suffix_records", context),
            decoy_name_suffix=_required_string(data, "decoy_name_suffix", context),
            source_prefix=_required_string(data, "source_prefix", context),
            scope=_required_string(data, "scope", context),
            preserve_sequence_bytes=data["preserve_sequence_bytes"],
            preserve_record_order=data["preserve_record_order"],
            mask_differences=_required_string(data, "mask_differences", context),
            mask_sentinels=tuple(
                MaskSentinel.from_mapping(item) for item in sentinel_data
            ),
        )
        expected = (
            transformation.schema
            == "genosapiens.masked-v2-plus-hs38d1-construction"
            and transformation.masked_component_records == 195
            and transformation.masked_component_terminal == "chrEBV"
            and transformation.decoy_suffix_records == 2385
            and transformation.decoy_name_suffix == "_decoy"
            and transformation.source_prefix == "chr"
            and transformation.scope == "first-fasta-header-token"
            and transformation.preserve_sequence_bytes
            and transformation.preserve_record_order
            and transformation.mask_differences == "unmasked-base-to-N-only"
        )
        if not expected:
            raise ReferenceConcordanceError(
                "reference identity requires the pinned masked-v2 construction"
            )
        return transformation


@dataclass(frozen=True, slots=True)
class AcceptanceRequirements:
    public_construction_validation_required: bool
    controlled_complete_concordance_required: bool
    operator_override_can_bypass: bool

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "AcceptanceRequirements":
        expected = {
            "public_construction_validation_required",
            "controlled_complete_concordance_required",
            "operator_override_can_bypass",
        }
        _reject_unknown(data, expected, "acceptance")
        if any(not isinstance(data.get(key), bool) for key in expected):
            raise ReferenceConcordanceError("acceptance fields must be booleans")
        acceptance = cls(**{key: data[key] for key in expected})
        if (
            not acceptance.public_construction_validation_required
            or not acceptance.controlled_complete_concordance_required
            or acceptance.operator_override_can_bypass
        ):
            raise ReferenceConcordanceError(
                "reference acceptance prerequisites must fail closed"
            )
        return acceptance


@dataclass(frozen=True, slots=True)
class ReferenceIdentityContract:
    contract_id: str
    resource_id: str
    release: str
    assembly: str
    challenge_reference_basename: str
    primary_assembly_accession: str
    decoy_assembly_accession: str
    source_artifacts: tuple[IdentitySourceArtifact, ...]
    reference_roles: tuple[ReferenceRole, ...]
    transformation: IdentityTransformation
    acceptance: AcceptanceRequirements
    expected_outputs: tuple[IdentityOutput, ...]
    sequence_count: int
    canonical_contig_lengths: Mapping[str, int]
    schema: str = REFERENCE_IDENTITY_SCHEMA

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "ReferenceIdentityContract":
        allowed = {
            "schema",
            "contract_id",
            "resource_id",
            "release",
            "assembly",
            "challenge_reference_basename",
            "primary_assembly_accession",
            "decoy_assembly_accession",
            "source_artifacts",
            "reference_roles",
            "transformation",
            "acceptance",
            "expected_outputs",
            "sequence_count",
            "canonical_contig_lengths",
        }
        _reject_unknown(data, allowed, "reference identity")
        if data.get("schema") != REFERENCE_IDENTITY_SCHEMA:
            raise ReferenceConcordanceError(
                f"reference identity schema must be {REFERENCE_IDENTITY_SCHEMA!r}"
            )
        source_data = data.get("source_artifacts")
        role_data = data.get("reference_roles")
        output_data = data.get("expected_outputs")
        transformation_data = data.get("transformation")
        acceptance_data = data.get("acceptance")
        contig_data = data.get("canonical_contig_lengths")
        if not isinstance(source_data, list) or not source_data:
            raise ReferenceConcordanceError("source_artifacts must be a non-empty array")
        if not isinstance(output_data, list) or not output_data:
            raise ReferenceConcordanceError("expected_outputs must be a non-empty array")
        if not isinstance(role_data, list) or not role_data or not all(
            isinstance(item, Mapping) for item in role_data
        ):
            raise ReferenceConcordanceError("reference_roles must be an object array")
        if not isinstance(transformation_data, Mapping):
            raise ReferenceConcordanceError("transformation must be an object")
        if not isinstance(acceptance_data, Mapping):
            raise ReferenceConcordanceError("acceptance must be an object")
        if not isinstance(contig_data, Mapping) or not contig_data:
            raise ReferenceConcordanceError(
                "canonical_contig_lengths must be a non-empty object"
            )
        contigs: dict[str, int] = {}
        for name, length in contig_data.items():
            if (
                not isinstance(name, str)
                or not name
                or not isinstance(length, int)
                or isinstance(length, bool)
                or length <= 0
            ):
                raise ReferenceConcordanceError(
                    "canonical contig names and lengths must be valid"
                )
            contigs[name] = length
        contract = cls(
            schema=data["schema"],
            contract_id=_required_string(data, "contract_id", "reference identity"),
            resource_id=_required_string(data, "resource_id", "reference identity"),
            release=_required_string(data, "release", "reference identity"),
            assembly=_required_string(data, "assembly", "reference identity"),
            challenge_reference_basename=_required_string(
                data, "challenge_reference_basename", "reference identity"
            ),
            primary_assembly_accession=_required_string(
                data, "primary_assembly_accession", "reference identity"
            ),
            decoy_assembly_accession=_required_string(
                data, "decoy_assembly_accession", "reference identity"
            ),
            source_artifacts=tuple(
                IdentitySourceArtifact.from_mapping(item) for item in source_data
            ),
            reference_roles=tuple(
                ReferenceRole.from_mapping(item) for item in role_data
            ),
            transformation=IdentityTransformation.from_mapping(transformation_data),
            acceptance=AcceptanceRequirements.from_mapping(acceptance_data),
            expected_outputs=tuple(
                IdentityOutput.from_mapping(item) for item in output_data
            ),
            sequence_count=_required_integer(
                data, "sequence_count", "reference identity"
            ),
            canonical_contig_lengths=contigs,
        )
        if not IDENTIFIER.fullmatch(contract.contract_id):
            raise ReferenceConcordanceError("reference contract_id is invalid")
        if not IDENTIFIER.fullmatch(contract.resource_id):
            raise ReferenceConcordanceError("reference resource_id is invalid")
        source_roles = [artifact.role for artifact in contract.source_artifacts]
        if set(source_roles) != EXPECTED_SOURCE_ROLES or len(source_roles) != len(
            set(source_roles)
        ):
            raise ReferenceConcordanceError(
                "source artifacts must define every pinned construction source once"
            )
        source_urls = [artifact.url for artifact in contract.source_artifacts]
        if len(source_urls) != len(set(source_urls)):
            raise ReferenceConcordanceError("source artifact URLs must be unique")
        source_paths = [artifact.path for artifact in contract.source_artifacts]
        if len(source_paths) != len(set(source_paths)):
            raise ReferenceConcordanceError("source artifact paths must be unique")
        roles = [role.role for role in contract.reference_roles]
        if set(roles) != EXPECTED_REFERENCE_ROLES or len(roles) != len(set(roles)):
            raise ReferenceConcordanceError(
                "reference_roles must define the three typed roles once"
            )
        role_ids = [role.identity_id for role in contract.reference_roles]
        if len(role_ids) != len(set(role_ids)):
            raise ReferenceConcordanceError("typed reference identity IDs must be unique")
        output_roles = {output.role for output in contract.expected_outputs}
        if output_roles != EXPECTED_OUTPUT_ROLES:
            raise ReferenceConcordanceError(
                "expected outputs do not define the complete candidate reference bundle"
            )
        output_paths = [output.path for output in contract.expected_outputs]
        if len(output_paths) != len(set(output_paths)):
            raise ReferenceConcordanceError("expected output paths must be unique")
        if contract.sequence_count <= 0:
            raise ReferenceConcordanceError("sequence_count must be positive")
        if (
            contract.sequence_count
            != contract.transformation.masked_component_records
            + contract.transformation.decoy_suffix_records
        ):
            raise ReferenceConcordanceError(
                "sequence_count does not match component composition"
            )
        if (
            contract.challenge_reference_basename
            != "GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_hs38d1_"
            "maskedGRC_exclusions_v2_no_chr.fasta"
        ):
            raise ReferenceConcordanceError(
                "challenge reference basename is not the accepted construction descriptor"
            )
        return contract

    @property
    def identity_id(self) -> str:
        return self.reference_role("challenge_alignment").identity_id

    @property
    def canonical_source(self) -> IdentitySourceArtifact:
        return self.source_artifact("masked_component_fasta_gzip")

    def source_artifact(self, role: str) -> IdentitySourceArtifact:
        for artifact in self.source_artifacts:
            if artifact.role == role:
                return artifact
        raise ReferenceConcordanceError(
            f"reference identity lacks {role} source"
        )

    def expected_output(self, role: str) -> IdentityOutput:
        aliases = {
            "canonical_fasta": "public_native_fasta",
            "canonical_fai": "public_native_fai",
            "challenge_fasta": "challenge_alignment_fasta",
            "challenge_fai": "challenge_alignment_fai",
        }
        role = aliases.get(role, role)
        for output in self.expected_outputs:
            if output.role == role:
                return output
        raise ReferenceConcordanceError(
            f"reference identity lacks {role} expected output"
        )

    def reference_role(self, role: str) -> ReferenceRole:
        for value in self.reference_roles:
            if value.role == role:
                return value
        raise ReferenceConcordanceError(f"reference identity lacks typed role {role}")


def load_reference_identity(path: str | Path) -> ReferenceIdentityContract:
    contract_path = Path(path)
    try:
        data = json.loads(contract_path.read_text(encoding="utf-8"))
    except OSError as error:
        raise ReferenceConcordanceError(
            f"cannot read reference identity {contract_path}: {error}"
        ) from error
    except json.JSONDecodeError as error:
        raise ReferenceConcordanceError(
            f"invalid reference identity JSON {contract_path}: {error}"
        ) from error
    if not isinstance(data, Mapping):
        raise ReferenceConcordanceError("reference identity root must be an object")
    return ReferenceIdentityContract.from_mapping(data)


@dataclass(frozen=True, slots=True)
class ReferenceConstructionReport:
    contract_id: str
    contract_sha256: str
    challenge_alignment_identity_id: str
    challenge_tool_view_identity_id: str
    public_annotation_benchmark_identity_id: str
    masked_component_records: int
    decoy_suffix_records: int
    total_records: int
    mask_delta_intervals: int
    masked_bases: int
    mask_delta_sha256: str
    source_sha256: Mapping[str, str]
    output_sha256: Mapping[str, str]
    accepted: bool
    schema: str = CONSTRUCTION_REPORT_SCHEMA

    def __post_init__(self) -> None:
        if self.schema != CONSTRUCTION_REPORT_SCHEMA:
            raise ReferenceConcordanceError(
                f"construction report schema must be {CONSTRUCTION_REPORT_SCHEMA!r}"
            )
        for label, value in (
            ("contract_id", self.contract_id),
            ("challenge_alignment_identity_id", self.challenge_alignment_identity_id),
            ("challenge_tool_view_identity_id", self.challenge_tool_view_identity_id),
            (
                "public_annotation_benchmark_identity_id",
                self.public_annotation_benchmark_identity_id,
            ),
        ):
            if not IDENTIFIER.fullmatch(value):
                raise ReferenceConcordanceError(
                    f"construction report {label} is invalid"
                )
        for label, digest in (
            ("contract_sha256", self.contract_sha256),
            ("mask_delta_sha256", self.mask_delta_sha256),
        ):
            _require_digest(digest, SHA256, f"construction report {label}")
        for label, value in (
            ("masked_component_records", self.masked_component_records),
            ("decoy_suffix_records", self.decoy_suffix_records),
            ("total_records", self.total_records),
            ("mask_delta_intervals", self.mask_delta_intervals),
            ("masked_bases", self.masked_bases),
        ):
            if not isinstance(value, int) or isinstance(value, bool) or value < 0:
                raise ReferenceConcordanceError(
                    f"construction report {label} must be a non-negative integer"
                )
        for label, values in (
            ("source_sha256", self.source_sha256),
            ("output_sha256", self.output_sha256),
        ):
            if not isinstance(values, Mapping) or not values:
                raise ReferenceConcordanceError(
                    f"construction report {label} must be a non-empty object"
                )
            for role, digest in values.items():
                if not isinstance(role, str) or not role:
                    raise ReferenceConcordanceError(
                        f"construction report {label} role is invalid"
                    )
                if not isinstance(digest, str):
                    raise ReferenceConcordanceError(
                        f"construction report {label} digest is invalid"
                    )
                _require_digest(digest, SHA256, f"construction report {label} digest")
        if not isinstance(self.accepted, bool):
            raise ReferenceConcordanceError(
                "construction report accepted must be a boolean"
            )
        expected_acceptance = (
            self.masked_component_records > 0
            and self.decoy_suffix_records > 0
            and self.total_records
            == self.masked_component_records + self.decoy_suffix_records
            and self.mask_delta_intervals > 0
            and self.masked_bases > 0
        )
        if self.accepted != expected_acceptance:
            raise ReferenceConcordanceError(
                "construction report acceptance conflicts with validation counts"
            )

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "ReferenceConstructionReport":
        expected = {
            "schema",
            "contract_id",
            "contract_sha256",
            "challenge_alignment_identity_id",
            "challenge_tool_view_identity_id",
            "public_annotation_benchmark_identity_id",
            "masked_component_records",
            "decoy_suffix_records",
            "total_records",
            "mask_delta_intervals",
            "masked_bases",
            "mask_delta_sha256",
            "source_sha256",
            "output_sha256",
            "accepted",
        }
        _reject_unknown(data, expected, "construction report")
        if set(data) != expected:
            missing = ", ".join(sorted(expected - set(data)))
            raise ReferenceConcordanceError(
                f"construction report is missing fields: {missing}"
            )
        try:
            return cls(**data)
        except TypeError as error:
            raise ReferenceConcordanceError(
                f"invalid construction report: {error}"
            ) from error


def load_construction_report(path: str | Path) -> ReferenceConstructionReport:
    report_path = Path(path)
    try:
        data = json.loads(report_path.read_text(encoding="utf-8"))
    except OSError as error:
        raise ReferenceConcordanceError(
            f"cannot read construction report {report_path}: {error}"
        ) from error
    except json.JSONDecodeError as error:
        raise ReferenceConcordanceError(
            f"invalid construction report JSON {report_path}: {error}"
        ) from error
    if not isinstance(data, Mapping):
        raise ReferenceConcordanceError("construction report root must be an object")
    return ReferenceConstructionReport.from_mapping(data)


@dataclass(frozen=True, slots=True)
class ReferenceIdentification:
    assembly: str
    reference_role: str
    reference_identity_id: str
    reference_basename: str
    contig_names_match: bool
    contig_lengths_match: bool
    contig_order_match: bool


@dataclass(frozen=True, slots=True)
class _VcfScan:
    reference_declarations: tuple[str, ...]
    header_valid: bool
    reference_declaration: str | None
    contigs: Mapping[str, int]
    input_sha256: str
    record_payload_sha256: str
    order_insensitive_header_sha256: str
    total_records: int
    checked_records: int
    matching_records: int
    mismatching_records: int
    unknown_contig_records: int
    out_of_bounds_records: int
    malformed_ref_records: int
    otherwise_unresolvable_records: int


def _mapping_sha256(values: Mapping[str, int]) -> str:
    digest = hashlib.sha256()
    for name, length in values.items():
        digest.update(name.encode("utf-8"))
        digest.update(b"\t")
        digest.update(str(length).encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


def _update_framed_digest(
    digest: Any, record_type: bytes, value: bytes
) -> None:
    digest.update(record_type)
    digest.update(len(value).to_bytes(8, "big"))
    digest.update(value)


def _canonical_json_sha256(value: Any) -> str:
    encoded = (
        json.dumps(value, indent=2, sort_keys=True, ensure_ascii=True) + "\n"
    ).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def _reference_basename(declaration: str) -> str:
    value = declaration.strip()
    if len(value) >= 2 and value[0] == value[-1] == '"':
        value = value[1:-1]
    parsed = urlparse(value)
    path = unquote(parsed.path) if parsed.scheme else value
    return PurePosixPath(path.replace("\\", "/")).name


def _structured_attributes(value: str) -> Mapping[str, str] | None:
    attributes: dict[str, str] = {}
    parts: list[str] = []
    current: list[str] = []
    quoted = False
    escaped = False
    for character in value:
        if escaped:
            current.append(character)
            escaped = False
        elif character == "\\" and quoted:
            current.append(character)
            escaped = True
        elif character == '"':
            current.append(character)
            quoted = not quoted
        elif character == "," and not quoted:
            parts.append("".join(current))
            current = []
        else:
            current.append(character)
    if quoted or escaped:
        return None
    parts.append("".join(current))
    for part in parts:
        if "=" not in part:
            return None
        key, item = part.split("=", 1)
        if not key or key in attributes:
            return None
        if len(item) >= 2 and item[0] == item[-1] == '"':
            item = item[1:-1]
        attributes[key] = item
    return attributes


class _DigestingReader:
    def __init__(self, source) -> None:
        self.source = source
        self.digest = hashlib.sha256()

    def read(self, size: int = -1) -> bytes:
        value = self.source.read(size)
        self.digest.update(value)
        return value

    def readline(self, size: int = -1) -> bytes:
        value = self.source.readline(size)
        self.digest.update(value)
        return value

    def tell(self) -> int:
        return self.source.tell()

    def __iter__(self):
        return self

    def __next__(self) -> bytes:
        value = self.readline()
        if not value:
            raise StopIteration
        return value

    def hexdigest(self) -> str:
        return self.digest.hexdigest()


def _stat_identity(value: os.stat_result) -> tuple[int, int, int, int, int]:
    return (
        value.st_dev,
        value.st_ino,
        value.st_size,
        value.st_mtime_ns,
        value.st_ctime_ns,
    )


@contextmanager
def _open_stable_binary(path: Path, label: str):
    try:
        if path.is_symlink():
            raise ReferenceConcordanceError(f"{label} must be a regular file")
        with path.open("rb") as source:
            initial_fd_identity = _stat_identity(os.fstat(source.fileno()))
            initial_path_identity = _stat_identity(path.stat())
            if initial_fd_identity != initial_path_identity:
                raise ReferenceConcordanceError(f"{label} changed while being opened")
            try:
                yield source
            finally:
                try:
                    final_fd_identity = _stat_identity(os.fstat(source.fileno()))
                    final_path_identity = _stat_identity(path.stat())
                except OSError as error:
                    raise ReferenceConcordanceError(
                        f"{label} changed while being used"
                    ) from error
                if (
                    initial_fd_identity != final_fd_identity
                    or initial_path_identity != final_path_identity
                ):
                    raise ReferenceConcordanceError(
                        f"{label} changed while being used"
                    )
    except OSError as error:
        raise ReferenceConcordanceError(f"cannot read {label}") from error


def _open_file_sha256(source) -> str:
    digest = hashlib.sha256()
    source.seek(0)
    while chunk := source.read(1024 * 1024):
        digest.update(chunk)
    source.seek(0)
    return digest.hexdigest()


@dataclass(frozen=True, slots=True)
class _OpenFastaIndexEntry:
    length: int
    offset: int
    line_bases: int
    line_width: int


def _parse_fai_bytes(content: bytes) -> Mapping[str, _OpenFastaIndexEntry]:
    entries: dict[str, _OpenFastaIndexEntry] = {}
    try:
        lines = content.decode("ascii").splitlines()
    except UnicodeDecodeError as error:
        raise ReferenceConcordanceError("reference FAI is not ASCII") from error
    for line in lines:
        fields = line.split("\t")
        if len(fields) < 5:
            raise ReferenceConcordanceError("reference FAI has an invalid record")
        try:
            entry = _OpenFastaIndexEntry(*(int(value) for value in fields[1:5]))
        except ValueError as error:
            raise ReferenceConcordanceError(
                "reference FAI has a non-integer field"
            ) from error
        name = fields[0]
        if (
            not name
            or entry.length <= 0
            or entry.offset < 0
            or entry.line_bases <= 0
            or entry.line_width < entry.line_bases
            or name in entries
        ):
            raise ReferenceConcordanceError("reference FAI has an invalid entry")
        entries[name] = entry
    if not entries:
        raise ReferenceConcordanceError("reference FAI has no sequences")
    return entries


class _OpenIndexedFasta:
    def __init__(self, source, index: Mapping[str, _OpenFastaIndexEntry]) -> None:
        self.source = source
        self.index = index

    def fetch(self, contig: str, position: int, length: int) -> str:
        if position <= 0 or length <= 0:
            raise ValueError("reference lookup requires positive position and length")
        try:
            entry = self.index[contig]
        except KeyError as error:
            raise KeyError("unknown reference contig") from error
        zero_based = position - 1
        if zero_based + length > entry.length:
            raise IndexError("reference lookup exceeds contig")
        chunks: list[bytes] = []
        remaining = length
        try:
            while remaining:
                line_number, within_line = divmod(zero_based, entry.line_bases)
                count = min(remaining, entry.line_bases - within_line)
                byte_offset = (
                    entry.offset + line_number * entry.line_width + within_line
                )
                self.source.seek(byte_offset)
                chunk = self.source.read(count)
                try:
                    decoded = chunk.decode("ascii")
                except UnicodeDecodeError as error:
                    raise ReferenceConcordanceError(
                        "invalid indexed reference sequence bytes"
                    ) from error
                if len(chunk) != count or not REFERENCE_ALLELE.fullmatch(decoded):
                    raise ReferenceConcordanceError(
                        "invalid indexed reference sequence bytes"
                    )
                chunks.append(chunk)
                zero_based += count
                remaining -= count
        except OSError as error:
            raise ReferenceConcordanceError(
                "cannot read indexed reference FASTA"
            ) from error
        return b"".join(chunks).decode("ascii")


@contextmanager
def _open_vcf(path: Path):
    try:
        with path.open("rb") as raw_source:
            initial_fd_identity = _stat_identity(os.fstat(raw_source.fileno()))
            initial_path_identity = _stat_identity(path.stat())
            compressed = raw_source.read(2) == b"\x1f\x8b"
            raw_source.seek(0)
            digesting_source = _DigestingReader(raw_source)
            source = (
                gzip.GzipFile(fileobj=digesting_source, mode="rb")
                if compressed
                else digesting_source
            )
            completed = False
            try:
                yield source, digesting_source
                completed = True
            finally:
                if source is not digesting_source:
                    source.close()
                if completed:
                    while digesting_source.read(1024 * 1024):
                        pass
                final_fd_identity = _stat_identity(os.fstat(raw_source.fileno()))
                final_path_identity = _stat_identity(path.stat())
                if (
                    initial_fd_identity != final_fd_identity
                    or initial_path_identity != final_path_identity
                    or initial_fd_identity != initial_path_identity
                ):
                    raise ReferenceConcordanceError(
                        "VCF changed while reference concordance was being evaluated"
                    )
    except OSError as error:
        raise ReferenceConcordanceError("cannot read VCF") from error


def _scan_vcf(
    path: Path,
    reference_contigs: Mapping[str, int],
    fetch_reference: Callable[[str, int, int], str],
    *,
    challenge_reference_basename: str | None = None,
    read_limits: VcfReadLimits = VcfReadLimits(),
) -> _VcfScan:
    record_digest = hashlib.sha256()
    header_digest = hashlib.sha256()
    contig_header_lines: list[bytes] = []
    reference_declarations: list[str] = []
    selected_reference = None
    contigs: dict[str, int] = {}
    header_valid = True
    fileformat_seen = False
    column_header_seen = False
    total_records = 0
    checked_records = 0
    matching_records = 0
    mismatching_records = 0
    unknown_contig_records = 0
    out_of_bounds_records = 0
    malformed_ref_records = 0
    otherwise_unresolvable_records = 0
    line_number = 0

    try:
        with _open_vcf(path) as (source, raw_source):
            for raw_line in bounded_vcf_lines(source, read_limits, error_type=ReferenceConcordanceError):
                line_number += 1
                line = raw_line.rstrip(b"\r\n")
                if not column_header_seen and line.startswith(b"##"):
                    if line.startswith(b"##contig="):
                        contig_header_lines.append(raw_line)
                    else:
                        _update_framed_digest(header_digest, b"H", raw_line)
                    if line.startswith(b"##fileformat=VCFv4."):
                        if line_number != 1 or fileformat_seen:
                            header_valid = False
                        fileformat_seen = True
                    if line.startswith(b"##reference="):
                        try:
                            reference_declarations.append(
                                line.removeprefix(b"##reference=").decode("ascii")
                            )
                        except UnicodeDecodeError:
                            header_valid = False
                    elif line.startswith(b"##contig="):
                        try:
                            text = line.decode("ascii")
                        except UnicodeDecodeError:
                            header_valid = False
                            continue
                        if not text.endswith(">") or not text.startswith("##contig=<"):
                            header_valid = False
                            continue
                        attributes = _structured_attributes(text[10:-1])
                        if attributes is None:
                            header_valid = False
                            continue
                        name = attributes.get("ID")
                        length_text = attributes.get("length")
                        try:
                            length = int(length_text) if length_text is not None else 0
                        except ValueError:
                            length = 0
                        if (
                            not name
                            or any(character.isspace() for character in name)
                            or length <= 0
                            or name in contigs
                        ):
                            header_valid = False
                            continue
                        contigs[name] = length
                    continue
                if not column_header_seen and line.startswith(b"#CHROM"):
                    _update_framed_digest(header_digest, b"H", raw_line)
                    columns = line.split(b"\t")
                    if columns[:8] != [
                        b"#CHROM",
                        b"POS",
                        b"ID",
                        b"REF",
                        b"ALT",
                        b"QUAL",
                        b"FILTER",
                        b"INFO",
                    ]:
                        header_valid = False
                    if len(reference_declarations) == 1:
                        selected_reference = reference_declarations[0]
                    elif (len(reference_declarations) > 1
                            and challenge_reference_basename is not None
                            and any(_reference_basename(value) == challenge_reference_basename
                                    for value in reference_declarations)):
                        # The challenge's original header may retain more than
                        # one pipeline reference. Resolve its declared challenge
                        # basename without rewriting bytes or bypassing REF checks.
                        selected_reference = challenge_reference_basename
                        warnings.warn(
                            'The VCF header contains additional reference declarations. '
                            'Unmatched declarations remain unresolved metadata. '
                            'The complete REF check uses the installed challenge reference.',
                            OperatorInputWarning, stacklevel=2,
                        )
                    else:
                        header_valid = False
                    column_header_seen = True
                    continue
                if not column_header_seen and line.startswith(b"#"):
                    _update_framed_digest(header_digest, b"H", raw_line)
                    header_valid = False
                    continue

                if not column_header_seen:
                    header_valid = False
                    column_header_seen = True
                record_digest.update(raw_line)
                total_records += 1
                try:
                    fields = line.decode("ascii").split("\t")
                except UnicodeDecodeError:
                    otherwise_unresolvable_records += 1
                    continue
                if len(fields) < 8:
                    otherwise_unresolvable_records += 1
                    continue
                contig, position_text, ref = fields[0], fields[1], fields[3]
                try:
                    position = int(position_text)
                except ValueError:
                    otherwise_unresolvable_records += 1
                    continue
                if not REFERENCE_ALLELE.fullmatch(ref):
                    malformed_ref_records += 1
                    continue
                if contig not in reference_contigs:
                    unknown_contig_records += 1
                    continue
                if position <= 0 or position + len(ref) - 1 > reference_contigs[contig]:
                    out_of_bounds_records += 1
                    continue
                if len(ref) > read_limits.max_reference_bytes:
                    raise ReferenceConcordanceError("VCF reference request exceeds its bounded size (capacity failure)")
                try:
                    observed = fetch_reference(contig, position, len(ref))
                except (IndexError, KeyError, OSError, TypeError, ValueError):
                    otherwise_unresolvable_records += 1
                    continue
                if not isinstance(observed, str) or len(observed) != len(ref):
                    otherwise_unresolvable_records += 1
                    continue
                checked_records += 1
                if observed == ref:
                    matching_records += 1
                else:
                    mismatching_records += 1
    except (EOFError, gzip.BadGzipFile, OSError, zlib.error) as error:
        raise ReferenceConcordanceError("cannot parse VCF") from error

    if line_number == 0 or not column_header_seen or not fileformat_seen:
        header_valid = False
    for raw_line in sorted(contig_header_lines):
        _update_framed_digest(header_digest, b"C", raw_line)
    return _VcfScan(
        reference_declarations=tuple(reference_declarations),
        header_valid=header_valid,
        reference_declaration=selected_reference,
        contigs=contigs,
        input_sha256=raw_source.hexdigest(),
        record_payload_sha256=record_digest.hexdigest(),
        order_insensitive_header_sha256=header_digest.hexdigest(),
        total_records=total_records,
        checked_records=checked_records,
        matching_records=matching_records,
        mismatching_records=mismatching_records,
        unknown_contig_records=unknown_contig_records,
        out_of_bounds_records=out_of_bounds_records,
        malformed_ref_records=malformed_ref_records,
        otherwise_unresolvable_records=otherwise_unresolvable_records,
    )


def _reference_metadata(scan, basename):
    return tuple(ReferenceDeclaration(value,
        "matched-operational-basename" if _reference_basename(value) == basename else "unresolved-metadata")
        for value in scan.reference_declarations)


def identify_vcf_reference(
    reference_declaration: str,
    declared_contigs: Mapping[str, int],
    challenge_contigs: Mapping[str, int],
    contract: ReferenceIdentityContract,
) -> ReferenceIdentification:
    """Resolve a VCF reference declaration without conflating typed roles."""
    basename = _reference_basename(reference_declaration)
    if (
        "hs37d5" in basename.lower()
        or declared_contigs.get("1") == 249250621
        or declared_contigs.get("chr1") == 249250621
    ):
        raise ReferenceConcordanceError("assembly-mismatch: VCF is GRCh37/b37")

    if basename == contract.challenge_reference_basename:
        declared_names = tuple(declared_contigs)
        challenge_names = tuple(challenge_contigs)
        names_match = set(declared_names) == set(challenge_names)
        lengths_match = names_match and all(
            declared_contigs[name] == challenge_contigs[name]
            for name in challenge_names
        )
        return ReferenceIdentification(
            assembly=contract.assembly,
            reference_role="challenge_alignment",
            reference_identity_id=contract.reference_role(
                "challenge_alignment"
            ).identity_id,
            reference_basename=basename,
            contig_names_match=names_match,
            contig_lengths_match=lengths_match,
            contig_order_match=lengths_match and declared_names == challenge_names,
        )

    if basename in PUBLIC_GIAB_REFERENCE_BASENAMES:
        expected_items = tuple(challenge_contigs.items())[
            : contract.transformation.masked_component_records
        ]
        expected_contigs = {f"chr{name}": length for name, length in expected_items}
        declared_names = tuple(declared_contigs)
        expected_names = tuple(expected_contigs)
        names_match = set(declared_names) == set(expected_names)
        lengths_match = names_match and all(
            declared_contigs[name] == expected_contigs[name]
            for name in expected_names
        )
        role = contract.reference_role("public_annotation_benchmark")
        return ReferenceIdentification(
            assembly="GRCh38",
            reference_role=role.role,
            reference_identity_id=role.identity_id,
            reference_basename=basename,
            contig_names_match=names_match,
            contig_lengths_match=lengths_match,
            contig_order_match=lengths_match and declared_names == expected_names,
        )

    raise ReferenceConcordanceError(
        "assembly-mismatch: VCF reference declaration is not an accepted GRCh38 role"
    )


@dataclass(frozen=True, slots=True)
class ContigOrderAdjudication:
    resolution: str
    source_report_sha256: str
    source_input_sha256: str
    revalidated_input_sha256: str
    source_reference_manifest_sha256: str
    revalidated_reference_manifest_sha256: str
    record_payload_sha256: str
    accepted: bool
    schema: str = CONTIG_ORDER_ADJUDICATION_SCHEMA

    def __post_init__(self) -> None:
        if self.schema != CONTIG_ORDER_ADJUDICATION_SCHEMA:
            raise ReferenceConcordanceError(
                "contig-order adjudication schema is invalid"
            )
        if self.resolution not in ORDER_ADJUDICATION_RESOLUTIONS:
            raise ReferenceConcordanceError(
                "contig-order adjudication resolution is invalid"
            )
        for label, value in (
            ("source_report_sha256", self.source_report_sha256),
            ("source_input_sha256", self.source_input_sha256),
            ("revalidated_input_sha256", self.revalidated_input_sha256),
            (
                "source_reference_manifest_sha256",
                self.source_reference_manifest_sha256,
            ),
            (
                "revalidated_reference_manifest_sha256",
                self.revalidated_reference_manifest_sha256,
            ),
            ("record_payload_sha256", self.record_payload_sha256),
        ):
            _require_digest(value, SHA256, f"adjudication {label}")
        if not isinstance(self.accepted, bool) or not self.accepted:
            raise ReferenceConcordanceError(
                "contig-order adjudication must be accepted"
            )
        if self.resolution == "vcf-header-reordered-without-record-change":
            if (
                self.source_input_sha256 == self.revalidated_input_sha256
                or self.source_reference_manifest_sha256
                != self.revalidated_reference_manifest_sha256
            ):
                raise ReferenceConcordanceError(
                    "VCF-header adjudication bindings are inconsistent"
                )
        elif (
            self.source_input_sha256 != self.revalidated_input_sha256
            or self.source_reference_manifest_sha256
            == self.revalidated_reference_manifest_sha256
        ):
            raise ReferenceConcordanceError(
                "reference-order adjudication bindings are inconsistent"
            )

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "ContigOrderAdjudication":
        expected = {
            "schema",
            "resolution",
            "source_report_sha256",
            "source_input_sha256",
            "revalidated_input_sha256",
            "source_reference_manifest_sha256",
            "revalidated_reference_manifest_sha256",
            "record_payload_sha256",
            "accepted",
        }
        _reject_unknown(data, expected, "contig-order adjudication")
        if set(data) != expected:
            raise ReferenceConcordanceError(
                "contig-order adjudication is missing required fields"
            )
        try:
            return cls(**data)
        except TypeError as error:
            raise ReferenceConcordanceError(
                f"invalid contig-order adjudication: {error}"
            ) from error


@dataclass(frozen=True, slots=True)
class ReferenceConcordanceReport:
    reference_identity_id: str
    reference_role: str
    reference_manifest_sha256: str
    input_sha256: str
    validated_input_sha256: str
    reference_basename: str
    vcf_contig_dictionary_sha256: str
    reference_contig_dictionary_sha256: str
    record_payload_sha256: str
    header_valid: bool
    reference_basename_match: bool
    contig_names_match: bool
    contig_lengths_match: bool
    contig_order_match: bool
    order_revalidated: bool
    order_adjudication_sha256: str | None
    declared_contigs: int
    reference_contigs: int
    total_records: int
    checked_records: int
    matching_records: int
    mismatching_records: int
    unknown_contig_records: int
    out_of_bounds_records: int
    malformed_ref_records: int
    otherwise_unresolvable_records: int
    unresolved_records: int
    accepted: bool
    schema: str = CONCORDANCE_REPORT_SCHEMA
    read_limits: VcfReadLimits = VcfReadLimits()
    reference_declarations: tuple[ReferenceDeclaration, ...] = ()

    def __post_init__(self) -> None:
        if (type(self.reference_declarations) is not tuple
                or any(type(item) is not ReferenceDeclaration for item in self.reference_declarations)):
            raise ReferenceConcordanceError("invalid reference metadata evidence")
        if type(self.read_limits) is not VcfReadLimits:
            raise ReferenceConcordanceError("invalid concordance read limits")
        if self.schema != CONCORDANCE_REPORT_SCHEMA:
            raise ReferenceConcordanceError(
                f"concordance report schema must be {CONCORDANCE_REPORT_SCHEMA!r}"
            )
        if not IDENTIFIER.fullmatch(self.reference_identity_id):
            raise ReferenceConcordanceError("report reference_identity_id is invalid")
        if self.reference_role != "challenge_alignment":
            raise ReferenceConcordanceError(
                "complete concordance report must use challenge_alignment"
            )
        if not self.reference_basename:
            raise ReferenceConcordanceError("report reference_basename is invalid")
        for label, digest in (
            ("reference manifest SHA-256", self.reference_manifest_sha256),
            ("input SHA-256", self.input_sha256),
            ("validated input SHA-256", self.validated_input_sha256),
            ("VCF contig dictionary SHA-256", self.vcf_contig_dictionary_sha256),
            (
                "reference contig dictionary SHA-256",
                self.reference_contig_dictionary_sha256,
            ),
            ("record payload SHA-256", self.record_payload_sha256),
        ):
            _require_digest(digest, SHA256, label)
        if self.order_adjudication_sha256 is not None:
            _require_digest(
                self.order_adjudication_sha256,
                SHA256,
                "order adjudication SHA-256",
            )
        for label, value in (
            ("declared_contigs", self.declared_contigs),
            ("reference_contigs", self.reference_contigs),
            ("total_records", self.total_records),
            ("checked_records", self.checked_records),
            ("matching_records", self.matching_records),
            ("mismatching_records", self.mismatching_records),
            ("unknown_contig_records", self.unknown_contig_records),
            ("out_of_bounds_records", self.out_of_bounds_records),
            ("malformed_ref_records", self.malformed_ref_records),
            (
                "otherwise_unresolvable_records",
                self.otherwise_unresolvable_records,
            ),
            ("unresolved_records", self.unresolved_records),
        ):
            if not isinstance(value, int) or isinstance(value, bool) or value < 0:
                raise ReferenceConcordanceError(
                    f"report {label} must be a non-negative integer"
                )
        boolean_fields = (
            self.header_valid,
            self.reference_basename_match,
            self.contig_names_match,
            self.contig_lengths_match,
            self.contig_order_match,
            self.order_revalidated,
            self.accepted,
        )
        if any(not isinstance(value, bool) for value in boolean_fields):
            raise ReferenceConcordanceError("report status fields must be booleans")
        expected_unresolved = (
            self.unknown_contig_records
            + self.out_of_bounds_records
            + self.malformed_ref_records
            + self.otherwise_unresolvable_records
        )
        if self.unresolved_records != expected_unresolved:
            raise ReferenceConcordanceError(
                "report unresolved_records conflicts with reason counts"
            )
        if self.checked_records != self.matching_records + self.mismatching_records:
            raise ReferenceConcordanceError(
                "report checked_records conflicts with match counts"
            )
        if self.total_records != self.checked_records + self.unresolved_records:
            raise ReferenceConcordanceError(
                "report total_records is not completely accounted for"
            )
        if self.contig_lengths_match and not self.contig_names_match:
            raise ReferenceConcordanceError(
                "report length match requires a complete name match"
            )
        if self.contig_order_match and not self.contig_lengths_match:
            raise ReferenceConcordanceError(
                "report order match requires complete names and lengths"
            )
        if self.contig_names_match and self.declared_contigs != self.reference_contigs:
            raise ReferenceConcordanceError(
                "report complete name match conflicts with contig counts"
            )
        if (
            self.contig_order_match
            and self.vcf_contig_dictionary_sha256
            != self.reference_contig_dictionary_sha256
        ):
            raise ReferenceConcordanceError(
                "report ordered dictionary match conflicts with dictionary digests"
            )
        if self.order_revalidated != (self.order_adjudication_sha256 is not None):
            raise ReferenceConcordanceError(
                "report order revalidation and adjudication binding conflict"
            )
        if not self.order_revalidated and self.input_sha256 != self.validated_input_sha256:
            raise ReferenceConcordanceError(
                "report changed input requires order adjudication and revalidation"
            )
        expected_acceptance = (
            self.header_valid
            and self.reference_basename_match
            and self.contig_names_match
            and self.contig_lengths_match
            and self.contig_order_match
            and self.total_records > 0
            and self.matching_records == self.total_records
            and self.mismatching_records == 0
            and self.unresolved_records == 0
        )
        if self.accepted != expected_acceptance:
            raise ReferenceConcordanceError(
                "report accepted status conflicts with complete validation"
            )

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "ReferenceConcordanceReport":
        expected = {
            "read_limits", "reference_declarations",
            "schema",
            "reference_identity_id",
            "reference_role",
            "reference_manifest_sha256",
            "input_sha256",
            "validated_input_sha256",
            "reference_basename",
            "vcf_contig_dictionary_sha256",
            "reference_contig_dictionary_sha256",
            "record_payload_sha256",
            "header_valid",
            "reference_basename_match",
            "contig_names_match",
            "contig_lengths_match",
            "contig_order_match",
            "order_revalidated",
            "order_adjudication_sha256",
            "declared_contigs",
            "reference_contigs",
            "total_records",
            "checked_records",
            "matching_records",
            "mismatching_records",
            "unknown_contig_records",
            "out_of_bounds_records",
            "malformed_ref_records",
            "otherwise_unresolvable_records",
            "unresolved_records",
            "accepted",
        }
        _reject_unknown(data, expected, "concordance report")
        if set(data) != expected:
            raise ReferenceConcordanceError(
                "concordance report is missing required fields"
            )
        try:
            values = dict(data)
            if type(data["read_limits"]) is not dict or set(data["read_limits"]) != set(VcfReadLimits.__dataclass_fields__):
                raise ValueError("incomplete concordance read policy")
            values["read_limits"] = VcfReadLimits(**data["read_limits"])
            values["reference_declarations"] = tuple(ReferenceDeclaration(**item) for item in data["reference_declarations"])
            return cls(**values)
        except (TypeError, ValueError) as error:
            raise ReferenceConcordanceError(
                f"invalid concordance report: {error}"
            ) from error


@dataclass(frozen=True, slots=True)
class PublicValidationConcordanceReport:
    """Record-level concordance for a bounded public VCF dictionary.

    Public fixtures declare only the contigs they use. This distinct type keeps
    that bounded dictionary from being accepted where the complete challenge
    dictionary is mandatory.
    """

    reference_identity_id: str
    reference_role: str
    reference_manifest_sha256: str
    input_sha256: str
    validated_input_sha256: str
    reference_basename: str
    alias_map_sha256: str
    vcf_contig_dictionary_sha256: str
    reference_contig_dictionary_sha256: str
    record_payload_sha256: str
    header_valid: bool
    reference_basename_match: bool
    contig_names_match: bool
    contig_lengths_match: bool
    contig_order_match: bool
    declared_contigs: int
    available_reference_contigs: int
    total_records: int
    checked_records: int
    matching_records: int
    mismatching_records: int
    unknown_contig_records: int
    out_of_bounds_records: int
    malformed_ref_records: int
    otherwise_unresolvable_records: int
    unresolved_records: int
    accepted: bool
    schema: str = PUBLIC_VALIDATION_CONCORDANCE_REPORT_SCHEMA
    read_limits: VcfReadLimits = VcfReadLimits()
    reference_declarations: tuple[ReferenceDeclaration, ...] = ()

    def __post_init__(self) -> None:
        if (type(self.reference_declarations) is not tuple
                or any(type(item) is not ReferenceDeclaration for item in self.reference_declarations)):
            raise ReferenceConcordanceError("invalid reference metadata evidence")
        if type(self.read_limits) is not VcfReadLimits:
            raise ReferenceConcordanceError("invalid concordance read limits")
        if (
            self.schema != PUBLIC_VALIDATION_CONCORDANCE_REPORT_SCHEMA
            or not IDENTIFIER.fullmatch(self.reference_identity_id)
            or self.reference_role != "public_annotation_benchmark"
            or not self.reference_basename
        ):
            raise ReferenceConcordanceError(
                "public validation concordance identity is invalid"
            )
        for value in (
            self.reference_manifest_sha256,
            self.input_sha256,
            self.validated_input_sha256,
            self.alias_map_sha256,
            self.vcf_contig_dictionary_sha256,
            self.reference_contig_dictionary_sha256,
            self.record_payload_sha256,
        ):
            _require_digest(value, SHA256, "public validation concordance digest")
        counts = (
            self.declared_contigs,
            self.available_reference_contigs,
            self.total_records,
            self.checked_records,
            self.matching_records,
            self.mismatching_records,
            self.unknown_contig_records,
            self.out_of_bounds_records,
            self.malformed_ref_records,
            self.otherwise_unresolvable_records,
            self.unresolved_records,
        )
        statuses = (
            self.header_valid,
            self.reference_basename_match,
            self.contig_names_match,
            self.contig_lengths_match,
            self.contig_order_match,
            self.accepted,
        )
        unresolved = (
            self.unknown_contig_records
            + self.out_of_bounds_records
            + self.malformed_ref_records
            + self.otherwise_unresolvable_records
        )
        expected_acceptance = (
            self.header_valid
            and self.reference_basename_match
            and self.contig_names_match
            and self.contig_lengths_match
            and self.contig_order_match
            and 0 < self.declared_contigs <= self.available_reference_contigs
            and self.total_records > 0
            and self.matching_records == self.total_records
            and self.mismatching_records == 0
            and self.unresolved_records == 0
        )
        if (
            self.input_sha256 != self.validated_input_sha256
            or any(type(value) is not int or value < 0 for value in counts)
            or any(type(value) is not bool for value in statuses)
            or self.unresolved_records != unresolved
            or self.checked_records != self.matching_records + self.mismatching_records
            or self.total_records != self.checked_records + self.unresolved_records
            or self.contig_order_match and not self.contig_lengths_match
            or self.contig_lengths_match and not self.contig_names_match
            or self.accepted != expected_acceptance
        ):
            raise ReferenceConcordanceError(
                "public validation concordance accounting is invalid"
            )

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "PublicValidationConcordanceReport":
        expected = set(cls.__dataclass_fields__)
        _reject_unknown(data, expected, "public validation concordance report")
        if set(data) != expected:
            raise ReferenceConcordanceError(
                "public validation concordance report is missing required fields"
            )
        try:
            values = dict(data)
            if type(data["read_limits"]) is not dict or set(data["read_limits"]) != set(VcfReadLimits.__dataclass_fields__):
                raise ValueError("incomplete concordance read policy")
            values["read_limits"] = VcfReadLimits(**data["read_limits"])
            values["reference_declarations"] = tuple(ReferenceDeclaration(**item) for item in data["reference_declarations"])
            return cls(**values)
        except (TypeError, ValueError) as error:
            raise ReferenceConcordanceError(
                "invalid public validation concordance report"
            ) from error


ReferenceFetcher = Callable[[str, int, int], str]


def _is_order_only_report(report: ReferenceConcordanceReport) -> bool:
    return (
        not report.accepted
        and report.header_valid
        and report.reference_basename_match
        and report.contig_names_match
        and report.contig_lengths_match
        and not report.contig_order_match
        and not report.order_revalidated
        and report.order_adjudication_sha256 is None
        and report.input_sha256 == report.validated_input_sha256
        and report.total_records > 0
        and report.matching_records == report.total_records
        and report.mismatching_records == 0
        and report.unresolved_records == 0
    )


def _validate_adjudication_binding(
    adjudication: ContigOrderAdjudication,
    *,
    source_report: ReferenceConcordanceReport,
    source_input_sha256: str,
    input_sha256: str,
    reference_manifest_sha256: str,
    record_payload_sha256: str,
    source_header_sha256: str,
    revalidated_header_sha256: str,
) -> str:
    if (
        adjudication.source_report_sha256
        != _canonical_json_sha256(asdict(source_report))
        or adjudication.source_input_sha256 != source_input_sha256
        or adjudication.revalidated_input_sha256 != input_sha256
        or adjudication.source_reference_manifest_sha256
        != source_report.reference_manifest_sha256
        or adjudication.revalidated_reference_manifest_sha256
        != reference_manifest_sha256
        or adjudication.record_payload_sha256 != record_payload_sha256
        or source_report.input_sha256 != source_input_sha256
        or source_report.validated_input_sha256 != source_input_sha256
        or source_report.record_payload_sha256 != record_payload_sha256
    ):
        raise ReferenceConcordanceError(
            "contig-order adjudication is stale for the revalidation inputs"
        )
    if not _is_order_only_report(source_report):
        raise ReferenceConcordanceError(
            "contig-order adjudication source report is not order-only"
        )
    if (
        adjudication.resolution == "vcf-header-reordered-without-record-change"
        and source_header_sha256 != revalidated_header_sha256
    ):
        raise ReferenceConcordanceError(
            "VCF-header adjudication changed content other than contig order"
        )
    return _canonical_json_sha256(asdict(adjudication))


def evaluate_vcf_reference_concordance(
    vcf_path: str | Path,
    fetch_reference: ReferenceFetcher,
    reference_contigs: Mapping[str, int],
    *,
    contract: ReferenceIdentityContract,
    reference_manifest_sha256: str,
    source_vcf_path: str | Path | None = None,
    source_report: ReferenceConcordanceReport | None = None,
    source_fetch_reference: ReferenceFetcher | None = None,
    source_reference_contigs: Mapping[str, int] | None = None,
    order_adjudication: ContigOrderAdjudication | None = None,
    read_limits: VcfReadLimits = VcfReadLimits(),
) -> ReferenceConcordanceReport:
    """Stream every VCF body record and account for each REF comparison outcome."""
    _require_digest(
        reference_manifest_sha256, SHA256, "reference manifest SHA-256"
    )
    if not reference_contigs:
        raise ReferenceConcordanceError("reference contig dictionary is empty")
    scan = _scan_vcf(Path(vcf_path), reference_contigs, fetch_reference,
        challenge_reference_basename=contract.challenge_reference_basename, read_limits=read_limits)
    if source_vcf_path is None:
        source_input_sha256 = scan.input_sha256
        source_payload_sha256 = scan.record_payload_sha256
        source_header_sha256 = scan.order_insensitive_header_sha256
    else:
        (
            source_input_sha256,
            source_payload_sha256,
            source_header_sha256,
        ) = _vcf_bindings(Path(source_vcf_path), read_limits)
        if source_payload_sha256 != scan.record_payload_sha256:
            raise ReferenceConcordanceError(
                "revalidated VCF record payload differs from the source VCF"
            )
    if scan.reference_declaration is None:
        raise ReferenceConcordanceError(
            "assembly-mismatch: VCF must declare one reference or include the pinned challenge reference"
        )
    identification = identify_vcf_reference(
        scan.reference_declaration,
        scan.contigs,
        reference_contigs,
        contract,
    )
    if identification.reference_role != "challenge_alignment":
        raise ReferenceConcordanceError(
            "controlled complete concordance requires the challenge alignment role"
        )
    adjudication_sha256 = None
    if source_input_sha256 != scan.input_sha256 and order_adjudication is None:
        raise ReferenceConcordanceError(
            "changed validation input requires contig-order adjudication"
        )
    if order_adjudication is not None:
        if source_vcf_path is None or source_report is None:
            raise ReferenceConcordanceError(
                "contig-order adjudication requires its source VCF and report"
            )
        if not identification.contig_order_match:
            raise ReferenceConcordanceError(
                "contig order remains unresolved after adjudication"
            )
        if source_fetch_reference is None and source_reference_contigs is None:
            if (
                order_adjudication.source_reference_manifest_sha256
                != reference_manifest_sha256
            ):
                raise ReferenceConcordanceError(
                    "reference-order adjudication requires the verified source reference"
                )
            source_fetch_reference = fetch_reference
            source_reference_contigs = reference_contigs
        elif source_fetch_reference is None or source_reference_contigs is None:
            raise ReferenceConcordanceError(
                "source reference fetcher and contigs must be provided together"
            )
        replayed_source_report = evaluate_vcf_reference_concordance(
            source_vcf_path,
            source_fetch_reference,
            source_reference_contigs,
            contract=contract,
            read_limits=read_limits,
            reference_manifest_sha256=(
                order_adjudication.source_reference_manifest_sha256
            ),
        )
        if replayed_source_report != source_report:
            raise ReferenceConcordanceError(
                "contig-order source report differs from complete validation replay"
            )
        adjudication_sha256 = _validate_adjudication_binding(
            order_adjudication,
            source_report=source_report,
            source_input_sha256=source_input_sha256,
            input_sha256=scan.input_sha256,
            reference_manifest_sha256=reference_manifest_sha256,
            record_payload_sha256=scan.record_payload_sha256,
            source_header_sha256=source_header_sha256,
            revalidated_header_sha256=scan.order_insensitive_header_sha256,
        )
    unresolved_records = (
        scan.unknown_contig_records
        + scan.out_of_bounds_records
        + scan.malformed_ref_records
        + scan.otherwise_unresolvable_records
    )
    accepted = (
        scan.header_valid
        and identification.reference_basename
        == contract.challenge_reference_basename
        and identification.contig_names_match
        and identification.contig_lengths_match
        and identification.contig_order_match
        and scan.total_records > 0
        and scan.matching_records == scan.total_records
        and scan.mismatching_records == 0
        and unresolved_records == 0
    )
    return ReferenceConcordanceReport(
        reference_declarations=_reference_metadata(scan, identification.reference_basename),
        read_limits=read_limits,
        reference_identity_id=identification.reference_identity_id,
        reference_role=identification.reference_role,
        reference_manifest_sha256=reference_manifest_sha256,
        input_sha256=source_input_sha256,
        validated_input_sha256=scan.input_sha256,
        reference_basename=identification.reference_basename,
        vcf_contig_dictionary_sha256=_mapping_sha256(scan.contigs),
        reference_contig_dictionary_sha256=_mapping_sha256(reference_contigs),
        record_payload_sha256=scan.record_payload_sha256,
        header_valid=scan.header_valid,
        reference_basename_match=(
            identification.reference_basename
            == contract.challenge_reference_basename
        ),
        contig_names_match=identification.contig_names_match,
        contig_lengths_match=identification.contig_lengths_match,
        contig_order_match=identification.contig_order_match,
        order_revalidated=adjudication_sha256 is not None,
        order_adjudication_sha256=adjudication_sha256,
        declared_contigs=len(scan.contigs),
        reference_contigs=len(reference_contigs),
        total_records=scan.total_records,
        checked_records=scan.checked_records,
        matching_records=scan.matching_records,
        mismatching_records=scan.mismatching_records,
        unknown_contig_records=scan.unknown_contig_records,
        out_of_bounds_records=scan.out_of_bounds_records,
        malformed_ref_records=scan.malformed_ref_records,
        otherwise_unresolvable_records=scan.otherwise_unresolvable_records,
        unresolved_records=unresolved_records,
        accepted=accepted,
    )


def evaluate_public_validation_vcf_concordance(
    vcf_path: str | Path,
    fetch_reference: ReferenceFetcher,
    reference_contigs: Mapping[str, int],
    source_to_reference_contigs: Mapping[str, str],
    *,
    contract: ReferenceIdentityContract,
    reference_manifest_sha256: str,
    read_limits: VcfReadLimits = VcfReadLimits(),
) -> PublicValidationConcordanceReport:
    """Validate every record in a bounded public fixture against public GRCh38.

    The fixture may declare a strict subset of the public reference dictionary,
    but every declared source name and length must map to the bound reference and
    preserve its order. This path cannot produce a challenge-alignment report.
    """
    _require_digest(
        reference_manifest_sha256, SHA256, "reference manifest SHA-256"
    )
    role = contract.reference_role("public_annotation_benchmark")
    if (
        not reference_contigs
        or not source_to_reference_contigs
        or any(
            type(source) is not str
            or not source
            or type(target) is not str
            or target not in reference_contigs
            for source, target in source_to_reference_contigs.items()
        )
    ):
        raise ReferenceConcordanceError(
            "public validation reference aliases are invalid"
        )
    compatible_contigs = {
        source: reference_contigs[target]
        for source, target in source_to_reference_contigs.items()
    }

    def mapped_fetch(contig: str, position: int, length: int) -> str:
        return fetch_reference(source_to_reference_contigs[contig], position, length)

    scan = _scan_vcf(Path(vcf_path), compatible_contigs, mapped_fetch, read_limits=read_limits)
    if scan.reference_declaration is None:
        raise ReferenceConcordanceError(
            "assembly-mismatch: public validation VCF must declare one reference"
        )
    basename = _reference_basename(scan.reference_declaration)
    names_match = bool(scan.contigs) and set(scan.contigs).issubset(compatible_contigs)
    lengths_match = names_match and all(
        scan.contigs[name] == compatible_contigs[name] for name in scan.contigs
    )
    order = {name: index for index, name in enumerate(reference_contigs)}
    projected = tuple(
        order[source_to_reference_contigs[name]]
        for name in scan.contigs
        if name in source_to_reference_contigs
    )
    order_match = (
        lengths_match
        and len(projected) == len(scan.contigs)
        and projected == tuple(sorted(projected))
    )
    unresolved_records = (
        scan.unknown_contig_records
        + scan.out_of_bounds_records
        + scan.malformed_ref_records
        + scan.otherwise_unresolvable_records
    )
    basename_match = basename in PUBLIC_GIAB_REFERENCE_BASENAMES
    accepted = (
        scan.header_valid
        and basename_match
        and names_match
        and lengths_match
        and order_match
        and scan.total_records > 0
        and scan.matching_records == scan.total_records
        and scan.mismatching_records == 0
        and unresolved_records == 0
    )
    return PublicValidationConcordanceReport(
        reference_declarations=_reference_metadata(scan, basename if basename_match else None),
        read_limits=read_limits,
        reference_identity_id=role.identity_id,
        reference_role=role.role,
        reference_manifest_sha256=reference_manifest_sha256,
        input_sha256=scan.input_sha256,
        validated_input_sha256=scan.input_sha256,
        reference_basename=basename,
        alias_map_sha256=_canonical_json_sha256(dict(source_to_reference_contigs)),
        vcf_contig_dictionary_sha256=_mapping_sha256(scan.contigs),
        reference_contig_dictionary_sha256=_mapping_sha256(reference_contigs),
        record_payload_sha256=scan.record_payload_sha256,
        header_valid=scan.header_valid,
        reference_basename_match=basename_match,
        contig_names_match=names_match,
        contig_lengths_match=lengths_match,
        contig_order_match=order_match,
        declared_contigs=len(scan.contigs),
        available_reference_contigs=len(reference_contigs),
        total_records=scan.total_records,
        checked_records=scan.checked_records,
        matching_records=scan.matching_records,
        mismatching_records=scan.mismatching_records,
        unknown_contig_records=scan.unknown_contig_records,
        out_of_bounds_records=scan.out_of_bounds_records,
        malformed_ref_records=scan.malformed_ref_records,
        otherwise_unresolvable_records=scan.otherwise_unresolvable_records,
        unresolved_records=unresolved_records,
        accepted=accepted,
    )


def _vcf_bindings(path: Path, read_limits: VcfReadLimits = VcfReadLimits()) -> tuple[str, str, str]:
    record_digest = hashlib.sha256()
    header_digest = hashlib.sha256()
    contig_header_lines: list[bytes] = []
    column_header_seen = False
    try:
        with _open_vcf(path) as (source, raw_source):
            for raw_line in bounded_vcf_lines(source, read_limits, error_type=ReferenceConcordanceError):
                line = raw_line.rstrip(b"\r\n")
                if not column_header_seen and line.startswith(b"#CHROM"):
                    _update_framed_digest(header_digest, b"H", raw_line)
                    column_header_seen = True
                    continue
                if not column_header_seen and line.startswith(b"#"):
                    if line.startswith(b"##contig="):
                        contig_header_lines.append(raw_line)
                    else:
                        _update_framed_digest(header_digest, b"H", raw_line)
                    continue
                column_header_seen = True
                record_digest.update(raw_line)
    except (EOFError, gzip.BadGzipFile, OSError, zlib.error) as error:
        raise ReferenceConcordanceError("cannot parse VCF") from error
    for raw_line in sorted(contig_header_lines):
        _update_framed_digest(header_digest, b"C", raw_line)
    return (
        raw_source.hexdigest(),
        record_digest.hexdigest(),
        header_digest.hexdigest(),
    )


def document_contig_order_adjudication(
    source_report: ReferenceConcordanceReport,
    *,
    source_input_path: str | Path,
    revalidated_input_path: str | Path,
    source_fetch_reference: ReferenceFetcher,
    source_reference_contigs: Mapping[str, int],
    contract: ReferenceIdentityContract,
    revalidated_reference_manifest_sha256: str,
    resolution: str,
) -> ContigOrderAdjudication:
    """Bind an order-only finding to a non-sequence-changing correction."""
    replayed_source_report = evaluate_vcf_reference_concordance(
        source_input_path,
        source_fetch_reference,
        source_reference_contigs,
        contract=contract,
        reference_manifest_sha256=source_report.reference_manifest_sha256,
    )
    if replayed_source_report != source_report:
        raise ReferenceConcordanceError(
            "contig-order source report differs from complete validation replay"
        )
    if not _is_order_only_report(source_report):
        raise ReferenceConcordanceError(
            "only an otherwise-valid order discrepancy can be adjudicated"
        )
    _require_digest(
        revalidated_reference_manifest_sha256,
        SHA256,
        "revalidated reference manifest SHA-256",
    )
    source_path = Path(source_input_path)
    revalidated_path = Path(revalidated_input_path)
    (
        source_input_sha256,
        source_payload_sha256,
        source_header_sha256,
    ) = _vcf_bindings(source_path)
    (
        revalidated_input_sha256,
        revalidated_payload_sha256,
        revalidated_header_sha256,
    ) = _vcf_bindings(revalidated_path)
    if source_input_sha256 != source_report.input_sha256:
        raise ReferenceConcordanceError("order adjudication source VCF is stale")
    if (
        source_payload_sha256 != source_report.record_payload_sha256
        or revalidated_payload_sha256 != source_payload_sha256
    ):
        raise ReferenceConcordanceError(
            "order adjudication changed or cannot bind the VCF record payload"
        )
    if (
        resolution == "vcf-header-reordered-without-record-change"
        and source_header_sha256 != revalidated_header_sha256
    ):
        raise ReferenceConcordanceError(
            "VCF-header adjudication changed content other than contig order"
        )
    return ContigOrderAdjudication(
        resolution=resolution,
        source_report_sha256=_canonical_json_sha256(asdict(source_report)),
        source_input_sha256=source_input_sha256,
        revalidated_input_sha256=revalidated_input_sha256,
        source_reference_manifest_sha256=source_report.reference_manifest_sha256,
        revalidated_reference_manifest_sha256=(
            revalidated_reference_manifest_sha256
        ),
        record_payload_sha256=source_payload_sha256,
        accepted=True,
    )


def concordance_report_as_json(
    report: ReferenceConcordanceReport | PublicValidationConcordanceReport,
) -> str:
    return json.dumps(asdict(report), indent=2, sort_keys=True) + "\n"


def contig_order_adjudication_as_json(
    adjudication: ContigOrderAdjudication,
) -> str:
    return json.dumps(asdict(adjudication), indent=2, sort_keys=True) + "\n"


def _parse_json_object(content: bytes, label: str) -> Mapping[str, Any]:
    try:
        data = json.loads(content.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ReferenceConcordanceError(f"invalid {label} JSON") from error
    if not isinstance(data, Mapping):
        raise ReferenceConcordanceError(f"{label} root must be an object")
    return data


@contextmanager
def _open_concordance_report(path: str | Path) -> Iterator[ReferenceConcordanceReport]:
    with _open_stable_binary(Path(path), "concordance report") as source:
        content = source.read()
        report = ReferenceConcordanceReport.from_mapping(
            _parse_json_object(content, "concordance report")
        )
        if hashlib.sha256(content).hexdigest() != _canonical_json_sha256(
            asdict(report)
        ):
            raise ReferenceConcordanceError("concordance report is not canonical")
        yield report


@contextmanager
def _open_contig_order_adjudication(
    path: str | Path,
) -> Iterator[ContigOrderAdjudication]:
    with _open_stable_binary(Path(path), "contig-order adjudication") as source:
        content = source.read()
        adjudication = ContigOrderAdjudication.from_mapping(
            _parse_json_object(content, "contig-order adjudication")
        )
        if hashlib.sha256(content).hexdigest() != _canonical_json_sha256(
            asdict(adjudication)
        ):
            raise ReferenceConcordanceError(
                "contig-order adjudication is not canonical"
            )
        yield adjudication


def load_concordance_report(path: str | Path) -> ReferenceConcordanceReport:
    with _open_concordance_report(path) as report:
        return report


def load_public_validation_concordance_report(
    path: str | Path,
) -> PublicValidationConcordanceReport:
    """Load only the bounded public-reference report type."""

    with _open_stable_binary(Path(path), "public validation concordance report") as source:
        content = source.read()
        report = PublicValidationConcordanceReport.from_mapping(
            _parse_json_object(content, "public validation concordance report")
        )
        if hashlib.sha256(content).hexdigest() != _canonical_json_sha256(
            asdict(report)
        ):
            raise ReferenceConcordanceError(
                "public validation concordance report is not canonical"
            )
        return report


def load_contig_order_adjudication(path: str | Path) -> ContigOrderAdjudication:
    with _open_contig_order_adjudication(path) as adjudication:
        return adjudication


def _write_private_json(
    path: Path,
    value: str,
    label: str,
    *,
    allowed_root: str | Path | None = None,
) -> str:
    path.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    temporary_name = f".{path.name}.{uuid.uuid4().hex}.partial"
    parent_fd = -1
    try:
        parent_fd = os.open(
            path.parent,
            os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
        )
        if allowed_root is not None:
            root = Path(allowed_root).resolve(strict=True)
            opened_parent = Path(f"/proc/self/fd/{parent_fd}").resolve(strict=True)
            if not opened_parent.is_relative_to(root):
                raise ReferenceConcordanceError(
                    f"{label} destination escapes controlled root"
                )
        descriptor = os.open(
            temporary_name,
            os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
            0o600,
            dir_fd=parent_fd,
        )
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as destination:
            destination.write(value)
            destination.flush()
            os.fsync(destination.fileno())
        os.replace(
            temporary_name,
            path.name,
            src_dir_fd=parent_fd,
            dst_dir_fd=parent_fd,
        )
        os.chmod(path.name, 0o600, dir_fd=parent_fd, follow_symlinks=False)
        os.fsync(parent_fd)
    except OSError as error:
        raise ReferenceConcordanceError(f"cannot write {label}") from error
    finally:
        if parent_fd >= 0:
            try:
                os.unlink(temporary_name, dir_fd=parent_fd)
            except FileNotFoundError:
                pass
            finally:
                os.close(parent_fd)
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def write_concordance_report(
    path: str | Path,
    report: ReferenceConcordanceReport,
    *,
    allowed_root: str | Path | None = None,
) -> str:
    return _write_private_json(
        Path(path),
        concordance_report_as_json(report),
        "concordance report",
        allowed_root=allowed_root,
    )


def write_contig_order_adjudication(
    path: str | Path,
    adjudication: ContigOrderAdjudication,
    *,
    allowed_root: str | Path | None = None,
) -> str:
    return _write_private_json(
        Path(path),
        contig_order_adjudication_as_json(adjudication),
        "contig-order adjudication",
        allowed_root=allowed_root,
    )


def verify_reference_manifest_contract(
    contract: ReferenceIdentityContract,
    manifest: PublicResourceManifest,
    *,
    contract_sha256: str | None = None,
) -> None:
    if manifest.resource_id != contract.resource_id:
        raise ReferenceConcordanceError("reference manifest resource_id mismatch")
    if manifest.release != contract.release:
        raise ReferenceConcordanceError("reference manifest release mismatch")
    if manifest.assembly != contract.assembly:
        raise ReferenceConcordanceError("reference manifest assembly mismatch")
    source = contract.canonical_source
    if manifest.source.url != source.url or manifest.source.sha256 != source.sha256:
        raise ReferenceConcordanceError("reference manifest source identity mismatch")
    if manifest.transformation.schema != "genosapiens.grch38-bootstrap":
        raise ReferenceConcordanceError("reference manifest transformation is obsolete")
    if contract_sha256 is not None:
        _require_digest(contract_sha256, SHA256, "reference contract SHA-256")
        if manifest.transformation.parameters_sha256 != contract_sha256:
            raise ReferenceConcordanceError("reference manifest contract is stale")
    expected_paths = {
        *(artifact.path for artifact in contract.source_artifacts),
        *(output.path for output in contract.expected_outputs),
    }
    manifest_paths = {output.path for output in manifest.outputs}
    if manifest_paths != expected_paths:
        raise ReferenceConcordanceError(
            "reference manifest does not bind the complete source and output set"
        )


def verify_reference_construction(
    contract: ReferenceIdentityContract,
    manifest: PublicResourceManifest,
    reference_root: str | Path,
    *,
    contract_sha256: str,
) -> ReferenceConstructionReport:
    """Require an accepted report whose source and output hashes match the manifest."""
    verify_reference_manifest_contract(
        contract, manifest, contract_sha256=contract_sha256
    )
    root = Path(reference_root).resolve()
    report_path = root / contract.expected_output("construction_report").path
    report = load_construction_report(report_path)
    expected_role_ids = {
        role.role: role.identity_id for role in contract.reference_roles
    }
    if (
        report.contract_id != contract.contract_id
        or report.contract_sha256 != contract_sha256
        or report.challenge_alignment_identity_id
        != expected_role_ids["challenge_alignment"]
        or report.challenge_tool_view_identity_id
        != expected_role_ids["challenge_tool_view"]
        or report.public_annotation_benchmark_identity_id
        != expected_role_ids["public_annotation_benchmark"]
        or report.masked_component_records
        != contract.transformation.masked_component_records
        or report.decoy_suffix_records != contract.transformation.decoy_suffix_records
        or report.total_records != contract.sequence_count
        or not report.accepted
    ):
        raise ReferenceConcordanceError(
            "reference construction report does not satisfy the pinned contract"
        )
    expected_sources = {
        artifact.role: artifact.sha256 for artifact in contract.source_artifacts
    }
    if dict(report.source_sha256) != expected_sources:
        raise ReferenceConcordanceError(
            "reference construction report source identities differ"
        )
    manifest_by_path = {output.path: output.sha256 for output in manifest.outputs}
    expected_outputs = {
        output.role: manifest_by_path.get(output.path)
        for output in contract.expected_outputs
        if output.role != "construction_report"
    }
    if None in expected_outputs.values() or dict(report.output_sha256) != expected_outputs:
        raise ReferenceConcordanceError(
            "reference construction report output identities differ"
        )
    mask_delta = contract.expected_output("mask_delta")
    if report.mask_delta_sha256 != manifest_by_path.get(mask_delta.path):
        raise ReferenceConcordanceError("reference mask-delta identity differs")
    return report


@contextmanager
def open_installed_reference(
    contract: ReferenceIdentityContract,
    *,
    identity_path: str | Path,
    reference_manifest_path: str | Path,
    reference_root: str | Path,
    allow_order_only_difference: bool = False,
    reuse_published_validation: bool = False,
) -> Iterator[tuple[str, ReferenceFetcher, Mapping[str, int]]]:
    """Open an installed reference with strict or manifest-trusted validation.

    Reuse trusts the completed public bundle validation and checks stable file
    metadata, FAI bytes and requested bases. It is not a new FASTA byte audit.
    """

    if type(reuse_published_validation) is not bool:
        raise ReferenceConcordanceError("installed reference validation mode is invalid")

    try:
        root = Path(reference_root).resolve(strict=True)
        manifest_argument = Path(reference_manifest_path)
        manifest_file = manifest_argument.resolve(strict=True)
        contract_sha256 = _file_sha256(Path(identity_path))
        manifest_directory = root / "manifests" / "grch38"
        allowed_manifests = {
            (manifest_directory / f"{contract.resource_id}.json").resolve(),
            (
                manifest_directory
                / f"{contract.resource_id}.{contract_sha256[:12]}.json"
            ).resolve(),
        }
    except OSError as error:
        raise ReferenceConcordanceError(
            "installed reference bundle is unavailable"
        ) from error
    if (
        manifest_argument.is_symlink()
        or manifest_file not in allowed_manifests
        or not manifest_file.is_file()
        or not manifest_file.is_relative_to(root)
    ):
        raise ReferenceConcordanceError(
            "reference manifest must be the canonical bundle manifest"
        )
    try:
        with _open_stable_binary(manifest_file, "reference manifest") as source:
            manifest_content = source.read()
            try:
                manifest_data = json.loads(manifest_content.decode("utf-8"))
            except (UnicodeDecodeError, json.JSONDecodeError) as error:
                raise ReferenceConcordanceError(
                    "invalid reference manifest JSON"
                ) from error
            if not isinstance(manifest_data, Mapping):
                raise ReferenceConcordanceError(
                    "reference manifest root must be an object"
                )
            manifest = PublicResourceManifest.from_mapping(manifest_data)
            verify_reference_manifest_contract(
                contract,
                manifest,
                contract_sha256=contract_sha256,
            )

            # Imported here to avoid the construction module's dependency on this module.
            from genosapiens.grch38_reference import (
                Grch38BootstrapError,
                validate_grch38_bundle,
            )

            if not reuse_published_validation:
                try:
                    validate_grch38_bundle(
                        contract,
                        root,
                        identity_path=identity_path,
                        allow_order_only_difference=allow_order_only_difference,
                    )
                except Grch38BootstrapError as error:
                    raise ReferenceConcordanceError(
                        "installed reference output verification failed"
                    ) from error
            role = contract.reference_role("challenge_alignment")
            fasta_relative = contract.expected_output(role.fasta_output_role).path
            fai_relative = contract.expected_output(role.fai_output_role).path
            manifest_outputs = {output.path: output for output in manifest.outputs}
            if reuse_published_validation:
                report_relative = contract.expected_output("construction_report").path
                report_output = manifest_outputs[report_relative]
                if report_output.size_bytes > 8 * 1024 * 1024:
                    raise ReferenceConcordanceError(
                        "installed reference construction evidence is invalid"
                    )
                report_path = root / report_relative
                with _open_stable_binary(report_path, "reference construction report") as report_source:
                    report_content = report_source.read(8 * 1024 * 1024 + 1)
                if (len(report_content) != report_output.size_bytes
                        or hashlib.sha256(report_content).hexdigest() != report_output.sha256):
                    raise ReferenceConcordanceError(
                        "installed reference construction evidence is invalid"
                    )
                report = load_construction_report(report_path)
                if not report.accepted or report.total_records != contract.sequence_count:
                    raise ReferenceConcordanceError(
                        "installed reference construction evidence is invalid"
                    )
            fasta_output = manifest_outputs[fasta_relative]
            fai_output = manifest_outputs[fai_relative]
            fasta_path = root / fasta_relative
            fai_path = root / fai_relative
            with (
                _open_stable_binary(fai_path, "reference FAI") as fai_source,
                _open_stable_binary(fasta_path, "reference FASTA") as fasta_source,
            ):
                fai_content = fai_source.read()
                if (
                    hashlib.sha256(fai_content).hexdigest() != fai_output.sha256
                    or len(fai_content) != fai_output.size_bytes
                    or not reuse_published_validation
                    and _open_file_sha256(fasta_source) != fasta_output.sha256
                    or os.fstat(fasta_source.fileno()).st_size
                    != fasta_output.size_bytes
                ):
                    raise ReferenceConcordanceError(
                        "open reference files differ from the accepted manifest"
                    )
                index = _parse_fai_bytes(fai_content)
                if (
                    len(index) != contract.sequence_count
                    or any(
                        name not in index or index[name].length != expected_length
                        for name, expected_length in contract.canonical_contig_lengths.items()
                    )
                ):
                    raise ReferenceConcordanceError(
                        "open reference dictionary differs from the pinned contract"
                    )
                reference = _OpenIndexedFasta(fasta_source, index)
                yield (
                    hashlib.sha256(manifest_content).hexdigest(),
                    reference.fetch,
                    {name: entry.length for name, entry in index.items()},
                )
    except ManifestError as error:
        raise ReferenceConcordanceError(
            "installed reference manifest validation failed"
        ) from error


def evaluate_installed_reference_concordance(
    *,
    identity_path: str | Path,
    input_path: str | Path,
    source_input_path: str | Path | None = None,
    source_report_path: str | Path | None = None,
    reference_manifest_path: str | Path,
    reference_root: str | Path,
    source_reference_manifest_path: str | Path | None = None,
    source_reference_root: str | Path | None = None,
    order_adjudication_path: str | Path | None = None,
    read_limits: VcfReadLimits = VcfReadLimits(),
) -> ReferenceConcordanceReport:
    """Validate an input VCF against a fully verified installed reference bundle."""
    identity_file = Path(identity_path)
    contract = load_reference_identity(identity_file)
    with ExitStack() as artifacts:
        adjudication = (
            artifacts.enter_context(
                _open_contig_order_adjudication(order_adjudication_path)
            )
            if order_adjudication_path is not None
            else None
        )
        source_report = (
            artifacts.enter_context(_open_concordance_report(source_report_path))
            if source_report_path is not None
            else None
        )
        effective_source_input_path = source_input_path
        if adjudication is not None:
            if source_report is None:
                raise ReferenceConcordanceError(
                    "contig-order adjudication requires its canonical source report"
                )
            if effective_source_input_path is None:
                effective_source_input_path = input_path
        with open_installed_reference(
            contract,
            identity_path=identity_file,
            reference_manifest_path=reference_manifest_path,
            reference_root=reference_root,
        ) as (manifest_sha256, fetch_reference, reference_contigs):
            if (
                adjudication is not None
                and adjudication.resolution
                == "reconstructed-reference-order-corrected"
            ):
                if source_reference_manifest_path is None:
                    raise ReferenceConcordanceError(
                        "reference-order adjudication requires a source reference manifest"
                    )
                with open_installed_reference(
                    contract,
                    identity_path=identity_file,
                    reference_manifest_path=source_reference_manifest_path,
                    reference_root=source_reference_root or reference_root,
                    allow_order_only_difference=True,
                ) as (
                    source_manifest_sha256,
                    source_fetch_reference,
                    source_reference_contigs,
                ):
                    if (
                        source_manifest_sha256
                        != adjudication.source_reference_manifest_sha256
                    ):
                        raise ReferenceConcordanceError(
                            "contig-order adjudication source reference manifest is stale"
                        )
                    return evaluate_vcf_reference_concordance(
                        input_path,
                        fetch_reference,
                        reference_contigs,
                        contract=contract, read_limits=read_limits,
                        reference_manifest_sha256=manifest_sha256,
                        source_vcf_path=effective_source_input_path,
                        source_report=source_report,
                        source_fetch_reference=source_fetch_reference,
                        source_reference_contigs=source_reference_contigs,
                        order_adjudication=adjudication,
                    )
            return evaluate_vcf_reference_concordance(
                input_path,
                fetch_reference,
                reference_contigs,
                contract=contract, read_limits=read_limits,
                reference_manifest_sha256=manifest_sha256,
                source_vcf_path=effective_source_input_path,
                source_report=source_report,
                order_adjudication=adjudication,
            )


def document_installed_contig_order_adjudication(
    *,
    identity_path: str | Path,
    source_report_path: str | Path,
    source_input_path: str | Path,
    revalidated_input_path: str | Path,
    reference_manifest_path: str | Path,
    reference_root: str | Path,
    source_reference_manifest_path: str | Path | None = None,
    source_reference_root: str | Path | None = None,
    resolution: str,
) -> ContigOrderAdjudication:
    identity_file = Path(identity_path)
    contract = load_reference_identity(identity_file)
    with (
        _open_concordance_report(source_report_path) as source_report,
        open_installed_reference(
            contract,
            identity_path=identity_file,
            reference_manifest_path=reference_manifest_path,
            reference_root=reference_root,
        ) as (
            revalidated_manifest_sha256,
            revalidated_fetch_reference,
            revalidated_reference_contigs,
        ),
    ):
        if resolution == "reconstructed-reference-order-corrected":
            if source_reference_manifest_path is None:
                raise ReferenceConcordanceError(
                    "reference-order adjudication requires a source reference manifest"
                )
            with open_installed_reference(
                contract,
                identity_path=identity_file,
                reference_manifest_path=source_reference_manifest_path,
                reference_root=source_reference_root or reference_root,
                allow_order_only_difference=True,
            ) as (
                source_manifest_sha256,
                source_fetch_reference,
                source_reference_contigs,
            ):
                if source_manifest_sha256 != source_report.reference_manifest_sha256:
                    raise ReferenceConcordanceError(
                        "contig-order source report reference manifest is stale"
                    )
                return document_contig_order_adjudication(
                    source_report,
                    source_input_path=source_input_path,
                    revalidated_input_path=revalidated_input_path,
                    source_fetch_reference=source_fetch_reference,
                    source_reference_contigs=source_reference_contigs,
                    contract=contract,
                    revalidated_reference_manifest_sha256=(
                        revalidated_manifest_sha256
                    ),
                    resolution=resolution,
                )
        return document_contig_order_adjudication(
            source_report,
            source_input_path=source_input_path,
            revalidated_input_path=revalidated_input_path,
            source_fetch_reference=revalidated_fetch_reference,
            source_reference_contigs=revalidated_reference_contigs,
            contract=contract,
            revalidated_reference_manifest_sha256=revalidated_manifest_sha256,
            resolution=resolution,
        )


@dataclass(frozen=True, slots=True)
class RescueReferenceAuthorization:
    reference_identity_id: str
    reference_manifest_sha256: str
    input_sha256: str
    validated_input_sha256: str
    concordance_report_sha256: str
    total_records: int
    matching_records: int
    authorized: bool = True


def require_rescue_reference_concordance(
    *,
    identity_path: str | Path,
    report_path: str | Path,
    input_path: str | Path,
    validated_input_path: str | Path | None = None,
    reference_manifest_path: str | Path,
    reference_root: str | Path,
    order_adjudication_path: str | Path | None = None,
    source_report_path: str | Path | None = None,
    source_reference_manifest_path: str | Path | None = None,
    source_reference_root: str | Path | None = None,
    operator_override: bool = False,
) -> RescueReferenceAuthorization:
    """Replay complete concordance before authorizing rescue."""
    try:
        report_file = Path(report_path)
        with _open_concordance_report(report_file) as report:
            validated_input = (
                Path(validated_input_path)
                if validated_input_path is not None
                else Path(input_path)
            )
            (
                source_input_sha256,
                _source_payload_sha256,
                _source_header_sha256,
            ) = _vcf_bindings(Path(input_path))
            if report.input_sha256 != source_input_sha256:
                raise ReferenceConcordanceError("concordance report input is stale")
            if validated_input_path is None:
                validated_input_sha256 = source_input_sha256
            else:
                (
                    validated_input_sha256,
                    _validated_payload_sha256,
                    _validated_header_sha256,
                ) = _vcf_bindings(validated_input)
            if report.validated_input_sha256 != validated_input_sha256:
                raise ReferenceConcordanceError(
                    "concordance report validated input is stale"
                )
            replayed = evaluate_installed_reference_concordance(
                identity_path=identity_path,
                input_path=validated_input,
                source_input_path=input_path,
                source_report_path=source_report_path,
                reference_manifest_path=reference_manifest_path,
                reference_root=reference_root,
                source_reference_manifest_path=source_reference_manifest_path,
                source_reference_root=source_reference_root,
                order_adjudication_path=order_adjudication_path,
            )
            if report.input_sha256 != replayed.input_sha256:
                raise ReferenceConcordanceError("concordance report input is stale")
            if report.validated_input_sha256 != replayed.validated_input_sha256:
                raise ReferenceConcordanceError(
                    "concordance report validated input is stale"
                )
            if report.reference_manifest_sha256 != replayed.reference_manifest_sha256:
                raise ReferenceConcordanceError(
                    "concordance report reference manifest is stale"
                )
            if report != replayed:
                raise ReferenceConcordanceError(
                    "concordance report differs from complete validation replay"
                )
            if (
                not report.accepted
                or report.reference_role != "challenge_alignment"
                or not report.reference_basename_match
                or not report.contig_names_match
                or not report.contig_lengths_match
                or not report.contig_order_match
                or report.total_records == 0
                or report.matching_records != report.total_records
                or report.mismatching_records != 0
                or report.unresolved_records != 0
            ):
                raise ReferenceConcordanceError(
                    "reference concordance report is not accepted with complete accounting"
                )
            return RescueReferenceAuthorization(
                reference_identity_id=report.reference_identity_id,
                reference_manifest_sha256=report.reference_manifest_sha256,
                input_sha256=report.input_sha256,
                validated_input_sha256=report.validated_input_sha256,
                concordance_report_sha256=_canonical_json_sha256(asdict(report)),
                total_records=report.total_records,
                matching_records=report.matching_records,
            )
    except (OSError, ReferenceConcordanceError) as error:
        reason = str(error)
        if operator_override:
            reason += "; operator override cannot bypass reference concordance"
        raise ReferenceConcordanceError(reason) from error


def rescue_authorization_as_json(authorization: RescueReferenceAuthorization) -> str:
    return json.dumps(asdict(authorization), indent=2, sort_keys=True) + "\n"
