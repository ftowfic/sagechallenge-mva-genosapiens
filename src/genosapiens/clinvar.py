"""Bootstrap and validate a dated ClinVar GRCh38 VCF and tabix index."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import asdict, dataclass
from datetime import date
import gzip
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
from typing import Any
import uuid

from genosapiens.download import DownloadRequest, DownloadResult, download_atomic
from genosapiens.grch38_reference import IndexedFasta
from genosapiens.manifest import (
    ManifestError,
    PublicResourceManifest,
    ResourceOutput,
    RetentionPolicy,
    SourceIdentity,
    Transformation,
    publish_manifest,
    verify_outputs,
)
from genosapiens.reference_concordance import (
    ReferenceConcordanceError,
    ReferenceIdentityContract,
    ReferenceRole,
    ReferenceRoleBinding,
    load_reference_identity,
)
from genosapiens.storage import DEFAULT_RESERVE_BYTES, CapacityPlan, require_capacity
from genosapiens.tabix import TabixError, TabixIndex


SCHEMA = "genosapiens.clinvar-grch38"
TRANSFORMATION_SCHEMA = "genosapiens.clinvar-source-preserving-bootstrap"
SHA256 = re.compile(r"^[0-9a-f]{64}$")
MD5 = re.compile(r"^[0-9a-f]{32}$")
DIGEST = re.compile(r"^sha256:[0-9a-f]{64}$")
IDENTIFIER = re.compile(r"^[a-z0-9][a-z0-9._-]*$")
INFO_HEADER = re.compile(r"^##INFO=<ID=([^,>]+),")
SEQUENCE = re.compile(r"^[ACGTN]+$")
EXPECTED_ROLES = {"vcf", "tbi"}
EXPECTED_SENTINEL_CATEGORIES = {"pathogenic", "benign", "conflicting"}


class ClinvarError(RuntimeError):
    """Raised when ClinVar acquisition, provenance, or validation fails."""


def _required_string(data: Mapping[str, Any], name: str, context: str) -> str:
    value = data.get(name)
    if not isinstance(value, str) or not value:
        raise ClinvarError(f"{context}.{name} must be a non-empty string")
    return value


def _required_positive_int(data: Mapping[str, Any], name: str, context: str) -> int:
    value = data.get(name)
    if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
        raise ClinvarError(f"{context}.{name} must be a positive integer")
    return value


def _optional_string(data: Mapping[str, Any], name: str, context: str) -> str | None:
    value = data.get(name)
    if value is not None and (not isinstance(value, str) or not value):
        raise ClinvarError(f"{context}.{name} must be null or a non-empty string")
    return value


def _safe_relative_path(value: str, context: str) -> PurePosixPath:
    path = PurePosixPath(value)
    if path.is_absolute() or path == PurePosixPath(".") or ".." in path.parts:
        raise ClinvarError(f"{context} must be a safe relative path")
    return path


def _string_tuple(data: Mapping[str, Any], name: str) -> tuple[str, ...]:
    value = data.get(name)
    if (
        not isinstance(value, list)
        or not value
        or any(not isinstance(item, str) or not item for item in value)
        or len(value) != len(set(value))
    ):
        raise ClinvarError(f"ClinVar {name} must be a non-empty unique string array")
    return tuple(value)


@dataclass(frozen=True, slots=True)
class ClinvarArtifact:
    role: str
    url: str
    path: str
    size_bytes: int
    sha256: str
    publisher_md5: str | None

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "ClinvarArtifact":
        expected = {"role", "url", "path", "size_bytes", "sha256", "publisher_md5"}
        if set(data) != expected:
            raise ClinvarError("ClinVar artifact fields are invalid")
        artifact = cls(
            role=_required_string(data, "role", "artifact"),
            url=_required_string(data, "url", "artifact"),
            path=_required_string(data, "path", "artifact"),
            size_bytes=_required_positive_int(data, "size_bytes", "artifact"),
            sha256=_required_string(data, "sha256", "artifact"),
            publisher_md5=_optional_string(data, "publisher_md5", "artifact"),
        )
        if (
            not artifact.url.startswith("https://")
            or not SHA256.fullmatch(artifact.sha256)
            or (
                artifact.publisher_md5 is not None
                and not MD5.fullmatch(artifact.publisher_md5)
            )
        ):
            raise ClinvarError(f"invalid ClinVar artifact identity: {artifact.role}")
        _safe_relative_path(artifact.path, "ClinVar artifact path")
        return artifact


@dataclass(frozen=True, slots=True)
class ClinvarSentinel:
    category: str
    chromosome: str
    position: int
    variation_id: str
    reference: str
    alternate: str
    allele_id: str
    clinical_significance: str
    conflicting_significance: str | None
    review_status: str
    condition_names: str
    condition_identifiers: str
    clinical_sources: str | None
    clinical_hgvs: str

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "ClinvarSentinel":
        expected = {
            "category",
            "chromosome",
            "position",
            "variation_id",
            "reference",
            "alternate",
            "allele_id",
            "clinical_significance",
            "conflicting_significance",
            "review_status",
            "condition_names",
            "condition_identifiers",
            "clinical_sources",
            "clinical_hgvs",
        }
        if set(data) != expected:
            raise ClinvarError("ClinVar sentinel fields are invalid")
        required = expected - {
            "position",
            "conflicting_significance",
            "clinical_sources",
        }
        sentinel = cls(
            **{name: _required_string(data, name, "sentinel") for name in required},
            position=_required_positive_int(data, "position", "sentinel"),
            conflicting_significance=_optional_string(
                data, "conflicting_significance", "sentinel"
            ),
            clinical_sources=_optional_string(data, "clinical_sources", "sentinel"),
        )
        if (
            sentinel.category not in EXPECTED_SENTINEL_CATEGORIES
            or not sentinel.variation_id.isdecimal()
            or not sentinel.allele_id.isdecimal()
            or not SEQUENCE.fullmatch(sentinel.reference)
            or not SEQUENCE.fullmatch(sentinel.alternate)
            or sentinel.reference == sentinel.alternate
            or (
                sentinel.category == "conflicting"
                and sentinel.conflicting_significance is None
            )
            or (
                sentinel.category != "conflicting"
                and sentinel.conflicting_significance is not None
            )
        ):
            raise ClinvarError("ClinVar sentinel identity or assertion is invalid")
        return sentinel


@dataclass(frozen=True, slots=True)
class ClinvarContract:
    resource_id: str
    release: str
    release_date: str
    assembly: str
    reference_binding: ReferenceRoleBinding
    license: str
    license_url: str
    published_checksum_source: str
    required_info_fields: tuple[str, ...]
    required_review_fields: tuple[str, ...]
    required_index_contigs: tuple[str, ...]
    artifacts: tuple[ClinvarArtifact, ...]
    sentinels: tuple[ClinvarSentinel, ...]
    schema: str = SCHEMA

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "ClinvarContract":
        expected = {
            "schema",
            "resource_id",
            "release",
            "release_date",
            "assembly",
            "reference_binding",
            "license",
            "license_url",
            "published_checksum_source",
            "required_info_fields",
            "required_review_fields",
            "required_index_contigs",
            "artifacts",
            "sentinels",
        }
        artifact_data = data.get("artifacts")
        sentinel_data = data.get("sentinels")
        reference_binding = data.get("reference_binding")
        if (
            set(data) != expected
            or data.get("schema") != SCHEMA
            or not isinstance(artifact_data, list)
            or not artifact_data
            or not all(isinstance(item, Mapping) for item in artifact_data)
            or not isinstance(sentinel_data, list)
            or not sentinel_data
            or not all(isinstance(item, Mapping) for item in sentinel_data)
            or not isinstance(reference_binding, Mapping)
        ):
            raise ClinvarError("ClinVar contract fields or schema are invalid")
        try:
            contract = cls(
                schema=SCHEMA,
                resource_id=_required_string(data, "resource_id", "contract"),
                release=_required_string(data, "release", "contract"),
                release_date=_required_string(data, "release_date", "contract"),
                assembly=_required_string(data, "assembly", "contract"),
                reference_binding=ReferenceRoleBinding.from_mapping(
                    reference_binding
                ),
                license=_required_string(data, "license", "contract"),
                license_url=_required_string(data, "license_url", "contract"),
                published_checksum_source=_required_string(
                    data, "published_checksum_source", "contract"
                ),
                required_info_fields=_string_tuple(data, "required_info_fields"),
                required_review_fields=_string_tuple(data, "required_review_fields"),
                required_index_contigs=_string_tuple(data, "required_index_contigs"),
                artifacts=tuple(
                    ClinvarArtifact.from_mapping(item) for item in artifact_data
                ),
                sentinels=tuple(
                    ClinvarSentinel.from_mapping(item) for item in sentinel_data
                ),
            )
        except ReferenceConcordanceError as error:
            raise ClinvarError(f"invalid ClinVar reference binding: {error}") from error
        try:
            parsed_date = date.fromisoformat(contract.release_date)
        except ValueError as error:
            raise ClinvarError("ClinVar release_date must be an ISO date") from error
        roles = [artifact.role for artifact in contract.artifacts]
        paths = [artifact.path for artifact in contract.artifacts]
        categories = [sentinel.category for sentinel in contract.sentinels]
        vcf_path = contract.artifact("vcf").path
        if (
            not IDENTIFIER.fullmatch(contract.resource_id)
            or contract.release != contract.release_date
            or parsed_date.isoformat() != contract.release_date
            or not contract.license_url.startswith("https://")
            or not contract.published_checksum_source.startswith("https://")
            or set(roles) != EXPECTED_ROLES
            or len(roles) != len(set(roles))
            or len(paths) != len(set(paths))
            or len({str(PurePosixPath(path).parent) for path in paths}) != 1
            or contract.artifact("tbi").path != f"{vcf_path}.tbi"
            or set(contract.required_review_fields)
            - set(contract.required_info_fields)
            or set(categories) != EXPECTED_SENTINEL_CATEGORIES
            or len(categories) != len(set(categories))
        ):
            raise ClinvarError("ClinVar contract identities are inconsistent")
        return contract

    def artifact(self, role: str) -> ClinvarArtifact:
        for artifact in self.artifacts:
            if artifact.role == role:
                return artifact
        raise ClinvarError(f"ClinVar contract lacks artifact role {role}")


def load_clinvar_contract(path: str | Path) -> ClinvarContract:
    contract_path = Path(path)
    try:
        data = json.loads(contract_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ClinvarError(f"cannot read ClinVar contract {contract_path}: {error}") from error
    if not isinstance(data, Mapping):
        raise ClinvarError("ClinVar contract root must be an object")
    return ClinvarContract.from_mapping(data)


@dataclass(frozen=True, slots=True)
class ClinvarAssertion:
    category: str
    chromosome: str
    position: int
    variation_id: str
    reference: str
    alternate: str
    allele_id: str
    clinical_significance: str
    conflicting_significance: str | None
    review_status: str
    condition_names: str
    condition_identifiers: str
    clinical_sources: str | None
    clinical_hgvs: str


@dataclass(frozen=True, slots=True)
class ClinvarValidation:
    release_date: str
    reference: str
    info_fields: tuple[str, ...]
    review_fields: tuple[str, ...]
    indexed_contigs: tuple[str, ...]
    assertions: tuple[ClinvarAssertion, ...]


@dataclass(frozen=True, slots=True)
class ClinvarBootstrapResult:
    resource_id: str
    release: str
    manifest_path: str
    manifest_sha256: str
    validation: ClinvarValidation
    reused: bool


Downloader = Callable[[DownloadRequest, str | Path], DownloadResult]
CapacityChecker = Callable[[str | Path, CapacityPlan], Any]


def _file_identity(path: Path) -> tuple[str, str, int]:
    sha256 = hashlib.sha256()
    md5 = hashlib.md5(usedforsecurity=False)
    size = 0
    try:
        with path.open("rb") as source:
            while chunk := source.read(1024 * 1024):
                sha256.update(chunk)
                md5.update(chunk)
                size += len(chunk)
    except OSError as error:
        raise ClinvarError(f"cannot read ClinVar artifact {path}: {error}") from error
    return sha256.hexdigest(), md5.hexdigest(), size


def _verify_artifact(path: Path, artifact: ClinvarArtifact) -> None:
    sha256, md5, size = _file_identity(path)
    if (
        size != artifact.size_bytes
        or sha256 != artifact.sha256
        or (artifact.publisher_md5 is not None and md5 != artifact.publisher_md5)
    ):
        raise ClinvarError(f"ClinVar artifact checksum mismatch: {artifact.role}")


def _read_vcf_header(
    path: Path, contract: ClinvarContract
) -> tuple[str, str, tuple[str, ...]]:
    release_date: str | None = None
    reference: str | None = None
    info_fields: list[str] = []
    columns: tuple[str, ...] | None = None
    try:
        with gzip.open(path, "rt", encoding="utf-8") as source:
            for line in source:
                value = line.rstrip("\r\n")
                if value.startswith("##fileDate="):
                    release_date = value.partition("=")[2]
                elif value.startswith("##reference="):
                    reference = value.partition("=")[2]
                else:
                    match = INFO_HEADER.match(value)
                    if match:
                        info_fields.append(match.group(1))
                if value.startswith("#CHROM\t"):
                    columns = tuple(value.split("\t"))
                    break
                if not value.startswith("#"):
                    raise ClinvarError("ClinVar VCF lacks a complete header")
    except (OSError, EOFError, gzip.BadGzipFile, UnicodeDecodeError) as error:
        raise ClinvarError(f"cannot read ClinVar VCF header: {error}") from error
    if release_date != contract.release_date or reference != contract.assembly:
        raise ClinvarError("ClinVar VCF release date or assembly differs from contract")
    if columns != ("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"):
        raise ClinvarError("ClinVar VCF must contain exactly eight site-only columns")
    if len(info_fields) != len(set(info_fields)):
        raise ClinvarError("ClinVar VCF has duplicate INFO declarations")
    missing = sorted(set(contract.required_info_fields) - set(info_fields))
    if missing:
        raise ClinvarError(f"ClinVar VCF lacks required INFO fields: {', '.join(missing)}")
    return release_date, reference, tuple(info_fields)


def _parse_info(value: str) -> Mapping[str, str | None]:
    parsed: dict[str, str | None] = {}
    for item in value.split(";"):
        key, separator, field_value = item.partition("=")
        if not key or key in parsed:
            raise ClinvarError("ClinVar record contains malformed or duplicate INFO fields")
        parsed[key] = field_value if separator else None
    return parsed


def _assertion_from_record(record: str, sentinel: ClinvarSentinel) -> ClinvarAssertion:
    fields = record.rstrip("\r\n").split("\t")
    if len(fields) != 8:
        raise ClinvarError("ClinVar indexed record does not have eight columns")
    chromosome, position, variation_id, reference, alternate = fields[:5]
    try:
        numeric_position = int(position)
    except ValueError as error:
        raise ClinvarError("ClinVar indexed record position is invalid") from error
    if (
        chromosome,
        numeric_position,
        variation_id,
        reference,
        alternate,
    ) != (
        sentinel.chromosome,
        sentinel.position,
        sentinel.variation_id,
        sentinel.reference,
        sentinel.alternate,
    ):
        raise ClinvarError(f"ClinVar sentinel identity mismatch: {sentinel.category}")
    info = _parse_info(fields[7])
    values = {
        "allele_id": info.get("ALLELEID"),
        "clinical_significance": info.get("CLNSIG"),
        "conflicting_significance": info.get("CLNSIGCONF"),
        "review_status": info.get("CLNREVSTAT"),
        "condition_names": info.get("CLNDN"),
        "condition_identifiers": info.get("CLNDISDB"),
        "clinical_sources": info.get("CLNVI"),
        "clinical_hgvs": info.get("CLNHGVS"),
    }
    expected = {
        "allele_id": sentinel.allele_id,
        "clinical_significance": sentinel.clinical_significance,
        "conflicting_significance": sentinel.conflicting_significance,
        "review_status": sentinel.review_status,
        "condition_names": sentinel.condition_names,
        "condition_identifiers": sentinel.condition_identifiers,
        "clinical_sources": sentinel.clinical_sources,
        "clinical_hgvs": sentinel.clinical_hgvs,
    }
    if values != expected:
        raise ClinvarError(f"ClinVar sentinel assertion mismatch: {sentinel.category}")
    return ClinvarAssertion(
        category=sentinel.category,
        chromosome=chromosome,
        position=numeric_position,
        variation_id=variation_id,
        reference=reference,
        alternate=alternate,
        **values,
    )


def _validate_paths(
    contract: ClinvarContract,
    reference_contract: ReferenceIdentityContract,
    reference_role: ReferenceRole,
    reference_root: Path,
    vcf_path: Path,
    tbi_path: Path,
) -> ClinvarValidation:
    _verify_artifact(vcf_path, contract.artifact("vcf"))
    _verify_artifact(tbi_path, contract.artifact("tbi"))
    release_date, reference, info_fields = _read_vcf_header(vcf_path, contract)
    try:
        index = TabixIndex.load(tbi_path)
        index.require_vcf_layout()
    except TabixError as error:
        raise ClinvarError(f"ClinVar tabix index is invalid: {error}") from error
    indexed_contigs = tuple(item.name for item in index.references)
    if indexed_contigs != contract.required_index_contigs:
        raise ClinvarError("ClinVar tabix contigs differ from the pinned contract")

    fasta_output = reference_contract.expected_output(reference_role.fasta_output_role)
    fasta_path = reference_root / fasta_output.path
    fai_path = reference_root / reference_contract.expected_output(
        reference_role.fai_output_role
    ).path
    reference_fasta = IndexedFasta(fasta_path, fai_path)
    try:
        aliases = json.loads(
            (
                reference_root
                / reference_contract.expected_output("contig_aliases").path
            ).read_text(encoding="utf-8")
        )["alignment_to_tool"]
    except (KeyError, OSError, json.JSONDecodeError, TypeError) as error:
        raise ClinvarError(f"cannot load reference contig aliases: {error}") from error
    if not isinstance(aliases, Mapping):
        raise ClinvarError("reference alignment_to_tool aliases must be an object")
    assertions: list[ClinvarAssertion] = []
    for sentinel in contract.sentinels:
        try:
            records = index.query(
                vcf_path,
                sentinel.chromosome,
                sentinel.position,
                sentinel.position,
            )
        except (KeyError, OSError, TabixError) as error:
            raise ClinvarError(
                f"ClinVar tabix query failed for {sentinel.category}: {error}"
            ) from error
        matches = [
            record
            for record in records
            if record.split("\t", 5)[:5]
            == [
                sentinel.chromosome,
                str(sentinel.position),
                sentinel.variation_id,
                sentinel.reference,
                sentinel.alternate,
            ]
        ]
        if len(matches) != 1:
            raise ClinvarError(
                f"ClinVar index did not recover exactly one {sentinel.category} sentinel"
            )
        reference_contig = aliases.get(sentinel.chromosome, sentinel.chromosome)
        if not isinstance(reference_contig, str) or (
            reference_fasta.fetch(
                reference_contig, sentinel.position, len(sentinel.reference)
            )
            != sentinel.reference
        ):
            raise ClinvarError(f"ClinVar sentinel REF mismatch: {sentinel.category}")
        assertions.append(_assertion_from_record(matches[0], sentinel))
    return ClinvarValidation(
        release_date=release_date,
        reference=reference,
        info_fields=tuple(
            field for field in contract.required_info_fields if field in info_fields
        ),
        review_fields=contract.required_review_fields,
        indexed_contigs=indexed_contigs,
        assertions=tuple(assertions),
    )


def validate_clinvar_bundle(
    contract: ClinvarContract,
    reference_contract: ReferenceIdentityContract,
    reference_root: str | Path,
) -> ClinvarValidation:
    root = Path(reference_root).resolve(strict=True)
    try:
        reference_role = contract.reference_binding.resolve(
            reference_contract, required_role="public_annotation_benchmark"
        )
    except ReferenceConcordanceError as error:
        raise ClinvarError(str(error)) from error
    return _validate_paths(
        contract,
        reference_contract,
        reference_role,
        root,
        root / contract.artifact("vcf").path,
        root / contract.artifact("tbi").path,
    )


def _manifest_path(root: Path, contract: ClinvarContract) -> Path:
    return root / f"manifests/clinvar/{contract.resource_id}.json"


def _load_manifest(path: Path) -> PublicResourceManifest:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ClinvarError(f"cannot read ClinVar manifest {path}: {error}") from error
    if not isinstance(data, Mapping):
        raise ClinvarError("ClinVar manifest root must be an object")
    try:
        return PublicResourceManifest.from_mapping(data)
    except ManifestError as error:
        raise ClinvarError(f"invalid ClinVar manifest: {error}") from error


def _copy_file(source: Path, destination: Path) -> None:
    try:
        with source.open("rb") as input_file, destination.open("xb") as output_file:
            shutil.copyfileobj(input_file, output_file, length=1024 * 1024)
            output_file.flush()
            os.fsync(output_file.fileno())
        os.chmod(destination, 0o444)
    except OSError as error:
        raise ClinvarError(f"cannot stage ClinVar artifact {source}: {error}") from error


def _require_within_root(root: Path, path: Path, context: str) -> None:
    if not path.is_absolute() or not path.is_relative_to(root):
        raise ClinvarError(f"{context} escapes reference_root: {path}")
    current = root
    for part in path.relative_to(root).parts:
        current /= part
        if current.is_symlink():
            raise ClinvarError(f"{context} contains a symlink: {current}")


def _verify_release_tree(release: Path, contract: ClinvarContract) -> None:
    if release.is_symlink() or not release.is_dir():
        raise ClinvarError(f"ClinVar release is not a regular directory: {release}")
    expected = {PurePosixPath(artifact.path).name for artifact in contract.artifacts}
    try:
        observed = {path.name for path in release.iterdir()}
    except OSError as error:
        raise ClinvarError(f"cannot inspect ClinVar release {release}: {error}") from error
    if observed != expected or any(not path.is_file() for path in release.iterdir()):
        raise ClinvarError("ClinVar release is incomplete or contaminated")


def _manifest(
    contract: ClinvarContract,
    contract_sha256: str,
    image_digest: str,
    retrieved_on: str,
) -> PublicResourceManifest:
    return PublicResourceManifest(
        resource_id=contract.resource_id,
        release=contract.release,
        source=SourceIdentity(
            url=contract.artifact("vcf").url,
            sha256=contract.artifact("vcf").sha256,
        ),
        license=contract.license,
        assembly=contract.assembly,
        retrieved_on=retrieved_on,
        transformation=Transformation(TRANSFORMATION_SCHEMA, contract_sha256),
        image_digest=image_digest,
        outputs=tuple(
            ResourceOutput(
                path=artifact.path,
                size_bytes=artifact.size_bytes,
                sha256=artifact.sha256,
                retention=RetentionPolicy("durable", None),
            )
            for artifact in contract.artifacts
        ),
    )


def _verify_manifest_contract(
    manifest: PublicResourceManifest,
    contract: ClinvarContract,
    contract_sha256: str,
) -> None:
    if (
        manifest.resource_id != contract.resource_id
        or manifest.release != contract.release
        or manifest.assembly != contract.assembly
        or manifest.license != contract.license
        or manifest.source.url != contract.artifact("vcf").url
        or manifest.source.sha256 != contract.artifact("vcf").sha256
        or manifest.transformation != Transformation(TRANSFORMATION_SCHEMA, contract_sha256)
        or {(output.path, output.size_bytes, output.sha256) for output in manifest.outputs}
        != {(artifact.path, artifact.size_bytes, artifact.sha256) for artifact in contract.artifacts}
    ):
        raise ClinvarError("ClinVar manifest identity differs from its contract")


def bootstrap_clinvar(
    *,
    contract_path: str | Path,
    reference_identity_path: str | Path,
    reference_root: str | Path,
    image_digest: str,
    manifest_path: str | Path | None = None,
    retrieved_on: str | None = None,
    scratch_root: str | Path | None = None,
    reserve_bytes: int = DEFAULT_RESERVE_BYTES,
    downloader: Downloader = download_atomic,
    capacity_checker: CapacityChecker = require_capacity,
) -> ClinvarBootstrapResult:
    contract_file = Path(contract_path)
    root = Path(reference_root).resolve(strict=True)
    contract = load_clinvar_contract(contract_file)
    try:
        reference_contract = load_reference_identity(reference_identity_path)
    except ReferenceConcordanceError as error:
        raise ClinvarError(str(error)) from error
    try:
        reference_role = contract.reference_binding.resolve(
            reference_contract, required_role="public_annotation_benchmark"
        )
    except ReferenceConcordanceError as error:
        raise ClinvarError(str(error)) from error
    if not DIGEST.fullmatch(image_digest):
        raise ClinvarError("ClinVar bootstrap image_digest must be a SHA-256 digest")
    contract_sha256 = _file_identity(contract_file)[0]
    release = root / str(PurePosixPath(contract.artifact("vcf").path).parent)
    if manifest_path is None:
        manifest_path = _manifest_path(root, contract)
    else:
        manifest_path = Path(manifest_path)
        if not manifest_path.is_absolute():
            manifest_path = root / manifest_path
        if (
            ".." in manifest_path.parts
            or not manifest_path.is_relative_to(root / "manifests/clinvar")
            or manifest_path.suffix != ".json"
        ):
            raise ClinvarError("ClinVar manifest destination must be under manifests/clinvar")
    _require_within_root(root, release, "ClinVar release")
    _require_within_root(root, manifest_path, "ClinVar manifest")

    if manifest_path.exists() or manifest_path.is_symlink():
        manifest = _load_manifest(manifest_path)
        _verify_manifest_contract(manifest, contract, contract_sha256)
        try:
            verify_outputs(manifest, root)
        except ManifestError as error:
            raise ClinvarError(f"completed ClinVar bundle is invalid: {error}") from error
        _verify_release_tree(release, contract)
        validation = validate_clinvar_bundle(contract, reference_contract, root)
        return ClinvarBootstrapResult(
            resource_id=contract.resource_id,
            release=contract.release,
            manifest_path=str(manifest_path),
            manifest_sha256=_file_identity(manifest_path)[0],
            validation=validation,
            reused=True,
        )

    retrieval_date = retrieved_on or date.today().isoformat()
    if release.exists() or release.is_symlink():
        _verify_release_tree(release, contract)
        validation = validate_clinvar_bundle(contract, reference_contract, root)
        publication = publish_manifest(
            _manifest(contract, contract_sha256, image_digest, retrieval_date),
            manifest_path,
            root,
        )
        return ClinvarBootstrapResult(
            resource_id=contract.resource_id,
            release=contract.release,
            manifest_path=str(publication.path),
            manifest_sha256=publication.manifest_sha256,
            validation=validation,
            reused=False,
        )

    scratch_base = (
        Path(scratch_root) if scratch_root is not None else root / "scratch/public-build"
    )
    if not scratch_base.is_absolute():
        raise ClinvarError("ClinVar scratch_root must be absolute")
    _require_within_root(root, scratch_base, "ClinVar scratch_root")
    scratch = scratch_base / "clinvar" / contract.release
    _require_within_root(root, scratch, "ClinVar scratch")
    scratch.mkdir(parents=True, exist_ok=True, mode=0o700)
    source_paths = {
        artifact.role: scratch / PurePosixPath(artifact.path).name
        for artifact in contract.artifacts
    }
    missing = [
        artifact
        for artifact in contract.artifacts
        if not source_paths[artifact.role].exists()
    ]
    capacity_checker(
        root,
        CapacityPlan(
            operation=f"bootstrap ClinVar GRCh38 {contract.release}",
            source_bytes=sum(artifact.size_bytes for artifact in missing),
            output_bytes=sum(artifact.size_bytes for artifact in contract.artifacts),
            reserve_bytes=reserve_bytes,
        ),
    )
    for artifact in contract.artifacts:
        downloader(
            DownloadRequest(
                url=artifact.url,
                destination=source_paths[artifact.role],
                expected_size_bytes=artifact.size_bytes,
                expected_sha256=artifact.sha256,
                expected_md5=artifact.publisher_md5,
                attempts=3,
                timeout_seconds=300,
            ),
            root,
        )

    release.parent.mkdir(parents=True, exist_ok=True)
    partial_root = release.parent / ".partial"
    _require_within_root(root, partial_root, "ClinVar release partial path")
    partial_root.mkdir(mode=0o700, exist_ok=True)
    staging = partial_root / f"{release.name}.{uuid.uuid4().hex}.partial"
    staging.mkdir(mode=0o700)
    try:
        for artifact in contract.artifacts:
            _copy_file(
                source_paths[artifact.role], staging / PurePosixPath(artifact.path).name
            )
        validation = _validate_paths(
            contract,
            reference_contract,
            reference_role,
            root,
            staging / PurePosixPath(contract.artifact("vcf").path).name,
            staging / PurePosixPath(contract.artifact("tbi").path).name,
        )
        os.chmod(staging, 0o755)
        try:
            os.replace(staging, release)
        except OSError:
            if not release.exists():
                raise
            shutil.rmtree(staging, ignore_errors=True)
            _verify_release_tree(release, contract)
            validation = validate_clinvar_bundle(contract, reference_contract, root)
        directory_fd = os.open(release.parent, os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    validation = validate_clinvar_bundle(contract, reference_contract, root)
    publication = publish_manifest(
        _manifest(contract, contract_sha256, image_digest, retrieval_date),
        manifest_path,
        root,
    )
    for path in source_paths.values():
        path.unlink(missing_ok=True)
    return ClinvarBootstrapResult(
        resource_id=contract.resource_id,
        release=contract.release,
        manifest_path=str(publication.path),
        manifest_sha256=publication.manifest_sha256,
        validation=validation,
        reused=publication.reused,
    )


def clinvar_result_as_json(result: ClinvarBootstrapResult) -> str:
    return json.dumps(asdict(result), indent=2, sort_keys=True) + "\n"
