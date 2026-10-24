"""Bootstrap and validate the public DeepVariant HG002 35x callset."""

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
from urllib.parse import urlparse
import uuid

from genosapiens.download import DownloadRequest, DownloadResult, download_atomic
from genosapiens.grch38_reference import IndexedFasta, load_fasta_index
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
from genosapiens.storage import CapacityPlan, require_capacity
from genosapiens.tabix import TabixError, TabixIndex


SCHEMA = "genosapiens.deepvariant-caller-vcf"
TRANSFORMATION_SCHEMA = "genosapiens.deepvariant-caller-vcf-bootstrap"
MD5 = re.compile(r"^[0-9a-f]{32}$")
GIT_SHA = re.compile(r"^[0-9a-f]{40}$")
SEQUENCE = re.compile(r"^[ACGTRYSWKMBDHVN]+$")
GENOTYPE = re.compile(r"^[0-9]+(?:[/|][0-9]+)+$")
CONTIG_HEADER = re.compile(r"^##contig=<ID=([^,>]+),length=([0-9]+)(?:,.*)?>$")
FORMAT_HEADER = re.compile(r"^##FORMAT=<ID=([^,>]+),")
EXPECTED_ROLES = {"caller_vcf", "caller_tbi", "execution_log"}
AUTOSOMES = {f"chr{number}" for number in range(1, 23)}


class DeepVariantError(RuntimeError):
    """Raised when the public DeepVariant callset is invalid."""


@dataclass(frozen=True, slots=True)
class DocumentationIdentity:
    url: str
    git_blob_sha: str

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "DocumentationIdentity":
        if set(data) != {"url", "git_blob_sha"}:
            raise DeepVariantError("DeepVariant documentation fields are invalid")
        try:
            value = cls(**data)
        except TypeError as error:
            raise DeepVariantError(f"invalid DeepVariant documentation: {error}") from error
        if not value.url.startswith("https://") or not GIT_SHA.fullmatch(
            value.git_blob_sha
        ):
            raise DeepVariantError("DeepVariant documentation identity is invalid")
        return value


@dataclass(frozen=True, slots=True)
class SourceInput:
    url: str
    size_bytes: int
    md5: str
    object_generation: str

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "SourceInput":
        if set(data) != {"url", "size_bytes", "md5", "object_generation"}:
            raise DeepVariantError("DeepVariant source input fields are invalid")
        try:
            value = cls(**data)
        except TypeError as error:
            raise DeepVariantError(f"invalid DeepVariant source input: {error}") from error
        if (
            not value.url.startswith("https://")
            or not isinstance(value.size_bytes, int)
            or isinstance(value.size_bytes, bool)
            or value.size_bytes <= 0
            or not MD5.fullmatch(value.md5)
            or not value.object_generation
        ):
            raise DeepVariantError("DeepVariant source input identity is invalid")
        return value


@dataclass(frozen=True, slots=True)
class DeepVariantArtifact:
    role: str
    url: str
    path: str
    size_bytes: int
    md5: str
    object_generation: str

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "DeepVariantArtifact":
        allowed = {"role", "url", "path", "size_bytes", "md5", "object_generation"}
        if set(data) != allowed:
            raise DeepVariantError("DeepVariant artifact fields are invalid")
        try:
            value = cls(**data)
        except TypeError as error:
            raise DeepVariantError(f"invalid DeepVariant artifact: {error}") from error
        path = PurePosixPath(value.path)
        if (
            not value.role
            or not value.url.startswith("https://")
            or path.is_absolute()
            or ".." in path.parts
            or path == PurePosixPath(".")
            or not isinstance(value.size_bytes, int)
            or isinstance(value.size_bytes, bool)
            or value.size_bytes <= 0
            or not MD5.fullmatch(value.md5)
            or not value.object_generation
        ):
            raise DeepVariantError(f"invalid DeepVariant artifact: {value.role}")
        return value


@dataclass(frozen=True, slots=True)
class DeepVariantContract:
    resource_id: str
    release: str
    sample: str
    assembly: str
    reference_binding: ReferenceRoleBinding
    caller: str
    caller_version: str
    model_type: str
    sequencing: str
    license: str
    documentation: DocumentationIdentity
    input_bam: SourceInput
    required_command_tokens: tuple[str, ...]
    artifacts: tuple[DeepVariantArtifact, ...]
    schema: str = SCHEMA

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "DeepVariantContract":
        allowed = {
            "schema",
            "resource_id",
            "release",
            "sample",
            "assembly",
            "reference_binding",
            "caller",
            "caller_version",
            "model_type",
            "sequencing",
            "license",
            "documentation",
            "input_bam",
            "required_command_tokens",
            "artifacts",
        }
        if set(data) != allowed or data.get("schema") != SCHEMA:
            raise DeepVariantError("DeepVariant contract fields or schema are invalid")
        documentation = data.get("documentation")
        input_bam = data.get("input_bam")
        reference_binding = data.get("reference_binding")
        tokens = data.get("required_command_tokens")
        artifacts = data.get("artifacts")
        if (
            not isinstance(documentation, Mapping)
            or not isinstance(input_bam, Mapping)
            or not isinstance(reference_binding, Mapping)
            or not isinstance(tokens, list)
            or not tokens
            or not all(isinstance(token, str) and token for token in tokens)
            or not isinstance(artifacts, list)
            or not all(isinstance(artifact, Mapping) for artifact in artifacts)
        ):
            raise DeepVariantError("DeepVariant contract nested fields are invalid")
        try:
            contract = cls(
                schema=data["schema"],
                resource_id=data["resource_id"],
                release=data["release"],
                sample=data["sample"],
                assembly=data["assembly"],
                reference_binding=ReferenceRoleBinding.from_mapping(
                    reference_binding
                ),
                caller=data["caller"],
                caller_version=data["caller_version"],
                model_type=data["model_type"],
                sequencing=data["sequencing"],
                license=data["license"],
                documentation=DocumentationIdentity.from_mapping(documentation),
                input_bam=SourceInput.from_mapping(input_bam),
                required_command_tokens=tuple(tokens),
                artifacts=tuple(
                    DeepVariantArtifact.from_mapping(artifact)
                    for artifact in artifacts
                ),
            )
        except (KeyError, TypeError, ReferenceConcordanceError) as error:
            raise DeepVariantError(f"invalid DeepVariant contract: {error}") from error
        strings = (
            contract.resource_id,
            contract.release,
            contract.sample,
            contract.assembly,
            contract.caller,
            contract.caller_version,
            contract.model_type,
            contract.sequencing,
            contract.license,
        )
        if not all(isinstance(value, str) and value for value in strings):
            raise DeepVariantError("DeepVariant contract strings must be non-empty")
        if (
            contract.sample != "HG002"
            or contract.caller != "DeepVariant"
            or contract.caller_version != "1.0.0"
            or contract.model_type != "WGS"
            or "35x" not in contract.sequencing
            or "--regions" in contract.required_command_tokens
        ):
            raise DeepVariantError("DeepVariant contract does not identify the 35x WGS callset")
        roles = {artifact.role for artifact in contract.artifacts}
        if roles != EXPECTED_ROLES or len(roles) != len(contract.artifacts):
            raise DeepVariantError("DeepVariant contract must define each artifact once")
        parents = {str(PurePosixPath(artifact.path).parent) for artifact in contract.artifacts}
        if len(parents) != 1:
            raise DeepVariantError("DeepVariant artifacts must share one release directory")
        return contract

    def artifact(self, role: str) -> DeepVariantArtifact:
        for artifact in self.artifacts:
            if artifact.role == role:
                return artifact
        raise DeepVariantError(f"DeepVariant contract lacks artifact role {role}")


def load_deepvariant_contract(path: str | Path) -> DeepVariantContract:
    contract_path = Path(path)
    try:
        data = json.loads(contract_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise DeepVariantError(
            f"cannot read DeepVariant contract {contract_path}: {error}"
        ) from error
    if not isinstance(data, Mapping):
        raise DeepVariantError("DeepVariant contract root must be an object")
    return DeepVariantContract.from_mapping(data)


@dataclass(frozen=True, slots=True)
class DeepVariantValidation:
    variant_records: int
    pass_records: int
    filtered_records: int
    indexed_contigs: int


@dataclass(frozen=True, slots=True)
class DeepVariantBootstrapResult:
    manifest_path: str
    manifest_sha256: str
    sample: str
    release: str
    validation: DeepVariantValidation
    reused: bool


@dataclass(frozen=True, slots=True)
class _Sentinel:
    contig: str
    position: int
    ref: str
    filter_value: str


Downloader = Callable[[DownloadRequest, str | Path], DownloadResult]


def _checksums(path: Path) -> tuple[str, str, int]:
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
        raise DeepVariantError(f"cannot read {path}: {error}") from error
    return sha256.hexdigest(), md5.hexdigest(), size


def _verify_artifact(path: Path, artifact: DeepVariantArtifact) -> str:
    sha256, md5, size = _checksums(path)
    if size != artifact.size_bytes or md5 != artifact.md5:
        raise DeepVariantError(
            f"DeepVariant published checksum mismatch: {artifact.role}"
        )
    return sha256


def _validate_execution_log(path: Path, contract: DeepVariantContract) -> None:
    found = False
    try:
        with path.open("rt", encoding="utf-8") as source:
            for line in source:
                if all(token in line for token in contract.required_command_tokens):
                    if "--regions" in line:
                        raise DeepVariantError(
                            "DeepVariant execution log describes a region-restricted callset"
                        )
                    found = True
                    break
    except (OSError, UnicodeDecodeError) as error:
        raise DeepVariantError(f"cannot read DeepVariant execution log: {error}") from error
    if not found:
        raise DeepVariantError(
            "DeepVariant execution log lacks the pinned full-WGS command"
        )


def _complete_genotype(format_keys: list[str], sample_values: list[str]) -> bool:
    if len(format_keys) != len(sample_values):
        return False
    values = dict(zip(format_keys, sample_values, strict=True))
    genotype = values.get("GT", "")
    if not GENOTYPE.fullmatch(genotype):
        return False
    if not values.get("GQ", "").isdigit() or not values.get("DP", "").isdigit():
        return False
    allelic_depths = values.get("AD", "").split(",")
    return len(allelic_depths) >= 2 and all(value.isdigit() for value in allelic_depths)


def _scan_vcf(
    path: Path,
    contract: DeepVariantContract,
    reference_contract: ReferenceIdentityContract,
    reference_role: ReferenceRole,
    reference_root: Path,
) -> tuple[DeepVariantValidation, _Sentinel, _Sentinel]:
    canonical_fasta = reference_root / reference_contract.expected_output(
        reference_role.fasta_output_role
    ).path
    canonical_fai = reference_root / reference_contract.expected_output(
        reference_role.fai_output_role
    ).path
    reference_index = load_fasta_index(canonical_fai)
    contig_lengths: dict[str, int] = {}
    format_definitions: set[str] = set()
    samples: tuple[str, ...] | None = None
    version_seen = False
    pass_count = 0
    filtered_count = 0
    record_count = 0
    pass_sentinel: _Sentinel | None = None
    filtered_sentinel: _Sentinel | None = None
    try:
        with gzip.open(path, "rt", encoding="utf-8") as source:
            for line in source:
                stripped = line.rstrip("\n")
                if stripped == f"##DeepVariant_version={contract.caller_version}":
                    version_seen = True
                    continue
                if line.startswith("##contig="):
                    match = CONTIG_HEADER.match(stripped)
                    if match:
                        contig_lengths[match.group(1)] = int(match.group(2))
                    continue
                if line.startswith("##FORMAT="):
                    match = FORMAT_HEADER.match(stripped)
                    if match:
                        format_definitions.add(match.group(1))
                    continue
                if line.startswith("#CHROM\t"):
                    samples = tuple(stripped.split("\t")[9:])
                    continue
                if line.startswith("#"):
                    continue
                fields = stripped.split("\t")
                if len(fields) != 10:
                    raise DeepVariantError(
                        "DeepVariant caller VCF must contain one sample column"
                    )
                contig, position_text, _id, ref, alt, _qual, filter_value = fields[:7]
                try:
                    position = int(position_text)
                except ValueError as error:
                    raise DeepVariantError("DeepVariant VCF position is invalid") from error
                if (
                    contig not in contig_lengths
                    or contig not in reference_index
                    or contig_lengths[contig] != reference_index[contig].length
                    or position <= 0
                    or position + len(ref) - 1 > contig_lengths[contig]
                    or not SEQUENCE.fullmatch(ref)
                    or not alt
                ):
                    raise DeepVariantError(
                        f"DeepVariant VCF record or contig is invalid: {contig}"
                    )
                format_keys = fields[8].split(":")
                sample_values = fields[9].split(":")
                complete = _complete_genotype(format_keys, sample_values)
                record_count += 1
                sentinel = _Sentinel(contig, position, ref, filter_value)
                if filter_value == "PASS":
                    pass_count += 1
                    if pass_sentinel is None and contig in AUTOSOMES and complete:
                        pass_sentinel = sentinel
                elif filter_value != ".":
                    filtered_count += 1
                    if filtered_sentinel is None and contig in AUTOSOMES and complete:
                        filtered_sentinel = sentinel
    except (OSError, EOFError, gzip.BadGzipFile, UnicodeDecodeError) as error:
        raise DeepVariantError(f"cannot read DeepVariant caller VCF: {error}") from error
    if samples != (contract.sample,):
        raise DeepVariantError(
            f"DeepVariant caller VCF sample is not exactly {contract.sample}"
        )
    if not version_seen:
        raise DeepVariantError("DeepVariant caller version header is missing")
    if not {"GT", "GQ", "DP", "AD"}.issubset(format_definitions):
        raise DeepVariantError("DeepVariant caller VCF lacks required genotype fields")
    for number in range(1, 23):
        contig = f"chr{number}"
        expected_length = reference_contract.canonical_contig_lengths[str(number)]
        if contig_lengths.get(contig) != expected_length:
            raise DeepVariantError(f"DeepVariant contig length mismatch: {contig}")
    if (
        record_count == 0
        or pass_count == 0
        or filtered_count == 0
        or pass_sentinel is None
        or filtered_sentinel is None
    ):
        raise DeepVariantError(
            "DeepVariant caller VCF lacks queryable PASS or filtered genotypes"
        )
    reference = IndexedFasta(canonical_fasta, canonical_fai)
    for sentinel in (pass_sentinel, filtered_sentinel):
        if reference.fetch(sentinel.contig, sentinel.position, len(sentinel.ref)) != sentinel.ref:
            raise DeepVariantError(
                f"DeepVariant caller REF differs from canonical reference: {sentinel.contig}"
            )
    return (
        DeepVariantValidation(record_count, pass_count, filtered_count, 0),
        pass_sentinel,
        filtered_sentinel,
    )


def _record_matches(record: str, sentinel: _Sentinel) -> bool:
    fields = record.split("\t")
    return (
        len(fields) >= 7
        and fields[0] == sentinel.contig
        and fields[1] == str(sentinel.position)
        and fields[3] == sentinel.ref
        and fields[6] == sentinel.filter_value
    )


def validate_deepvariant_bundle(
    contract: DeepVariantContract,
    reference_contract: ReferenceIdentityContract,
    reference_root: str | Path,
) -> DeepVariantValidation:
    root = Path(reference_root).resolve()
    try:
        reference_role = contract.reference_binding.resolve(
            reference_contract, required_role="public_annotation_benchmark"
        )
    except ReferenceConcordanceError as error:
        raise DeepVariantError(str(error)) from error
    for artifact in contract.artifacts:
        _verify_artifact(root / artifact.path, artifact)
    _validate_execution_log(root / contract.artifact("execution_log").path, contract)
    vcf_path = root / contract.artifact("caller_vcf").path
    validation, pass_sentinel, filtered_sentinel = _scan_vcf(
        vcf_path, contract, reference_contract, reference_role, root
    )
    try:
        index = TabixIndex.load(root / contract.artifact("caller_tbi").path)
        index.require_vcf_layout()
        for sentinel in (pass_sentinel, filtered_sentinel):
            records = index.query(
                vcf_path, sentinel.contig, sentinel.position, sentinel.position
            )
            if not any(_record_matches(record, sentinel) for record in records):
                raise DeepVariantError(
                    f"DeepVariant index did not recover {sentinel.filter_value} sentinel"
                )
    except (OSError, KeyError, TabixError) as error:
        raise DeepVariantError(f"DeepVariant tabix query failed: {error}") from error
    indexed_names = {reference.name for reference in index.references}
    if not AUTOSOMES.issubset(indexed_names):
        raise DeepVariantError("DeepVariant tabix index lacks one or more autosomes")
    return DeepVariantValidation(
        variant_records=validation.variant_records,
        pass_records=validation.pass_records,
        filtered_records=validation.filtered_records,
        indexed_contigs=len(index.references),
    )


def _manifest_path(contract: DeepVariantContract) -> str:
    return f"manifests/deepvariant/{contract.resource_id}.json"


def _load_manifest(path: Path) -> PublicResourceManifest:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise DeepVariantError(f"cannot read DeepVariant manifest: {error}") from error
    if not isinstance(data, Mapping):
        raise DeepVariantError("DeepVariant manifest root must be an object")
    try:
        return PublicResourceManifest.from_mapping(data)
    except ManifestError as error:
        raise DeepVariantError(f"invalid DeepVariant manifest: {error}") from error


def _verify_manifest_contract(
    manifest: PublicResourceManifest, contract: DeepVariantContract
) -> None:
    caller_vcf = contract.artifact("caller_vcf")
    if (
        manifest.resource_id != contract.resource_id
        or manifest.release != contract.release
        or manifest.assembly != contract.assembly
        or manifest.license != contract.license
        or manifest.source.url != caller_vcf.url
        or manifest.source.object_generation != caller_vcf.object_generation
        or {output.path for output in manifest.outputs}
        != {artifact.path for artifact in contract.artifacts}
    ):
        raise DeepVariantError("DeepVariant manifest identity differs from its contract")


def _copy_file(source: Path, destination: Path) -> None:
    try:
        with source.open("rb") as input_file, destination.open("xb") as output_file:
            shutil.copyfileobj(input_file, output_file, length=1024 * 1024)
            output_file.flush()
            os.fsync(output_file.fileno())
    except OSError as error:
        raise DeepVariantError(f"cannot stage DeepVariant artifact: {error}") from error


def _outputs(
    root: Path, contract: DeepVariantContract
) -> tuple[ResourceOutput, ...]:
    outputs = []
    for artifact in contract.artifacts:
        sha256, _md5, size = _checksums(root / artifact.path)
        outputs.append(
            ResourceOutput(
                path=artifact.path,
                size_bytes=size,
                sha256=sha256,
                retention=RetentionPolicy("durable", None),
            )
        )
    return tuple(outputs)


def bootstrap_deepvariant_hg002(
    *,
    contract_path: str | Path,
    reference_identity_path: str | Path,
    reference_root: str | Path,
    image_digest: str,
    retrieved_on: str | None = None,
    downloader: Downloader = download_atomic,
) -> DeepVariantBootstrapResult:
    contract_file = Path(contract_path)
    identity_file = Path(reference_identity_path)
    root = Path(reference_root).resolve(strict=True)
    contract = load_deepvariant_contract(contract_file)
    try:
        reference_contract = load_reference_identity(identity_file)
    except ReferenceConcordanceError as error:
        raise DeepVariantError(str(error)) from error
    manifest_path = root / _manifest_path(contract)
    if manifest_path.exists():
        manifest = _load_manifest(manifest_path)
        _verify_manifest_contract(manifest, contract)
        try:
            verify_outputs(manifest, root)
        except ManifestError as error:
            raise DeepVariantError(
                f"completed DeepVariant bundle is invalid: {error}"
            ) from error
        validation = validate_deepvariant_bundle(contract, reference_contract, root)
        return DeepVariantBootstrapResult(
            str(manifest_path),
            _checksums(manifest_path)[0],
            contract.sample,
            contract.release,
            validation,
            True,
        )
    total_size = sum(artifact.size_bytes for artifact in contract.artifacts)
    require_capacity(
        root,
        CapacityPlan(
            operation="bootstrap DeepVariant HG002 NovaSeq 35x caller VCF",
            source_bytes=total_size,
            output_bytes=total_size,
        ),
    )
    target_relative = PurePosixPath(contract.artifacts[0].path).parent
    target = root / str(target_relative)
    scratch = root / "scratch/public-build/deepvariant/HG002/1.0.0/GRCh38"
    if not target.exists():
        scratch.mkdir(parents=True, exist_ok=True, mode=0o700)
        downloaded: dict[str, Path] = {}
        for artifact in contract.artifacts:
            destination = scratch / Path(urlparse(artifact.url).path).name
            downloader(
                DownloadRequest(
                    url=artifact.url,
                    destination=destination,
                    expected_size_bytes=artifact.size_bytes,
                    expected_md5=artifact.md5,
                    expected_generation=artifact.object_generation,
                    attempts=3,
                    timeout_seconds=120,
                ),
                root,
            )
            downloaded[artifact.role] = destination
        target.parent.mkdir(parents=True, exist_ok=True)
        partial_root = target.parent / ".partial"
        partial_root.mkdir(mode=0o700, exist_ok=True)
        staging = partial_root / f"{target.name}.{uuid.uuid4().hex}.partial"
        staging.mkdir(mode=0o700)
        try:
            for artifact in contract.artifacts:
                staged = staging / PurePosixPath(artifact.path).name
                _copy_file(downloaded[artifact.role], staged)
                _verify_artifact(staged, artifact)
                os.chmod(staged, 0o444)
            os.replace(staging, target)
            directory_fd = os.open(target.parent, os.O_RDONLY)
            try:
                os.fsync(directory_fd)
            finally:
                os.close(directory_fd)
        except Exception:
            shutil.rmtree(staging, ignore_errors=True)
            raise
    elif not target.is_dir():
        raise DeepVariantError(
            f"DeepVariant release path is not a directory: {target}"
        )
    validation = validate_deepvariant_bundle(contract, reference_contract, root)
    caller_vcf = contract.artifact("caller_vcf")
    vcf_sha256 = _checksums(root / caller_vcf.path)[0]
    manifest = PublicResourceManifest(
        resource_id=contract.resource_id,
        release=contract.release,
        source=SourceIdentity(
            url=caller_vcf.url,
            sha256=vcf_sha256,
            object_generation=caller_vcf.object_generation,
        ),
        license=contract.license,
        assembly=contract.assembly,
        retrieved_on=retrieved_on or date.today().isoformat(),
        transformation=Transformation(
            schema=TRANSFORMATION_SCHEMA,
            parameters_sha256=_checksums(contract_file)[0],
        ),
        image_digest=image_digest,
        outputs=_outputs(root, contract),
    )
    try:
        publication = publish_manifest(manifest, manifest_path, root)
    except ManifestError as error:
        raise DeepVariantError(f"cannot publish DeepVariant manifest: {error}") from error
    if scratch.exists():
        for artifact in contract.artifacts:
            (scratch / Path(urlparse(artifact.url).path).name).unlink(missing_ok=True)
    return DeepVariantBootstrapResult(
        str(publication.path),
        publication.manifest_sha256,
        contract.sample,
        contract.release,
        validation,
        publication.reused,
    )


def deepvariant_result_as_json(result: DeepVariantBootstrapResult) -> str:
    return json.dumps(asdict(result), indent=2, sort_keys=True) + "\n"
