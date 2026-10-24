"""Bootstrap and validate the public GIAB HG002 NIST v4.2.1 bundle."""

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
import struct
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


GIAB_SCHEMA = "genosapiens.giab-hg002"
GIAB_TRANSFORMATION_SCHEMA = "genosapiens.giab-hg002-bootstrap"
MD5 = re.compile(r"^[0-9a-f]{32}$")
SEQUENCE = re.compile(rb"^[ACGTRYSWKMBDHVN]+$")
CONTIG_HEADER = re.compile(r"^##contig=<ID=([^,>]+),length=([0-9]+)(?:,.*)?>$")
EXPECTED_ROLES = {
    "truth_vcf",
    "truth_tbi",
    "confident_bed",
    "reference_fasta_bgzip",
    "reference_fai",
    "reference_gzi",
}


class GiabError(RuntimeError):
    """Raised when GIAB acquisition or validation fails."""


@dataclass(frozen=True, slots=True)
class GiabArtifact:
    role: str
    url: str
    path: str
    size_bytes: int
    md5: str

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "GiabArtifact":
        allowed = {"role", "url", "path", "size_bytes", "md5"}
        unknown = sorted(set(data) - allowed)
        if unknown:
            raise GiabError(f"unknown GIAB artifact fields: {', '.join(unknown)}")
        try:
            artifact = cls(**data)
        except TypeError as error:
            raise GiabError(f"invalid GIAB artifact: {error}") from error
        if (
            not artifact.role
            or not artifact.url.startswith("https://")
            or not isinstance(artifact.size_bytes, int)
            or isinstance(artifact.size_bytes, bool)
            or artifact.size_bytes <= 0
            or not MD5.fullmatch(artifact.md5)
        ):
            raise GiabError(f"invalid GIAB artifact identity: {artifact.role}")
        path = PurePosixPath(artifact.path)
        if path.is_absolute() or ".." in path.parts or path == PurePosixPath("."):
            raise GiabError("GIAB artifact path must be safe and relative")
        return artifact


@dataclass(frozen=True, slots=True)
class GiabResourceContract:
    resource_id: str
    release: str
    sample: str
    assembly: str
    reference_binding: ReferenceRoleBinding
    published_checksum_source: str
    autosomal_contigs: tuple[str, ...]
    artifacts: tuple[GiabArtifact, ...]
    schema: str = GIAB_SCHEMA

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "GiabResourceContract":
        allowed = {
            "schema",
            "resource_id",
            "release",
            "sample",
            "assembly",
            "reference_binding",
            "published_checksum_source",
            "autosomal_contigs",
            "artifacts",
        }
        unknown = sorted(set(data) - allowed)
        if unknown:
            raise GiabError(f"unknown GIAB contract fields: {', '.join(unknown)}")
        if data.get("schema") != GIAB_SCHEMA:
            raise GiabError(f"GIAB schema must be {GIAB_SCHEMA!r}")
        contigs = data.get("autosomal_contigs")
        artifact_data = data.get("artifacts")
        reference_binding = data.get("reference_binding")
        if not isinstance(contigs, list) or not all(
            isinstance(contig, str) and contig for contig in contigs
        ):
            raise GiabError("GIAB autosomal_contigs must be a string array")
        if not isinstance(reference_binding, Mapping):
            raise GiabError("GIAB reference_binding must be an object")
        if not isinstance(artifact_data, list) or not all(
            isinstance(item, Mapping) for item in artifact_data
        ):
            raise GiabError("GIAB artifacts must be an object array")
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
                published_checksum_source=data["published_checksum_source"],
                autosomal_contigs=tuple(contigs),
                artifacts=tuple(
                    GiabArtifact.from_mapping(item) for item in artifact_data
                ),
            )
        except (KeyError, TypeError, ReferenceConcordanceError) as error:
            raise GiabError(f"invalid GIAB contract: {error}") from error
        for value in (
            contract.resource_id,
            contract.release,
            contract.sample,
            contract.assembly,
            contract.published_checksum_source,
        ):
            if not isinstance(value, str) or not value:
                raise GiabError("GIAB contract string fields must be non-empty")
        expected_contigs = tuple(f"chr{number}" for number in range(1, 23))
        if contract.autosomal_contigs != expected_contigs:
            raise GiabError("GIAB contract must declare ordered chr1-chr22 scope")
        roles = {artifact.role for artifact in contract.artifacts}
        if roles != EXPECTED_ROLES or len(contract.artifacts) != len(roles):
            raise GiabError("GIAB contract must define each required artifact once")
        parents = {str(PurePosixPath(item.path).parent) for item in contract.artifacts}
        if len(parents) != 1:
            raise GiabError("GIAB artifacts must share one release directory")
        return contract

    def artifact(self, role: str) -> GiabArtifact:
        for artifact in self.artifacts:
            if artifact.role == role:
                return artifact
        raise GiabError(f"GIAB contract lacks artifact role {role}")


def load_giab_contract(path: str | Path) -> GiabResourceContract:
    contract_path = Path(path)
    try:
        data = json.loads(contract_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise GiabError(f"cannot read GIAB contract {contract_path}: {error}") from error
    if not isinstance(data, Mapping):
        raise GiabError("GIAB contract root must be an object")
    return GiabResourceContract.from_mapping(data)


@dataclass(frozen=True, slots=True)
class GiabValidation:
    variant_records: int
    confident_intervals: int
    reference_sequences: int
    indexed_contigs: int


@dataclass(frozen=True, slots=True)
class GiabBootstrapResult:
    manifest_path: str
    manifest_sha256: str
    sample: str
    release: str
    validation: GiabValidation
    reused: bool


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
        raise GiabError(f"cannot read {path}: {error}") from error
    return sha256.hexdigest(), md5.hexdigest(), size


def _verify_artifact(path: Path, artifact: GiabArtifact) -> str:
    sha256, md5, size = _checksums(path)
    if size != artifact.size_bytes or md5 != artifact.md5:
        raise GiabError(f"GIAB published checksum mismatch: {artifact.role}")
    return sha256


def _load_sequence_dictionary(path: Path) -> Mapping[str, tuple[int, str]]:
    entries: dict[str, tuple[int, str]] = {}
    try:
        with path.open("rt", encoding="ascii") as source:
            for line in source:
                if not line.startswith("@SQ\t"):
                    continue
                fields = dict(
                    field.split(":", 1) for field in line.rstrip("\n").split("\t")[1:]
                )
                name = fields["SN"]
                if name in entries:
                    raise GiabError(f"duplicate dictionary sequence: {name}")
                entries[name] = (int(fields["LN"]), fields["M5"])
    except (KeyError, OSError, ValueError) as error:
        raise GiabError(f"invalid sequence dictionary {path}: {error}") from error
    if not entries:
        raise GiabError(f"sequence dictionary has no entries: {path}")
    return entries


def _header_name_and_md5(line: bytes) -> tuple[str, str]:
    try:
        tokens = line[1:].decode("ascii").split()
    except UnicodeDecodeError as error:
        raise GiabError("GIAB reference header is not ASCII") from error
    if not tokens:
        raise GiabError("GIAB reference header has no sequence name")
    md5_values = [token[3:] for token in tokens[1:] if token.startswith("M5:")]
    if len(md5_values) != 1 or not MD5.fullmatch(md5_values[0]):
        raise GiabError(f"GIAB reference header lacks M5: {tokens[0]}")
    return tokens[0], md5_values[0]


def _validate_matching_reference(
    fasta_path: Path,
    fai_path: Path,
    gzi_path: Path,
    reference_contract: ReferenceIdentityContract,
    reference_role: ReferenceRole,
    reference_root: Path,
) -> int:
    dictionary_path = reference_root / reference_contract.expected_output(
        reference_role.dictionary_output_role
    ).path
    expected_dictionary = _load_sequence_dictionary(dictionary_path)
    giab_index = load_fasta_index(fai_path)
    observed: dict[str, tuple[int, str]] = {}
    current_name: str | None = None
    expected_md5: str | None = None
    digest = hashlib.md5(usedforsecurity=False)
    length = 0

    def finish() -> None:
        nonlocal current_name, expected_md5, digest, length
        if current_name is None or expected_md5 is None:
            return
        actual_md5 = digest.hexdigest()
        if actual_md5 != expected_md5:
            raise GiabError(f"GIAB reference sequence M5 mismatch: {current_name}")
        if current_name not in giab_index or giab_index[current_name].length != length:
            raise GiabError(f"GIAB reference FAI mismatch: {current_name}")
        if expected_dictionary.get(current_name) != (length, actual_md5):
            raise GiabError(
                f"GIAB reference differs from pinned canonical reference: {current_name}"
            )
        observed[current_name] = (length, actual_md5)

    try:
        with gzip.open(fasta_path, "rb") as source:
            while line := source.readline():
                if line.startswith(b">"):
                    finish()
                    current_name, expected_md5 = _header_name_and_md5(line)
                    if current_name in observed:
                        raise GiabError(f"duplicate GIAB reference sequence: {current_name}")
                    digest = hashlib.md5(usedforsecurity=False)
                    length = 0
                    continue
                sequence = line.rstrip(b"\r\n")
                if current_name is None or not sequence or not SEQUENCE.fullmatch(sequence):
                    raise GiabError("GIAB reference contains invalid sequence data")
                digest.update(sequence)
                length += len(sequence)
            finish()
    except (OSError, EOFError, gzip.BadGzipFile) as error:
        raise GiabError(f"cannot read GIAB reference FASTA: {error}") from error
    if set(observed) != set(giab_index):
        raise GiabError("GIAB reference FASTA and FAI names differ")
    for contig, expected_length in reference_contract.canonical_contig_lengths.items():
        if contig == "M":
            canonical_name = "chrM"
        else:
            canonical_name = f"chr{contig}"
        if canonical_name in {f"chr{number}" for number in range(1, 23)} and (
            observed.get(canonical_name, (None,))[0] != expected_length
        ):
            raise GiabError(f"GIAB canonical contig length mismatch: {canonical_name}")
    try:
        data = gzi_path.read_bytes()
    except OSError as error:
        raise GiabError(f"cannot read GIAB GZI: {error}") from error
    if len(data) < 8:
        raise GiabError("GIAB GZI is truncated")
    entry_count = struct.unpack("<Q", data[:8])[0]
    if len(data) != 8 + entry_count * 16:
        raise GiabError("GIAB GZI size does not match its entry count")
    previous_compressed = -1
    previous_uncompressed = -1
    for offset in range(8, len(data), 16):
        compressed, uncompressed = struct.unpack("<QQ", data[offset : offset + 16])
        if compressed <= previous_compressed or uncompressed <= previous_uncompressed:
            raise GiabError("GIAB GZI offsets are not strictly increasing")
        previous_compressed = compressed
        previous_uncompressed = uncompressed
    return len(observed)


def _validate_vcf(
    path: Path,
    contract: GiabResourceContract,
    reference_contract: ReferenceIdentityContract,
    reference_role: ReferenceRole,
    reference_root: Path,
) -> tuple[int, Mapping[str, tuple[int, str]]]:
    expected_contigs = set(contract.autosomal_contigs)
    header_lengths: dict[str, int] = {}
    samples: tuple[str, ...] | None = None
    first_records: dict[str, tuple[int, str]] = {}
    record_count = 0
    try:
        with gzip.open(path, "rt", encoding="utf-8") as source:
            for line in source:
                if line.startswith("##contig="):
                    match = CONTIG_HEADER.match(line.rstrip("\n"))
                    if match:
                        header_lengths[match.group(1)] = int(match.group(2))
                    continue
                if line.startswith("#CHROM\t"):
                    samples = tuple(line.rstrip("\n").split("\t")[9:])
                    continue
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 10:
                    raise GiabError("GIAB truth VCF record has fewer than 10 columns")
                contig = fields[0]
                if contig not in expected_contigs:
                    raise GiabError(f"GIAB truth VCF is not autosomal-only: {contig}")
                try:
                    position = int(fields[1])
                except ValueError as error:
                    raise GiabError("GIAB truth VCF position is invalid") from error
                ref = fields[3]
                if position <= 0 or not ref or not re.fullmatch(r"[ACGTRYSWKMBDHVN]+", ref):
                    raise GiabError("GIAB truth VCF REF is invalid")
                first_records.setdefault(contig, (position, ref))
                record_count += 1
    except (OSError, EOFError, gzip.BadGzipFile, UnicodeDecodeError) as error:
        raise GiabError(f"cannot read GIAB truth VCF: {error}") from error
    if samples != (contract.sample,):
        raise GiabError(f"GIAB truth VCF sample is not exactly {contract.sample}")
    if set(first_records) != expected_contigs or record_count == 0:
        raise GiabError("GIAB truth VCF does not cover every autosome")
    for number in range(1, 23):
        contig = f"chr{number}"
        expected_length = reference_contract.canonical_contig_lengths[str(number)]
        if header_lengths.get(contig) != expected_length:
            raise GiabError(f"GIAB VCF header length mismatch: {contig}")
    canonical_fasta = reference_root / reference_contract.expected_output(
        reference_role.fasta_output_role
    ).path
    canonical_fai = reference_root / reference_contract.expected_output(
        reference_role.fai_output_role
    ).path
    reference = IndexedFasta(canonical_fasta, canonical_fai)
    for contig, (position, ref) in first_records.items():
        if reference.fetch(contig, position, len(ref)) != ref:
            raise GiabError(f"GIAB truth REF mismatch: {contig}")
    return record_count, first_records


def _validate_bed(
    path: Path,
    contract: GiabResourceContract,
    reference_contract: ReferenceIdentityContract,
) -> int:
    contig_order = {name: index for index, name in enumerate(contract.autosomal_contigs)}
    seen: set[str] = set()
    previous_order = -1
    previous_start = -1
    intervals = 0
    try:
        with path.open("rt", encoding="utf-8") as source:
            for line in source:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3 or fields[0] not in contig_order:
                    raise GiabError("GIAB confident BED is not autosomal-only")
                try:
                    start, end = int(fields[1]), int(fields[2])
                except ValueError as error:
                    raise GiabError("GIAB confident BED coordinate is invalid") from error
                order = contig_order[fields[0]]
                length = reference_contract.canonical_contig_lengths[fields[0][3:]]
                if start < 0 or end <= start or end > length:
                    raise GiabError("GIAB confident BED interval is out of range")
                if order < previous_order or (order == previous_order and start < previous_start):
                    raise GiabError("GIAB confident BED is not coordinate-sorted")
                if order != previous_order:
                    previous_start = -1
                previous_order = order
                previous_start = start
                seen.add(fields[0])
                intervals += 1
    except (OSError, UnicodeDecodeError) as error:
        raise GiabError(f"cannot read GIAB confident BED: {error}") from error
    if seen != set(contract.autosomal_contigs) or intervals == 0:
        raise GiabError("GIAB confident BED does not cover every autosome")
    return intervals


def validate_giab_bundle(
    contract: GiabResourceContract,
    reference_contract: ReferenceIdentityContract,
    reference_root: str | Path,
) -> GiabValidation:
    root = Path(reference_root).resolve()
    try:
        reference_role = contract.reference_binding.resolve(
            reference_contract, required_role="public_annotation_benchmark"
        )
    except ReferenceConcordanceError as error:
        raise GiabError(str(error)) from error
    for artifact in contract.artifacts:
        _verify_artifact(root / artifact.path, artifact)
    reference_sequences = _validate_matching_reference(
        root / contract.artifact("reference_fasta_bgzip").path,
        root / contract.artifact("reference_fai").path,
        root / contract.artifact("reference_gzi").path,
        reference_contract,
        reference_role,
        root,
    )
    variant_count, first_records = _validate_vcf(
        root / contract.artifact("truth_vcf").path,
        contract,
        reference_contract,
        reference_role,
        root,
    )
    intervals = _validate_bed(
        root / contract.artifact("confident_bed").path,
        contract,
        reference_contract,
    )
    try:
        index = TabixIndex.load(root / contract.artifact("truth_tbi").path)
        index.require_vcf_layout()
    except TabixError as error:
        raise GiabError(f"GIAB tabix index is invalid: {error}") from error
    if tuple(reference.name for reference in index.references) != contract.autosomal_contigs:
        raise GiabError("GIAB tabix index is not ordered chr1-chr22")
    vcf_path = root / contract.artifact("truth_vcf").path
    for contig, (position, ref) in first_records.items():
        try:
            records = index.query(vcf_path, contig, position, position)
        except (OSError, TabixError) as error:
            raise GiabError(f"GIAB tabix query failed for {contig}: {error}") from error
        if not any(
            record.split("\t", 4)[:4] == [contig, str(position), ".", ref]
            or (
                record.split("\t", 4)[0] == contig
                and record.split("\t", 4)[1] == str(position)
                and record.split("\t", 4)[3] == ref
            )
            for record in records
        ):
            raise GiabError(f"GIAB tabix query did not recover sentinel: {contig}")
    return GiabValidation(
        variant_records=variant_count,
        confident_intervals=intervals,
        reference_sequences=reference_sequences,
        indexed_contigs=len(index.references),
    )


def _manifest_path(contract: GiabResourceContract) -> str:
    return f"manifests/giab/{contract.resource_id}.json"


def _load_manifest(path: Path) -> PublicResourceManifest:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise GiabError(f"cannot read GIAB manifest {path}: {error}") from error
    if not isinstance(data, Mapping):
        raise GiabError("GIAB manifest root must be an object")
    try:
        return PublicResourceManifest.from_mapping(data)
    except ManifestError as error:
        raise GiabError(f"invalid GIAB manifest: {error}") from error


def _verify_manifest_contract(
    manifest: PublicResourceManifest, contract: GiabResourceContract
) -> None:
    if (
        manifest.resource_id != contract.resource_id
        or manifest.release != contract.release
        or manifest.assembly != contract.assembly
        or manifest.source.url != contract.artifact("truth_vcf").url
    ):
        raise GiabError("GIAB manifest identity differs from its contract")
    if {output.path for output in manifest.outputs} != {
        artifact.path for artifact in contract.artifacts
    }:
        raise GiabError("GIAB manifest output paths differ from its contract")


def _copy_file(source: Path, destination: Path) -> None:
    try:
        with source.open("rb") as input_file, destination.open("xb") as output_file:
            shutil.copyfileobj(input_file, output_file, length=1024 * 1024)
            output_file.flush()
            os.fsync(output_file.fileno())
    except OSError as error:
        raise GiabError(f"cannot stage GIAB artifact {source}: {error}") from error


def _outputs(root: Path, contract: GiabResourceContract) -> tuple[ResourceOutput, ...]:
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


def bootstrap_giab_hg002(
    *,
    contract_path: str | Path,
    reference_identity_path: str | Path,
    reference_root: str | Path,
    image_digest: str,
    retrieved_on: str | None = None,
    downloader: Downloader = download_atomic,
) -> GiabBootstrapResult:
    resource_file = Path(contract_path)
    identity_file = Path(reference_identity_path)
    root = Path(reference_root).resolve(strict=True)
    contract = load_giab_contract(resource_file)
    try:
        reference_contract = load_reference_identity(identity_file)
    except ReferenceConcordanceError as error:
        raise GiabError(str(error)) from error
    manifest_path = root / _manifest_path(contract)
    if manifest_path.exists():
        manifest = _load_manifest(manifest_path)
        _verify_manifest_contract(manifest, contract)
        try:
            verify_outputs(manifest, root)
        except ManifestError as error:
            raise GiabError(f"completed GIAB bundle is invalid: {error}") from error
        validation = validate_giab_bundle(contract, reference_contract, root)
        return GiabBootstrapResult(
            manifest_path=str(manifest_path),
            manifest_sha256=_checksums(manifest_path)[0],
            sample=contract.sample,
            release=contract.release,
            validation=validation,
            reused=True,
        )
    total_size = sum(artifact.size_bytes for artifact in contract.artifacts)
    require_capacity(
        root,
        CapacityPlan(
            operation="bootstrap GIAB HG002 NIST v4.2.1",
            source_bytes=total_size,
            output_bytes=total_size,
        ),
    )
    target_relative = PurePosixPath(contract.artifacts[0].path).parent
    target = root / str(target_relative)
    if not target.exists():
        scratch = root / "scratch/public-build/giab/HG002/NISTv4.2.1/GRCh38"
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
        raise GiabError(f"GIAB release path is not a directory: {target}")
    validation = validate_giab_bundle(contract, reference_contract, root)
    truth = contract.artifact("truth_vcf")
    truth_sha256 = _checksums(root / truth.path)[0]
    contract_sha256 = _checksums(resource_file)[0]
    manifest = PublicResourceManifest(
        resource_id=contract.resource_id,
        release=contract.release,
        source=SourceIdentity(url=truth.url, sha256=truth_sha256),
        license="CC0-1.0",
        assembly=contract.assembly,
        retrieved_on=retrieved_on or date.today().isoformat(),
        transformation=Transformation(
            schema=GIAB_TRANSFORMATION_SCHEMA,
            parameters_sha256=contract_sha256,
        ),
        image_digest=image_digest,
        outputs=_outputs(root, contract),
    )
    publication = publish_manifest(manifest, manifest_path, root)
    scratch = root / "scratch/public-build/giab/HG002/NISTv4.2.1/GRCh38"
    if scratch.exists():
        for artifact in contract.artifacts:
            (scratch / Path(urlparse(artifact.url).path).name).unlink(missing_ok=True)
    return GiabBootstrapResult(
        manifest_path=str(publication.path),
        manifest_sha256=publication.manifest_sha256,
        sample=contract.sample,
        release=contract.release,
        validation=validation,
        reused=publication.reused,
    )


def giab_result_as_json(result: GiabBootstrapResult) -> str:
    return json.dumps(asdict(result), indent=2, sort_keys=True) + "\n"
