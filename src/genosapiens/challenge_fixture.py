"""Build a challenge-shaped fixture from public HG002 alleles only."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict, dataclass, replace
from datetime import date
import binascii
import gzip
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import struct
from typing import Any
import uuid
import zlib

from genosapiens.giab import GiabResourceContract, load_giab_contract
from genosapiens.governance import (
    ArtifactClassification,
    GovernanceError,
    validate_relative_classification,
)
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
from genosapiens.storage import CapacityPlan, require_capacity
from genosapiens.tabix import TabixError, TabixIndex


SCHEMA = "genosapiens.challenge-shaped-fixture"
PROVENANCE_SCHEMA = "genosapiens.fixture-allele-provenance"
TRANSFORMATION_SCHEMA = "genosapiens.challenge-shaped-fixture-build"
MD5 = re.compile(r"^[0-9a-f]{32}$")
SEQUENCE = re.compile(r"^[ACGTRYSWKMBDHVN]+$")
SELECTION_POLICY = (
    "first-autosomal-biallelic-snv",
    "second-autosomal-biallelic-snv",
    "first-autosomal-biallelic-insertion",
    "first-autosomal-biallelic-deletion",
    "first-autosomal-multiallelic-small-variant",
)
OUTPUT_ROLES = {"fixture_vcf", "fixture_tbi", "provenance"}
AUTOSOMES = tuple(f"chr{number}" for number in range(1, 23))
AUTOSOME_ORDER = {name: index for index, name in enumerate(AUTOSOMES)}


class ChallengeFixtureError(RuntimeError):
    """Raised when public fixture generation or validation fails."""


@dataclass(frozen=True, slots=True)
class ChallengeFixtureContract:
    resource_id: str
    release: str
    sample: str
    assembly: str
    reference_binding: ReferenceRoleBinding
    source_resource_id: str
    source_release: str
    source_vcf_role: str
    source_vcf_md5: str
    source_manifest_path: str
    selection_policy: tuple[str, ...]
    structural_basis: Mapping[str, Any]
    outputs: Mapping[str, str]
    schema: str = SCHEMA

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "ChallengeFixtureContract":
        allowed = {
            "schema",
            "resource_id",
            "release",
            "sample",
            "assembly",
            "reference_binding",
            "source_resource_id",
            "source_release",
            "source_vcf_role",
            "source_vcf_md5",
            "source_manifest_path",
            "selection_policy",
            "structural_basis",
            "outputs",
        }
        if set(data) != allowed or data.get("schema") != SCHEMA:
            raise ChallengeFixtureError("challenge fixture contract fields are invalid")
        policy = data.get("selection_policy")
        basis = data.get("structural_basis")
        outputs = data.get("outputs")
        reference_binding = data.get("reference_binding")
        if (
            not isinstance(policy, list)
            or not isinstance(basis, Mapping)
            or set(basis) != {"source", "facts"}
            or not isinstance(basis.get("source"), str)
            or not isinstance(basis.get("facts"), list)
            or not isinstance(outputs, Mapping)
            or set(outputs) != OUTPUT_ROLES
            or not isinstance(reference_binding, Mapping)
        ):
            raise ChallengeFixtureError("challenge fixture nested fields are invalid")
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
                source_resource_id=data["source_resource_id"],
                source_release=data["source_release"],
                source_vcf_role=data["source_vcf_role"],
                source_vcf_md5=data["source_vcf_md5"],
                source_manifest_path=data["source_manifest_path"],
                selection_policy=tuple(policy),
                structural_basis=dict(basis),
                outputs=dict(outputs),
            )
        except (KeyError, TypeError, ReferenceConcordanceError) as error:
            raise ChallengeFixtureError(f"invalid challenge fixture contract: {error}") from error
        strings = (
            contract.resource_id,
            contract.release,
            contract.sample,
            contract.assembly,
            contract.source_resource_id,
            contract.source_release,
            contract.source_vcf_role,
            contract.source_manifest_path,
        )
        if not all(isinstance(value, str) and value for value in strings):
            raise ChallengeFixtureError("challenge fixture strings must be non-empty")
        if (
            contract.sample != "SYNTHETIC"
            or not MD5.fullmatch(contract.source_vcf_md5)
            or contract.selection_policy != SELECTION_POLICY
        ):
            raise ChallengeFixtureError("challenge fixture identity or policy is invalid")
        parents = set()
        for output in contract.outputs.values():
            if not isinstance(output, str):
                raise ChallengeFixtureError("challenge fixture output paths must be strings")
            try:
                validate_relative_classification(output, ArtifactClassification.PUBLIC)
            except GovernanceError as error:
                raise ChallengeFixtureError(str(error)) from error
            parents.add(str(PurePosixPath(output).parent))
        if len(parents) != 1:
            raise ChallengeFixtureError("challenge fixture outputs must share one directory")
        return contract


def load_challenge_fixture_contract(path: str | Path) -> ChallengeFixtureContract:
    contract_path = Path(path)
    try:
        data = json.loads(contract_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ChallengeFixtureError(
            f"cannot read challenge fixture contract {contract_path}: {error}"
        ) from error
    if not isinstance(data, Mapping):
        raise ChallengeFixtureError("challenge fixture contract root must be an object")
    return ChallengeFixtureContract.from_mapping(data)


@dataclass(frozen=True, slots=True)
class ChallengeFixtureValidation:
    records: int
    public_alleles: int
    phased_records: int
    unphased_records: int
    filtered_records: int
    missing_quality_records: int
    indexed_contigs: int


@dataclass(frozen=True, slots=True)
class ChallengeFixtureResult:
    manifest_path: str
    manifest_sha256: str
    fixture_vcf: str
    provenance: str
    validation: ChallengeFixtureValidation
    reused: bool


@dataclass(frozen=True, slots=True)
class _PublicRecord:
    profile: str
    contig: str
    position: int
    ref: str
    alts: tuple[str, ...]
    source_record_sha256: str


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
        raise ChallengeFixtureError(f"cannot read {path}: {error}") from error
    return sha256.hexdigest(), md5.hexdigest(), size


def _load_manifest(path: Path) -> PublicResourceManifest:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ChallengeFixtureError(f"cannot read public manifest {path}: {error}") from error
    if not isinstance(data, Mapping):
        raise ChallengeFixtureError("public manifest root must be an object")
    try:
        return PublicResourceManifest.from_mapping(data)
    except ManifestError as error:
        raise ChallengeFixtureError(f"invalid public manifest: {error}") from error


def _source_identity(
    contract: ChallengeFixtureContract,
    giab_contract: GiabResourceContract,
    root: Path,
) -> tuple[Path, str]:
    if (
        contract.source_resource_id != giab_contract.resource_id
        or contract.source_release != giab_contract.release
    ):
        raise ChallengeFixtureError("fixture source does not match the GIAB contract")
    source = giab_contract.artifact(contract.source_vcf_role)
    if source.md5 != contract.source_vcf_md5:
        raise ChallengeFixtureError("fixture source MD5 differs from the GIAB contract")
    try:
        validate_relative_classification(source.path, ArtifactClassification.PUBLIC)
        validate_relative_classification(
            contract.source_manifest_path, ArtifactClassification.PUBLIC
        )
    except GovernanceError as error:
        raise ChallengeFixtureError(str(error)) from error
    source_path = root / source.path
    source_sha256, source_md5, source_size = _checksums(source_path)
    if source_md5 != source.md5 or source_size != source.size_bytes:
        raise ChallengeFixtureError("installed fixture source differs from published GIAB")
    manifest = _load_manifest(root / contract.source_manifest_path)
    if manifest.resource_id != contract.source_resource_id:
        raise ChallengeFixtureError("fixture source manifest identifies another resource")
    matching_outputs = [output for output in manifest.outputs if output.path == source.path]
    if (
        len(matching_outputs) != 1
        or matching_outputs[0].sha256 != source_sha256
        or matching_outputs[0].size_bytes != source_size
    ):
        raise ChallengeFixtureError("fixture source is not bound by the GIAB manifest")
    return source_path, source_sha256


def _select_public_records(path: Path, sample: str) -> tuple[_PublicRecord, ...]:
    snvs: list[_PublicRecord] = []
    insertion: _PublicRecord | None = None
    deletion: _PublicRecord | None = None
    multiallelic: _PublicRecord | None = None
    samples: tuple[str, ...] | None = None
    try:
        with gzip.open(path, "rt", encoding="utf-8") as source:
            for line in source:
                if line.startswith("#CHROM\t"):
                    samples = tuple(line.rstrip("\n").split("\t")[9:])
                    continue
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 10 or fields[0] not in AUTOSOME_ORDER:
                    continue
                try:
                    position = int(fields[1])
                except ValueError as error:
                    raise ChallengeFixtureError("public HG002 VCF position is invalid") from error
                ref = fields[3]
                alts = tuple(fields[4].split(","))
                if (
                    position <= 0
                    or not SEQUENCE.fullmatch(ref)
                    or not alts
                    or not all(SEQUENCE.fullmatch(alt) and len(alt) <= 50 for alt in alts)
                ):
                    continue
                record_sha256 = hashlib.sha256(line.rstrip("\n").encode()).hexdigest()
                if len(alts) >= 2 and multiallelic is None:
                    multiallelic = _PublicRecord(
                        "multiallelic",
                        fields[0],
                        position,
                        ref,
                        alts[:2],
                        record_sha256,
                    )
                if len(alts) == 1 and len(ref) == 1 and len(alts[0]) == 1 and len(snvs) < 2:
                    snvs.append(
                        _PublicRecord(
                            "phased_snv" if not snvs else "missing_quality_snv",
                            fields[0],
                            position,
                            ref,
                            alts,
                            record_sha256,
                        )
                    )
                elif len(alts) == 1 and len(alts[0]) > len(ref) and insertion is None:
                    insertion = _PublicRecord(
                        "unphased_insertion",
                        fields[0],
                        position,
                        ref,
                        alts,
                        record_sha256,
                    )
                elif len(alts) == 1 and len(ref) > len(alts[0]) and deletion is None:
                    deletion = _PublicRecord(
                        "filtered_deletion",
                        fields[0],
                        position,
                        ref,
                        alts,
                        record_sha256,
                    )
                if len(snvs) == 2 and insertion and deletion and multiallelic:
                    break
    except (OSError, EOFError, gzip.BadGzipFile, UnicodeDecodeError) as error:
        raise ChallengeFixtureError(f"cannot read public HG002 VCF: {error}") from error
    if samples != (sample,):
        raise ChallengeFixtureError(f"public fixture source sample is not exactly {sample}")
    if len(snvs) != 2 or insertion is None or deletion is None or multiallelic is None:
        raise ChallengeFixtureError("public HG002 VCF lacks a required fixture allele class")
    selected = (*snvs, insertion, deletion, multiallelic)
    return tuple(
        sorted(
            selected,
            key=lambda record: (
                AUTOSOME_ORDER[record.contig],
                record.position,
                record.profile,
            ),
        )
    )


def _load_aliases(
    reference_contract: ReferenceIdentityContract, root: Path
) -> tuple[Path, Mapping[str, str]]:
    alias_relative = reference_contract.expected_output("contig_aliases").path
    try:
        validate_relative_classification(alias_relative, ArtifactClassification.PUBLIC)
        data = json.loads((root / alias_relative).read_text(encoding="utf-8"))
        aliases = data["tool_to_alignment"]
    except (GovernanceError, KeyError, OSError, json.JSONDecodeError, TypeError) as error:
        raise ChallengeFixtureError(f"cannot load public contig aliases: {error}") from error
    if not isinstance(aliases, Mapping):
        raise ChallengeFixtureError("tool_to_alignment aliases must be an object")
    return root / alias_relative, aliases


def _profile_fields(record: _PublicRecord) -> tuple[str, str, str, str]:
    if record.profile == "phased_snv":
        return "60", "PASS", "GT:DP:AD:GQ:PL:PS", "0|1:35:18,17:50:50,0,55:1001"
    if record.profile == "missing_quality_snv":
        return ".", "PASS", "GT:DP:AD:GQ:PL", "0/1:.:.,.:.:."
    if record.profile == "unphased_insertion":
        return "45", "PASS", "GT:DP:AD:GQ:PL", "0/1:31:16,15:42:42,0,48"
    if record.profile == "filtered_deletion":
        return "8", "LowQual", "GT:DP:AD:GQ:PL", "0/1:9:7,2:7:7,0,35"
    if record.profile == "multiallelic":
        return "50", "PASS", "GT:DP:AD:GQ:PL", "1/2:40:10,16,14:45:80,50,45,55,0,60"
    raise ChallengeFixtureError(f"unknown fixture profile: {record.profile}")


def _bgzf_block(payload: bytes) -> bytes:
    compressor = zlib.compressobj(level=6, wbits=-15)
    compressed = compressor.compress(payload) + compressor.flush()
    total_size = 12 + 6 + len(compressed) + 8
    if len(payload) > 65536 or total_size > 65536:
        raise ChallengeFixtureError("fixture payload exceeds one BGZF block")
    header = (
        b"\x1f\x8b\x08\x04"
        + b"\x00\x00\x00\x00"
        + b"\x00\xff"
        + struct.pack("<H", 6)
    )
    extra = b"BC" + struct.pack("<H", 2) + struct.pack("<H", total_size - 1)
    footer = struct.pack("<II", binascii.crc32(payload) & 0xFFFFFFFF, len(payload))
    return header + extra + compressed + footer


def _tabix_bytes(record_offsets: Mapping[str, tuple[int, int]]) -> bytes:
    names = b"".join(name.encode("ascii") + b"\x00" for name in record_offsets)
    raw = b"TBI\x01" + struct.pack(
        "<8i", len(record_offsets), 2, 1, 2, 0, ord("#"), 0, len(names)
    )
    raw += names
    for begin, end in record_offsets.values():
        raw += struct.pack("<iIiQQi", 1, 0, 1, begin, end, 0)
    return _bgzf_block(raw) + _bgzf_block(b"")


def _build_fixture_bytes(
    records: tuple[_PublicRecord, ...],
    aliases: Mapping[str, str],
    source_path: Path,
    source_sha256: str,
    alias_path: Path,
    contract: ChallengeFixtureContract,
    reference_contract: ReferenceIdentityContract,
    reference_role: ReferenceRole,
    root: Path,
) -> tuple[bytes, bytes, bytes]:
    used_contigs = []
    for record in records:
        challenge_name = aliases.get(record.contig)
        if challenge_name != record.contig.removeprefix("chr"):
            raise ChallengeFixtureError(f"invalid unprefixed alias for {record.contig}")
        if challenge_name not in used_contigs:
            used_contigs.append(challenge_name)
    header = (
        "##fileformat=VCFv4.2\n"
        "##source=Genosapiens-public-challenge-shaped-fixture-v2\n"
        "##reference=GCA_000001405.15+GCA_000786075.2\n"
        "##FILTER=<ID=LowQual,Description=Deterministic synthetic low-quality profile>\n"
        "##INFO=<ID=PUBLIC_SOURCE,Number=1,Type=String,Description=Public allele source>\n"
        "##INFO=<ID=PROFILE,Number=1,Type=String,Description=Deterministic fixture profile>\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Read depth>\n"
        "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allelic depths>\n"
        "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>\n"
        "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Phred genotype likelihoods>\n"
        "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=Phase set>\n"
        "##FORMAT=<ID=PGT,Number=1,Type=String,Description=Physical phasing genotype>\n"
        "##FORMAT=<ID=PID,Number=1,Type=String,Description=Physical phasing ID>\n"
        + "".join(
            f"##contig=<ID={name},length={reference_contract.canonical_contig_lengths[name]}>\n"
            for name in used_contigs
        )
        + f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{contract.sample}\n"
    ).encode()
    payload = bytearray(header)
    contig_offsets: dict[str, tuple[int, int]] = {}
    provenance_records = []
    challenge_fasta_output = reference_contract.expected_output(
        reference_role.fasta_output_role
    )
    challenge_fasta = root / challenge_fasta_output.path
    challenge_fai = root / reference_contract.expected_output(
        reference_role.fai_output_role
    ).path
    reference = IndexedFasta(challenge_fasta, challenge_fai)
    for index, record in enumerate(records, start=1):
        challenge_name = str(aliases[record.contig])
        if reference.fetch(challenge_name, record.position, len(record.ref)) != record.ref:
            raise ChallengeFixtureError(f"public fixture REF mismatch: {record.contig}")
        quality, filter_value, format_value, sample_value = _profile_fields(record)
        begin = len(payload)
        payload.extend(
            (
                f"{challenge_name}\t{record.position}\t.\t{record.ref}\t"
                f"{','.join(record.alts)}\t{quality}\t{filter_value}\t"
                f"PUBLIC_SOURCE=GIAB_HG002_NISTv4.2.1;PROFILE={record.profile}\t"
                f"{format_value}\t{sample_value}\n"
            ).encode()
        )
        previous = contig_offsets.get(challenge_name)
        contig_offsets[challenge_name] = (
            previous[0] if previous else begin,
            len(payload),
        )
        provenance_records.append(
            {
                "fixture_record_id": f"public-{index:02d}",
                "chrom": challenge_name,
                "position": record.position,
                "ref": record.ref,
                "alts": list(record.alts),
                "profile": record.profile,
                "synthetic_fields": ["QUAL", "FILTER", "INFO", "FORMAT", "SYNTHETIC"],
                "allele_sources": [
                    {
                        "alt": alt,
                        "kind": "public_giab_hg002",
                        "resource_id": contract.source_resource_id,
                        "release": contract.source_release,
                        "source_chrom": record.contig,
                        "source_position": record.position,
                        "source_ref": record.ref,
                        "source_alt": alt,
                        "source_record_sha256": record.source_record_sha256,
                    }
                    for alt in record.alts
                ],
            }
        )
    source_relative = source_path.relative_to(root).as_posix()
    alias_relative = alias_path.relative_to(root).as_posix()
    provenance = {
        "schema": PROVENANCE_SCHEMA,
        "fixture_id": contract.resource_id,
        "release": contract.release,
        "classification": "public",
        "sample": contract.sample,
        "assembly": contract.assembly,
        "reference_binding": asdict(contract.reference_binding),
        "controlled_inputs_used": False,
        "generation_inputs": [
            {
                "path": source_relative,
                "classification": "public",
                "sha256": source_sha256,
            },
            {
                "path": alias_relative,
                "classification": "public",
                "sha256": _checksums(alias_path)[0],
            },
            {
                "path": challenge_fasta_output.path,
                "classification": "public",
                "sha256": _checksums(challenge_fasta)[0],
            },
        ],
        "structural_basis": contract.structural_basis,
        "records": provenance_records,
    }
    return (
        _bgzf_block(bytes(payload)) + _bgzf_block(b""),
        _tabix_bytes(contig_offsets),
        (json.dumps(provenance, indent=2, sort_keys=True) + "\n").encode(),
    )


def _validate_provenance(
    path: Path, contract: ChallengeFixtureContract
) -> Mapping[tuple[str, int, str, tuple[str, ...]], Mapping[str, Any]]:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ChallengeFixtureError(f"cannot read fixture provenance: {error}") from error
    if (
        not isinstance(data, Mapping)
        or data.get("schema") != PROVENANCE_SCHEMA
        or data.get("fixture_id") != contract.resource_id
        or data.get("classification") != "public"
        or data.get("controlled_inputs_used") is not False
        or data.get("sample") != contract.sample
        or data.get("reference_binding") != asdict(contract.reference_binding)
    ):
        raise ChallengeFixtureError("fixture provenance identity or data boundary is invalid")
    inputs = data.get("generation_inputs")
    records = data.get("records")
    if not isinstance(inputs, list) or not isinstance(records, list) or not records:
        raise ChallengeFixtureError("fixture provenance inputs or records are invalid")
    for item in inputs:
        if not isinstance(item, Mapping) or item.get("classification") != "public":
            raise ChallengeFixtureError("fixture provenance contains a non-public input")
        try:
            validate_relative_classification(
                str(item["path"]), ArtifactClassification.PUBLIC
            )
        except (KeyError, GovernanceError) as error:
            raise ChallengeFixtureError("fixture provenance input path is invalid") from error
    by_variant = {}
    for record in records:
        if not isinstance(record, Mapping):
            raise ChallengeFixtureError("fixture provenance record must be an object")
        try:
            key = (
                str(record["chrom"]),
                int(record["position"]),
                str(record["ref"]),
                tuple(record["alts"]),
            )
            sources = record["allele_sources"]
        except (KeyError, TypeError, ValueError) as error:
            raise ChallengeFixtureError("fixture provenance record is invalid") from error
        if (
            key in by_variant
            or not isinstance(sources, list)
            or len(sources) != len(key[3])
            or {source.get("alt") for source in sources if isinstance(source, Mapping)}
            != set(key[3])
            or any(
                not isinstance(source, Mapping)
                or source.get("kind") != "public_giab_hg002"
                or source.get("resource_id") != contract.source_resource_id
                for source in sources
            )
        ):
            raise ChallengeFixtureError("fixture allele provenance is incomplete")
        by_variant[key] = record
    return by_variant


def validate_challenge_fixture(
    contract: ChallengeFixtureContract,
    reference_contract: ReferenceIdentityContract,
    reference_root: str | Path,
) -> ChallengeFixtureValidation:
    root = Path(reference_root).resolve()
    try:
        reference_role = contract.reference_binding.resolve(
            reference_contract, required_role="challenge_alignment"
        )
    except ReferenceConcordanceError as error:
        raise ChallengeFixtureError(str(error)) from error
    provenance = _validate_provenance(root / contract.outputs["provenance"], contract)
    challenge_fasta = root / reference_contract.expected_output(
        reference_role.fasta_output_role
    ).path
    challenge_fai = root / reference_contract.expected_output(
        reference_role.fai_output_role
    ).path
    reference = IndexedFasta(challenge_fasta, challenge_fai)
    records = []
    samples = None
    phase_headers = set()
    contig_lengths = {}
    previous_key = (-1, -1)
    phased = unphased = filtered = missing = alleles = 0
    try:
        with gzip.open(root / contract.outputs["fixture_vcf"], "rt", encoding="utf-8") as source:
            for line in source:
                if line.startswith("##FORMAT=<ID="):
                    phase_headers.add(line.split("=", 2)[2].split(",", 1)[0])
                    continue
                if line.startswith("##contig=<ID="):
                    content = line.split("<", 1)[1].split(">", 1)[0]
                    fields = dict(item.split("=", 1) for item in content.split(","))
                    contig_lengths[fields["ID"]] = int(fields["length"])
                    continue
                if line.startswith("#CHROM\t"):
                    samples = tuple(line.rstrip("\n").split("\t")[9:])
                    continue
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 10 or fields[0].startswith("chr"):
                    raise ChallengeFixtureError("fixture is not a one-sample unprefixed VCF")
                position = int(fields[1])
                key_order = (int(fields[0]), position)
                if key_order < previous_key:
                    raise ChallengeFixtureError("fixture VCF is not coordinate-sorted")
                previous_key = key_order
                alts = tuple(fields[4].split(","))
                key = (fields[0], position, fields[3], alts)
                if key not in provenance:
                    raise ChallengeFixtureError("fixture record lacks allele provenance")
                if reference.fetch(fields[0], position, len(fields[3])) != fields[3]:
                    raise ChallengeFixtureError("fixture record REF mismatch")
                format_keys = fields[8].split(":")
                sample_values = fields[9].split(":")
                if len(format_keys) != len(sample_values) or format_keys[0] != "GT":
                    raise ChallengeFixtureError("fixture FORMAT/sample fields are inconsistent")
                genotype = sample_values[0]
                phased += int("|" in genotype)
                unphased += int("/" in genotype)
                filtered += int(fields[6] not in {"PASS", "."})
                missing += int("." in sample_values[1:] or fields[5] == ".")
                alleles += len(alts)
                records.append((fields[0], position, fields[3], alts, fields[6]))
    except (KeyError, OSError, ValueError, EOFError, gzip.BadGzipFile) as error:
        raise ChallengeFixtureError(f"invalid challenge-shaped fixture VCF: {error}") from error
    if samples != (contract.sample,) or {"GT", "PS", "PGT", "PID"} - phase_headers:
        raise ChallengeFixtureError("fixture sample or phasing-capable headers are invalid")
    if set(provenance) != {record[:4] for record in records}:
        raise ChallengeFixtureError("fixture provenance and VCF records differ")
    if phased == 0 or unphased == 0 or filtered == 0 or missing == 0:
        raise ChallengeFixtureError("fixture does not exercise required synthetic edge cases")
    for contig, length in contig_lengths.items():
        if length != reference_contract.canonical_contig_lengths[contig]:
            raise ChallengeFixtureError(f"fixture contig length mismatch: {contig}")
    try:
        index = TabixIndex.load(root / contract.outputs["fixture_tbi"])
        index.require_vcf_layout()
        for contig, position, ref, alts, filter_value in records:
            queried = index.query(
                root / contract.outputs["fixture_vcf"], contig, position, position
            )
            if not any(
                line.split("\t")[0] == contig
                and line.split("\t")[1] == str(position)
                and line.split("\t")[3] == ref
                and tuple(line.split("\t")[4].split(",")) == alts
                and line.split("\t")[6] == filter_value
                for line in queried
            ):
                raise ChallengeFixtureError("fixture tabix query missed a record")
    except (OSError, KeyError, TabixError) as error:
        raise ChallengeFixtureError(f"fixture tabix index is invalid: {error}") from error
    return ChallengeFixtureValidation(
        records=len(records),
        public_alleles=alleles,
        phased_records=phased,
        unphased_records=unphased,
        filtered_records=filtered,
        missing_quality_records=missing,
        indexed_contigs=len(index.references),
    )


def _resource_outputs(
    root: Path, contract: ChallengeFixtureContract
) -> tuple[ResourceOutput, ...]:
    outputs = []
    for path in contract.outputs.values():
        sha256, _md5, size = _checksums(root / path)
        outputs.append(
            ResourceOutput(
                path=path,
                size_bytes=size,
                sha256=sha256,
                retention=RetentionPolicy("durable", None),
            )
        )
    return tuple(outputs)


def _verify_manifest_contract(
    manifest: PublicResourceManifest,
    contract: ChallengeFixtureContract,
    source_url: str,
) -> None:
    if (
        manifest.resource_id != contract.resource_id
        or manifest.release != contract.release
        or manifest.assembly != contract.assembly
        or manifest.source.url != source_url
        or {output.path for output in manifest.outputs} != set(contract.outputs.values())
    ):
        raise ChallengeFixtureError("fixture manifest differs from its contract")


def bootstrap_challenge_fixture(
    *,
    contract_path: str | Path,
    giab_contract_path: str | Path,
    reference_identity_path: str | Path,
    reference_root: str | Path,
    image_digest: str,
    generated_on: str | None = None,
) -> ChallengeFixtureResult:
    contract_file = Path(contract_path)
    contract = load_challenge_fixture_contract(contract_file)
    giab_contract = load_giab_contract(giab_contract_path)
    try:
        reference_contract = load_reference_identity(reference_identity_path)
    except ReferenceConcordanceError as error:
        raise ChallengeFixtureError(str(error)) from error
    root = Path(reference_root).resolve(strict=True)
    try:
        reference_role = contract.reference_binding.resolve(
            reference_contract, required_role="challenge_alignment"
        )
    except ReferenceConcordanceError as error:
        raise ChallengeFixtureError(str(error)) from error
    source_path, source_sha256 = _source_identity(contract, giab_contract, root)
    source_url = giab_contract.artifact(contract.source_vcf_role).url
    manifest_path = root / f"manifests/fixtures/{contract.resource_id}.json"
    if manifest_path.exists():
        manifest = _load_manifest(manifest_path)
        _verify_manifest_contract(manifest, contract, source_url)
        try:
            verify_outputs(manifest, root)
        except ManifestError as error:
            raise ChallengeFixtureError(f"completed fixture is invalid: {error}") from error
        validation = validate_challenge_fixture(contract, reference_contract, root)
        return ChallengeFixtureResult(
            str(manifest_path),
            _checksums(manifest_path)[0],
            str(root / contract.outputs["fixture_vcf"]),
            str(root / contract.outputs["provenance"]),
            validation,
            True,
        )
    target = root / str(PurePosixPath(contract.outputs["fixture_vcf"]).parent)
    if target.exists():
        if not target.is_dir():
            raise ChallengeFixtureError("fixture release path is not a directory")
        validation = validate_challenge_fixture(contract, reference_contract, root)
    else:
        require_capacity(
            root,
            CapacityPlan(
                operation="build public challenge-shaped small-variant fixture",
                output_bytes=10_000_000,
                scratch_bytes=10_000_000,
            ),
        )
        records = _select_public_records(source_path, "HG002")
        alias_path, aliases = _load_aliases(reference_contract, root)
        vcf_bytes, tbi_bytes, provenance_bytes = _build_fixture_bytes(
            records,
            aliases,
            source_path,
            source_sha256,
            alias_path,
            contract,
            reference_contract,
            reference_role,
            root,
        )
        target.parent.mkdir(parents=True, exist_ok=True)
        partial_root = target.parent / ".partial"
        partial_root.mkdir(mode=0o700, exist_ok=True)
        staging = partial_root / f"{target.name}.{uuid.uuid4().hex}.partial"
        staging.mkdir(mode=0o700)
        try:
            content_by_role = {
                "fixture_vcf": vcf_bytes,
                "fixture_tbi": tbi_bytes,
                "provenance": provenance_bytes,
            }
            staging_outputs = {}
            for role, content in content_by_role.items():
                output = staging / PurePosixPath(contract.outputs[role]).name
                with output.open("xb") as destination:
                    destination.write(content)
                    destination.flush()
                    os.fsync(destination.fileno())
                os.chmod(output, 0o444)
                staging_outputs[role] = output.relative_to(root).as_posix()
            validation = validate_challenge_fixture(
                replace(contract, outputs=staging_outputs), reference_contract, root
            )
            os.replace(staging, target)
            directory_fd = os.open(target.parent, os.O_RDONLY)
            try:
                os.fsync(directory_fd)
            finally:
                os.close(directory_fd)
        except Exception:
            shutil.rmtree(staging, ignore_errors=True)
            raise
    manifest = PublicResourceManifest(
        resource_id=contract.resource_id,
        release=contract.release,
        source=SourceIdentity(url=source_url, sha256=source_sha256),
        license="CC0-1.0",
        assembly=contract.assembly,
        retrieved_on=generated_on or date.today().isoformat(),
        transformation=Transformation(
            schema=TRANSFORMATION_SCHEMA,
            parameters_sha256=_checksums(contract_file)[0],
        ),
        image_digest=image_digest,
        outputs=_resource_outputs(root, contract),
    )
    try:
        publication = publish_manifest(manifest, manifest_path, root)
    except ManifestError as error:
        raise ChallengeFixtureError(f"cannot publish fixture manifest: {error}") from error
    return ChallengeFixtureResult(
        str(publication.path),
        publication.manifest_sha256,
        str(root / contract.outputs["fixture_vcf"]),
        str(root / contract.outputs["provenance"]),
        validation,
        publication.reused,
    )


def challenge_fixture_result_as_json(result: ChallengeFixtureResult) -> str:
    return json.dumps(asdict(result), indent=2, sort_keys=True) + "\n"
