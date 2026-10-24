"""Deterministically materialize the public Track 1 validation seed package."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict, dataclass
import hashlib
import json
import os
from pathlib import Path
import re
import stat

from genosapiens.bounded_json import BoundedJSONError, load_bounded_json


GENERATOR_SCHEMA = "genosapiens.track1-synthetic-case-generator"
RECIPE_SCHEMA = "genosapiens.track1-synthetic-case-recipe"
PERTURBATION_SCHEMA = "genosapiens.track1-synthetic-case-perturbations"
BACKGROUND_SCHEMA = "genosapiens.track1-synthetic-public-background"
TRUTH_SCHEMA = "genosapiens.track1-synthetic-causal-truth"
NATIVE_SCHEMA = "genosapiens.track1-synthetic-native-input"
SOURCE_FORMAT_SCHEMA = "genosapiens.supplied-vcf-format"
MATERIALIZED_ROOT = Path("tests/fixtures/track1/materialized")
REFERENCE_ROLE = "public_annotation_benchmark"
REFERENCE_IDENTITY = "grch38-no-alt-hs38d1-public-native-v1"
REFERENCE_MANIFEST_SHA256 = (
    "27eaf00fde35e5fa83f945c63e5dee871960374d90488cc0797b5ebea2fb1c25"
)
REPEAT_IDENTITY = "longtr-fd73e0d-vcf4.1/v1"
SMALL_VARIANT_REPEAT_IDENTITY = "repid-ru-repcn-pos-end"
PUBLIC_BACKGROUND_CONTRACT = Path(
    "config/public-resources/challenge-shaped-fixture.json"
)
PUBLIC_BACKGROUND_SOURCE = {
    "contract_sha256": "79a1a5c64a587bcdf81924f6eff530ae7d8d2b1206a9d72aaa5fd31df6e2d21e",
    "manifest_sha256": "3a7901ba0814842fc799781916e3faf1213ac15ab45567f6b4758c90dcd403cb",
    "vcf_sha256": "758456587bdcf53e7aa56ce8abfbbeb65460d1f3c80ca5f10c9472810cd49e7a",
    "provenance_sha256": "3e0820dba8474c87c572ce786a46d00d4419d09ea5552041f423542824f0e966",
}
PUBLIC_BACKGROUND_RECORDS = {
    "unphased_insertion": (
        "public-01", "1", 602113, "T", ("TGCCCA",),
        "15f598aa7240d04822262d8d0a56b4372b4e12653de482e989d0a065bbffe461",
    ),
    "phased_snv": (
        "public-02", "1", 604358, "G", ("C",),
        "838676330d71289fd571b558d4dd412e0b93e403a4a34d0dd7e779bcfe00ff7d",
    ),
    "missing_quality_snv": (
        "public-03", "1", 604360, "T", ("C",),
        "bd4c4b0e7202a2a849363a523bff0d2821a1fbda32659c89f92728a9bd2865c5",
    ),
    "filtered_deletion": (
        "public-04", "1", 631859, "CG", ("C",),
        "5f437304b94030541a9c44c4d487a07d4c59bba2144cadd1a097749fe2e9e8c9",
    ),
    "multiallelic": (
        "public-05", "1", 1735008, "G", ("GGTTT", "GTTTT"),
        "797d2fc60e07d0860c9725cb9d8bceba7db03d1f5b059c5c76107feac381638c",
    ),
}
SHA256 = re.compile(r"^[0-9a-f]{64}$")
SEQUENCE = re.compile(r"^[ACGTN]+$")
CONTIG_ORDER = {
    name: index
    for index, name in enumerate(("1", "19", "21", "22", "X", "M"))
}
REQUIRED_SCOPES = (
    "compound_heterozygous",
    "dominant",
    "homozygous_recessive",
    "incomplete_recessive",
    "mitochondrial",
    "mixed_class",
    "mosaic",
    "regulatory",
    "repeat",
    "sv_cnv",
    "x_linked",
)
EXPECTED_PERTURBATIONS = {
    "development": {
        "source-filtered-call",
        "missing-quality-fields",
        "low-allele-balance-call",
        "unphased-recessive-members",
    },
    "held_out": {
        "bounded-structural-uncertainty",
        "conflicting-repeat-correspondence",
        "missing-quality-fields",
        "partial-repeat-projection",
    },
}
GENERATED_FILENAMES = (
    "causal-truth.manifest.json",
    "development.hpo.txt",
    "development.source-format.json",
    "development.vcf",
    "held-out-longtr.source-format.json",
    "held-out-longtr.vcf",
    "held-out-native.json",
    "held-out-small-native.vcf",
    "held-out-small.source-format.json",
    "held-out.hpo.txt",
)


class ValidationSeedError(RuntimeError):
    """Raised when public seed inputs or generated outputs are invalid."""


@dataclass(frozen=True, slots=True)
class GeneratedArtifact:
    filename: str
    size_bytes: int
    sha256: str


@dataclass(frozen=True, slots=True)
class ValidationSeedResult:
    generator_id: str
    generator_version: str
    output_root: str
    artifacts: tuple[GeneratedArtifact, ...]
    controlled_sources_used: bool = False
    validation_only: bool = True


@dataclass(frozen=True, slots=True)
class _VcfRecord:
    contig: str
    position: int
    record_id: str
    ref: str
    alt: str
    quality: str
    filter_value: str
    info: str
    format_value: str
    sample_value: str

    @property
    def identity(self) -> str:
        return f"{self.contig}:{self.position}:{self.ref}:{self.alt}"

    def render(self) -> str:
        return "\t".join(
            (
                self.contig,
                str(self.position),
                self.record_id,
                self.ref,
                self.alt,
                self.quality,
                self.filter_value,
                self.info,
                self.format_value,
                self.sample_value,
            )
        )


def _json_bytes(document: Mapping[str, object]) -> bytes:
    return (json.dumps(document, indent=2, ensure_ascii=True) + "\n").encode("ascii")


def _sha256(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def _load_public_json(
    repository_root: Path,
    relative_path: Path,
    *,
    fields: set[str],
    schema: str,
) -> Mapping[str, object]:
    try:
        document, _content = load_bounded_json(
            repository_root / relative_path,
            reference_root=repository_root,
        )
    except BoundedJSONError as error:
        raise ValidationSeedError("validation seed input cannot be read") from error
    if (
        not isinstance(document, Mapping)
        or set(document) != fields
        or document.get("schema") != schema
        or document.get("classification") != "public"
        or document.get("validation_only") is not True
    ):
        raise ValidationSeedError("validation seed input contract is invalid")
    return document


def _load_recipe(repository_root: Path, split: str) -> Mapping[str, object]:
    name = "development" if split == "development" else "held-out"
    document = _load_public_json(
        repository_root,
        MATERIALIZED_ROOT / f"{name}-recipe.manifest.json",
        fields={
            "schema",
            "classification",
            "validation_only",
            "case_id",
            "split",
            "recipe_group_id",
            "generator_id",
            "seed",
            "background_id",
            "input_roles",
            "hpo_terms",
            "output_profile",
            "scope",
        },
        schema=RECIPE_SCHEMA,
    )
    hpo_terms = document.get("hpo_terms")
    seed = document.get("seed")
    expected_roles = (
        ["plain-hpo-text", "synthetic-overlay-vcf"]
        if split == "development"
        else ["plain-hpo-text", "native-input-manifest"]
    )
    if (
        document.get("split") != split
        or document.get("case_id") != f"{name}-case"
        or document.get("recipe_group_id") != f"{name}-recipe"
        or document.get("generator_id") != "track1-synthetic-case-generator"
        or document.get("background_id") != "hg002-public-background"
        or document.get("input_roles") != expected_roles
        or document.get("output_profile") != "all-class-seed-v2"
        or document.get("scope") != "materialization-seed-only"
        or type(seed) is not int
        or seed < 0
        or not isinstance(hpo_terms, list)
        or len(hpo_terms) != len(set(hpo_terms))
        or not hpo_terms
        or any(
            type(term) is not str or re.fullmatch(r"HP:[0-9]{7}", term) is None
            for term in hpo_terms
        )
    ):
        raise ValidationSeedError("validation seed recipe is invalid")
    return document


def _load_perturbations(repository_root: Path, split: str) -> Mapping[str, object]:
    name = "development" if split == "development" else "held-out"
    document = _load_public_json(
        repository_root,
        MATERIALIZED_ROOT / f"{name}-perturbations.manifest.json",
        fields={
            "schema",
            "classification",
            "validation_only",
            "case_id",
            "split",
            "perturbation_group_id",
            "operations",
        },
        schema=PERTURBATION_SCHEMA,
    )
    operations = document.get("operations")
    if (
        document.get("split") != split
        or document.get("case_id") != f"{name}-case"
        or document.get("perturbation_group_id") != f"{name}-perturbations"
        or not isinstance(operations, list)
        or set(operations) != EXPECTED_PERTURBATIONS[split]
        or len(operations) != len(set(operations))
    ):
        raise ValidationSeedError("validation seed perturbations are invalid")
    return document


def _load_background(repository_root: Path) -> tuple[_VcfRecord, ...]:
    try:
        source_contract, source_contract_bytes = load_bounded_json(
            repository_root / PUBLIC_BACKGROUND_CONTRACT,
            reference_root=repository_root,
        )
    except BoundedJSONError as error:
        raise ValidationSeedError(
            "validation seed public background contract cannot be read"
        ) from error
    if (
        not isinstance(source_contract, Mapping)
        or source_contract.get("schema")
        != "genosapiens.challenge-shaped-fixture"
        or source_contract.get("resource_id")
        != "challenge-shaped-public-small-variant-fixture-v2"
        or source_contract.get("release") != "v2"
        or _sha256(source_contract_bytes)
        != PUBLIC_BACKGROUND_SOURCE["contract_sha256"]
    ):
        raise ValidationSeedError(
            "validation seed public background contract is invalid"
        )
    document = _load_public_json(
        repository_root,
        MATERIALIZED_ROOT / "hg002-public-background.json",
        fields={
            "schema",
            "classification",
            "validation_only",
            "background_id",
            "subject_id",
            "controlled_inputs_used",
            "source_fixture",
            "records",
        },
        schema=BACKGROUND_SCHEMA,
    )
    source = document.get("source_fixture")
    records = document.get("records")
    if (
        document.get("background_id") != "hg002-public-background"
        or document.get("subject_id") != "HG002"
        or document.get("controlled_inputs_used") is not False
        or not isinstance(source, Mapping)
        or set(source)
        != {
            "resource_id",
            "release",
            "contract_sha256",
            "manifest_sha256",
            "vcf_sha256",
            "provenance_sha256",
        }
        or source.get("resource_id")
        != "challenge-shaped-public-small-variant-fixture-v2"
        or source.get("release") != "v2"
        or any(source.get(field) != value for field, value in PUBLIC_BACKGROUND_SOURCE.items())
        or not isinstance(records, list)
        or len(records) != 5
    ):
        raise ValidationSeedError("validation seed public background is invalid")

    profile_values = {
        "unphased_insertion": (
            "45",
            "PASS",
            "GT:DP:AD:GQ:PL",
            "0/1:31:16,15:42:42,0,48",
        ),
        "phased_snv": (
            "60",
            "PASS",
            "GT:DP:AD:GQ:PL:PS",
            "0|1:35:18,17:50:50,0,55:1001",
        ),
        "missing_quality_snv": (
            ".",
            "PASS",
            "GT:DP:AD:GQ:PL",
            "0/1:.:.,.:.:.",
        ),
        "filtered_deletion": (
            "8",
            "LowQual",
            "GT:DP:AD:GQ:PL",
            "0/1:9:7,2:7:7,0,35",
        ),
        "multiallelic": (
            "50",
            "PASS",
            "GT:DP:AD:GQ:PL",
            "1/2:40:10,16,14:45:80,50,45,55,0,60",
        ),
    }
    loaded = []
    seen_ids = set()
    seen_identities = set()
    seen_profiles = set()
    for item in records:
        if not isinstance(item, Mapping) or set(item) != {
            "fixture_record_id",
            "chrom",
            "position",
            "ref",
            "alts",
            "profile",
            "source_record_sha256",
        }:
            raise ValidationSeedError("validation seed background record is invalid")
        record_id = item.get("fixture_record_id")
        contig = item.get("chrom")
        position = item.get("position")
        ref = item.get("ref")
        alts = item.get("alts")
        profile = item.get("profile")
        source_sha256 = item.get("source_record_sha256")
        if (
            type(record_id) is not str
            or record_id in seen_ids
            or type(contig) is not str
            or contig not in CONTIG_ORDER
            or type(position) is not int
            or position <= 0
            or type(ref) is not str
            or SEQUENCE.fullmatch(ref) is None
            or not isinstance(alts, list)
            or not alts
            or any(
                type(alt) is not str or SEQUENCE.fullmatch(alt) is None
                for alt in alts
            )
            or type(profile) is not str
            or profile not in profile_values
            or profile in seen_profiles
            or (
                record_id,
                contig,
                position,
                ref,
                tuple(alts),
                source_sha256,
            )
            != PUBLIC_BACKGROUND_RECORDS.get(profile)
        ):
            raise ValidationSeedError("validation seed background record is invalid")
        identity = (contig, position, ref, tuple(alts))
        if identity in seen_identities:
            raise ValidationSeedError("validation seed background record is duplicated")
        quality, filter_value, format_value, sample_value = profile_values[profile]
        loaded.append(
            _VcfRecord(
                contig=str(contig),
                position=position,
                record_id=record_id,
                ref=ref,
                alt=",".join(alts),
                quality=quality,
                filter_value=filter_value,
                info=".",
                format_value=format_value,
                sample_value=sample_value,
            )
        )
        seen_ids.add(record_id)
        seen_identities.add(identity)
        seen_profiles.add(profile)
    if seen_profiles != set(profile_values):
        raise ValidationSeedError("validation seed background profiles are incomplete")
    return tuple(sorted(loaded, key=_record_order))


def _record_order(record: _VcfRecord) -> tuple[int, int, str]:
    return CONTIG_ORDER[record.contig], record.position, record.record_id


def _opaque_id(*parts: object) -> str:
    content = "\0".join(str(part) for part in parts).encode("ascii")
    return f"seed-{hashlib.sha256(content).hexdigest()[:16]}"


def _assign_opaque_ids(
    records: tuple[_VcfRecord, ...], *, split: str, seed: int
) -> tuple[_VcfRecord, ...]:
    assigned = []
    ids = set()
    for record in records:
        record_id = _opaque_id(split, seed, record.identity)
        if record_id in ids:
            raise ValidationSeedError("validation seed record identifier collision")
        assigned.append(
            _VcfRecord(
                record.contig,
                record.position,
                record_id,
                record.ref,
                record.alt,
                record.quality,
                record.filter_value,
                record.info,
                record.format_value,
                record.sample_value,
            )
        )
        ids.add(record_id)
    return tuple(assigned)


def _synthetic_records(
    split: str, seed: int
) -> tuple[tuple[_VcfRecord, ...], Mapping[str, str]]:
    if split == "development":
        positions = {
            "structural": 1_510_100,
            "regulatory": 1_510_400,
            "homozygous": 1_510_404,
            "incomplete": 1_510_800,
            "repeat": 1_511_000,
            "dominant": 44_101_000,
            "compound_a": 20_000_000,
            "compound_b": 20_000_200,
            "mosaic": 20_000_000,
            "x_linked": 10_101_000,
            "mitochondrial": 201,
            "decoys": (1_511_500, 1_511_600, 1_511_700),
        }
        refs = {
            "structural": "A",
            "regulatory": "C",
            "homozygous": "T",
            "incomplete": "TTC",
            "repeat": "G",
            "dominant": "C",
            "compound_a": "ACA",
            "compound_b": "A",
            "mosaic": "G",
            "x_linked": "A",
            "mitochondrial": "A",
            "decoys": ("G", "G", "G"),
        }
    else:
        positions = {
            "structural": 2_520_100,
            "regulatory": 2_520_400,
            "homozygous": 2_520_404,
            "incomplete": 2_520_800,
            "dominant": 44_202_000,
            "compound_a": 21_000_000,
            "compound_b": 21_000_200,
            "mosaic": 21_000_000,
            "x_linked": 10_202_000,
            "mitochondrial": 302,
            "decoys": (2_521_500, 2_521_600, 2_521_700),
        }
        refs = {
            "structural": "A",
            "regulatory": "G",
            "homozygous": "T",
            "incomplete": "ACT",
            "dominant": "A",
            "compound_a": "TAA",
            "compound_b": "A",
            "mosaic": "C",
            "x_linked": "A",
            "mitochondrial": "A",
            "decoys": ("G", "T", "G"),
        }
    structural_position = int(positions["structural"])
    structural_info = f"SVTYPE=DEL;END={structural_position + 600}"
    if split == "held_out":
        structural_info += ";CIPOS=-3,5;CIEND=-4,7"
    records = [
        _VcfRecord(
            "1",
            structural_position,
            ".",
            str(refs["structural"]),
            "<DEL>",
            "45" if split == "development" else "20",
            "PASS",
            structural_info,
            "GT:GQ:DP:AD:CN",
            (
                "0/1:35:22:12,10:1"
                if split == "development"
                else "0/1:18:19:11,8:1.2"
            ),
        ),
        _VcfRecord(
            "1",
            int(positions["regulatory"]),
            ".",
            str(refs["regulatory"]),
            "A",
            "12" if split == "development" else ".",
            "LowQual" if split == "development" else ".",
            ".",
            "GT:GQ:DP:AD",
            "0/1:8:18:12,6" if split == "development" else "0/1:.:.:.",
        ),
        _VcfRecord(
            "1",
            int(positions["homozygous"]),
            ".",
            str(refs["homozygous"]),
            "A",
            "50",
            "PASS",
            ".",
            "GT:GQ:DP:AD",
            "1/1:45:24:0,24",
        ),
        _VcfRecord(
            "1",
            int(positions["incomplete"]),
            ".",
            str(refs["incomplete"]),
            str(refs["incomplete"])[0],
            "." if split == "development" else "22",
            "." if split == "development" else "LowQual",
            ".",
            "GT:GQ:DP:AD",
            "0/1:.:.:." if split == "development" else "0/1:9:16:10,6",
        ),
        _VcfRecord(
            "19",
            int(positions["dominant"]),
            ".",
            str(refs["dominant"]),
            "T" if split == "development" else "C",
            "60" if split == "development" else "55",
            "PASS",
            ".",
            "GT:GQ:DP:AD",
            (
                "0/1:50:30:15,15"
                if split == "development"
                else "0/1:48:32:16,16"
            ),
        ),
        _VcfRecord(
            "21",
            int(positions["compound_a"]),
            ".",
            str(refs["compound_a"]),
            str(refs["compound_a"])[0],
            "40",
            "PASS",
            ".",
            "GT:GQ:DP:AD",
            "0/1:35:20:10,10",
        ),
        _VcfRecord(
            "21",
            int(positions["compound_b"]),
            ".",
            str(refs["compound_b"]),
            f"{refs['compound_b']}G",
            "38",
            "PASS",
            ".",
            "GT:GQ:DP:AD",
            "0/1:32:18:9,9",
        ),
        _VcfRecord(
            "22",
            int(positions["mosaic"]),
            ".",
            str(refs["mosaic"]),
            "A" if split == "development" else "T",
            "35" if split == "development" else "40",
            "PASS",
            ".",
            "GT:GQ:DP:AD",
            (
                "0/0:40:100:90,10"
                if split == "development"
                else "0/1:35:100:88,12"
            ),
        ),
        _VcfRecord(
            "X",
            int(positions["x_linked"]),
            ".",
            str(refs["x_linked"]),
            "G" if split == "development" else "C",
            "25" if split == "development" else "28",
            "PASS",
            ".",
            "GT:GQ:DP:AD",
            "1:20:16:0,16" if split == "development" else "1:24:18:0,18",
        ),
        _VcfRecord(
            "M",
            int(positions["mitochondrial"]),
            ".",
            str(refs["mitochondrial"]),
            "C",
            "42" if split == "development" else "39",
            "PASS",
            ".",
            "GT:GQ:DP:AD",
            (
                "0/1:30:80:48,32"
                if split == "development"
                else "0/1:28:100:65,35"
            ),
        ),
    ]
    repeat_id = _opaque_id(split, seed, "small-repeat")
    if split == "development":
        records.append(
            _VcfRecord(
                "1",
                int(positions["repeat"]),
                ".",
                str(refs["repeat"]),
                "<STR>",
                "30",
                "PASS",
                f"SVTYPE=STR;END={int(positions['repeat'])};REPID={repeat_id};RU=A",
                "GT:GQ:DP:AD:REPCN",
                "0/1:25:20:10,10:4/9",
            )
        )
    for index, (position, ref) in enumerate(
        zip(positions["decoys"], refs["decoys"], strict=True)
    ):
        records.append(
            _VcfRecord(
                "1",
                int(position),
                ".",
                str(ref),
                "A" if ref != "A" else "C",
                "40",
                "PASS",
                ".",
                "GT:GQ:DP:AD",
                "0/1:35:20:10,10" if index != 1 else "0/0:45:22:21,1",
            )
        )
    by_tag = {
        "structural": records[0].identity,
        "regulatory": records[1].identity,
        "homozygous": records[2].identity,
        "incomplete": records[3].identity,
        "dominant": records[4].identity,
        "compound_a": records[5].identity,
        "compound_b": records[6].identity,
        "mosaic": records[7].identity,
        "x_linked": records[8].identity,
        "mitochondrial": records[9].identity,
        "repeat": repeat_id,
    }
    return tuple(sorted(records, key=_record_order)), by_tag


def _vcf_bytes(records: tuple[_VcfRecord, ...], *, split: str) -> bytes:
    used_contigs = tuple(
        name
        for name in CONTIG_ORDER
        if any(record.contig == name for record in records)
    )
    lengths = {
        "1": 248956422,
        "19": 58617616,
        "21": 46709983,
        "22": 50818468,
        "X": 156040895,
        "M": 16569,
    }
    lines = [
        "##fileformat=VCFv4.2",
        "##reference=GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta",
        f"##source=Genosapiens-track1-{split}-validation-seed-v2",
        *(f"##contig=<ID={contig},length={lengths[contig]}>" for contig in used_contigs),
        '##FILTER=<ID=LowQual,Description="Deterministic synthetic low-quality profile">',
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Synthetic structural variant type">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="Synthetic end coordinate">',
        '##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Synthetic POS confidence interval">',
        '##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Synthetic END confidence interval">',
        '##INFO=<ID=REPID,Number=1,Type=String,Description="Synthetic repeat locus identifier">',
        '##INFO=<ID=RU,Number=1,Type=String,Description="Synthetic repeat unit">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
        '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred genotype likelihoods">',
        '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">',
        '##FORMAT=<ID=CN,Number=1,Type=Float,Description="Synthetic copy number">',
        '##FORMAT=<ID=REPCN,Number=1,Type=String,Description="Synthetic repeat copy counts">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC",
        *(record.render() for record in records),
    ]
    return ("\n".join(lines) + "\n").encode("ascii")


def _longtr_bytes(seed: int) -> tuple[bytes, str]:
    starts = (1_401_903, 1_402_106, 1_402_260)
    seed_digest = hashlib.sha256(f"held_out\0{seed}".encode("ascii")).hexdigest()
    lines = [
        "##fileformat=VCFv4.1",
        "##reference=GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta",
        "##source=Genosapiens-track1-held-out-longtr-validation-seed-v2",
        f"##genosapiens_seed_sha256={seed_digest}",
        "##contig=<ID=1,length=248956422>",
        '##INFO=<ID=START,Number=1,Type=Integer,Description="Synthetic LongTR field">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="Synthetic LongTR field">',
        '##INFO=<ID=MOTIF,Number=.,Type=String,Description="Synthetic LongTR field">',
        '##INFO=<ID=PERIOD,Number=.,Type=Integer,Description="Synthetic LongTR field">',
        '##INFO=<ID=INEXACT_ALLELE,Number=A,Type=Integer,Description="Synthetic LongTR field">',
        '##INFO=<ID=BPDIFFS,Number=A,Type=Integer,Description="Synthetic LongTR field">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=GB,Number=1,Type=String,Description="Synthetic LongTR field">',
        '##FORMAT=<ID=Q,Number=1,Type=Float,Description="Synthetic LongTR field">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC",
        f"1\t{starts[0]}\t.\tAAAA\tAA,AAAAA\t.\t.\tSTART={starts[0]};END={starts[0] + 3};MOTIF=A;PERIOD=1;INEXACT_ALLELE=0,0;BPDIFFS=-2,1\tGT:GB:Q:DP\t2|1:1|-2:0.99:24",
        f"1\t{starts[1]}\t.\tAAAA\tAA,AAAAA\t.\t.\tSTART={starts[1]};END={starts[1] + 3};MOTIF=A;PERIOD=1;INEXACT_ALLELE=0,1;BPDIFFS=-2,1\tGT:GB:Q:DP\t1|2:-2|1:0.80:18",
        f"1\t{starts[2]}\t.\tAAAA\tAAAAA\t.\t.\tSTART={starts[2]};END={starts[2] + 3};MOTIF=A;PERIOD=1;INEXACT_ALLELE=0;BPDIFFS=1\tGT:GB:Q:DP\t1|1:2|2:0.95:20",
    ]
    return (
        ("\n".join(lines) + "\n").encode("ascii"),
        f"longtr-region:1:{starts[0]}-{starts[0] + 3}",
    )


def _truth_case(case_id: str, identities: Mapping[str, str]) -> Mapping[str, object]:
    members = {
        "compound_heterozygous": [identities["compound_a"], identities["compound_b"]],
        "dominant": [identities["dominant"]],
        "homozygous_recessive": [identities["homozygous"]],
        "incomplete_recessive": [identities["incomplete"]],
        "mitochondrial": [identities["mitochondrial"]],
        "mixed_class": [identities["structural"], identities["regulatory"]],
        "mosaic": [identities["mosaic"]],
        "regulatory": [identities["regulatory"]],
        "repeat": [identities["repeat"]],
        "sv_cnv": [identities["structural"]],
        "x_linked": [identities["x_linked"]],
    }
    return {
        "case_id": case_id,
        "hypotheses": [
            {"scope": scope, "members": members[scope]} for scope in REQUIRED_SCOPES
        ],
    }


def _source_format_bytes(
    input_bytes: bytes,
    *,
    mitochondrial_contigs: list[str],
    repeat_identity: str | None,
) -> bytes:
    return _json_bytes(
        {
            "schema": SOURCE_FORMAT_SCHEMA,
            "input_sha256": _sha256(input_bytes),
            "reference_role": REFERENCE_ROLE,
            "reference_identity_id": REFERENCE_IDENTITY,
            "reference_manifest_sha256": REFERENCE_MANIFEST_SHA256,
            "mitochondrial_contigs": mitochondrial_contigs,
            "repeat_identity": repeat_identity,
        }
    )


def _build_outputs(
    development: Mapping[str, object],
    held_out: Mapping[str, object],
    background: tuple[_VcfRecord, ...],
) -> Mapping[str, bytes]:
    development_seed = int(development["seed"])
    held_out_seed = int(held_out["seed"])
    if development_seed == held_out_seed:
        raise ValidationSeedError("validation seed split seeds must differ")
    development_records, development_truth = _synthetic_records(
        "development", development_seed
    )
    held_out_records, held_out_truth = _synthetic_records("held_out", held_out_seed)
    development_vcf = _vcf_bytes(
        tuple(
            sorted(
                _assign_opaque_ids(
                    (*background, *development_records),
                    split="development",
                    seed=development_seed,
                ),
                key=_record_order,
            )
        ),
        split="development",
    )
    held_out_vcf = _vcf_bytes(
        tuple(
            sorted(
                _assign_opaque_ids(
                    (*background, *held_out_records),
                    split="held_out",
                    seed=held_out_seed,
                ),
                key=_record_order,
            )
        ),
        split="held-out",
    )
    development_source_format = _source_format_bytes(
        development_vcf,
        mitochondrial_contigs=["M"],
        repeat_identity=SMALL_VARIANT_REPEAT_IDENTITY,
    )
    held_out_small_source_format = _source_format_bytes(
        held_out_vcf,
        mitochondrial_contigs=["M"],
        repeat_identity=None,
    )
    longtr_vcf, held_out_repeat_identity = _longtr_bytes(held_out_seed)
    held_out_truth = dict(held_out_truth)
    held_out_truth["repeat"] = held_out_repeat_identity
    longtr_source_format = _source_format_bytes(
        longtr_vcf,
        mitochondrial_contigs=[],
        repeat_identity=REPEAT_IDENTITY,
    )
    held_out_native = _json_bytes(
        {
            "schema": NATIVE_SCHEMA,
            "classification": "public",
            "validation_only": True,
            "case_id": "held-out-case",
            "sample_id": "proband",
            "reference_role": REFERENCE_ROLE,
            "reference_identity_id": REFERENCE_IDENTITY,
            "reference_manifest_sha256": REFERENCE_MANIFEST_SHA256,
            "members": [
                {
                    "role": "supplied-small-and-native-vcf",
                    "path": (MATERIALIZED_ROOT / "held-out-small-native.vcf").as_posix(),
                    "size_bytes": len(held_out_vcf),
                    "sha256": _sha256(held_out_vcf),
                    "sample_name": "SYNTHETIC",
                    "source_format_binding": {
                        "path": (
                            MATERIALIZED_ROOT
                            / "held-out-small.source-format.json"
                        ).as_posix(),
                        "size_bytes": len(held_out_small_source_format),
                        "sha256": _sha256(held_out_small_source_format),
                        "repeat_identity": None,
                    },
                },
                {
                    "role": "supplied-longtr-vcf",
                    "path": (MATERIALIZED_ROOT / "held-out-longtr.vcf").as_posix(),
                    "size_bytes": len(longtr_vcf),
                    "sha256": _sha256(longtr_vcf),
                    "sample_name": "SYNTHETIC",
                    "source_format_binding": {
                        "path": (
                            MATERIALIZED_ROOT / "held-out-longtr.source-format.json"
                        ).as_posix(),
                        "size_bytes": len(longtr_source_format),
                        "sha256": _sha256(longtr_source_format),
                        "repeat_identity": REPEAT_IDENTITY,
                    },
                },
            ],
            "limitations": [
                "This manifest preserves raw native members but is not yet a completed SuppliedSourceSetInputs replay artifact.",
                "Held-out execution remains blocked until development policy derivation and authenticated candidate-stack publication.",
            ],
        }
    )
    truth = _json_bytes(
        {
            "schema": TRUTH_SCHEMA,
            "classification": "public",
            "validation_only": True,
            "audience": "evaluator_only",
            "cases": [
                _truth_case("development-case", development_truth),
                _truth_case("held-out-case", held_out_truth),
            ],
            "limitations": [
                "These seed truths are synthetic recovery targets, not clinical labels.",
                "They cover the aggregate seed inputs but do not yet define model-specific cases or evaluator expectations for task 12.3.",
            ],
        }
    )
    return {
        "causal-truth.manifest.json": truth,
        "development.hpo.txt": (
            "\n".join(str(term) for term in development["hpo_terms"]) + "\n"
        ).encode("ascii"),
        "development.source-format.json": development_source_format,
        "development.vcf": development_vcf,
        "held-out-longtr.source-format.json": longtr_source_format,
        "held-out-longtr.vcf": longtr_vcf,
        "held-out-native.json": held_out_native,
        "held-out-small-native.vcf": held_out_vcf,
        "held-out-small.source-format.json": held_out_small_source_format,
        "held-out.hpo.txt": (
            "\n".join(str(term) for term in held_out["hpo_terms"]) + "\n"
        ).encode("ascii"),
    }


def _open_trusted_output_parent(parent: Path) -> int:
    """Open a caller-owned parent without following path-component symlinks."""

    flags = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW
    effective_uid = os.geteuid()
    descriptor = os.open("/", flags)
    try:
        for component in parent.parts[1:]:
            ancestor = os.fstat(descriptor)
            if ancestor.st_uid not in {0, effective_uid}:
                raise ValidationSeedError(
                    "validation seed output path has an untrusted owner"
                )
            writable = ancestor.st_mode & (stat.S_IWGRP | stat.S_IWOTH)
            if writable and not ancestor.st_mode & stat.S_ISVTX:
                raise ValidationSeedError(
                    "validation seed output path has an unsafe writable ancestor"
                )
            child_descriptor = os.open(component, flags, dir_fd=descriptor)
            child = os.fstat(child_descriptor)
            if writable and child.st_uid not in {0, effective_uid}:
                os.close(child_descriptor)
                raise ValidationSeedError(
                    "validation seed output path leaves a sticky trusted boundary"
                )
            os.close(descriptor)
            descriptor = child_descriptor
        metadata = os.fstat(descriptor)
        if metadata.st_uid != effective_uid or metadata.st_mode & (
            stat.S_IWGRP | stat.S_IWOTH
        ):
            raise ValidationSeedError(
                "validation seed output parent must be caller-owned and not writable by others"
            )
        return descriptor
    except Exception:
        os.close(descriptor)
        raise


def materialize_validation_seed(
    *,
    repository_root: str | Path,
    output_root: str | Path,
) -> ValidationSeedResult:
    """Generate the canonical public seed artifacts into a new directory."""

    try:
        root = Path(repository_root).resolve(strict=True)
        requested_output = Path(output_root)
        if (
            requested_output.name in {"", ".", ".."}
            or ".." in requested_output.parts
        ):
            raise ValidationSeedError("validation seed output path is invalid")
        output = (
            requested_output
            if requested_output.is_absolute()
            else Path.cwd() / requested_output
        )
    except OSError as error:
        raise ValidationSeedError("validation seed root is unavailable") from error

    development = _load_recipe(root, "development")
    held_out = _load_recipe(root, "held_out")
    _load_perturbations(root, "development")
    _load_perturbations(root, "held_out")
    background = _load_background(root)
    outputs = _build_outputs(development, held_out, background)
    if tuple(sorted(outputs)) != GENERATED_FILENAMES:
        raise ValidationSeedError("validation seed output contract is incomplete")

    directory_flags = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW
    try:
        parent_descriptor = _open_trusted_output_parent(output.parent)
    except ValidationSeedError:
        raise
    except OSError as error:
        raise ValidationSeedError(
            "validation seed output parent cannot be opened safely"
        ) from error
    try:
        try:
            os.mkdir(output.name, mode=0o700, dir_fd=parent_descriptor)
        except FileExistsError as error:
            raise ValidationSeedError(
                "validation seed output must be a new directory"
            ) from error
        try:
            output_descriptor = os.open(
                output.name,
                directory_flags,
                dir_fd=parent_descriptor,
            )
        except OSError as error:
            raise ValidationSeedError(
                "validation seed output directory cannot be opened safely"
            ) from error
        try:
            opened = os.fstat(output_descriptor)
            current = os.stat(
                output.name,
                dir_fd=parent_descriptor,
                follow_symlinks=False,
            )
            if (current.st_dev, current.st_ino) != (opened.st_dev, opened.st_ino):
                raise ValidationSeedError(
                    "validation seed output directory changed during creation"
                )
            file_flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW
            for filename in GENERATED_FILENAMES:
                file_descriptor = os.open(
                    filename,
                    file_flags,
                    0o600,
                    dir_fd=output_descriptor,
                )
                with os.fdopen(file_descriptor, "wb") as destination:
                    destination.write(outputs[filename])
                    destination.flush()
                    os.fsync(destination.fileno())
            os.fsync(output_descriptor)
            current = os.stat(
                output.name,
                dir_fd=parent_descriptor,
                follow_symlinks=False,
            )
            if (current.st_dev, current.st_ino) != (opened.st_dev, opened.st_ino):
                raise ValidationSeedError(
                    "validation seed output directory changed during publication"
                )
            os.fsync(parent_descriptor)
        except OSError as error:
            # Preserve the reserved path on failure rather than risk deleting a
            # concurrently replaced entry during cleanup.
            raise ValidationSeedError(
                "validation seed output could not be published without clobbering"
            ) from error
        finally:
            os.close(output_descriptor)
    finally:
        os.close(parent_descriptor)

    artifacts = tuple(
        GeneratedArtifact(
            filename=filename,
            size_bytes=len(outputs[filename]),
            sha256=_sha256(outputs[filename]),
        )
        for filename in GENERATED_FILENAMES
    )
    return ValidationSeedResult(
        generator_id="track1-synthetic-case-generator",
        generator_version="v2",
        output_root=str(output),
        artifacts=artifacts,
    )


def validation_seed_result_as_json(result: ValidationSeedResult) -> str:
    return json.dumps(asdict(result), indent=2, sort_keys=True) + "\n"


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser(
        description="Materialize the public Track 1 validation seed package"
    )
    parser.add_argument("repository_root")
    parser.add_argument("output_root")
    arguments = parser.parse_args()
    result = materialize_validation_seed(
        repository_root=arguments.repository_root,
        output_root=arguments.output_root,
    )
    print(validation_seed_result_as_json(result), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
