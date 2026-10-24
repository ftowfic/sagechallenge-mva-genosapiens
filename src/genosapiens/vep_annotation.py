"""Offline, all-transcript VEP evidence. No acquisition, scoring, or mechanisms.

The caller supplies inventoried evidence and an immutable container image pin.
Only genotype-free sequence/interval transports reach VEP; identities never change.
Public verification uses the complete task 4.2 cache inventory, not just info.txt.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict, dataclass, replace
from functools import lru_cache
import hashlib
import json
import os
from pathlib import Path
import re
import socket
import stat
import subprocess
import tempfile
from types import MappingProxyType
from typing import ClassVar
from urllib.parse import unquote

from genosapiens.allele_identity import (
    MitochondrialAllele, MosaicAllele, ReferenceRole, SequenceAllele, StructuralAllele, RepeatAllele,
)
from genosapiens.evidence_model import (
    AlleleEvidence, IdentityMapping, IdentityStage, SourcedRecord, TranscriptConsequence,
)
from genosapiens.evidence_serialization import DUCKDB_VERSION, _plain, to_json
from genosapiens.governance import ArtifactClassification
from genosapiens.grch38_reference import IndexedFasta
from genosapiens.llm_guard import LLMPolicy
from genosapiens import open_targets as ot
from genosapiens.reference_concordance import load_reference_identity, ReferenceRoleBinding
from genosapiens.submission_rendering import load_submission_contig_renderer
from genosapiens.variant_normalization import (
    SmallVariantNormalizationResult, VARIANT_NORMALIZATION_SCHEMA,
    variant_normalization_as_json,
)
from genosapiens.vep_bundle import (
    VepSentinelResult, _bundle_data, _bundle_manifest_path, _file_identity,
    _manifest_path,
    _require_within_root, _validate_release, _verify_resource_manifests,
    load_vep_bundle_contract,
)


SCHEMA = "genosapiens.vep-annotation"
CONSTRAINT_SCHEMA = "genosapiens.open-targets-constraint-attachment"
EXECUTABLE = Path("/opt/vep/src/ensembl-vep/vep")
FLAGS = (
    "--cache", "--offline", "--cache_version", "116", "--species", "homo_sapiens",
    "--assembly", "GRCh38", "--format", "vcf", "--vcf", "--no_stats",
    "--hgvs", "--hgvsg", "--protein", "--protein_version", "--gene_version",
    "--transcript_version", "--symbol", "--numbers", "--domains", "--canonical",
    "--mane", "--biotype", "--uniprot", "--tsl", "--appris", "--variant_class",
    "--sift", "b", "--polyphen", "b", "--allele_number", "--check_ref",
    "--use_given_ref", "--terms", "SO", "--no_progress",
)
CORE_FIELDS = frozenset(("Allele", "Consequence", "IMPACT", "Gene", "Feature_type", "Feature", "ALLELE_NUM"))
RICH_FIELDS = CORE_FIELDS | frozenset((
    "SYMBOL", "SYMBOL_SOURCE", "HGNC_ID", "BIOTYPE", "EXON", "INTRON", "HGVSc",
    "HGVSp", "HGVSg", "Protein_position", "Amino_acids", "FLAGS", "VARIANT_CLASS",
    "CANONICAL", "MANE", "MANE_SELECT", "MANE_PLUS_CLINICAL", "TSL", "APPRIS",
    "ENSP", "SWISSPROT", "TREMBL", "UNIPARC", "UNIPROT_ISOFORM", "SIFT", "PolyPhen", "DOMAINS",
))
_SHA = re.compile(r"[0-9a-f]{64}")
_NAME = re.compile(r"[A-Za-z_][A-Za-z0-9_.]*")
_COLUMNS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"


class VepAnnotationError(RuntimeError):
    """Static diagnostics only. Never include an input, path, or subprocess log."""


class _Controlled:
    __slots__ = ()
    classification: ClassVar[ArtifactClassification] = ArtifactClassification.CONTROLLED
    llm_policy: ClassVar[LLMPolicy] = LLMPolicy.NO_LLM

    def __repr__(self):
        return "<controlled VEP annotation; no-llm>"


def _digest(content: str) -> str:
    return hashlib.sha256(content.encode("utf-8")).hexdigest()


def _json(data) -> str:
    return json.dumps(data, sort_keys=True, separators=(",", ":"), ensure_ascii=True) + "\n"


def _sequence(identity):
    if type(identity) is SequenceAllele:
        return identity
    if type(identity) in (MitochondrialAllele, MosaicAllele):
        return identity.sequence
    # A structural/repeat exact_sequence is only nested evidence, NOT equivalence.
    return None


@dataclass(frozen=True, slots=True, repr=False)
class VepRequest(_Controlled):
    allele: AlleleEvidence

    @property
    def identity_id(self):
        return self.allele.native_identity_id

    def __post_init__(self):
        try:
            if type(self.allele) is not AlleleEvidence:
                raise ValueError
            identities = {item.identity_id: item for item in self.allele.identities}
            selected = identities[self.identity_id]
            if selected.stage is not IdentityStage.NATIVE:
                raise ValueError
        except (ValueError, TypeError, KeyError, AttributeError):
            raise VepAnnotationError("invalid VEP identity selection") from None

    @classmethod
    def from_normalization(cls, result: SmallVariantNormalizationResult, allele_id: str):
        """Select an existing normalized allele, retaining every contributing source.

        Does not process quarantines or invent class designations, GT, or rescue
        authorization. For caller-designated classes use existing AlleleEvidence.
        """
        try:
            if type(result) is not SmallVariantNormalizationResult:
                raise ValueError
            result.__post_init__()
            normalized = next(a for a in result.normalized_alleles if a.normalized_allele_id == allele_id)
            mappings = tuple(m for m in result.source_to_normalized if m.normalized_allele_id == allele_id)
            source_sha = _digest(variant_normalization_as_json(result))
            role = ReferenceRole(normalized.reference_role)
            sequence = SequenceAllele(role, normalized.reference_identity_id, normalized.contig,
                                      normalized.position, normalized.reference, normalized.alternate)
            source = SourcedRecord(
                "genosapiens.variant-normalization", "current",
                allele_id, source_sha, None,
                {"normalized": {**asdict(normalized), "kind": normalized.kind.value},
                 "source_to_normalized": [asdict(m) for m in mappings],
                 "input_sha256": result.input_sha256,
                 "reference_concordance_sha256": result.reference_concordance_sha256},
            )
            identities = []
            for mapping in mappings:
                original = mapping.source
                record = original.record
                identities.append(IdentityMapping(
                    original.source_allele_id, IdentityStage.SOURCE,
                    SequenceAllele(role, normalized.reference_identity_id, record.contig,
                                   record.position, record.reference, original.alternate),
                    result.reference_manifest_sha256,
                    replace(source, source_record_id=original.source_allele_id,
                            source_fields=asdict(original)),
                ))
            identities.append(IdentityMapping(
                allele_id, IdentityStage.NORMALIZED, sequence, result.reference_manifest_sha256,
                source, tuple(i.identity_id for i in identities), VARIANT_NORMALIZATION_SCHEMA,
            ))
            native_id = "native-" + allele_id
            identities.append(IdentityMapping(native_id, IdentityStage.NATIVE, sequence,
                                              result.reference_manifest_sha256, source, (allele_id,), SCHEMA))
            return cls(AlleleEvidence(allele_id, native_id, tuple(identities)))
        except (ValueError, TypeError, RuntimeError, AttributeError, StopIteration):
            raise VepAnnotationError("invalid normalized VEP input") from None


@dataclass(frozen=True, slots=True, repr=False)
class VepVcfRecord(_Controlled):
    """Expected genotype-free transport row; ALT indexes are local, not source GT."""

    record_id: str
    contig: str
    position: int
    reference: str
    alternates: tuple[str, ...]

    def __post_init__(self):
        if (type(self.record_id) is not str or re.fullmatch(r"[A-Za-z0-9_.-]+", self.record_id) is None
                or type(self.contig) is not str or re.fullmatch(r"[A-Za-z0-9_.-]+", self.contig) is None
                or type(self.position) is not int or self.position < 1
                or type(self.reference) is not str or re.fullmatch(r"[ACGTN]+", self.reference) is None
                or type(self.alternates) is not tuple or not self.alternates
                or any(type(a) is not str or re.fullmatch(r"[ACGTN]+", a) is None or a == self.reference
                       for a in self.alternates) or len(set(self.alternates)) != len(self.alternates)):
            raise VepAnnotationError("invalid VEP transport identity")

    def line(self):
        return f"{self.contig}\t{self.position}\t{self.record_id}\t{self.reference}\t{','.join(self.alternates)}\t.\t.\t.\n"

    def csq_allele(self, number: int):
        # Exact pinned Parser/VCF.pm default indel conversion, without --minimal.
        alleles = (self.reference, *self.alternates)
        trim = any(len(a) != len(self.reference) for a in self.alternates) and len({a[0] for a in alleles}) == 1
        return (self.alternates[number - 1][1:] or "-") if trim else self.alternates[number - 1]


@dataclass(frozen=True, slots=True, repr=False)
class ParsedVep(_Controlled):
    fields: tuple[str, ...]
    transcripts: Mapping[tuple[str, int], tuple[TranscriptConsequence, ...]]
    non_transcripts: Mapping[tuple[str, int], tuple[SourcedRecord, ...]]


@dataclass(frozen=True, slots=True)
class ParsedVepHeader:
    text: str
    fields: tuple[str, ...]
    info_headers: frozenset[str]


def parse_vep_header(text: str, *, required_fields: frozenset[str] = CORE_FIELDS) -> ParsedVepHeader:
    """Validate a complete VEP header once for bounded row parsing."""
    try:
        lines = text.splitlines()
        if (not text.endswith("\n") or not lines
                or lines[0] not in ("##fileformat=VCFv4.2", "##fileformat=VCFv4.3")
                or lines[-1] != _COLUMNS):
            raise ValueError
        fields = None
        info_headers = set()
        for line in lines[1:-1]:
            if not line or any(ord(c) < 32 and c != "\t" for c in line):
                raise ValueError
            if not line.startswith("##") or line.startswith(("##fileformat=", "##FORMAT=")):
                raise ValueError
            if line.startswith("##INFO="):
                header = re.fullmatch(r'##INFO=<ID=([A-Za-z_][A-Za-z0-9_.]*),Number=([^,]+),Type=([^,]+),Description="((?:[^"\\]|\\.)*)">', line)
                if header is None or header[1] in info_headers:
                    raise ValueError
                tag, number, dtype, description = header.groups()
                info_headers.add(tag)
                if tag == "CSQ":
                    schema = re.fullmatch(r".*?Format: ([A-Za-z0-9_.|]+)", description)
                    if number != "." or dtype != "String" or schema is None:
                        raise ValueError
                    fields = tuple(schema[1].split("|"))
                    if (len(set(fields)) != len(fields)
                            or any(_NAME.fullmatch(field) is None for field in fields)
                            or not required_fields <= set(fields)):
                        raise ValueError
        if fields is None:
            raise ValueError
        return ParsedVepHeader(text, fields, frozenset(info_headers))
    except (TypeError, ValueError, AttributeError):
        raise VepAnnotationError("invalid or incomplete VEP CSQ output") from None


def _value(value):
    if value in (None, "", "."):
        return None
    if re.search(r"%(?![0-9A-Fa-f]{2})", value):
        raise ValueError
    decoded = unquote(value, encoding="utf-8", errors="strict")
    if any(ord(c) < 32 or ord(c) == 127 for c in decoded):
        raise ValueError
    return decoded


def _transcript(raw, source, identifier):
    values = {key: _value(value) for key, value in raw.items()}
    def versioned(key, prefix):
        value = values.get(key)
        if value is None:
            return None
        if re.fullmatch(prefix + r"[0-9]{11}(?:\.[1-9][0-9]*)?", value) is None:
            raise ValueError
        return value if "." in value else None

    position = values.get("Protein_position")
    start = end = None
    if position is not None:
        match = re.fullmatch(r"([1-9][0-9]*|\?)(?:-([1-9][0-9]*|\?))?(?:/([1-9][0-9]*|\?))?", position)
        if match is None:
            raise ValueError
        left, right, total = match.groups()
        right = right or left
        start, end = (int(v) if v != "?" else None for v in (left, right))
        if start is not None and end is not None and start > end:
            raise ValueError
        if total not in (None, "?") and any(v is not None and v > int(total) for v in (start, end)):
            raise ValueError
    canonical = None
    if "CANONICAL" in raw:
        if raw["CANONICAL"] not in ("YES", "", "."):
            raise ValueError
        canonical = None if raw["CANONICAL"] == "." else raw["CANONICAL"] == "YES"
    for key in ("EXON", "INTRON"):
        if values.get(key) and re.fullmatch(r"[1-9][0-9]*(?:-[1-9][0-9]*)?/[1-9][0-9]*", values[key]) is None:
            raise ValueError
    for key in ("SIFT", "PolyPhen"):
        if values.get(key):
            prediction = re.fullmatch(r"[a-z_]+\(([0-9]+(?:\.[0-9]+)?(?:e-[0-9]+)?)\)", values[key])
            if prediction is None or not 0 <= float(prediction[1]) <= 1:
                raise ValueError
    status = []
    if canonical:
        status.append("CANONICAL")
    for key in ("MANE", "TSL", "APPRIS", "FLAGS"):
        if values.get(key):
            status.extend(v if key == "MANE" else key + ":" + v for v in values[key].split("&"))
    mane = tuple(v for key in ("MANE_SELECT", "MANE_PLUS_CLINICAL") if (v := values.get(key)))
    for key, label in (("MANE_SELECT", "MANE_Select"), ("MANE_PLUS_CLINICAL", "MANE_Plus_Clinical")):
        if values.get(key) and label not in status:
            status.append(label)
    def sourced(key, value, index=0):
        return replace(source, source_record_id=f"{identifier}:{key}:{index}",
                       source_fields={"csq_field": key, "raw_value": value})
    domains = values.get("DOMAINS")
    predictions = tuple(sourced(key, raw[key]) for key in ("IMPACT", "SIFT", "PolyPhen") if values.get(key))
    accessions = tuple(value for key in ("SWISSPROT", "TREMBL")
                       if values.get(key) for value in values[key].split("&"))
    # VEP emits accession.entry-version; preserve that original in raw CSQ.
    if any(re.fullmatch(r"[A-Z0-9]+(?:\.[0-9]+)?", v) is None for v in accessions):
        raise ValueError
    accessions = tuple(v.split(".", 1)[0] for v in accessions)
    isoforms = values.get("UNIPROT_ISOFORM")
    hgnc_id = values.get("HGNC_ID")
    if hgnc_id is not None and re.fullmatch(r"(?:HGNC:)?[1-9][0-9]*", hgnc_id) is None:
        raise ValueError
    return TranscriptConsequence(
        identifier, source, ensembl_gene_id=versioned("Gene", "ENSG"),
        ensembl_transcript_id=versioned("Feature", "ENST"), ensembl_protein_id=versioned("ENSP", "ENSP"),
        hgnc_symbol=values.get("SYMBOL") if values.get("SYMBOL_SOURCE") == "HGNC" else None,
        hgnc_id=(hgnc_id if hgnc_id.startswith("HGNC:") else "HGNC:" + hgnc_id) if hgnc_id else None,
        mane_refseq=mane[0] if mane and len(set(mane)) == 1 else None,
        uniprot_accessions=accessions or None,
        uniprot_isoforms=tuple(isoforms.split("&")) if isoforms else None,
        hgvs_genomic=values.get("HGVSg"), hgvs_coding=values.get("HGVSc"), hgvs_protein=values.get("HGVSp"),
        exon=values.get("EXON"), intron=values.get("INTRON"), protein_position=position,
        protein_start=start, protein_end=end, amino_acids=values.get("Amino_acids"),
        domain_overlaps=tuple(sourced("DOMAINS", v, i) for i, v in enumerate(domains.split("&"))) if domains else None,
        transcript_status=tuple(status) if status else None, canonical=canonical,
        variant_class=values.get("VARIANT_CLASS"), consequences=tuple(values["Consequence"].split("&")) if values.get("Consequence") else None,
        biotype=values.get("BIOTYPE"), predictions=predictions or None,
        conservation=None, constraint=None,
    )


def parse_vep_vcf(
    content: str, *, expected: tuple[VepVcfRecord, ...], source: SourcedRecord,
    required_fields: frozenset[str] = CORE_FIELDS,
    parsed_header: ParsedVepHeader | None = None,
) -> ParsedVep:
    """Strict CSQ parser for unfiltered, non-minimal, genotype-free VEP VCF.

    Unknown declared CSQ fields survive verbatim. Missing required schema fields,
    absent ALT assessments, malformed rows/tags, changed alleles and partial CSQs
    fail closed. No padding, ordinal guessing, worst-transcript selection or GT.
    The supplied source SHA-256 must bind the exact UTF-8 output bytes.
    """
    try:
        if (type(content) is not str or type(expected) is not tuple or type(source) is not SourcedRecord
                or source.source_sha256 != _digest(content) or not CORE_FIELDS <= required_fields):
            raise ValueError
        records = {r.record_id: r for r in expected}
        if len(records) != len(expected) or any(type(r) is not VepVcfRecord for r in expected):
            raise ValueError
        transcripts = {(r.record_id, i): [] for r in expected for i in range(1, len(r.alternates) + 1)}
        non_transcripts = {key: [] for key in transcripts}
        fields = parsed_header.fields if parsed_header is not None else None
        seen_records, seen_keys = set(), set()
        info_headers = set(parsed_header.info_headers) if parsed_header is not None else set()
        columns = parsed_header is not None
        if parsed_header is not None:
            if (type(parsed_header) is not ParsedVepHeader
                    or not content.startswith(parsed_header.text)
                    or not required_fields <= set(parsed_header.fields)):
                raise ValueError
            lines = content[len(parsed_header.text):].splitlines()
        else:
            lines = content.splitlines()
            if not lines or lines[0] not in ("##fileformat=VCFv4.2", "##fileformat=VCFv4.3"):
                raise ValueError
            lines = lines[1:]
        for line in lines:
            if not line or any(ord(c) < 32 and c != "\t" for c in line):
                raise ValueError
            if line.startswith("##"):
                if columns or line.startswith(("##fileformat=", "##FORMAT=")):
                    raise ValueError
                if line.startswith("##INFO="):
                    header = re.fullmatch(r'##INFO=<ID=([A-Za-z_][A-Za-z0-9_.]*),Number=([^,]+),Type=([^,]+),Description="((?:[^"\\]|\\.)*)">', line)
                    if header is None or header[1] in info_headers:
                        raise ValueError
                    tag, number, dtype, description = header.groups()
                    info_headers.add(tag)
                    if tag == "CSQ":
                        schema = re.fullmatch(r".*?Format: ([A-Za-z0-9_.|]+)", description)
                        if number != "." or dtype != "String" or schema is None:
                            raise ValueError
                        fields = tuple(schema[1].split("|"))
                        if (len(set(fields)) != len(fields) or any(_NAME.fullmatch(f) is None for f in fields)
                                or not required_fields <= set(fields)):
                            raise ValueError
                continue
            if line.startswith("#"):
                if columns or line != _COLUMNS or fields is None:
                    raise ValueError
                columns = True
                continue
            if not columns:
                raise ValueError
            row = line.split("\t")
            if len(row) != 8 or row[2] not in records or row[2] in seen_records:
                raise ValueError
            expected_row = records[row[2]]
            if row[:7] != expected_row.line().rstrip("\n").split("\t")[:7]:
                raise ValueError
            seen_records.add(row[2])
            tags = {}
            for item in row[7].split(";"):
                pair = item.split("=", 1)
                if len(pair) != 2 or pair[0] in tags or pair[0] not in info_headers or not pair[1]:
                    raise ValueError
                tags[pair[0]] = pair[1]
            if set(tags) != {"CSQ"}:
                raise ValueError
            for ordinal, raw_csq in enumerate(tags["CSQ"].split(",")):
                parts = raw_csq.split("|")
                if len(parts) != len(fields):
                    raise ValueError
                raw = dict(zip(fields, parts, strict=True))
                # VEP only emits CHECK_REF when checking failed. Never turn a
                # failed assessment into ordinary evidence via raw-field retention.
                if _value(raw.get("CHECK_REF")) is not None:
                    raise ValueError
                if re.fullmatch(r"[1-9][0-9]*", raw["ALLELE_NUM"]) is None:
                    raise ValueError
                number = int(raw["ALLELE_NUM"])
                key = (row[2], number)
                if key not in transcripts or raw["Allele"] != expected_row.csq_allele(number):
                    raise ValueError
                if (not raw["Consequence"] or any(re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_]*", v) is None
                                                  for v in raw["Consequence"].split("&"))
                        or raw["IMPACT"] not in ("HIGH", "MODERATE", "LOW", "MODIFIER")):
                    raise ValueError
                feature_type = raw["Feature_type"]
                if feature_type not in ("Transcript", "RegulatoryFeature", "MotifFeature", "", "."):
                    raise ValueError
                if feature_type == "Transcript" and not _value(raw["Feature"]):
                    raise ValueError
                if feature_type in ("", ".") and (raw["Consequence"] != "intergenic_variant"
                                                     or _value(raw["Feature"]) or _value(raw["Gene"])):
                    raise ValueError
                for value in raw.values():
                    _value(value)
                seen_keys.add(key)
                identifier = "vep-csq-" + _digest(_json([
                    source.source_dataset, source.source_release, source.public_manifest_sha256,
                    asdict(expected_row), ordinal, raw,
                ]))
                bound_source = replace(source, source_record_id=identifier, source_fields={
                    "annotation_binding": source.source_fields, "vcf_record_id": row[2],
                    "allele_number": number, "csq_ordinal": ordinal, "csq_fields": list(fields),
                    "raw_csq": raw_csq, "csq": raw,
                })
                if feature_type == "Transcript":
                    transcripts[key].append(_transcript(raw, bound_source, identifier))
                else:
                    non_transcripts[key].append(bound_source)
        if not columns or seen_records != set(records) or seen_keys != set(transcripts):
            raise ValueError
        return ParsedVep(fields,
                         MappingProxyType({key: tuple(v) for key, v in transcripts.items()}),
                         MappingProxyType({key: tuple(v) for key, v in non_transcripts.items()}))
    except (ValueError, TypeError, KeyError, IndexError, AttributeError, RuntimeError, OverflowError):
        raise VepAnnotationError("invalid or incomplete VEP CSQ output") from None


@dataclass(frozen=True, slots=True, repr=False)
class VepAssessment(_Controlled):
    allele: AlleleEvidence
    requested_identity_id: str
    resource_allele: SequenceAllele | None
    status: str
    reason: str | None
    non_transcript_consequences: tuple[SourcedRecord, ...] = ()
    native_geometry_assessments: tuple[SourcedRecord, ...] = ()


@dataclass(frozen=True, slots=True, repr=False)
class VepAnnotationResult(_Controlled):
    assessments: tuple[VepAssessment, ...]
    provenance: SourcedRecord

    @property
    def alleles(self):
        return tuple(item.allele for item in self.assessments)


def _readonly(path: Path):
    if not os.statvfs(path).f_flag & os.ST_RDONLY:
        raise VepAnnotationError("VEP resources require read-only mounts")


def _local(path):
    path = Path(path)
    if not path.is_absolute() or ".." in path.parts:
        raise VepAnnotationError("VEP paths require absolute non-symlink local paths")
    # Inspect ancestors in root-to-leaf order: even lstat of a leaf would follow
    # an unchecked intermediate directory symlink.
    for component in (*reversed(path.parents), path):
        try:
            mode = component.lstat().st_mode
        except FileNotFoundError:
            break  # New output paths are permitted; input preflight requires existence.
        if stat.S_ISLNK(mode) or component != path and not stat.S_ISDIR(mode):
            raise VepAnnotationError("VEP paths require absolute non-symlink local paths")
    return path


def _preflight_resources(files=(), trees=()):
    """Check topology without following links, then actual mount flags, before I/O.

    Never hash/read an inventory to discover whether the surrounding tree is
    safe. Directory traversal itself uses no-follow entry metadata.
    """
    checked = set()
    pending = []
    for path in trees:
        path = _local(path)
        if not stat.S_ISDIR(path.lstat().st_mode):
            raise VepAnnotationError("invalid annotation resource tree")
        pending.append(path)
    for path in files:
        path = _local(path)
        if not stat.S_ISREG(path.lstat().st_mode):
            raise VepAnnotationError("invalid annotation resource file")
        checked.add(path)
    while pending:
        directory = pending.pop()
        _local(directory)
        checked.add(directory)
        with os.scandir(directory) as entries:
            for entry in entries:
                mode = entry.stat(follow_symlinks=False).st_mode
                path = Path(entry.path)
                if stat.S_ISDIR(mode):
                    pending.append(path)
                elif not stat.S_ISREG(mode):
                    raise VepAnnotationError("invalid annotation resource tree entry")
                checked.add(path)
    for path in sorted(checked):
        _local(path)
        _readonly(path)


@dataclass(frozen=True, slots=True, repr=False)
class VepReferenceCompatibility(_Controlled):
    """Pinned coordinate contract, not a case/reference acceptance assertion.

    The optional reconstructed bundle supplies verified naming metadata for its
    three roles. Without it, those source identities remain unassessed; the
    independently verified VEP-native reference remains usable.
    """

    identity_path: Path
    identity_sha256: str
    reference_root: Path | None = None
    manifest_path: Path | None = None
    manifest_sha256: str | None = None

    def __post_init__(self):
        try:
            object.__setattr__(self, "identity_path", _local(self.identity_path))
            if type(self.identity_sha256) is not str or _SHA.fullmatch(self.identity_sha256) is None:
                raise ValueError
            supplied = (self.reference_root, self.manifest_path, self.manifest_sha256)
            if any(v is not None for v in supplied):
                if any(v is None for v in supplied) or _SHA.fullmatch(self.manifest_sha256) is None:
                    raise ValueError
                object.__setattr__(self, "reference_root", _local(self.reference_root))
                object.__setattr__(self, "manifest_path", _local(self.manifest_path))
        except (ValueError, TypeError, OSError):
            raise VepAnnotationError("invalid source reference compatibility binding") from None


def _source_coordinate_maps(
    binding, vep_contract, vep_contract_sha256, vep_root, reference
):
    if type(binding) is not VepReferenceCompatibility:
        raise VepAnnotationError("source reference compatibility binding is required")
    binding.__post_init__()
    _preflight_resources(files=(binding.identity_path,))
    if _file_identity(binding.identity_path)[0] != binding.identity_sha256:
        raise VepAnnotationError("source reference identity contract pin mismatch")
    contract = load_reference_identity(binding.identity_path)
    if (contract.primary_assembly_accession != "GCA_000001405.15"
            or contract.assembly != "GRCh38-masked-v2-no-alt-plus-hs38d1"):
        raise VepAnnotationError("source reference coordinate assembly mismatch")
    for canonical, length in contract.canonical_contig_lengths.items():
        native_name = "MT" if canonical == "M" else canonical
        if native_name not in reference.index or reference.index[native_name].length != length:
            raise VepAnnotationError("source coordinate contract and VEP dictionary disagree")
    fasta_manifest_path = _manifest_path(
        Path(vep_root),
        vep_contract.fasta.source.resource_id,
        vep_contract_sha256,
    )
    fasta_manifest_sha = _file_identity(fasta_manifest_path)[0]
    maps = {(ReferenceRole.PUBLIC_ANNOTATION_BENCHMARK.value, vep_contract.fasta.source.resource_id):
            (fasta_manifest_sha, dict(resource_native_aliases(reference)))}
    proof = {"identity_contract_sha256": binding.identity_sha256, "identity_contract_id": contract.contract_id,
             "primary_assembly_accession": contract.primary_assembly_accession,
             "vep_fasta_manifest_sha256": fasta_manifest_sha, "source_manifest_sha256": binding.manifest_sha256,
             "scope": "coordinate-compatibility-only; not case acceptance or sequence identity", "source_naming": []}
    if binding.manifest_path is not None:
        metadata_paths = tuple(binding.reference_root / contract.expected_output(role).path for role in (
            "contig_aliases", *(r.fai_output_role for r in contract.reference_roles),
        ))
        _preflight_resources(files=(binding.manifest_path, *metadata_paths))
        if _file_identity(binding.manifest_path)[0] != binding.manifest_sha256:
            raise VepAnnotationError("source reference manifest pin mismatch")
        for role in contract.reference_roles:
            # Reuse the existing verified name-map loader, not its rendering or
            # submission logic. This authenticates manifest, aliases and role FAI.
            names = load_submission_contig_renderer(
                identity_path=binding.identity_path, reference_manifest_path=binding.manifest_path,
                reference_root=binding.reference_root,
                reference_binding=ReferenceRoleBinding(role.role, role.identity_id),
            )
            aliases = {}
            for source_name, (canonical, length) in names.accepted_names.items():
                native_name = "MT" if canonical == "M" else canonical
                if native_name in reference.index and reference.index[native_name].length == length:
                    aliases[source_name] = native_name
            maps[(role.role, role.identity_id)] = (binding.manifest_sha256, aliases)
            proof["source_naming"].append({"role": role.role, "identity_id": role.identity_id,
                "aliases_sha256": names.contig_aliases_sha256, "aliases": aliases})
        _preflight_resources(files=(binding.identity_path, binding.manifest_path, *metadata_paths))
        if (_file_identity(binding.manifest_path)[0] != binding.manifest_sha256
                or _file_identity(binding.identity_path)[0] != binding.identity_sha256):
            raise VepAnnotationError("source reference compatibility inputs changed")
    return maps, proof


@dataclass(frozen=True, slots=True, repr=False)
class OpenTargetsConstraintReference(_Controlled):
    """Explicit local metadata/target view and the frozen full task 4.4 manifest.

    release_directory contains the three metadata files and output/target only;
    other datasets need not be mounted and are never opened by this adapter.
    """

    contract_path: Path
    release_directory: Path
    manifest_path: Path
    manifest_sha256: str

    def __post_init__(self):
        try:
            for name in ("contract_path", "release_directory", "manifest_path"):
                object.__setattr__(self, name, _local(getattr(self, name)))
            if type(self.manifest_sha256) is not str or _SHA.fullmatch(self.manifest_sha256) is None:
                raise ValueError
        except (ValueError, TypeError, OSError):
            raise VepAnnotationError("invalid Open Targets constraint reference") from None


def attach_open_targets_constraints(
    transcripts: tuple[TranscriptConsequence, ...], *, reference: OpenTargetsConstraintReference,
) -> tuple[tuple[TranscriptConsequence, ...], SourcedRecord]:
    """Verify all canonical target partitions and append source-separated evidence.

    No scoring, row deduplication, constraint-type selection, or source-version
    inference. Only the explicit Ensembl stable ID is joined; an existing VEP
    identifier version remains untouched. Missing IDs/matches and null/empty
    source collections are audited without inventing a constraint measurement.
    """
    try:
        import duckdb

        if (type(reference) is not OpenTargetsConstraintReference or type(transcripts) is not tuple
                or any(type(t) is not TranscriptConsequence for t in transcripts)
                or duckdb.__version__ != DUCKDB_VERSION
                or any(name != "lo" for _, name in socket.if_nameindex())):
            raise ValueError
        reference.__post_init__()
        release = reference.release_directory
        _preflight_resources(files=(reference.contract_path, reference.manifest_path,
                                    *(release / name for name in ot.METADATA_NAMES)),
                             trees=(release / "output/target",))
        if _file_identity(reference.manifest_path)[0] != reference.manifest_sha256:
            raise ValueError
        contract = ot.load_open_targets_contract(reference.contract_path)
        if contract.resource != "open-targets-platform" or contract.release != "26.06" or contract.assembly != "GRCh38":
            raise ValueError
        contract_sha = _file_identity(reference.contract_path)[0]
        metadata = tuple(ot.ReleaseArtifact(None, m.name, m.size_bytes, m.sha1, m.sha256) for m in contract.metadata)
        for artifact in metadata:
            path = _local(release / artifact.relative_path)
            if not path.is_file():
                raise ValueError
            _readonly(path)
        ot._verify_artifacts(release, metadata)
        # The full contract/model authenticates the full manifest's identities.
        # Only the target view below is opened, never any other dataset partition.
        artifacts, models, _joins = ot._release_model(release, contract)
        manifest = ot._load_manifest(reference.manifest_path)
        ot._verify_manifest_contract(manifest, contract, contract_sha, artifacts)
        outputs = {o.path: o for o in manifest.outputs}
        selected = []
        for artifact in artifacts:
            if artifact.dataset not in (None, "target"):
                continue
            output = outputs[f"{contract.output_prefix}/{artifact.relative_path}"]
            if output.size_bytes != artifact.size_bytes or artifact.sha256 not in (None, output.sha256):
                raise ValueError
            _readonly(_local(release / artifact.relative_path))
            selected.append(replace(artifact, sha256=output.sha256))
        partitions = tuple(a for a in selected if a.dataset == "target")
        target_contract = next(d for d in contract.datasets if d.name == "target")
        target_model = next(m for m in models if m.name == "target")
        if len(partitions) != target_contract.partition_count or target_model.key_fields != ("target/id",):
            raise ValueError
        model_fields = {f.field_id: f for f in target_model.fields}
        if (model_fields["target/id"].data_type != "sc:Text" or model_fields["target/id"].repeated_depth != 0
                or model_fields["target/constraint"].repeated_depth != 1):
            raise ValueError
        ot._verify_release_tree(release / "output/target", tuple(
            replace(a, relative_path=str(Path(a.relative_path).relative_to("output/target"))) for a in partitions
        ))
        verified = ot._verify_artifacts(release, selected)
        validation = ot._validate_parquet(
            release, replace(contract, datasets=(target_contract,), semantic_joins=()), verified, (target_model,),
        )[0]
        croissant = ot._json_object(release / "croissant.json", "Open Targets metadata")
        target_definition = next(r for r in croissant["recordSet"] if r["name"] == "target")
        definition = next(f for f in target_definition["field"] if f["name"] == "constraint")
        definition_sha = _digest(_json(definition))
        gene_ids = []
        for transcript in transcripts:
            raw_fields = transcript.source.source_fields or {}
            csq = raw_fields.get("csq", {})
            if isinstance(csq, Mapping) and _value(csq.get("CHECK_REF")) is not None:
                raise ValueError
            raw_gene = _value(csq.get("Gene")) if isinstance(csq, Mapping) else None
            if raw_fields.get("annotation_scope") in ("structural_effect", "transcript_context"):
                native_raw = raw_fields.get("raw_record", {})
                native_gene = native_raw.get("gene_id") if isinstance(native_raw, Mapping) else None
                if native_gene in ("", ".", "-"):
                    native_gene = None
                if native_gene is not None:
                    if type(native_gene) is not str or raw_gene is not None and raw_gene != native_gene:
                        raise ValueError
                    raw_gene = native_gene
            typed_gene = transcript.ensembl_gene_id
            if raw_gene and typed_gene and raw_gene != typed_gene:
                # No reconciliation of conflicting gene or version identities.
                raise ValueError
            gene = raw_gene or typed_gene
            match = re.fullmatch(r"(ENSG[0-9]{11})(?:\.[1-9][0-9]*)?", gene) if gene else None
            if gene and match is None:
                raise ValueError
            gene_ids.append(match[1] if match else None)
        rows_by_gene = {gene: [] for gene in sorted({g for g in gene_ids if g is not None})}
        connection = duckdb.connect(database=":memory:", config={
            "threads": "1", "memory_limit": "1GB", "preserve_insertion_order": "true",
            "autoinstall_known_extensions": "false", "autoload_known_extensions": "false",
        })
        try:
            for artifact in partitions:
                # file_row_number is the physical zero-based Parquet row ordinal,
                # computed before filtering; equal id/value rows stay distinct.
                rows = connection.execute(
                    'SELECT file_row_number, id, "constraint" FROM read_parquet(?, file_row_number=true) '
                    'WHERE id IN (SELECT unnest(?::VARCHAR[])) ORDER BY file_row_number',
                    [str(release / artifact.relative_path), list(rows_by_gene)],
                ).fetchall()
                for ordinal, gene, values in rows:
                    if values is not None and type(values) is not list:
                        raise ValueError
                    record = SourcedRecord(
                        "Open Targets target", contract.release, f"{artifact.relative_path}#row={ordinal}",
                        artifact.sha256, reference.manifest_sha256,
                        {"id": gene, "constraint": values, "file_row_number": ordinal,
                         "artifact_path": f"{contract.output_prefix}/{artifact.relative_path}",
                         "contract_sha256": contract_sha, "constraint_field_schema_sha256": definition_sha},
                    )
                    rows_by_gene[gene].append(record)
        finally:
            connection.close()
        # Read-only mounts prevent writes from this process. Rehash to reject
        # externally replaced artifacts during schema validation or lookup too.
        _preflight_resources(files=(reference.contract_path, reference.manifest_path,
                                    *(release / a.relative_path for a in verified)),
                             trees=(release / "output/target",))
        ot._verify_artifacts(release, verified)
        if (_file_identity(reference.manifest_path)[0] != reference.manifest_sha256
                or _file_identity(reference.contract_path)[0] != contract_sha):
            raise ValueError
        attached = []
        for transcript, gene in zip(transcripts, gene_ids, strict=True):
            records = rows_by_gene.get(gene, ())
            populated = tuple(r for r in records if r.source_fields["constraint"])
            constraint = transcript.constraint
            if populated:
                constraint = (constraint or ()) + populated
            elif records and all(r.source_fields["constraint"] == () for r in records) and constraint is None:
                constraint = ()
            attached.append(replace(transcript, constraint=constraint))
        binding = {
            "schema": CONSTRAINT_SCHEMA, "contract_sha256": contract_sha,
            "manifest_sha256": reference.manifest_sha256, "release": contract.release,
            "published_date": contract.published_date, "croissant_date_published": contract.croissant_date_published,
            "retrieved_on": manifest.retrieved_on, "resource_image_digest": manifest.image_digest,
            "validation_scope": "canonical-output-target-all-partitions", "validation": asdict(validation),
            "artifacts": [asdict(a) for a in verified], "constraint_field_definition": definition,
            "constraint_field_schema_sha256": definition_sha, "duckdb_version": duckdb.__version__,
            "join": "exact-Ensembl-stable-ID; explicit-version-retained-not-invented",
            "lookups": {gene: {
                "status": "matched" if rows else "no-target-row",
                "rows": [_plain(r) for r in rows],
            } for gene, rows in rows_by_gene.items()},
            "unresolved_consequence_ids": [t.consequence_id for t, g in zip(transcripts, gene_ids, strict=True) if g is None],
        }
        provenance = SourcedRecord("Open Targets target", contract.release, "target-constraint-lookup",
                                  reference.manifest_sha256, reference.manifest_sha256, binding)
        return tuple(attached), provenance
    except Exception:
        # DuckDB and upstream verifiers can include query values and local paths.
        raise VepAnnotationError("offline Open Targets constraint attachment failed") from None


def _verify_bundle(contract_path, reference_root, manifest_sha256):
    root, contract_path = _local(reference_root), _local(contract_path)
    _preflight_resources(files=(contract_path,))
    contract = load_vep_bundle_contract(contract_path)
    if contract.tool_version != "116.2" or contract.versions.cache != "116" or contract.versions.assembly != "GRCh38":
        raise VepAnnotationError("unsupported VEP release")
    contract_sha = _file_identity(contract_path)[0]
    release = root / contract.cache.output_root
    bundle_path = _bundle_manifest_path(root, contract, contract_sha)
    resource_manifest_paths = tuple(
        _manifest_path(root, resource.resource_id, contract_sha)
        for resource in (contract.cache.source, contract.fasta.source)
    )
    _preflight_resources(
        files=(bundle_path, release / "cache.inventory.json", root / contract.cache.info_path,
               root / contract.fasta.output_path, root / contract.fasta.fai_path,
               *resource_manifest_paths),
        trees=(release,),
    )
    for path in (release, bundle_path):
        _require_within_root(root, path, "VEP resource")
        _readonly(path)
    for path in resource_manifest_paths:
        _require_within_root(root, path, "VEP resource manifest")
        _readonly(path)
    if not isinstance(manifest_sha256, str) or _SHA.fullmatch(manifest_sha256) is None or _file_identity(bundle_path)[0] != manifest_sha256:
        raise VepAnnotationError("VEP bundle manifest binding mismatch")
    bundle = json.loads(bundle_path.read_text(encoding="utf-8"))
    publications = _verify_resource_manifests(root, contract, contract_sha)
    if any(p is None for p in publications):
        raise VepAnnotationError("VEP resource manifests are incomplete")
    # This hashes EVERY extracted cache file and verifies the exact directory set.
    validation = _validate_release(release, contract)
    expected = _bundle_data(contract, contract_sha, publications, validation,
                            VepSentinelResult(**bundle["offline_sentinel"]), root)
    if bundle != expected:
        raise VepAnnotationError("VEP bundle contents differ from verified resources")
    fasta_path = root / contract.fasta.output_path
    reference = IndexedFasta(fasta_path, Path(str(fasta_path) + ".fai"))
    return contract, contract_sha, release, reference, validation


@lru_cache(maxsize=4)
def _verify_readonly_bundle_for_process(
    contract_path: str,
    reference_root: str,
    manifest_sha256: str,
):
    """Reuse exact bundle verification only inside one isolated execution.

    ``_verify_bundle`` proves the entire cache tree and requires it to be
    read-only.  Canonical paths plus the expected manifest digest form the
    process-local key; nothing is persisted or trusted in a later process.
    """
    return _verify_bundle(contract_path, reference_root, manifest_sha256)


def resource_native_aliases(reference: IndexedFasta) -> Mapping[str, str]:
    """Exact FASTA-native names plus explicit GRCh38 canonical aliases only.

    No generic chr stripping, accession guessing, decoy conversion or liftover.
    MT, M and chrM map to resource-native MT, not a new reference role.
    """
    aliases = {name: name for name in reference.index}
    for name in (*(str(i) for i in range(1, 23)), "X", "Y"):
        if name in reference.index:
            aliases["chr" + name] = name
    if "MT" in reference.index:
        aliases.update({"M": "MT", "chrM": "MT"})
    return MappingProxyType(aliases)


def _annotation_child(command, *, work, env, timeout, name, execution_config=None,
                      output_files=()):
    """Use bounded pipes and watched file outputs for bounded WGS callers."""
    if execution_config is None:
        return subprocess.run(command, cwd=work, env=env, capture_output=True,
                              timeout=timeout, check=True)
    from types import SimpleNamespace
    from genosapiens.bounded_json import read_bounded_bytes
    from genosapiens.wgs_annotation import run_bounded_child
    limit = execution_config.annotation.max_shard_bytes
    config = replace(execution_config, child_process=replace(execution_config.child_process,
        wall_timeout_seconds=min(timeout, execution_config.child_process.wall_timeout_seconds),
        max_total_output_bytes=min(limit, execution_config.child_process.max_total_output_bytes),
        max_stderr_bytes=min(limit, execution_config.child_process.max_stderr_bytes)))
    stdout, stderr = work / (name + ".capture.stdout"), work / (name + ".capture.stderr")
    run_bounded_child(tuple(command), cwd=work, stdout_path=stdout, stderr_path=stderr,
        execution_config=config, environment=env, output_files=tuple(output_files))
    return SimpleNamespace(stdout=read_bounded_bytes(stdout, maximum_bytes=limit),
                           stderr=read_bounded_bytes(stderr, maximum_bytes=limit))


def _annotation_content(path, execution_config=None):
    if execution_config is None:
        return path.read_bytes()
    from genosapiens.bounded_json import read_bounded_bytes
    return read_bounded_bytes(path, maximum_bytes=execution_config.annotation.max_shard_bytes)


def run_vep_annotation(
    requests: tuple[VepRequest, ...], *, contract_path: str | Path,
    reference_root: str | Path, bundle_manifest_sha256: str, runtime_image: str,
    upstream_image_digest: str, output_directory: str | Path,
    constraint_reference: OpenTargetsConstraintReference,
    reference_compatibility: VepReferenceCompatibility, timeout: int = 600,
    buffer_size: int = 1000,
    execution_config=None,
) -> VepAnnotationResult:
    """Annotate a caller-bounded batch in the current offline pinned container.

    output_directory must be a new child of existing inventoried writable storage.
    Emits transport.vcf, annotated.vcf, alleles.json and binding.json (controlled).
    runtime_image is an externally verified immutable OCI reference, not a tag or
    self-attestation of container identity. No Docker socket/daemon is accessed.
    A failed call publishes no directory. Reference failures are never repaired.
    """
    try:
        from genosapiens.vep_native_annotation import configuration as native_configuration, run_native_annotation

        if (type(constraint_reference) is not OpenTargetsConstraintReference
                or type(reference_compatibility) is not VepReferenceCompatibility
                or type(requests) is not tuple or any(type(r) is not VepRequest for r in requests)
                or len({r.allele.allele_id for r in requests}) != len(requests)
                or type(runtime_image) is not str
                or re.fullmatch(r"[A-Za-z0-9_./:-]+@sha256:[0-9a-f]{64}", runtime_image) is None
                or type(timeout) is not int or timeout <= 0 or type(buffer_size) is not int or buffer_size <= 0):
            raise VepAnnotationError("invalid VEP run configuration")
        if any(name != "lo" for _, name in socket.if_nameindex()):
            raise VepAnnotationError("VEP runtime requires an isolated network namespace")
        destination = _local(output_directory)
        if destination.exists() or not destination.parent.is_dir():
            raise VepAnnotationError("VEP output requires a new directory in existing storage")
        canonical_contract = str(Path(contract_path).resolve(strict=True))
        canonical_root = str(Path(reference_root).resolve(strict=True))
        contract, contract_sha, release, reference, validation = (
            _verify_readonly_bundle_for_process(
                canonical_contract,
                canonical_root,
                bundle_manifest_sha256,
            )
        )
        if upstream_image_digest != contract.tool_image_digest:
            raise VepAnnotationError("VEP upstream image provenance mismatch")
        _preflight_resources(files=(EXECUTABLE,))
        executable_sha = _file_identity(EXECUTABLE)[0]
        coordinate_maps, compatibility_proof = _source_coordinate_maps(
            reference_compatibility,
            contract,
            contract_sha,
            reference_root,
            reference,
        )
        configuration = {"schema": SCHEMA, "flags": list(FLAGS), "buffer_size": buffer_size,
                          "timeout_seconds": timeout, "reference_compatibility": compatibility_proof,
                          "constraint_adapter": CONSTRAINT_SCHEMA,
                          "constraint_manifest_sha256": constraint_reference.manifest_sha256,
                          "native_backend": native_configuration(),
                         "resource_reference_role": ReferenceRole.PUBLIC_ANNOTATION_BENCHMARK.value}
        binding = {
            "schema": SCHEMA, "classification": "controlled", "llm_policy": "no-llm",
            "runtime_image": runtime_image, "upstream_image_digest": upstream_image_digest,
            "executable_sha256": executable_sha, "contract_sha256": contract_sha,
            "bundle_manifest_sha256": bundle_manifest_sha256,
            "cache_inventory_sha256": validation.cache_inventory_sha256,
            "fasta_sha256": contract.fasta.output_sha256, "fai_sha256": contract.fasta.fai_sha256,
            "resource_reference_identity_id": contract.fasta.source.resource_id,
            "configuration": configuration, "configuration_sha256": _digest(_json(configuration)),
            "input_evidence_sha256": _digest(_json([
                {"identity_id": r.identity_id, "allele": json.loads(to_json(r.allele))} for r in requests
            ])), "assessments": [],
        }
        if execution_config is not None:
            from genosapiens.bounded_wgs_config import bounded_wgs_execution_config_as_json
            configuration["bounded_execution_sha256"] = hashlib.sha256(
                bounded_wgs_execution_config_as_json(execution_config)).hexdigest()
            binding["configuration_sha256"] = _digest(_json(configuration))
        assessments, records, record_to_index, native_inputs = [], [], {}, []
        for request in requests:
            request.__post_init__()
            identities = {i.identity_id: i for i in request.allele.identities}
            native = identities[request.allele.native_identity_id]
            if type(native.identity) in (StructuralAllele, RepeatAllele):
                native_inputs.append((len(assessments), request))
                assessments.append(VepAssessment(request.allele, request.identity_id, None, "unassessed", None))
                binding["assessments"].append({
                    "allele_id": request.allele.allele_id, "native_identity_id": native.identity_id,
                    "requested_identity_id": request.identity_id,
                    "source_reference_manifest_sha256": native.reference_manifest_sha256,
                    "transport_record_id": None,
                })
                continue
            sequence = _sequence(native.identity)
            resource_allele, reason = None, None
            coordinate_map = None if sequence is None else coordinate_maps.get(
                (sequence.reference_role.value, sequence.reference_identity_id))
            if sequence is None:
                reason = "native-identity-has-no-established-sequence-equivalence"
            elif coordinate_map is None:
                reason = "source-reference-not-coordinate-bound"
            elif native.reference_manifest_sha256 != coordinate_map[0]:
                reason = "source-reference-manifest-mismatch"
            elif any(re.fullmatch(r"[ACGT]+", v) is None for v in (sequence.reference, sequence.alternate)):
                reason = "ambiguous-sequence-bases"
            elif sequence.contig not in coordinate_map[1]:
                reason = "contig-outside-resource-alias-coverage"
            else:
                contig = coordinate_map[1][sequence.contig]
                end = sequence.position + len(sequence.reference) - 1
                if end > reference.index[contig].length:
                    reason = "reference-span-out-of-bounds"
                elif reference.fetch(contig, sequence.position, len(sequence.reference)).upper() != sequence.reference:
                    reason = "resource-reference-mismatch"
                else:
                    resource_allele = SequenceAllele(
                        ReferenceRole.PUBLIC_ANNOTATION_BENCHMARK, contract.fasta.source.resource_id,
                        contig, sequence.position, sequence.reference, sequence.alternate,
                    )
            index = len(assessments)
            assessments.append(VepAssessment(request.allele, request.identity_id, resource_allele,
                                             "unassessed" if reason else "annotated", reason))
            if reason is None:
                key = f"vep-{index + 1:09d}"
                records.append(VepVcfRecord(key, resource_allele.contig, resource_allele.position,
                                            resource_allele.reference, (resource_allele.alternate,)))
                record_to_index[key] = index
            binding["assessments"].append({
                "allele_id": request.allele.allele_id, "native_identity_id": native.identity_id,
                "requested_identity_id": request.identity_id,
                "source_reference_manifest_sha256": native.reference_manifest_sha256,
                "status": "unassessed" if reason else "annotated", "reason": reason,
                "transport_record_id": None if reason else key,
            })
        order = {name: i for i, name in enumerate(reference.index)}
        records.sort(key=lambda r: (order[r.contig], r.position, r.record_id))
        transport = "##fileformat=VCFv4.2\n" + "".join(
            f"##contig=<ID={name},length={reference.index[name].length}>\n"
            for name in dict.fromkeys(r.contig for r in records)
        ) + _COLUMNS + "\n" + "".join(r.line() for r in records)
        binding["transport_sha256"] = _digest(transport)
        with tempfile.TemporaryDirectory(prefix=".vep-", dir=destination.parent) as directory:
            work = Path(directory)
            os.chmod(work, 0o700)
            env = {"PATH": "/usr/local/bin:/usr/bin:/bin", "HOME": str(work), "TMPDIR": str(work),
                   "LANG": "C", "LC_ALL": "C", "PERL5LIB": "/opt/vep/src/ensembl-vep:/opt/vep/src/ensembl-vep/modules",
                   "PERL_HASH_SEED": "0", "PERL_PERTURB_KEYS": "0"}
            if execution_config is None:
                version = subprocess.run([str(EXECUTABLE), "--help"], cwd=work, env=env,
                                         capture_output=True, text=True, timeout=120, check=True)
            else:
                version = _annotation_child([str(EXECUTABLE), "--help"], work=work, env=env,
                    timeout=120, name="version", execution_config=execution_config)
                version.stdout, version.stderr = version.stdout.decode(), version.stderr.decode()
            versions = re.findall(r"ensembl-vep\s*:\s*([0-9]+(?:\.[0-9]+)?)", version.stdout + version.stderr)
            if versions != [contract.tool_version]:
                raise VepAnnotationError("VEP executable version mismatch")
            binding["tool_version"] = versions[0]
            binding["runtime_version_probe_sha256"] = _digest(version.stdout + version.stderr)
            input_path, output_path = work / "transport.vcf", work / "annotated.vcf"
            input_path.write_text(transport, encoding="ascii")
            os.chmod(input_path, 0o600)
            publication_hashes = {input_path: binding["transport_sha256"]}
            log_paths = []
            if records:
                command = [str(EXECUTABLE), *FLAGS, "--buffer_size", str(buffer_size),
                           "--dir", str(work), "--dir_plugins", str(work), "--dir_cache", str(release),
                           "--fasta", str(reference.fasta_path), "--input_file", str(input_path),
                           "--output_file", str(output_path)]
                completed = _annotation_child(command, work=work, env=env, timeout=timeout,
                    name="vep", execution_config=execution_config, output_files=(output_path,))
                for name, data in (("vep.stdout", completed.stdout), ("vep.stderr", completed.stderr)):
                    path = work / name
                    path.write_bytes(data)
                    os.chmod(path, 0o600)
                    log_paths.append(path)
                    binding[name + "_sha256"] = hashlib.sha256(data).hexdigest()
                    publication_hashes[path] = binding[name + "_sha256"]
                content = _annotation_content(output_path, execution_config).decode("utf-8")
                binding["output_sha256"] = _file_identity(output_path)[0]
                publication_hashes[output_path] = binding["output_sha256"]
                source = SourcedRecord("Ensembl VEP", contract.versions.cache, "vep-output",
                    binding["output_sha256"], bundle_manifest_sha256,
                    {key: binding[key] for key in (
                        "schema", "runtime_image", "upstream_image_digest", "executable_sha256",
                        "contract_sha256", "cache_inventory_sha256", "fasta_sha256", "fai_sha256",
                        "resource_reference_identity_id", "configuration_sha256", "input_evidence_sha256",
                        "transport_sha256", "tool_version",
                    )})
                parsed = parse_vep_vcf(content, expected=tuple(records), source=source, required_fields=RICH_FIELDS)
                for record in records:
                    index = record_to_index[record.record_id]
                    assessment = assessments[index]
                    previous = assessment.allele.transcript_consequences or ()
                    new = parsed.transcripts[(record.record_id, 1)]
                    assessments[index] = replace(assessment,
                        allele=replace(assessment.allele, transcript_consequences=previous + new),
                        non_transcript_consequences=parsed.non_transcripts[(record.record_id, 1)])
                binding["csq_fields"] = list(parsed.fields)
                os.chmod(output_path, 0o600)
            else:
                binding["output_sha256"] = None
                binding["csq_fields"] = None
            native_paths = ()
            if native_inputs:
                outcomes, proof, native_paths = run_native_annotation(
                    native_inputs, coordinate_maps=coordinate_maps, reference=reference, release=release,
                    work=work, env=env, executable=EXECUTABLE, timeout=timeout, buffer_size=buffer_size,
                    binding=binding, versions=asdict(contract.versions),
                    execution_config=execution_config,
                )
                binding["native_annotation"] = proof
                publication_hashes.update({p: proof["artifacts"][p.name] for p in native_paths})
                for index, (status, reason, transcripts, other, geometries) in outcomes.items():
                    previous = assessments[index]
                    allele = previous.allele
                    if transcripts:
                        allele = replace(allele, transcript_consequences=(allele.transcript_consequences or ()) + transcripts)
                    assessments[index] = replace(previous, allele=allele, status=status, reason=reason,
                        non_transcript_consequences=other, native_geometry_assessments=geometries)
                    binding["assessments"][index].update(status=status, reason=reason,
                        native_geometry_assessments=[_plain(g) for g in geometries])
            if _file_identity(input_path)[0] != binding["transport_sha256"] or _file_identity(EXECUTABLE)[0] != executable_sha:
                raise VepAnnotationError("VEP execution identity changed")
            enriched, constraint_source = attach_open_targets_constraints(
                tuple(t for a in assessments for t in a.allele.transcript_consequences or ()),
                reference=constraint_reference,
            )
            remaining = iter(enriched)
            for index, assessment in enumerate(assessments):
                previous = assessment.allele.transcript_consequences
                if previous is not None:
                    assessments[index] = replace(assessment, allele=replace(assessment.allele,
                        transcript_consequences=tuple(next(remaining) for _ in previous)))
            binding["constraint_source"] = _plain(constraint_source)
            alleles_path = work / "alleles.json"
            alleles_path.write_text(_json([json.loads(to_json(a.allele)) for a in assessments]), encoding="ascii")
            os.chmod(alleles_path, 0o600)
            binding["alleles_sha256"] = _file_identity(alleles_path)[0]
            publication_hashes[alleles_path] = binding["alleles_sha256"]
            binding_path = work / "binding.json"
            binding_path.write_text(_json(binding), encoding="ascii")
            os.chmod(binding_path, 0o600)
            publication_hashes[binding_path] = _file_identity(binding_path)[0]
            # Every previously bound byte stream must still be present and exact,
            # including sequence outputs/logs produced before native execution.
            for path, digest in publication_hashes.items():
                _local(path)
                if not stat.S_ISREG(path.lstat().st_mode) or _file_identity(path)[0] != digest:
                    raise VepAnnotationError("annotation artifacts changed before publication")
            # Publication is no-clobber; caller must serialize writers to its storage.
            destination.mkdir(mode=0o700)
            published = []
            try:
                for path in publication_hashes:
                    target = destination / path.name
                    os.link(path, target)
                    published.append(target)
            except OSError:
                for path in published:
                    path.unlink()
                destination.rmdir()
                raise
            provenance = SourcedRecord("Ensembl VEP", contract.versions.cache, "vep-run",
                                      _file_identity(binding_path)[0], bundle_manifest_sha256, binding)
            return VepAnnotationResult(tuple(assessments), provenance)
    except VepAnnotationError:
        raise
    except (OSError, ValueError, TypeError, RuntimeError, AttributeError, KeyError, subprocess.SubprocessError):
        raise VepAnnotationError("offline VEP annotation failed") from None
