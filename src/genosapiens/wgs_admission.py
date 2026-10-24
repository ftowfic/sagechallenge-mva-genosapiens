"""Ledger-backed retained evidence for bounded annotation and admission."""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from contextlib import closing
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sqlite3
from typing import Iterator, Literal

from genosapiens.allele_identity import (
    Breakpoint,
    CopyNumberEvidence,
    GenomicInterval,
    MitochondrialAllele,
    MosaicAllele,
    ReadSupport,
    ReferenceRole,
    RepeatAllele,
    RepeatSize,
    RepeatSizeUnit,
    SequenceAllele,
    StructuralAllele,
    StructuralVariantType,
)
from genosapiens.evidence_model import (
    AlleleEvidence,
    CallEvidence,
    IdentityMapping,
    IdentityStage,
    SourcedRecord,
    _typed,
)
from genosapiens.evidence_serialization import from_json, to_json
from genosapiens.variant_inventory import (
    ScreenOutcome,
    VariantRecordInventoryFacts,
    WgsClass,
)
from genosapiens.variant_normalization import (
    SourceAllele,
    SourceAlleleNormalization,
    SourceSmallVariantRecord,
    VARIANT_NORMALIZATION_SCHEMA,
)
from genosapiens.vcf_contracts import POLICY_SHA256 as VCF_POLICY_SHA256
from genosapiens.vcf_qc import FilterStatus, VcfQcStreamHeader, VcfRecordQc, technical_flags_for
from genosapiens.wgs_annotation import (
    RetainedAllele,
    VerifiedVepRecord,
    WgsAnnotationError,
    _require_prepared_projection_parent,
    _working_binds_source,
    iter_retained_alleles,
    open_retained_projection,
)
from genosapiens.wgs_prepared_case import PreparedIndex
from genosapiens.wgs_ledger import (
    MAX_ROW_BYTES,
    LedgerManifest,
    WorkingSetManifest,
    canonical,
    iter_ledger,
    read_artifact,
)
from genosapiens.bounded_wgs_config import BoundedWgsExecutionConfig
from genosapiens.annotation_admission import (
    AnnotationAssociationProof,
    sequence_annotation_associations,
)
from genosapiens.longtr_admission import FORMAT as LONGTR_FORMAT
from genosapiens.wgs_source_routing import (
    NativeSourceContract,
    REPEAT_IDENTITY_FORMAT,
    SourceOccurrenceFacts as _Occurrence,
    SourceRoutingContext,
    WgsSourceRoutingError,
    classify_source_occurrence,
    decode_source_occurrence,
)
from genosapiens.wgs_failures import failure_boundary


LEDGER_ADMISSION_INPUT_SCHEMA = "genosapiens.wgs-ledger-admission-input"
ADMISSION_INVENTORY_SCHEMA = "genosapiens.admission-inventory"
RULE_IDS = (
    "sequence-transcript-overlap",
    "native-structural-overlap",
    "regulatory-link",
    "repeat-locus-support",
    "explicit-clinical-link",
)
_DIGEST = re.compile(r"^[0-9a-f]{64}$")
_SOURCE_ALLELE = re.compile(r"^source-allele-[0-9a-f]{24}$")
_IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:/-]{0,127}$")


class WgsAdmissionError(RuntimeError):
    """Static bounded admission-input failure."""


def _digest(value: object) -> bool:
    return isinstance(value, str) and _DIGEST.fullmatch(value) is not None


def _count(value: object, minimum: int = 0) -> bool:
    return isinstance(value, int) and not isinstance(value, bool) and value >= minimum


@dataclass(frozen=True, slots=True, repr=False)
class SourceTechnicalObligation:
    source_allele_id: str
    sample_name: str
    technical_flags: tuple[str, ...]
    read_validation_required: bool
    source_filtered: bool

    def __post_init__(self) -> None:
        if (
            _SOURCE_ALLELE.fullmatch(self.source_allele_id) is None
            or not isinstance(self.sample_name, str)
            or not self.sample_name
            or not isinstance(self.technical_flags, tuple)
            or any(not isinstance(value, str) or not value for value in self.technical_flags)
            or type(self.read_validation_required) is not bool
            or type(self.source_filtered) is not bool
        ):
            raise WgsAdmissionError("source technical obligation is invalid")


@dataclass(frozen=True, slots=True, repr=False)
class LedgerAdmissionInput:
    retained_ordinal: int
    identity_id: str
    source_allele_ids: tuple[str, ...]
    transport: Literal[
        "sequence_vep", "native_vep", "repeat_catalog", "unresolved"
    ]
    allele: AlleleEvidence | None
    unresolved_source: SourceAllele | None
    required_rule_ids: tuple[str, ...]
    technical_obligations: tuple[SourceTechnicalObligation, ...]
    parent_source_manifest_sha256: str
    parent_working_set_manifest_sha256: str
    native_source_contract_sha256: str
    annotation_status: Literal["pending", "annotated", "unresolved", "unassessed"]
    unresolved_reasons: tuple[str, ...]

    vcf_evidence_policy_sha256: str = VCF_POLICY_SHA256

    def __post_init__(self) -> None:
        if (
            self.vcf_evidence_policy_sha256 != VCF_POLICY_SHA256
            or not _count(self.retained_ordinal, 1)
            or _IDENTIFIER.fullmatch(self.identity_id) is None
            or not isinstance(self.source_allele_ids, tuple)
            or not self.source_allele_ids
            or tuple(item.source_allele_id for item in self.technical_obligations)
            != self.source_allele_ids
            or any(_SOURCE_ALLELE.fullmatch(value) is None for value in self.source_allele_ids)
            or len(set(self.source_allele_ids)) != len(self.source_allele_ids)
            or self.transport not in {
                "sequence_vep", "native_vep", "repeat_catalog", "unresolved",
            }
            or self.required_rule_ids != tuple(
                rule for rule in RULE_IDS if rule in set(self.required_rule_ids)
            )
            or any(not _digest(value) for value in (
                self.parent_source_manifest_sha256,
                self.parent_working_set_manifest_sha256,
                self.native_source_contract_sha256,
            ))
            or self.annotation_status not in {"pending", "annotated", "unresolved", "unassessed"}
            or not isinstance(self.unresolved_reasons, tuple)
            or any(not isinstance(value, str) or not value for value in self.unresolved_reasons)
        ):
            raise WgsAdmissionError("ledger admission input is invalid")
        if self.transport == "unresolved":
            if (
                self.allele is not None
                or not isinstance(self.unresolved_source, SourceAllele)
                or self.annotation_status != "unresolved"
                or not self.unresolved_reasons
            ):
                raise WgsAdmissionError("ledger unresolved input is invalid")
        elif (
            not isinstance(self.allele, AlleleEvidence)
            or self.unresolved_source is not None
            or self.annotation_status == "unresolved"
            or self.allele.allele_id != self.identity_id
            and self.identity_id not in {
                mapping.identity_id for mapping in self.allele.identities
            }
        ):
            raise WgsAdmissionError("ledger derived input is invalid")
        if self.annotation_status == "annotated" and self.transport not in {
            "sequence_vep", "native_vep",
        }:
            raise WgsAdmissionError("ledger annotation state is invalid")
        if self.annotation_status == 'unassessed':
            from genosapiens.wgs_annotation_coverage import REASON_POLICIES
            if (self.transport != 'sequence_vep' or len(self.unresolved_reasons) != 1
                    or self.unresolved_reasons[0] not in REASON_POLICIES):
                raise WgsAdmissionError('ledger unassessed annotation state is invalid')

    @property
    def subject_id(self) -> str:
        return (
            self.unresolved_source.source_allele_id
            if self.unresolved_source is not None else self.allele.allele_id
        )

    def as_bytes(self) -> bytes:
        from genosapiens.evidence_model import _validation_session
        with _validation_session(isolated=True):
            return _owned_admission_input_bytes(_typed(self, LedgerAdmissionInput))


def _owned_admission_input_bytes(value):
    """Encode the record owner's immutable copy without another codec walk."""
    from genosapiens.evidence_serialization import _owned_json_envelope,_owned_canonical_bytes
    content = _owned_canonical_bytes({
        "schema": LEDGER_ADMISSION_INPUT_SCHEMA,
        "version": 1,
        "classification": "controlled",
        "llm_policy": "no-llm",
        "data": {
            "vcf_evidence_policy_sha256": value.vcf_evidence_policy_sha256,
            "retained_ordinal": value.retained_ordinal,
            "identity_id": value.identity_id,
            "source_allele_ids": value.source_allele_ids,
            "transport": value.transport,
            "allele": None if value.allele is None else _owned_json_envelope(value.allele),
            "unresolved_source": (
                None if value.unresolved_source is None else asdict(value.unresolved_source)
            ),
            "required_rule_ids": value.required_rule_ids,
            "technical_obligations": tuple(
                asdict(obligation) for obligation in value.technical_obligations
            ),
            "parent_source_manifest_sha256": value.parent_source_manifest_sha256,
            "parent_working_set_manifest_sha256":
                value.parent_working_set_manifest_sha256,
            "native_source_contract_sha256": value.native_source_contract_sha256,
            "annotation_status": value.annotation_status,
            "unresolved_reasons": value.unresolved_reasons,
        },
    })
    if len(content) > MAX_ROW_BYTES:
        raise WgsAdmissionError("ledger admission row exceeds its bounded size")
    return content


def load_ledger_admission_input(content: bytes) -> LedgerAdmissionInput:
    """Strictly decode one canonical bounded admission row."""

    from genosapiens.evidence_serialization import _from_mapping, _plain
    try:
        if not isinstance(content, bytes) or len(content) > MAX_ROW_BYTES:
            raise ValueError
        value = json.loads(content)
        if (
            canonical(value) != content
            or set(value) != {
                "schema", "version", "classification", "llm_policy", "data",
            }
            or value["schema"] != LEDGER_ADMISSION_INPUT_SCHEMA
            or type(value["version"]) is not int or value["version"] != 1
            or value["classification"] != "controlled"
            or value["llm_policy"] != "no-llm"
            or set(value["data"]) != {
                "vcf_evidence_policy_sha256",
                "retained_ordinal", "identity_id", "source_allele_ids", "transport",
                "allele", "unresolved_source", "required_rule_ids",
                "technical_obligations", "parent_source_manifest_sha256",
                "parent_working_set_manifest_sha256", "native_source_contract_sha256",
                "annotation_status", "unresolved_reasons",
            }
        ):
            raise ValueError
        data = value["data"]
        allele = (
            None if data["allele"] is None
            else _from_mapping(data["allele"], AlleleEvidence)
        )
        unresolved = (
            None if data["unresolved_source"] is None
            else _typed(data["unresolved_source"], SourceAllele, decode=True)
        )
        result = LedgerAdmissionInput(
            data["retained_ordinal"],
            data["identity_id"],
            tuple(data["source_allele_ids"]),
            data["transport"],
            allele,
            unresolved,
            tuple(data["required_rule_ids"]),
            tuple(
                _typed(item, SourceTechnicalObligation, decode=True)
                for item in data["technical_obligations"]
            ),
            data["parent_source_manifest_sha256"],
            data["parent_working_set_manifest_sha256"],
            data["native_source_contract_sha256"],
            data["annotation_status"],
            tuple(data["unresolved_reasons"]),
            vcf_evidence_policy_sha256=data["vcf_evidence_policy_sha256"],
        )
        decoded = _plain(result)
        if decoded['allele'] is not None:
            decoded['allele'] = {**data['allele'], 'data': decoded['allele']}
        if canonical(decoded) != canonical(data):
            raise ValueError
        return result
    except (
        KeyError, TypeError, ValueError, RuntimeError, json.JSONDecodeError,
    ):
        raise WgsAdmissionError("ledger admission row is invalid") from None


def _open_index(path: Path, config: BoundedWgsExecutionConfig) -> sqlite3.Connection:
    descriptor = -1
    try:
        descriptor = os.open(
            path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600,
        )
        os.fchmod(descriptor, 0o600)
        os.close(descriptor)
        descriptor = -1
        connection = sqlite3.connect(path)
        connection.enable_load_extension(False)
        connection.execute("PRAGMA foreign_keys=ON")
        connection.execute(f"PRAGMA cache_size=-{max(1, config.sqlite.cache_bytes // 1024)}")
        connection.execute("PRAGMA temp_store=FILE")
        connection.execute("PRAGMA mmap_size=0")
        connection.execute("PRAGMA trusted_schema=OFF")
        return connection
    except (OSError, sqlite3.Error):
        if descriptor >= 0:
            os.close(descriptor)
        raise WgsAdmissionError("cannot create ledger admission index") from None


def _rows(cursor, size: int):
    while True:
        batch = cursor.fetchmany(size)
        if not batch:
            return
        yield from batch


def _source_header(root: Path, manifest: LedgerManifest) -> VcfQcStreamHeader:
    try:
        value = json.loads(read_artifact(root, manifest.header))
        if (
            type(value) is not dict
            or set(value) != {
                "schema", "version", "classification", "llm_policy",
                "selected_sample_name", "header",
            }
            or value["selected_sample_name"] not in value["header"]["sample_names"]
        ):
            raise ValueError
        return _typed(value["header"], VcfQcStreamHeader, decode=True)
    except (KeyError, TypeError, ValueError, json.JSONDecodeError):
        raise WgsAdmissionError("source header admission binding is invalid") from None


def _occurrence(
    payload_json: str,
    *,
    source_allele_id: str,
    source_manifest: LedgerManifest,
    contract: NativeSourceContract,
) -> _Occurrence:
    with failure_boundary(
        "source_decode", "source_payload_decode",
        known_errors=(WgsSourceRoutingError,),
    ):
        return decode_source_occurrence(
            payload_json,
            source_allele_id=source_allele_id,
            source_manifest=source_manifest,
            contract=contract,
        )


def _technical(value: _Occurrence) -> SourceTechnicalObligation:
    flags = technical_flags_for(
        value.record, value.sample, value.source.allele_number,
    )
    screen = value.screen
    read_required = bool(flags) or (
        screen is not None and screen.outcome is ScreenOutcome.CANDIDATE
    )
    return SourceTechnicalObligation(
        value.source.source_allele_id,
        value.selected_sample_name,
        flags,
        read_required,
        value.record.filter_status is FilterStatus.FILTERED,
    )


def _source_record(value: _Occurrence, obligation: SourceTechnicalObligation) -> SourcedRecord:
    from genosapiens.vcf_qc import quality_acceptance_for
    return SourcedRecord(
        "genosapiens.wgs-source-ledger",
        "1",
        value.source.source_allele_id,
        value.payload_sha256,
        None,
        {
            "source_manifest_sha256": value.parent_source_manifest_sha256,
            "source_input_sha256": value.source.input_sha256,
            "source_record_id": value.source.source_record_id,
            "source_allele_id": value.source.source_allele_id,
            "record_number": value.record.record_number,
            "allele_number": value.source.allele_number,
            "raw_line_sha256": value.record.raw_line_sha256,
            "columns": tuple(value.record.raw_line.rstrip("\r\n").split("\t")),
            "filter_value": value.record.filter_value,
            "technical_flags": obligation.technical_flags,
            "read_validation_required": obligation.read_validation_required,
            "class_info_fields": value.record.class_info_fields,
            "quality": value.record.quality,
            "quality_acceptance": quality_acceptance_for(value.record, value.sample, value.source.allele_number),
            "storage_envelope": "bounded-ledger-v1",
        },
    )


def _call(
    value: _Occurrence,
    source: SourcedRecord,
) -> CallEvidence:
    return CallEvidence(
        source,
        value.sample,
        None if value.record.filter_status is FilterStatus.NOT_APPLIED
        else value.record.filters,
        value.record.quality,
    )


def _normalized_sequence(value: _Occurrence) -> SequenceAllele:
    normalized = value.normalization.normalized_allele
    if normalized is None:
        raise WgsAdmissionError("native sequence normalization is unavailable")
    return SequenceAllele(
        ReferenceRole(normalized.reference_role),
        normalized.reference_identity_id,
        normalized.contig,
        normalized.position,
        normalized.reference,
        normalized.alternate,
    )


def _source_sequence(value: _Occurrence) -> SequenceAllele:
    normalized = value.normalization.normalized_allele
    if normalized is None:
        raise WgsAdmissionError("source sequence normalization is unavailable")
    return SequenceAllele(
        ReferenceRole(normalized.reference_role),
        normalized.reference_identity_id,
        value.record.contig,
        value.record.position,
        value.record.reference,
        value.source.alternate,
    )


def _sequence_input(
    retained: RetainedAllele,
    occurrences: tuple[_Occurrence, ...],
    *,
    source_manifest: LedgerManifest,
    working_manifest: WorkingSetManifest,
    contract: NativeSourceContract,
) -> LedgerAdmissionInput:
    obligations = tuple(_technical(value) for value in occurrences)
    sequence = _normalized_sequence(occurrences[0])
    if (
        retained.identity_id != occurrences[0].normalization.normalized_allele.normalized_allele_id
        or any(_normalized_sequence(value) != sequence for value in occurrences)
    ):
        raise WgsAdmissionError("retained sequence identity is inconsistent")
    identities = []
    calls = []
    for value, obligation in zip(occurrences, obligations, strict=True):
        source = _source_record(value, obligation)
        identities.append(IdentityMapping(
            value.source.source_allele_id,
            IdentityStage.SOURCE,
            _source_sequence(value),
            source_manifest.source.reference_manifest_sha256,
            source,
        ))
        calls.append(_call(value, source))
    normalized_source = SourcedRecord(
        "genosapiens.wgs-normalization-ledger",
        "1",
        retained.identity_id,
        hashlib.sha256(canonical({
            "source_manifest_sha256": source_manifest.sha256,
            "working_set_manifest_sha256": working_manifest.sha256,
            "identity_id": retained.identity_id,
            "source_allele_ids": retained.source_allele_ids,
        })).hexdigest(),
        None,
        {
            "source_manifest_sha256": source_manifest.sha256,
            "working_set_manifest_sha256": working_manifest.sha256,
            "source_allele_ids": retained.source_allele_ids,
        },
    )
    identities.append(IdentityMapping(
        retained.identity_id,
        IdentityStage.NORMALIZED,
        sequence,
        source_manifest.source.reference_manifest_sha256,
        normalized_source,
        retained.source_allele_ids,
        VARIANT_NORMALIZATION_SCHEMA,
    ))
    native_id = "native-" + retained.identity_id
    identities.append(IdentityMapping(
        native_id,
        IdentityStage.NATIVE,
        sequence,
        source_manifest.source.reference_manifest_sha256,
        normalized_source,
        (retained.identity_id,),
        ADMISSION_INVENTORY_SCHEMA,
    ))
    allele = AlleleEvidence(
        retained.identity_id, native_id, tuple(identities), tuple(calls),
    )
    return LedgerAdmissionInput(
        retained.ordinal,
        retained.identity_id,
        retained.source_allele_ids,
        "sequence_vep",
        allele,
        None,
        (RULE_IDS[0], RULE_IDS[2], RULE_IDS[4]),
        obligations,
        source_manifest.sha256,
        working_manifest.sha256,
        contract.sha256,
        "pending",
        (),
    )


def _breakpoints(
    info: dict[str, str | None], definitions, contig: str, start: int, end: int,
) -> tuple[Breakpoint, ...]:
    values = []
    for name, position in (("CIPOS", start), ("CIEND", end)):
        raw = info.get(name)
        if raw in (None, "."):
            continue
        if (
            definitions.get(("INFO", name)) != ("2", "Integer")
            or not isinstance(raw, str)
            or re.fullmatch(r"-?[0-9]+,-?[0-9]+", raw) is None
        ):
            raise WgsAdmissionError("native breakpoint contract is invalid")
        values.append(Breakpoint(
            contig, position, confidence_interval=tuple(int(item) for item in raw.split(",")),
        ))
    return tuple(values)


def _structural_identity(
    value: _Occurrence, definitions, contig_lengths, *, reference_role, reference_identity_id,
) -> StructuralAllele:
    info = dict(value.record.class_info_fields)
    symbol = re.fullmatch(r"<(DEL|DUP|CNV)>", value.source.alternate)
    if (
        symbol is None
        or definitions.get(("INFO", "SVTYPE")) != ("1", "String")
        or definitions.get(("INFO", "END")) != ("1", "Integer")
        or info.get("SVTYPE") != symbol[1]
        or not isinstance(info.get("END"), str)
        or not info["END"].isdigit()
    ):
        raise WgsAdmissionError("native structural source format is unsupported")
    end = int(info["END"])
    if not value.record.position <= end <= contig_lengths[value.record.contig]:
        raise WgsAdmissionError("native structural interval is invalid")
    fields = dict(value.sample.class_format_fields)
    count = None
    if "CN" in fields and fields["CN"] not in {None, "."}:
        contract = definitions.get(("FORMAT", "CN"))
        raw = fields["CN"]
        if contract not in {("1", "Integer"), ("1", "Float")}:
            raise WgsAdmissionError("native copy-number contract is invalid")
        try:
            count = int(raw) if contract[1] == "Integer" else float(raw)
        except (TypeError, ValueError):
            raise WgsAdmissionError("native copy-number value is invalid") from None
        if not math.isfinite(count) or count < 0:
            raise WgsAdmissionError("native copy-number value is invalid")
    return StructuralAllele(
        ReferenceRole(reference_role),
        reference_identity_id,
        StructuralVariantType(symbol[1]),
        GenomicInterval(value.record.contig, value.record.position, end),
        _breakpoints(info, definitions, value.record.contig, value.record.position, end),
        CopyNumberEvidence(count, None, value.sample.genotype, value.sample.ploidy),
    )


def _repeat_identity(
    value: _Occurrence, definitions, contig_lengths, *, reference_role, reference_identity_id,
) -> RepeatAllele:
    info = dict(value.record.class_info_fields)
    if (
        definitions.get(("INFO", "REPID")) != ("1", "String")
        or definitions.get(("INFO", "RU")) != ("1", "String")
        or definitions.get(("INFO", "END")) != ("1", "Integer")
        or definitions.get(("FORMAT", "REPCN")) != ("1", "String")
        or not isinstance(info.get("REPID"), str)
        or info["REPID"] in {"", "."}
        or not isinstance(info.get("RU"), str)
        or re.fullmatch(r"[ACGT]+", info["RU"]) is None
        or not isinstance(info.get("END"), str)
        or not info["END"].isdigit()
    ):
        raise WgsAdmissionError("native repeat source format is unsupported")
    end = int(info["END"])
    if not value.record.position <= end <= contig_lengths[value.record.contig]:
        raise WgsAdmissionError("native repeat interval is invalid")
    raw_sizes = dict(value.sample.class_format_fields).get("REPCN")
    sizes = None
    if raw_sizes not in {None, "."}:
        if not isinstance(raw_sizes, str) or re.fullmatch(
            r"(?:[0-9]+|\.)(?:/(?:[0-9]+|\.))*", raw_sizes,
        ) is None:
            raise WgsAdmissionError("native repeat sizes are invalid")
        sizes = tuple(
            RepeatSize(None if item == "." else int(item), RepeatSizeUnit.REPEAT_UNITS)
            for item in raw_sizes.split("/")
        )
    return RepeatAllele(
        ReferenceRole(reference_role),
        reference_identity_id,
        info["REPID"],
        GenomicInterval(value.record.contig, value.record.position, end),
        info["RU"],
        sizes,
        value.sample.genotype,
        value.sample.ploidy,
    )


def _native_input(
    retained: RetainedAllele,
    value: _Occurrence,
    *,
    source_manifest: LedgerManifest,
    working_manifest: WorkingSetManifest,
    routing: SourceRoutingContext,
    definitions,
    contig_lengths,
) -> LedgerAdmissionInput:
    contract = routing.contract
    obligation = _technical(value)
    source_record = _source_record(value, obligation)
    represented = set(value.inventory.represented_classes)
    try:
        route = classify_source_occurrence(
            value,
            contract=contract,
            renderer=routing.renderer,
        )
    except WgsSourceRoutingError:
        raise WgsAdmissionError("native admission routing is invalid") from None
    if (
        retained.route != "native"
        or retained.identity_id != value.source.source_allele_id
        or route.preparation_disposition != "native_class"
        or route.normalized_allele_id is not None
    ):
        raise WgsAdmissionError("native admission routing is inconsistent")
    native_tags = set(route.native_class_tags)
    reason = None
    identity = None
    transport = "unresolved"
    required = RULE_IDS
    try:
        normalized = value.normalization.normalized_allele
        if "mitochondrial" in native_tags:
            sequence = _normalized_sequence(value)
            identity = MitochondrialAllele(
                sequence,
                None,
                ReadSupport(value.sample.depth, value.sample.allelic_depths, value.source.allele_number),
                value.sample.genotype,
                value.sample.ploidy,
            )
            transport = "sequence_vep"
            required = (RULE_IDS[0], RULE_IDS[2], RULE_IDS[4])
        elif "mosaic" in native_tags:
            sequence = _normalized_sequence(value)
            identity = MosaicAllele(
                sequence,
                value.screen.allele_balance,
                ReadSupport(value.sample.depth, value.sample.allelic_depths, value.source.allele_number),
                value.sample.genotype,
                value.sample.ploidy,
            )
            transport = "sequence_vep"
            required = (RULE_IDS[0], RULE_IDS[2], RULE_IDS[4])
        elif "repeat" in native_tags:
            if contract.repeat_identity_format != REPEAT_IDENTITY_FORMAT:
                reason = "native-repeat-source-contract-unavailable"
                raise WgsAdmissionError(reason)
            identity = _repeat_identity(
                value, definitions, contig_lengths,
                reference_role=contract.reference_role,
                reference_identity_id=contract.reference_identity_id,
            )
            transport = "repeat_catalog"
            required = (RULE_IDS[2], RULE_IDS[3], RULE_IDS[4])
        elif "structural" in native_tags:
            identity = _structural_identity(
                value, definitions, contig_lengths,
                reference_role=contract.reference_role,
                reference_identity_id=contract.reference_identity_id,
            )
            transport = "native_vep"
            required = (RULE_IDS[1], RULE_IDS[2], RULE_IDS[4])
        elif native_tags.intersection({"x_linked", "y_linked"}) and normalized is not None:
            identity = _normalized_sequence(value)
            transport = "sequence_vep"
            required = (RULE_IDS[0], RULE_IDS[2], RULE_IDS[4])
        else:
            reason = "unsupported-native-source-semantics"
            raise WgsAdmissionError(reason)
    except WgsAdmissionError:
        if reason is None:
            reason = (
                "unsupported-or-unbound-repeat-source-format"
                if WgsClass.REPEAT in represented else
                "unsupported-structural-source-format"
                if represented.intersection({WgsClass.STRUCTURAL, WgsClass.COPY_NUMBER}) else
                "unsupported-native-source-semantics"
            )
    if identity is None:
        return LedgerAdmissionInput(
            retained.ordinal,
            retained.identity_id,
            retained.source_allele_ids,
            "unresolved",
            None,
            value.source,
            tuple(required),
            (obligation,),
            source_manifest.sha256,
            working_manifest.sha256,
            contract.sha256,
            "unresolved",
            (reason,),
        )

    source_identity = (
        _source_sequence(value)
        if value.normalization.normalized_allele is not None else identity
    )
    identities = [IdentityMapping(
        value.source.source_allele_id,
        IdentityStage.SOURCE,
        source_identity,
        source_manifest.source.reference_manifest_sha256,
        source_record,
    )]
    parents = (value.source.source_allele_id,)
    if value.normalization.normalized_allele is not None:
        normalized_id = value.normalization.normalized_allele.normalized_allele_id
        identities.append(IdentityMapping(
            normalized_id,
            IdentityStage.NORMALIZED,
            _normalized_sequence(value),
            source_manifest.source.reference_manifest_sha256,
            source_record,
            parents,
            VARIANT_NORMALIZATION_SCHEMA,
        ))
        parents = (normalized_id,)
    allele_id = "allele-" + hashlib.sha256(canonical([
        ADMISSION_INVENTORY_SCHEMA,
        value.source.source_allele_id,
        value.selected_sample_index + 1,
        identity.kind.value,
    ])).hexdigest()[:24]
    native_id = "native-" + allele_id
    identities.append(IdentityMapping(
        native_id,
        IdentityStage.NATIVE,
        identity,
        source_manifest.source.reference_manifest_sha256,
        source_record,
        parents,
        ADMISSION_INVENTORY_SCHEMA,
    ))
    allele = AlleleEvidence(
        allele_id,
        native_id,
        tuple(identities),
        (_call(value, source_record),),
    )
    return LedgerAdmissionInput(
        retained.ordinal,
        retained.identity_id,
        retained.source_allele_ids,
        transport,
        allele,
        None,
        tuple(required),
        (obligation,),
        source_manifest.sha256,
        working_manifest.sha256,
        contract.sha256,
        "pending",
        (),
    )


def iter_admission_source_rows(connection, *, first_ordinal=None, last_ordinal=None,
                               ordinals=None, execution_config, from_retained_projection=True):
    """Read a bounded range or exact ordinal set through existing indexes.

    Each range uses one identity query and one ordered occurrence/source join.
    Only one identity's bounded payloads are retained while the cursors advance.
    Source occurrences remain distinct and in their original occurrence order.
    """
    if ordinals is None:
        if (not _count(first_ordinal, 1) or not _count(last_ordinal, first_ordinal)
                or last_ordinal - first_ordinal + 1 > execution_config.sqlite.cursor_batch_rows):
            raise WgsAdmissionError("admission source range exceeds capacity")
        expected = range(first_ordinal, last_ordinal+1)
        parameters = (first_ordinal, last_ordinal)
        clause = 'BETWEEN ? AND ?'
    else:
        if (first_ordinal is not None or last_ordinal is not None
                or type(ordinals) is not tuple or not ordinals
                or len(ordinals) > execution_config.sqlite.cursor_batch_rows
                or any(not _count(v, 1) for v in ordinals)
                or ordinals != tuple(sorted(set(ordinals)))):
            raise WgsAdmissionError("admission source ordinals exceed capacity")
        expected = parameters = ordinals
        clause = 'IN (' + ','.join('?' for _ in ordinals) + ')'
    identities_sql = (
        "SELECT i.ordinal,i.identity_id,i.route,f.contig,f.position,f.reference,"
        "f.alternate,f.kind,f.payload_json FROM identities i JOIN facts f "
        "USING(identity_id) WHERE i.ordinal BETWEEN ? AND ? ORDER BY i.ordinal"
        if from_retained_projection else
        "SELECT ordinal,identity_id,route,contig,position,reference,alternate,kind,"
        "payload_json FROM identities WHERE ordinal BETWEEN ? AND ? ORDER BY ordinal"
    )
    occurrences_sql = (
        "SELECT i.ordinal,o.source_allele_id,p.payload_json,o.occurrence_ordinal "
        "FROM identities i CROSS JOIN occurrences o ON o.identity_id=i.identity_id "
        "JOIN source_payloads p USING(record_number) WHERE i.ordinal BETWEEN ? AND ? "
        "ORDER BY i.ordinal,o.occurrence_ordinal"
        if from_retained_projection else
        "SELECT i.ordinal,o.source_allele_id,o.payload_json,o.occurrence_ordinal "
        "FROM identities i CROSS JOIN occurrences o ON o.identity_id=i.identity_id "
        "WHERE i.ordinal BETWEEN ? AND ? ORDER BY i.ordinal,o.occurrence_ordinal"
    )
    # CROSS JOIN fixes the outer indexed ordinal range; SQLite must not choose a
    # whole occurrence/source scan for each small worker job.
    identities_sql = identities_sql.replace('BETWEEN ? AND ?', clause)
    occurrences_sql = occurrences_sql.replace('BETWEEN ? AND ?', clause)
    with closing(connection.execute(identities_sql, parameters)) as identities, \
            closing(connection.execute(occurrences_sql, parameters)) as occurrences:
        pending = occurrences.fetchone()
        count = 0
        for identity in identities:
            ordinal = identity[0]
            if count >= len(expected) or ordinal != expected[count]:
                raise WgsAdmissionError("admission source range is incomplete")
            payloads = []
            total = 0
            prior_occurrence = 0
            while pending is not None and pending[0] == ordinal:
                _, source_id, payload, occurrence_ordinal = pending
                if (type(payload) is not str or type(occurrence_ordinal) is not int
                        or occurrence_ordinal <= prior_occurrence):
                    raise WgsAdmissionError("admission source occurrence is invalid")
                total += len(payload.encode("ascii"))
                if total > execution_config.annotation.max_response_row_bytes:
                    raise WgsAdmissionError("ledger admission identity exceeds capacity")
                payloads.append((source_id, payload))
                prior_occurrence = occurrence_ordinal
                pending = occurrences.fetchone()
            if not payloads or pending is not None and pending[0] < ordinal:
                raise WgsAdmissionError("admission source occurrence coverage is incomplete")
            yield identity, tuple(payloads)
            count += 1
        if count != len(expected) or pending is not None:
            raise WgsAdmissionError("admission source range is incomplete")


def _iter_admission_projection_rows(connection, *, source_root, source_manifest,
                                    working_manifest, source_routing_context, header,
                                    definitions, contig_lengths, execution_config,
                                    from_retained_projection=False,
                                    first_ordinal=1, last_ordinal=None):
    native_source_contract = source_routing_context.contract
    if native_source_contract.repeat_identity_format == LONGTR_FORMAT:
        if first_ordinal != 1 or last_ordinal is not None:
            raise WgsAdmissionError("LongTR projection requires complete locus routing")
        from genosapiens.wgs_longtr import iter_longtr_admission_inputs
        yield from iter_longtr_admission_inputs(
            connection, source_root=source_root, source_manifest=source_manifest,
            working_manifest=working_manifest, routing=source_routing_context,
            header=header, execution_config=execution_config,
            from_retained_projection=from_retained_projection,
        )
        return
    def ranges():
        total = (working_manifest.summary.retained_unique_identities
                 if last_ordinal is None else last_ordinal)
        size = execution_config.sqlite.cursor_batch_rows
        for first in range(first_ordinal, total + 1, size):
            yield from iter_admission_source_rows(connection, first_ordinal=first,
                last_ordinal=min(total, first + size - 1), execution_config=execution_config,
                from_retained_projection=from_retained_projection)

    previous = first_ordinal - 1
    occurrence_count = 0
    with closing(ranges()) as rows:
        for identity_row, payloads in rows:
            (
                retained_ordinal, identity_id, route, contig, position,
                reference, alternate, kind, payload_json,
            ) = identity_row
            if retained_ordinal != previous + 1:
                raise WgsAdmissionError("ledger admission projection order is invalid")
            previous = retained_ordinal
            values = []
            source_ids = []
            for source_id, payload in payloads:
                source_ids.append(source_id)
                values.append(_occurrence(
                    payload,
                    source_allele_id=source_id,
                    source_manifest=source_manifest,
                    contract=native_source_contract,
                ))
            occurrences = tuple(values)
            occurrence_count += len(occurrences)
            retained = RetainedAllele(
                retained_ordinal,
                identity_id,
                route,
                tuple(source_ids),
                source_manifest.sha256,
                working_manifest.sha256,
                contig,
                position,
                reference,
                alternate,
                kind,
                payload_json,
            )
            if tuple(value.source.source_allele_id for value in occurrences) != (
                retained.source_allele_ids
            ):
                raise WgsAdmissionError("ledger admission occurrence order is invalid")
            if retained.route == "sequence":
                result = _sequence_input(
                    retained,
                    occurrences,
                    source_manifest=source_manifest,
                    working_manifest=working_manifest,
                    contract=native_source_contract,
                )
            elif len(occurrences) == 1:
                result = _native_input(
                    retained,
                    occurrences[0],
                    source_manifest=source_manifest,
                    working_manifest=working_manifest,
                    routing=source_routing_context,
                    definitions=definitions,
                    contig_lengths=contig_lengths,
                )
            else:
                raise WgsAdmissionError("native retained identity is duplicated")
            # The final admission encoder owns canonical validation and capacity.
            yield result
    expected_last = (working_manifest.summary.retained_unique_identities
                     if last_ordinal is None else last_ordinal)
    if (previous != expected_last or first_ordinal == 1 and last_ordinal is None
            and occurrence_count != working_manifest.summary.retained_occurrences):
        raise WgsAdmissionError("ledger admission projection coverage is incomplete")


def iter_ledger_admission_inputs(
    *,
    source_root: str | Path,
    source_manifest: LedgerManifest,
    working_root: str | Path,
    working_manifest: WorkingSetManifest,
    source_routing_context: SourceRoutingContext,
    scratch_path: str | Path,
    execution_config: BoundedWgsExecutionConfig,
    prepared_index: PreparedIndex | None = None,
    prepared_state_root: str | Path | None = None,
    retained_projection_path: str | Path | None = None,
) -> Iterator[LedgerAdmissionInput]:
    """Project retained ledger rows without reconstructing supplied-call inputs."""

    source_root, working_root, scratch_path = (
        Path(source_root), Path(working_root), Path(scratch_path)
    )
    retained_scratch = scratch_path.with_name(scratch_path.name + ".retained")
    if (
        not source_root.is_absolute()
        or not working_root.is_absolute()
        or not scratch_path.is_absolute()
        or not isinstance(source_manifest, LedgerManifest)
        or not isinstance(working_manifest, WorkingSetManifest)
        or not _working_binds_source(working_manifest, source_manifest)
        or not isinstance(source_routing_context, SourceRoutingContext)
        or source_routing_context.contract.parent_source_manifest_sha256
        != source_manifest.sha256
        or source_routing_context.contract.reference_manifest_sha256
        != source_manifest.source.reference_manifest_sha256
        or not isinstance(execution_config, BoundedWgsExecutionConfig)
    ):
        raise WgsAdmissionError("ledger admission cursor binding is invalid")
    header = _source_header(source_root, source_manifest)
    definitions = {
        (namespace, name): (number, kind)
        for namespace, name, number, kind in header.class_field_definitions
    }
    contig_lengths = {value.contig: value.length for value in header.declared_contigs}
    if prepared_index is not None:
        if prepared_state_root is None or retained_projection_path is None:
            raise WgsAdmissionError("retained prepared projection binding is missing")
        try:
            _require_prepared_projection_parent(prepared_index, Path(prepared_state_root),
                source_manifest, working_manifest)
            with closing(open_retained_projection(Path(retained_projection_path),
                    prepared_index=prepared_index, routing=source_routing_context,
                    execution_config=execution_config,
                    expected_identities=working_manifest.summary.retained_unique_identities,
                    expected_occurrences=working_manifest.summary.retained_occurrences)) as retained:
                yield from _iter_admission_projection_rows(retained,
                    source_root=source_root, source_manifest=source_manifest,
                    working_manifest=working_manifest,
                    source_routing_context=source_routing_context, header=header,
                    definitions=definitions, contig_lengths=contig_lengths,
                    execution_config=execution_config, from_retained_projection=True)
            return
        except (WgsAnnotationError, sqlite3.Error):
            raise WgsAdmissionError("retained prepared projection is invalid") from None
    if prepared_state_root is not None or retained_projection_path is not None:
        raise WgsAdmissionError("unbound retained projection path")
    connection = _open_index(scratch_path, execution_config)
    sealed = False
    try:
        connection.executescript("""
            CREATE TABLE identities(
                ordinal INTEGER PRIMARY KEY,
                identity_id TEXT NOT NULL UNIQUE,
                route TEXT NOT NULL,
                contig TEXT,
                position INTEGER,
                reference TEXT,
                alternate TEXT,
                kind TEXT,
                payload_json TEXT NOT NULL
            );
            CREATE TABLE occurrences(
                source_allele_id TEXT PRIMARY KEY,
                identity_id TEXT NOT NULL REFERENCES identities(identity_id),
                occurrence_ordinal INTEGER NOT NULL,
                payload_json TEXT,
                UNIQUE(identity_id,occurrence_ordinal)
            );
        """)
        cursor = iter_retained_alleles(
            source_root=source_root,
            source_manifest=source_manifest,
            working_root=working_root,
            working_manifest=working_manifest,
            source_routing_context=source_routing_context,
            scratch_path=retained_scratch,
            execution_config=execution_config,
        )
        for retained in cursor:
            connection.execute(
                "INSERT INTO identities VALUES (?,?,?,?,?,?,?,?,?)",
                (
                    retained.ordinal, retained.identity_id, retained.route,
                    retained.contig, retained.position, retained.reference,
                    retained.alternate, retained.kind, retained.payload_json,
                ),
            )
            connection.executemany(
                "INSERT INTO occurrences VALUES (?,?,?,NULL)",
                (
                    (source_id, retained.identity_id, ordinal)
                    for ordinal, source_id in enumerate(
                        retained.source_allele_ids, start=1,
                    )
                ),
            )
            if retained.ordinal % execution_config.sqlite.cursor_batch_rows == 0:
                connection.commit()
        connection.commit()
        for row in iter_ledger(source_root, source_manifest):
            for source_id in row.source_allele_ids:
                updated = connection.execute(
                    "UPDATE occurrences SET payload_json=? "
                    "WHERE source_allele_id=? AND payload_json IS NULL",
                    (row.payload_json, source_id),
                )
                if updated.rowcount not in {0, 1}:
                    raise WgsAdmissionError("ledger admission source join is invalid")
        connection.commit()
        if connection.execute(
            "SELECT count(*) FROM occurrences WHERE payload_json IS NULL"
        ).fetchone()[0]:
            raise WgsAdmissionError("ledger admission source coverage is incomplete")
        sealed = True
        yield from _iter_admission_projection_rows(connection,
            source_root=source_root, source_manifest=source_manifest,
            working_manifest=working_manifest, source_routing_context=source_routing_context,
            header=header, definitions=definitions, contig_lengths=contig_lengths,
            execution_config=execution_config)
    except sqlite3.Error:
        raise WgsAdmissionError("ledger admission index failed") from None
    except WgsAnnotationError:
        raise WgsAdmissionError("retained admission source is invalid") from None
    finally:
        connection.close()
        for path in (scratch_path, retained_scratch):
            try:
                path.unlink()
            except (FileNotFoundError, OSError):
                pass
        if not sealed:
            pass


def enrich_sequence_admission_input(
    value: LedgerAdmissionInput,
    verified: VerifiedVepRecord,
) -> LedgerAdmissionInput:
    """Attach a fully reconciled VEP row without changing native identity or calls."""

    if (
        not isinstance(value, LedgerAdmissionInput)
        or value.transport != "sequence_vep"
        or value.annotation_status != "pending"
        or not isinstance(verified, VerifiedVepRecord)
        or verified.request.identity_id != value.identity_id
        or verified.request.source_allele_ids != value.source_allele_ids
    ):
        raise WgsAdmissionError("sequence admission response binding is invalid")
    key = (verified.request.request_id, 1)
    transcripts = verified.parsed.transcripts[key]
    allele = replace(
        value.allele,
        transcript_consequences=(value.allele.transcript_consequences or ()) + transcripts,
    )
    return replace(value, allele=allele, annotation_status="annotated")


def validated_sequence_annotation_intake(
    value: LedgerAdmissionInput,
    verified: VerifiedVepRecord,
    *,
    annotation_stage_manifest_sha256: str,
    bundle_manifest_sha256: str,
) -> tuple[LedgerAdmissionInput, SourcedRecord]:
    """Attach authenticated original CSQ evidence before constructing proofs."""
    provenance = sequence_annotation_intake_provenance(verified,
        annotation_stage_manifest_sha256=annotation_stage_manifest_sha256,
        bundle_manifest_sha256=bundle_manifest_sha256)
    return enrich_sequence_admission_input(value, verified), provenance


def sequence_annotation_intake_provenance(
    verified: VerifiedVepRecord,
    *,
    annotation_stage_manifest_sha256: str,
    bundle_manifest_sha256: str,
) -> SourcedRecord:
    """Preserve the original complete annotation/source binding, without proofs."""

    if (
        not _digest(annotation_stage_manifest_sha256)
        or not _digest(bundle_manifest_sha256)
        or not isinstance(verified, VerifiedVepRecord)
        or any(
            source.public_manifest_sha256 != bundle_manifest_sha256
            for values in verified.parsed.transcripts.values()
            for source in (item.source for item in values)
        )
        or any(
            source.public_manifest_sha256 != bundle_manifest_sha256
            for values in verified.parsed.non_transcripts.values()
            for source in values
        )
    ):
        raise WgsAdmissionError("sequence admission replay binding is invalid")
    return SourcedRecord(
        "genosapiens.wgs-annotation-replay",
        "1",
        "annotation-replay-" + verified.request.request_id,
        annotation_stage_manifest_sha256,
        bundle_manifest_sha256,
        {
            "annotation_stage_manifest_sha256": annotation_stage_manifest_sha256,
            "partition_output_sha256": verified.parent_output_sha256,
            "raw_line_sha256": verified.raw_line_sha256,
            "transport_request_id": verified.request.request_id,
            "working_set_manifest_sha256":
                verified.request.parent_working_set_manifest_sha256,
        },
    )


def sequence_annotation_intake_proofs(
    enriched: LedgerAdmissionInput,
    verified: VerifiedVepRecord,
    provenance: SourcedRecord,
) -> tuple[AnnotationAssociationProof, ...]:
    """The existing scientific interpretation, only after an input advances."""
    if type(verified) is not VerifiedVepRecord or type(provenance) is not SourcedRecord:
        raise WgsAdmissionError('sequence annotation intake proof source is invalid')
    key = (verified.request.request_id, 1)
    if (type(enriched) is not LedgerAdmissionInput or enriched.transport != 'sequence_vep'
            or enriched.annotation_status != 'annotated'
            or enriched.identity_id != verified.request.identity_id
            or enriched.source_allele_ids != verified.request.source_allele_ids
            or provenance.source_fields['transport_request_id'] != verified.request.request_id
            or provenance.source_fields['working_set_manifest_sha256'] != verified.request.parent_working_set_manifest_sha256
            or provenance.source_fields['raw_line_sha256'] != verified.raw_line_sha256
            or provenance.source_fields['partition_output_sha256'] != verified.parent_output_sha256):
        raise WgsAdmissionError('sequence annotation intake proof binding is invalid')
    transcripts = verified.parsed.transcripts[key]
    if transcripts and (enriched.allele.transcript_consequences or ())[-len(transcripts):] != transcripts:
        raise WgsAdmissionError('sequence annotation intake transcript binding is invalid')
    proofs = sequence_annotation_associations(
        enriched.allele,
        assessment_status="annotated",
        transcripts=verified.parsed.transcripts[key],
        non_transcripts=verified.parsed.non_transcripts[key],
        provenance=provenance,
    )
    if {proof.rule_key for proof in proofs} != {
        (enriched.allele.allele_id, RULE_IDS[0])
    }:
        raise WgsAdmissionError("sequence admission rule coverage is incomplete")
    return proofs


def replay_sequence_annotation_admission(
    value: LedgerAdmissionInput,
    verified: VerifiedVepRecord,
    *,
    annotation_stage_manifest_sha256: str,
    bundle_manifest_sha256: str,
) -> tuple[LedgerAdmissionInput, tuple[AnnotationAssociationProof, ...]]:
    """Shared full replay for advancing or original-format producer inputs."""
    enriched, provenance = validated_sequence_annotation_intake(value, verified,
        annotation_stage_manifest_sha256=annotation_stage_manifest_sha256,
        bundle_manifest_sha256=bundle_manifest_sha256)
    return enriched, sequence_annotation_intake_proofs(enriched, verified, provenance)
