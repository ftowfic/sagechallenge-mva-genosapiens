"""Authenticated source-occurrence facts and pure preparation routing."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import re
from typing import Literal

from genosapiens.allele_identity import ReferenceRole
from genosapiens.evidence_model import _typed
from genosapiens.longtr_admission import FORMAT as LONGTR_FORMAT
from genosapiens.submission_rendering import SubmissionContigRenderer
from genosapiens.variant_inventory import (
    AlleleBalanceScreenPolicy,
    ScreenOutcome,
    VariantRecordInventoryFacts,
    WgsClass,
    _screen_call,
)
from genosapiens.variant_normalization import (
    QuarantineReason,
    SourceAllele,
    SourceAlleleNormalization,
    SourceRecordNormalization,
    SourceSmallVariantRecord,
)
from genosapiens.vcf_contracts import SOURCE_PAYLOAD_VERSION
from genosapiens.vcf_qc import VariantClass, VcfRecordQc
from genosapiens.wgs_ledger import LedgerManifest, canonical


NATIVE_SOURCE_CONTRACT_SCHEMA = "genosapiens.wgs-native-source-contract"
REPEAT_IDENTITY_FORMAT = "repid-ru-repcn-pos-end"
_DIGEST = re.compile(r"^[0-9a-f]{64}$")
_IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:/-]{0,127}$")
_NATIVE_CLASSES = frozenset({
    WgsClass.STRUCTURAL,
    WgsClass.COPY_NUMBER,
    WgsClass.MITOCHONDRIAL,
    WgsClass.MOSAIC,
    WgsClass.REPEAT,
    WgsClass.REGULATORY,
})
_NATIVE_VARIANT_CLASSES = frozenset({
    VariantClass.SYMBOLIC,
    VariantClass.BREAKEND,
    VariantClass.SPANNING_DELETION,
    VariantClass.OTHER,
})
_LONGTR_FIELDS = frozenset({"START", "MOTIF", "INEXACT_ALLELE"})


class WgsSourceRoutingError(RuntimeError):
    """Invalid source routing facts, without controlled values in the message."""


def _digest(value: object) -> bool:
    return isinstance(value, str) and _DIGEST.fullmatch(value) is not None


@dataclass(frozen=True, slots=True, repr=False)
class NativeSourceContract:
    """Explicit caller-format facts; contig spelling alone is not class authority."""

    parent_source_manifest_sha256: str
    reference_role: ReferenceRole
    reference_identity_id: str
    reference_manifest_sha256: str
    mitochondrial_contigs: tuple[str, ...]
    repeat_identity_format: str | None

    def __post_init__(self) -> None:
        if (
            not _digest(self.parent_source_manifest_sha256)
            or not isinstance(self.reference_role, ReferenceRole)
            or _IDENTIFIER.fullmatch(self.reference_identity_id) is None
            or not _digest(self.reference_manifest_sha256)
            or not isinstance(self.mitochondrial_contigs, tuple)
            or self.mitochondrial_contigs
            != tuple(sorted(set(self.mitochondrial_contigs)))
            or any(
                not isinstance(value, str) or not value
                for value in self.mitochondrial_contigs
            )
            or self.repeat_identity_format
            not in {None, REPEAT_IDENTITY_FORMAT, LONGTR_FORMAT}
        ):
            raise WgsSourceRoutingError("native source contract is invalid")

    def as_bytes(self) -> bytes:
        return canonical({
            "schema": NATIVE_SOURCE_CONTRACT_SCHEMA,
            "version": 1,
            "classification": "controlled",
            "llm_policy": "no-llm",
            "data": asdict(self),
        })

    @property
    def sha256(self) -> str:
        return hashlib.sha256(self.as_bytes()).hexdigest()


@dataclass(frozen=True, slots=True)
class SourceRoutingContext:
    """Authenticated reference inputs shared by production and replay."""

    contract: NativeSourceContract
    renderer: SubmissionContigRenderer

    def __post_init__(self) -> None:
        if (
            not isinstance(self.contract, NativeSourceContract)
            or not isinstance(self.renderer, SubmissionContigRenderer)
            or self.renderer.reference_binding.role
            != self.contract.reference_role.value
            or self.renderer.reference_binding.identity_id
            != self.contract.reference_identity_id
            or self.renderer.reference_manifest_sha256
            != self.contract.reference_manifest_sha256
        ):
            raise WgsSourceRoutingError("source routing context is invalid")


@dataclass(frozen=True, slots=True)
class SourceOccurrenceFacts:
    """One fully validated physical ALT and its exact selected sample call."""

    source: SourceAllele
    normalization: SourceAlleleNormalization
    record: VcfRecordQc
    inventory: VariantRecordInventoryFacts
    selected_sample_index: int
    selected_sample_name: str
    payload_sha256: str
    parent_source_manifest_sha256: str
    reference_role: ReferenceRole
    reference_identity_id: str
    reference_manifest_sha256: str

    def __post_init__(self) -> None:
        normalized = self.normalization.normalized_allele
        if (
            not isinstance(self.source, SourceAllele)
            or not isinstance(self.normalization, SourceAlleleNormalization)
            or not isinstance(self.record, VcfRecordQc)
            or not isinstance(self.inventory, VariantRecordInventoryFacts)
            or self.source.record != SourceSmallVariantRecord(
                self.record.record_number,
                self.record.contig,
                self.record.position,
                self.record.reference,
                self.record.alternates,
            )
            or self.source.source_record_id != self.record.source_record_id
            or self.source.allele_number != self.normalization.allele_number
            or self.source.source_allele_id != self.normalization.source_allele_id
            or self.source.alternate != self.normalization.alternate
            or self.inventory.record_number != self.record.record_number
            or self.inventory.source_record_id != self.record.source_record_id
            or type(self.selected_sample_index) is not int
            or not 0 <= self.selected_sample_index < len(self.record.sample_calls)
            or self.record.sample_calls[self.selected_sample_index].sample_name
            != self.selected_sample_name
            or not _digest(self.payload_sha256)
            or not _digest(self.parent_source_manifest_sha256)
            or not isinstance(self.reference_role, ReferenceRole)
            or _IDENTIFIER.fullmatch(self.reference_identity_id) is None
            or not _digest(self.reference_manifest_sha256)
            or normalized is not None
            and (
                normalized.reference_role != self.reference_role.value
                or normalized.reference_identity_id != self.reference_identity_id
            )
        ):
            raise WgsSourceRoutingError("source occurrence facts are invalid")
        _validate_screens(self.record, self.inventory)

    @property
    def sample(self):
        return self.record.sample_calls[self.selected_sample_index]

    @property
    def screen(self):
        matches = tuple(
            value
            for value in self.inventory.screen_calls
            if value.sample_call.sample_name == self.selected_sample_name
            and value.source_allele_id == self.source.source_allele_id
        )
        if len(matches) > 1:
            raise WgsSourceRoutingError("source mosaic screen is ambiguous")
        return matches[0] if matches else None


@dataclass(frozen=True, slots=True)
class SourceOccurrenceRoute:
    preparation_disposition: Literal[
        "normalized", "native_class", "non_variant", "quarantined"
    ]
    normalized_allele_id: str | None
    contig_class: Literal[
        "autosomal", "x_linked", "y_linked", "mitochondrial", "other"
    ]
    native_class_tags: tuple[str, ...]

    def __post_init__(self) -> None:
        normalized = self.preparation_disposition == "normalized"
        if (
            self.preparation_disposition
            not in {"normalized", "native_class", "non_variant", "quarantined"}
            or normalized != (self.normalized_allele_id is not None)
            or self.contig_class
            not in {"autosomal", "x_linked", "y_linked", "mitochondrial", "other"}
            or self.native_class_tags
            != tuple(sorted(set(self.native_class_tags)))
            or any(
                tag not in {"x_linked", "y_linked", "mitochondrial", "structural", "repeat", "mosaic"}
                for tag in self.native_class_tags
            )
        ):
            raise WgsSourceRoutingError("source occurrence route is invalid")


def _validate_screens(
    record: VcfRecordQc,
    inventory: VariantRecordInventoryFacts,
) -> None:
    calls = inventory.screen_calls
    if not calls:
        return
    if (
        len(calls) != len(record.sample_calls)
        or tuple(value.sample_call.sample_name for value in calls)
        != tuple(value.sample_name for value in record.sample_calls)
    ):
        raise WgsSourceRoutingError("source mosaic screen coverage is invalid")
    represented = set(inventory.represented_classes)
    unsupported = set(inventory.unsupported_classes)
    policy = AlleleBalanceScreenPolicy()
    for actual, sample in zip(calls, record.sample_calls, strict=True):
        expected = _screen_call(record, sample, represented, unsupported, policy)
        if (
            actual != expected
            or actual.sample_call != sample
            or actual.read_validation_required
            != (actual.outcome is ScreenOutcome.CANDIDATE)
        ):
            raise WgsSourceRoutingError("source mosaic screen binding is invalid")


def decode_source_occurrence(
    payload_json: str,
    *,
    source_allele_id: str,
    source_manifest: LedgerManifest,
    contract: NativeSourceContract,
    expected_record_number: int | None = None,
    expected_source_record_id: str | None = None,
    expected_source_allele_ids: tuple[str, ...] | None = None,
) -> SourceOccurrenceFacts:
    """Decode one source ALT and reject any inconsistent record-level facts."""

    try:
        if (
            not isinstance(source_manifest, LedgerManifest)
            or not isinstance(contract, NativeSourceContract)
            or contract.parent_source_manifest_sha256 != source_manifest.sha256
            or contract.reference_manifest_sha256
            != source_manifest.source.reference_manifest_sha256
        ):
            raise ValueError
        content = payload_json.encode("ascii")
        value = json.loads(content)
        if (
            canonical(value) != content
            or set(value) != {
                "schema", "version", "classification", "llm_policy",
                "selected_sample_name", "selected_sample_index", "raw_line_sha256",
                "qc", "normalization", "inventory",
            }
            or value["schema"] != "genosapiens.wgs-source-payload"
            or value["version"] != SOURCE_PAYLOAD_VERSION
            or value["classification"] != "controlled"
            or value["llm_policy"] != "no-llm"
        ):
            raise ValueError
        record = _typed(value["qc"], VcfRecordQc, decode=True)
        inventory = _typed(
            value["inventory"], VariantRecordInventoryFacts, decode=True,
        )
        source_record = SourceSmallVariantRecord(
            record.record_number,
            record.contig,
            record.position,
            record.reference,
            record.alternates,
        )
        normalization_data = value["normalization"]
        if set(normalization_data) != {
            "reference_identity_id", "reference_role", "alleles",
        }:
            raise ValueError
        normalization = SourceRecordNormalization(
            source_manifest.source.input_sha256,
            normalization_data["reference_identity_id"],
            normalization_data["reference_role"],
            record.source_record_id,
            source_record,
            tuple(
                _typed(item, SourceAlleleNormalization, decode=True)
                for item in normalization_data["alleles"]
            ),
        )
        source_ids = tuple(item.source_allele_id for item in normalization.alleles)
        matches = tuple(
            item
            for item in normalization.alleles
            if item.source_allele_id == source_allele_id
        )
        selected_index = value["selected_sample_index"]
        if (
            len(matches) != 1
            or inventory.record_number != record.record_number
            or inventory.source_record_id != record.source_record_id
            or record.raw_line_sha256 != value["raw_line_sha256"]
            or record.raw_line is None
            or type(selected_index) is not int
            or not 0 <= selected_index < len(record.sample_calls)
            or record.sample_calls[selected_index].sample_name
            != value["selected_sample_name"]
            or normalization.reference_identity_id != contract.reference_identity_id
            or normalization.reference_role != contract.reference_role.value
            or expected_record_number is not None
            and record.record_number != expected_record_number
            or expected_source_record_id is not None
            and record.source_record_id != expected_source_record_id
            or expected_source_allele_ids is not None
            and source_ids != expected_source_allele_ids
        ):
            raise ValueError
        _validate_screens(record, inventory)
        source = normalization.source_allele(matches[0].allele_number)
        return SourceOccurrenceFacts(
            source,
            matches[0],
            record,
            inventory,
            selected_index,
            value["selected_sample_name"],
            hashlib.sha256(content).hexdigest(),
            source_manifest.sha256,
            contract.reference_role,
            contract.reference_identity_id,
            contract.reference_manifest_sha256,
        )
    except (
        KeyError, TypeError, ValueError, RuntimeError, UnicodeError,
        json.JSONDecodeError,
    ):
        raise WgsSourceRoutingError(
            "source occurrence routing replay is invalid"
        ) from None


def classify_source_occurrence(
    facts: SourceOccurrenceFacts,
    *,
    contract: NativeSourceContract,
    renderer: SubmissionContigRenderer,
) -> SourceOccurrenceRoute:
    """Apply the source preparation precedence without query or policy state."""

    try:
        if (
            not isinstance(facts, SourceOccurrenceFacts)
            or not isinstance(contract, NativeSourceContract)
            or not isinstance(renderer, SubmissionContigRenderer)
            or facts.parent_source_manifest_sha256
            != contract.parent_source_manifest_sha256
            or facts.reference_role is not contract.reference_role
            or facts.reference_identity_id != contract.reference_identity_id
            or facts.reference_manifest_sha256 != contract.reference_manifest_sha256
            or renderer.reference_binding.role != contract.reference_role.value
            or renderer.reference_binding.identity_id != contract.reference_identity_id
            or renderer.reference_manifest_sha256 != contract.reference_manifest_sha256
        ):
            raise ValueError
        source_name = renderer.source_names.get(facts.record.contig)
        if (
            source_name is None
            or facts.record.position + len(facts.record.reference) - 1 > source_name[1]
        ):
            raise ValueError
        canonical_contig = source_name[0]
        contig_class = (
            "x_linked" if canonical_contig == "X" else
            "y_linked" if canonical_contig == "Y" else
            "mitochondrial" if canonical_contig in {"M", "MT"} else
            "autosomal"
            if canonical_contig.isdecimal() and 1 <= int(canonical_contig) <= 22 else
            "other"
        )
        if (
            contig_class == "mitochondrial"
            and facts.record.contig not in contract.mitochondrial_contigs
        ):
            raise ValueError
        normalized = facts.normalization.normalized_allele
        if normalized is not None:
            normalized_name = renderer.source_names.get(normalized.contig)
            if (
                normalized_name is None
                or normalized.position + len(normalized.reference) - 1
                > normalized_name[1]
            ):
                raise ValueError

        variant_class = facts.record.variant_classes[facts.source.allele_number - 1]
        represented = set(facts.inventory.represented_classes)
        represented_native = represented.intersection(_NATIVE_CLASSES)
        represented_native_for_alt = bool(
            represented_native
            and (
                WgsClass.SMALL_VARIANT not in represented
                or variant_class in _NATIVE_VARIANT_CLASSES
            )
        )
        screen = facts.screen
        selected_mosaic = (
            screen is not None and screen.outcome is ScreenOutcome.CANDIDATE
        )
        native_tags = set()
        if contig_class in {"x_linked", "y_linked", "mitochondrial"}:
            native_tags.add(contig_class)
        if represented_native_for_alt and represented.intersection({
            WgsClass.STRUCTURAL, WgsClass.COPY_NUMBER,
        }):
            native_tags.add("structural")
        if represented_native_for_alt and WgsClass.REPEAT in represented:
            native_tags.add("repeat")
        if selected_mosaic:
            native_tags.add("mosaic")

        complete_repeat_locus = _LONGTR_FIELDS.issubset(
            dict(facts.record.class_info_fields)
        )
        non_variant = (
            variant_class is VariantClass.NON_VARIANT
            or facts.normalization.quarantine_reason
            is QuarantineReason.NON_VARIANT_ALTERNATE
        )
        declared_native = (
            contig_class in {"x_linked", "y_linked", "mitochondrial"}
            or selected_mosaic
            or represented_native_for_alt
        )
        if complete_repeat_locus:
            disposition = "native_class"
            normalized_id = None
            native_tags.add("repeat")
        elif non_variant:
            disposition = "non_variant"
            normalized_id = None
        elif declared_native:
            disposition = "native_class"
            normalized_id = None
        elif normalized is not None:
            disposition = "normalized"
            normalized_id = normalized.normalized_allele_id
        else:
            disposition = "quarantined"
            normalized_id = None
        return SourceOccurrenceRoute(
            disposition,
            normalized_id,
            contig_class,
            tuple(sorted(native_tags)),
        )
    except (KeyError, TypeError, ValueError, RuntimeError, OverflowError):
        raise WgsSourceRoutingError(
            "source occurrence classification is invalid"
        ) from None
