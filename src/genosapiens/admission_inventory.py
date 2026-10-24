"""Pure supplied-VCF accounting for Decision 6a, not an admission evaluator.

Build from SuppliedCallInputs, annotate the returned ``alleles``, then rebuild
with that *complete* enriched tuple. Consumers must call
``validate_admission_accounting`` with keys projected from their replayed rule
results, not a caller's coverage declaration. Its return value is the unresolved
source entries that must reach context/partner-enumeration limits. It does not
authenticate rule outcomes or establish completed assessment/partner search.

The closed INVENTORY_POLICY_BYTES pins applicability only; the admission policy
must additionally bind its executable adapters, annotations and public resources.
Neither pin is empirical acceptance or controlled-use authorization. This module
checks supplied reference bindings, not FASTA contents or reference acceptance.

All source bytes (plain VCF or gzip/BGZF), headers, records and ALT slots are
accounted for, including ALT='.' (one slot, zero actual alternate alleles).
The bounded standard symbolic SV adapter produces nominal interval context, not
precise affected bases. Optional fields never veto an established sequence
identity. Unsupported native conversion remains explicit; no proxy SNVs.
There is no filesystem access, acquisition, read processing or model callback.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping
from dataclasses import asdict, dataclass, replace
import gzip
import hashlib
import json
import math
import re
from types import MappingProxyType
from typing import Literal
import zlib

from genosapiens.allele_identity import (
    Breakpoint, CopyNumberEvidence, GenomicInterval, MitochondrialAllele,
    MosaicAllele, ReadSupport, ReferenceRole, RepeatAllele, RepeatSize,
    RepeatSizeUnit, SequenceAllele, StructuralAllele, StructuralVariantType,
)
from genosapiens.evidence_model import (
    AlleleEvidence, CallEvidence, IdentityMapping, IdentityStage, SourcedRecord,
    _Controlled, _ids, _sha, _typed,
)
from genosapiens.evidence_serialization import _loads, _plain
from genosapiens.variant_inventory import (
    BREAKEND, CN_SYMBOL, REPEAT_SYMBOL, SEQUENCE_CLASSES, SV_TYPES,
    ScreenOutcome, VariantClassInventory, WgsClass,
    _record_classes, build_variant_class_inventory, variant_inventory_as_json,
)
from genosapiens.variant_normalization import (
    QuarantineReason, SmallVariantNormalizationResult, SourceAllele,
    SourceSmallVariantRecord, VARIANT_NORMALIZATION_SCHEMA,
    source_allele_id, source_record_id, variant_normalization_as_json,
)
from genosapiens import vcf_qc as qc
from genosapiens import longtr_admission as longtr


SCHEMA = "genosapiens.admission-inventory"
RULE_IDS = (
    "sequence-transcript-overlap", "native-structural-overlap", "regulatory-link",
    "repeat-locus-support", "explicit-clinical-link",
)
# Unknown native identity semantics require all rules; optional uninterpreted
# annotations are retained but have no authority over identity or applicability.
_CLASS_RULES = MappingProxyType({
    "sequence": (RULE_IDS[0], RULE_IDS[2], RULE_IDS[4]),
    "structural": (RULE_IDS[1], RULE_IDS[2], RULE_IDS[4]),
    "repeat": (RULE_IDS[2], RULE_IDS[3], RULE_IDS[4]),
    "unknown": RULE_IDS,
    "non-variant": (),
})
_BASE_INFO = frozenset({"AN", "AF"})
_BASE_FORMAT = frozenset(qc.FORMAT_CONTRACTS)
SOURCE_FORMAT_SCHEMA = "genosapiens.supplied-vcf-format"
SOURCE_SET_SCHEMA = "genosapiens.supplied-source-set"
REPEAT_IDENTITY_FORMAT = "repid-ru-repcn-pos-end"
_STRUCTURAL_INFO = frozenset({"SVTYPE", "END", "CIPOS", "CIEND", "IMPRECISE"})
INVENTORY_POLICY_BYTES = (json.dumps({
    "schema": "genosapiens.admission-inventory-policy",
    "class_rules": dict(_CLASS_RULES),
    "base_info": sorted(_BASE_INFO), "base_format": sorted(_BASE_FORMAT),
    "native_classifier": "genosapiens.variant-class-inventory",
    "unknown_fields": "retain-verbatim; no-identity-or-applicability-veto",
    "identity_hints": "SVTYPE,CN,REPCN,CIPOS,CIEND,IMPRECISE; RU alone cannot replace concrete sequence",
    "mixed_source_classes": "ALT-specific; conflicting native hints unresolved",
    "non_variant": "dot-or-ref-equal-without-native-hints; placeholders-unresolved",
    "sequence_projection": "retain supported sources independently of unresolved siblings",
    "sequence_repeat_context": "STR/RU/RPA retained as context; explicit repeat SVTYPE/REPCN adds source-level unresolved repeat rule, never removes sequence projection",
    "structural": {
        "format": "VCFv4.2-or-v4.3; single-symbolic-DEL-DUP-CNV",
        "required": {"INFO/SVTYPE": ["1", "String"], "INFO/END": ["1", "Integer"]},
        "type": "exact-symbol-and-SVTYPE-agreement; no-CN-direction-inference",
        "interval": "unchanged-POS-END-1-based-closed-nominal-context; not affected-bases",
        "precision": "unknown; no-PRECISE-or-absence-of-IMPRECISE-inference",
        "confidence": "declared-CIPOS-CIEND-Number2-Integer; source-anchors-and-signed-offsets",
        "copy_number": "optional-declared-FORMAT/CN-Number1-Integer-or-Float; baseline-unknown",
        "samples": "one-native-projection-per-source-ALT-and-sample; no-dosage-aggregation",
    },
    "mitochondrial": "explicit-input-and-reference-bound-source-format-contigs; heteroplasmy-unknown",
    "mosaic": "replayed-opt-in-screen-candidate-only; measured-AD/DP; mandatory-unperformed-read-validation",
    "source_format_schema": SOURCE_FORMAT_SCHEMA,
    "source_set": {"schema": SOURCE_SET_SCHEMA, "samples": "explicit-case-and-per-input-sample-assignment",
                   "reference": "common-role-identity-manifest; no-remapping",
                   "roster": "independently-pinned-full-manifest; preserve-per-file-artifacts",
                   "views": "source-addressed-sequence-views; native-parent-lineage-unchanged"},
    "repeat": {"format": REPEAT_IDENTITY_FORMAT, "locus": "explicit-INFO/REPID-not-VCF-ID",
               "required": {"INFO/REPID": ["1", "String"], "INFO/RU": ["1", "String"],
                            "INFO/END": ["1", "Integer"], "FORMAT/REPCN": ["1", "String"]},
                "interval": "POS-END-1-based-closed-repeat-region-defined-by-source-format",
                "sizes": "source-ordered-slash-separated-integer-repeat-units-or-dot; no-ALT-size-join",
                "representation": "authoritative-sample-locus-call; no-individual-allele-projection",
                "missing_sizes": "retain-locus; whole-dot-or-absent-is-unknown; cardinality-independent-of-GT-ploidy",
                "malformed_sizes": "retain-locus-and-original-fields; unresolved-measurement-list; never-repair-or-drop-other-samples",
                "scope": "single-STR-or-VNTR-symbol; no-catalog-or-gene-inference"},
    "qc": qc.VCF_QC_SCHEMA,
    "vcf_evidence_policy_sha256": qc.POLICY_SHA256,
    "longtr": longtr.POLICY,
    "scope": "supplied-records-only; no-partner-search-or-reference-acceptance",
}, sort_keys=True, separators=(",", ":")) + "\n").encode("ascii")
INVENTORY_POLICY_SHA256 = hashlib.sha256(INVENTORY_POLICY_BYTES).hexdigest()
RuleEvaluationKey = tuple[str, str]


class AdmissionInventoryError(ValueError):
    """Static validation error, distinct from a retained unsupported source entry."""


@dataclass(frozen=True, slots=True, repr=False)
class SuppliedCallInputs(_Controlled):
    """Required original bytes and existing artifacts; never reconstructed reports.

    ``class_inventory`` is mandatory even when its existing opt-in screen was
    disabled. It is replayed from QC, not trusted as a class/coverage assertion.
    Bytes and their externally selected input pin must both be supplied.

    Optional ``source_format_bytes`` and its SHA-256 are an explicit source
    interpretation, NOT a precision/coverage flag or acceptance certificate.
    The closed JSON object requires exactly: schema=SOURCE_FORMAT_SCHEMA,
    input_sha256, reference_role, reference_identity_id, reference_manifest_sha256,
    mitochondrial_contigs (an array of declared contigs), and repeat_identity
    (null, REPEAT_IDENTITY_FORMAT, or longtr_admission.FORMAT). All reference/input fields must match the
    supplied artifacts. The repeat format defines REPID as a native locus ID,
    RU as its exact motif, POS/END as its closed repeat region and REPCN as
    source-ordered sizes in repeat units. It is not an inferred caller dialect.
    The pinned LongTR profile admits original VCF 4.1 only and validates its
    interpreted header contracts; it retains authoritative locus calls and
    conditional source-bound member views. It does not authorize other 4.1
    dialects, rewrite headers, or replace reference/concordance verification.
    Standard symbolic SV and opt-in mosaic conversion need no such extra input.
    """

    source_bytes: bytes
    input_sha256: str
    qc_report: qc.VcfQcReport
    normalization: SmallVariantNormalizationResult
    class_inventory: VariantClassInventory
    policy_bytes: bytes = INVENTORY_POLICY_BYTES
    policy_sha256: str = INVENTORY_POLICY_SHA256
    source_format_bytes: bytes | None = None
    source_format_sha256: str | None = None

    def __post_init__(self):
        try:
            _Controlled.__post_init__(self)
            if (not self.source_bytes
                    or self.input_sha256 != hashlib.sha256(self.source_bytes).hexdigest()
                    or self.policy_bytes != INVENTORY_POLICY_BYTES
                    or self.policy_sha256 != hashlib.sha256(self.policy_bytes).hexdigest()
                    or self.qc_report.input_sha256 != self.input_sha256
                    or self.normalization.input_sha256 != self.input_sha256
                     or self.qc_report.reference_concordance_sha256
                    != self.normalization.reference_concordance_sha256
                     or (self.source_format_bytes is None) != (self.source_format_sha256 is None)
                    or self.source_format_bytes is not None and (
                        not self.source_format_bytes or len(self.source_format_bytes) > 65536
                        or hashlib.sha256(self.source_format_bytes).hexdigest() != self.source_format_sha256)):
                raise ValueError
        except (ValueError, TypeError, RuntimeError, AttributeError, OverflowError):
            raise AdmissionInventoryError("invalid or missing supplied-call input binding") from None


def _source_set_bytes(value):
    return (json.dumps(_plain(value), sort_keys=True, separators=(",", ":"), ensure_ascii=True,
                       allow_nan=False) + "\n").encode("ascii")


def _source_set_hash(value):
    return hashlib.sha256(_source_set_bytes(value)).hexdigest()


@dataclass(frozen=True, slots=True, repr=False)
class SuppliedSourceBinding(_Controlled):
    """One original file and its explicitly selected source sample; no relabeling."""

    input_sha256: str
    input_size_bytes: int
    sample_name: str
    reference_role: ReferenceRole
    reference_identity_id: str
    reference_manifest_sha256: str
    reference_concordance_sha256: str
    vcf_qc_sha256: str
    normalization_sha256: str
    class_inventory_sha256: str
    source_format_sha256: str | None
    inventory_sha256: str

    def __post_init__(self):
        _Controlled.__post_init__(self)
        _ids((self.sample_name,))
        _ids((self.reference_identity_id,))
        for value in (self.input_sha256, self.reference_manifest_sha256, self.reference_concordance_sha256,
                      self.vcf_qc_sha256, self.normalization_sha256, self.class_inventory_sha256, self.inventory_sha256):
            _sha(value)
        if self.source_format_sha256 is not None:
            _sha(self.source_format_sha256)
        if self.input_size_bytes <= 0:
            raise AdmissionInventoryError("invalid source-set file binding")


@dataclass(frozen=True, slots=True, repr=False)
class SuppliedSourceSetManifest(_Controlled):
    """Explicit input assignment, not biological identity or reference acceptance."""

    case_id: str
    sample_id: str
    members: tuple[SuppliedSourceBinding, ...]
    schema: Literal["genosapiens.supplied-source-set"] = SOURCE_SET_SCHEMA

    def __post_init__(self):
        _Controlled.__post_init__(self)
        _ids((self.case_id,))
        _ids((self.sample_id,))
        keys = tuple(m.input_sha256 for m in self.members)
        _ids(keys)
        if (not self.members or keys != tuple(sorted(keys)) or len({
            (m.reference_role, m.reference_identity_id, m.reference_manifest_sha256) for m in self.members
        }) != 1):
            raise AdmissionInventoryError("invalid source-set roster or reference bindings")

    @property
    def sha256(self):
        return _source_set_hash(self)

    def sample_for(self, input_sha256):
        return next(m.sample_name for m in self.members if m.input_sha256 == input_sha256)


@dataclass(frozen=True, slots=True, repr=False)
class SuppliedSourceSetInputs(_Controlled):
    """Original per-file inputs plus externally selected canonical manifest bytes/hash.

    File order is immaterial. Manifest members are canonically sorted by input
    hash. Every original sample label and file artifact remains unchanged. The
    caller's independent manifest pin defines completeness; a result cannot supply
    its own shortened roster. This does not accept read-derived producer batches.
    """

    sources: tuple[SuppliedCallInputs, ...]
    manifest_bytes: bytes
    manifest_sha256: str

    def __post_init__(self):
        try:
            _Controlled.__post_init__(self)
            _sha(self.manifest_sha256)
            if (not self.manifest_bytes or hashlib.sha256(self.manifest_bytes).hexdigest() != self.manifest_sha256
                    or self.manifest_bytes != _source_set_bytes(self.manifest)):
                raise ValueError
            keys = tuple(s.input_sha256 for s in self.sources)
            _ids(keys)
            if set(keys) != {m.input_sha256 for m in self.manifest.members}:
                raise ValueError
        except (ValueError, TypeError, RuntimeError, AttributeError, UnicodeError):
            raise AdmissionInventoryError("invalid or incomplete source-set input binding") from None

    @property
    def manifest(self):
        return _typed(_loads(self.manifest_bytes), SuppliedSourceSetManifest, decode=True)

    @property
    def input_sha256(self):
        """Collection-manifest identity, explicitly not a raw VCF checksum."""
        return self.manifest_sha256


SuppliedInputs = SuppliedCallInputs | SuppliedSourceSetInputs


@dataclass(frozen=True, slots=True, repr=False)
class AdmissionSourceRecord(_Controlled):
    source: SourcedRecord
    qc: qc.VcfRecordQc
    source_allele_ids: tuple[str, ...]
    uninterpreted_fields: tuple[tuple[str, str], ...]


@dataclass(frozen=True, slots=True, repr=False)
class SourceAlleleAccounting(_Controlled):
    """Projection accounting; additional native limits live in unresolved_entries.

    A derived sequence may coexist with an unresolved source-level repeat rule.
    Such a rule is not claimed non-applicable to the source, but is not applied
    to the concrete SequenceAllele as if it were a native RepeatAllele either.
    """

    source: SourceAllele
    source_classes: tuple[WgsClass, ...]
    normalized_allele_id: str | None
    derived_allele_ids: tuple[str, ...]
    status: Literal["derived", "unresolved", "non_variant"]
    required_rule_ids: tuple[str, ...]
    non_applicable_rules: tuple[tuple[str, str], ...]
    reasons: tuple[str, ...]
    read_validation_required: bool


@dataclass(frozen=True, slots=True, repr=False)
class UnresolvedSourceEntry(_Controlled):
    """No target is invented, even when a normalized sequence view is available."""

    entry_id: str
    source_record_id: str
    source_allele_id: str
    normalized_allele_id: str | None
    reasons: tuple[str, ...]
    required_rule_ids: tuple[str, ...]
    target_id: None = None
    status: Literal["unresolved"] = "unresolved"


@dataclass(frozen=True, slots=True, repr=False)
class AdmissionInventory(_Controlled):
    """Complete source and derived membership, not a candidate shortlist.

    ``records`` owns the full raw fields and bound QC; ``source_alleles`` owns
    every ALT-slot mapping and source-backed applicability reason. ``alleles``
    contains established sequence/native projections; ``unresolved_entries``
    must be carried forward even without a target. Nominal structural context
    does not establish endpoint precision, genotype validity or molecular effect.

    With ``source_set``, input_sha256 identifies its manifest; stage hashes bind
    ordered per-file digest collections, not synthesized reports. Raw per-file
    hashes live in source_set.members and original headers in source_headers.
    """

    input_sha256: str
    input_size_bytes: int
    policy_sha256: str
    reference_role: ReferenceRole
    reference_identity_id: str
    reference_manifest_sha256: str
    reference_concordance_sha256: str
    vcf_qc_sha256: str
    normalization_sha256: str
    class_inventory_sha256: str
    source_format_sha256: str | None
    headers: tuple[str, ...]
    records: tuple[AdmissionSourceRecord, ...]
    source_alleles: tuple[SourceAlleleAccounting, ...]
    alleles: tuple[AlleleEvidence, ...]
    unresolved_entries: tuple[UnresolvedSourceEntry, ...]
    required_rule_keys: tuple[RuleEvaluationKey, ...]
    schema: Literal["genosapiens.admission-inventory"] = SCHEMA
    source_set: SuppliedSourceSetManifest | None = None
    source_headers: tuple[tuple[str, tuple[str, ...]], ...] = ()

    @property
    def total_source_records(self) -> int:
        return len(self.records)

    @property
    def total_source_alt_slots(self) -> int:
        return len(self.source_alleles)

    @property
    def total_alternate_alleles(self) -> int:
        return sum(qc._alternate_allele_count(record.qc.alternates) for record in self.records)

    @property
    def source_input_sha256s(self):
        return (self.input_sha256,) if self.source_set is None else tuple(m.input_sha256 for m in self.source_set.members)


def _member_binding(inventory, sample_name):
    if any(sample_name not in {s.sample_name for s in r.qc.sample_calls} for r in inventory.records):
        raise AdmissionInventoryError("source-set selected sample unavailable")
    return SuppliedSourceBinding(inventory.input_sha256, inventory.input_size_bytes, sample_name,
        inventory.reference_role, inventory.reference_identity_id, inventory.reference_manifest_sha256,
        inventory.reference_concordance_sha256, inventory.vcf_qc_sha256, inventory.normalization_sha256,
        inventory.class_inventory_sha256, inventory.source_format_sha256,
        hashlib.sha256(admission_inventory_as_json(inventory).encode("ascii")).hexdigest())


def source_set_manifest_bytes(sources: tuple[SuppliedCallInputs, ...], *, case_id: str, sample_id: str,
                              sample_names: Mapping[str, str]) -> bytes:
    """Prepare canonical manifest bytes for independent selection, not self-approval.

    ``sample_names`` must explicitly cover each original input hash exactly once.
    The operator/case-input manifest owns these assignments. This helper replays
    file artifacts but does not infer that similarly named samples are one case.
    """
    try:
        sources = _typed(sources, tuple[SuppliedCallInputs, ...])
        keys = tuple(s.input_sha256 for s in sources)
        _ids(keys)
        if set(sample_names) != set(keys):
            raise ValueError
        members = tuple(_member_binding(build_admission_inventory(s), sample_names[s.input_sha256])
                        for s in sorted(sources, key=lambda s: s.input_sha256))
        return _source_set_bytes(SuppliedSourceSetManifest(case_id, sample_id, members))
    except (ValueError, TypeError, KeyError, AttributeError):
        raise AdmissionInventoryError("cannot bind supplied source-set manifest") from None


def _retain_enrichment(derived, alleles):
    if alleles is None:
        return tuple(derived)
    enriched = _typed(alleles, tuple[AlleleEvidence, ...])
    by_id = {a.allele_id: a for a in enriched}
    if len(by_id) != len(enriched) or set(by_id) != {a.allele_id for a in derived}:
        raise ValueError
    for baseline in derived:
        a = by_id[baseline.allele_id]
        if (a.native_identity_id != baseline.native_identity_id or a.identities != baseline.identities
                or a.calls != baseline.calls):
            raise ValueError
    return tuple(by_id[a.allele_id] for a in derived)


def _build_source_set(inputs, alleles):
    inputs = _typed(inputs, SuppliedSourceSetInputs)
    manifest = inputs.manifest
    originals = {s.input_sha256: s for s in inputs.sources}
    inventories, derived, slots, keys = [], [], [], set()
    for member in manifest.members:
        inventory = build_admission_inventory(originals[member.input_sha256])
        if _member_binding(inventory, member.sample_name) != member:
            raise ValueError
        inventories.append(inventory)
        # Normalized sequence IDs are reference/event-addressed and can collide
        # across files with different source calls. Retain each source view and
        # its unchanged identity graph; genotype evaluation relates exact events.
        names = {a.allele_id: "source-view-" + _source_set_hash((member.input_sha256, a.allele_id))
                 if isinstance(next(i.identity for i in a.identities if i.identity_id == a.native_identity_id), SequenceAllele)
                 else a.allele_id for a in inventory.alleles}
        derived.extend(replace(a, allele_id=names[a.allele_id]) for a in inventory.alleles)
        slots.extend(replace(s, derived_allele_ids=tuple(names[k] for k in s.derived_allele_ids))
                     for s in inventory.source_alleles)
        keys.update((names.get(subject, subject), rule) for subject, rule in inventory.required_rule_keys)
    # Per-file complete reference checks remain upstream, and the source-set
    # manifest requires one common reference identity.  VCF headers may declare
    # only the contigs used by that file, so a native-call supplement is allowed
    # to carry a strict subset of the primary VCF's declarations.  Still reject
    # incompatible lengths for any source-named contig declared by both files.
    declared_lengths = {}
    for source in inputs.sources:
        for contig in source.qc_report.declared_contigs:
            if contig.contig in declared_lengths and declared_lengths[contig.contig] != contig.length:
                raise ValueError
            declared_lengths[contig.contig] = contig.length
    _ids(tuple(a.allele_id for a in derived))
    first = inventories[0]

    def collection_hash(field):
        return _source_set_hash({"schema": SOURCE_SET_SCHEMA, "artifact": field,
                                 "members": tuple((i.input_sha256, getattr(i, field)) for i in inventories)})

    return AdmissionInventory(inputs.input_sha256, sum(i.input_size_bytes for i in inventories), INVENTORY_POLICY_SHA256,
        first.reference_role, first.reference_identity_id, first.reference_manifest_sha256,
        collection_hash("reference_concordance_sha256"), collection_hash("vcf_qc_sha256"),
        collection_hash("normalization_sha256"), collection_hash("class_inventory_sha256"),
        collection_hash("source_format_sha256"), (), tuple(r for i in inventories for r in i.records), tuple(slots),
        _retain_enrichment(derived, alleles), tuple(u for i in inventories for u in i.unresolved_entries),
        tuple(sorted(keys)), source_set=manifest,
        source_headers=tuple((i.input_sha256, i.headers) for i in inventories))


def _source_format(inputs):
    if inputs.source_format_bytes is None:
        return None
    profile = _loads(inputs.source_format_bytes.decode("utf-8"))
    normal = inputs.normalization
    binding = {"schema": SOURCE_FORMAT_SCHEMA, "input_sha256": inputs.input_sha256,
        "reference_role": normal.reference_role, "reference_identity_id": normal.reference_identity_id,
        "reference_manifest_sha256": normal.reference_manifest_sha256}
    if (type(profile) is not dict or set(profile) != set(binding) | {"mitochondrial_contigs", "repeat_identity"}
            or any(profile[k] != v for k, v in binding.items())
            or type(profile["mitochondrial_contigs"]) is not list
            or any(type(c) is not str for c in profile["mitochondrial_contigs"])
            or profile["repeat_identity"] not in (None, REPEAT_IDENTITY_FORMAT, longtr.FORMAT)):
        raise ValueError
    return profile


def _bound_records(inputs: SuppliedCallInputs):
    """Replay the existing QC primitives in memory against *every* original row."""
    content = (gzip.decompress(inputs.source_bytes)
               if inputs.source_bytes.startswith(b"\x1f\x8b") else inputs.source_bytes)
    text = content.decode("ascii")
    if any(ord(c) < 32 and c not in "\t\r\n" or ord(c) == 127 for c in text):
        raise ValueError
    lines = text.split("\n")
    if lines[-1] == "":
        lines.pop()
    lines = [line.removesuffix("\r") for line in lines]
    profile = _source_format(inputs)
    is_longtr = profile is not None and profile["repeat_identity"] == longtr.FORMAT
    versions = ("##fileformat=VCFv4.1",) if is_longtr else ("##fileformat=VCFv4.2", "##fileformat=VCFv4.3")
    if not lines or lines[0] not in versions:
        raise ValueError
    headers, records = [], []
    definitions, formats, contigs = {}, set(), {}
    samples = None
    for line in lines:
        if not line or any(ord(c) < 32 and c != "\t" or ord(c) == 127 for c in line):
            raise ValueError
        if line.startswith("#"):
            if samples is not None:
                raise ValueError
            headers.append(line)
            if line.startswith(("##INFO=", "##FORMAT=")):
                namespace = line[2:].split("=", 1)[0]
                match = re.fullmatch(r"##(?:INFO|FORMAT)=<(.+)>", line)
                attributes = qc._structured_attributes(match[1]) if match else None
                if attributes is None:
                    raise ValueError
                identifier = attributes.get("ID", "")
                number, dtype = attributes.get("Number"), attributes.get("Type")
                if (not identifier or any(c.isspace() for c in identifier)
                        or number is None or re.fullmatch(r"[0-9]+|[ARG.]", number) is None
                        or dtype not in {"Integer", "Float", "Flag", "Character", "String"}
                        or dtype == "Flag" and (namespace != "INFO" or number != "0")
                        or namespace == "INFO" and identifier in _BASE_INFO
                        and (number, dtype) != {"AN": ("1", "Integer"), "AF": ("A", "Float")}[identifier]
                        or (namespace, identifier) in definitions):
                    raise ValueError
                definitions[namespace, identifier] = (number, dtype)
                if namespace == "FORMAT":
                    formats.add(qc._parse_format_header(line))
            elif line.startswith("##contig="):
                contig, length = qc._parse_contig_header(line)
                if contig in contigs:
                    raise ValueError
                contigs[contig] = length
            elif line.startswith("#CHROM\t"):
                columns = line.split("\t")
                if columns[:9] != ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]:
                    raise ValueError
                samples = tuple(columns[9:])
                if not samples or len(set(samples)) != len(samples):
                    raise ValueError
            elif line.startswith("##fileformat=") and len(headers) != 1:
                raise ValueError
            elif not line.startswith("##"):
                raise ValueError
            continue
        if samples is None:
            raise ValueError
        columns = tuple(line.split("\t"))
        if len(columns) != 9 + len(samples) or any(not value for value in columns):
            raise ValueError
        original = SourceSmallVariantRecord(
            len(records) + 1, columns[0], int(columns[1]), columns[3], tuple(columns[4].split(",")),
        )
        record_id = source_record_id(inputs.input_sha256, original)
        info = qc._parse_info(columns[7])
        keys = () if columns[8] == "." else tuple(columns[8].split(":"))
        if (len(set(keys)) != len(keys) or not set(keys) <= formats
                or "GT" in keys and keys[0] != "GT"):
            raise ValueError
        calls = tuple(qc._sample_call(name, keys, value, qc._alternate_allele_count(original.alternates),
                                      {"PGT", "PID"} <= formats, inputs.qc_report.pid_phase_contract)
                      for name, value in zip(samples, columns[9:], strict=True))
        filter_status, filters = qc._filter(columns[6])
        technical_scopes = qc._technical_scopes(
            filter_status, calls, len(original.alternates)
        )
        required = qc._technical_scope_required(technical_scopes)
        quality = None if columns[5] == "." else float(columns[5])
        if quality is not None and (not math.isfinite(quality) or quality < 0):
            raise ValueError
        record = qc.VcfRecordQc(
            original.record_number, record_id, original.contig, original.position, original.reference,
            original.alternates, tuple(qc._variant_class(original.reference, a) for a in original.alternates),
            columns[6], filter_status, filters, qc._parse_optional_integer(info.get("AN")),
            qc._parse_callset_frequencies(info.get("AF")), calls, required, required,
            *technical_scopes,
            tuple((key, info[key]) for key in sorted(qc.CLASS_INFO_FIELDS.intersection(info))),
            quality,
            hashlib.sha256((line + "\n").encode("ascii")).hexdigest(),
            line + "\n",
            raw_quality=columns[5],
        )
        unknown = tuple(sorted(
            [("INFO", key) for key in info if key not in _BASE_INFO | qc.CLASS_INFO_FIELDS | _STRUCTURAL_INFO | {"REPID"}]
            + [("FORMAT", key) for key in keys if key not in _BASE_FORMAT | qc.CLASS_FORMAT_FIELDS]
        ))
        source = SourcedRecord("supplied-vcf", lines[0].split("=", 1)[1], record_id,
                               inputs.input_sha256, None, {
            "record_number": original.record_number, "columns": columns, "info": info,
            "format_keys": keys, "sample_values": dict(zip(samples, columns[9:], strict=True)),
            "quality": quality,
        })
        records.append(AdmissionSourceRecord(source, record,
            tuple(source_allele_id(record_id, n) for n in range(1, len(original.alternates) + 1)), unknown))
    if samples is None:
        raise ValueError
    rows = tuple(item.qc for item in records)
    class_definitions = tuple(sorted((ns, key, *contract) for (ns, key), contract in definitions.items()
                                    if key in (qc.CLASS_INFO_FIELDS if ns == "INFO" else qc.CLASS_FORMAT_FIELDS)))
    expected = qc.VcfQcReport(
        inputs.input_sha256, inputs.qc_report.reference_concordance_sha256, samples, len(samples),
        tuple(sorted(formats.intersection({"PGT", "PID", "PS"}))),
        tuple(qc.DeclaredContig(name, length) for name, length in contigs.items()),
        len(rows), sum(qc._alternate_allele_count(r.alternates) for r in rows), rows,
        qc._contig_summaries(rows, contigs), qc._sample_summaries(rows, samples),
        qc._named_counts(c.value for r in rows for c in r.variant_classes),
        qc._named_counts(r.filter_value for r in rows), class_definitions,
        inputs.qc_report.classification, inputs.qc_report.llm_policy,
        pid_phase_contract=inputs.qc_report.pid_phase_contract,
    )
    if expected != inputs.qc_report:
        raise ValueError
    if is_longtr:
        longtr.validate_header(headers, definitions)
    return tuple(headers), tuple(records), definitions


def _call_source(row, source, sample, *, read_required=False, interpretation=None):
    flags = qc.technical_flags_for(row.qc, sample, source.allele_number)
    return replace(row.source, source_fields={
        **row.source.source_fields, "source_allele_id": source.source_allele_id,
        "allele_number": source.allele_number, "vcf_qc_schema": qc.VCF_QC_SCHEMA,
        "technical_flags": flags,
        "quality_acceptance": qc.quality_acceptance_for(row.qc, sample, source.allele_number),
        "read_validation_required": read_required or bool(flags),
        "native_interpretation": interpretation,
    })


def supplied_structural_envelope(allele: AlleleEvidence) -> GenomicInterval | None:
    """Bound a replayed supplied SV's locus, without asserting affected bases.

    Both original endpoint confidence intervals are required from every identity
    and call contributor. Upstream inventory replay authenticates headers/bytes;
    this pure helper checks their retained correspondence. It never fills an
    absent interval from POS/END, absence of IMPRECISE, or dosage.
    """
    allele = _typed(allele, AlleleEvidence)
    native = next(i for i in allele.identities if i.identity_id == allele.native_identity_id)
    identity = native.identity
    if type(identity) is not StructuralAllele or identity.interval is None or not allele.calls:
        return None
    sources = [i.source for i in allele.identities] + [c.source for c in allele.calls]
    envelopes = []
    for source in sources:
        fields = source.source_fields or {}
        interpretation = fields.get("native_interpretation") or {}
        info = fields.get("info") or {}
        if (source.source_dataset != "supplied-vcf" or source.source_release not in ("VCFv4.2", "VCFv4.3")
                or interpretation.get("adapter") != SCHEMA
                or interpretation.get("policy_sha256") != INVENTORY_POLICY_SHA256
                or interpretation.get("interval_semantics") != "nominal-vcf-pos-end-context"
                or any(info.get(k) in (None, ".") for k in ("CIPOS", "CIEND"))):
            return None
        try:
            columns = fields["columns"]
            interval = identity.interval
            if (columns[0], int(columns[1]), columns[4], int(info["END"]), info["SVTYPE"]) != (
                interval.contig, interval.start, "<" + identity.sv_type.value + ">", interval.end, identity.sv_type.value
            ) or identity.sv_type not in (StructuralVariantType.DEL, StructuralVariantType.DUP, StructuralVariantType.CNV):
                raise ValueError
            points = []
            for name, position in (("CIPOS", interval.start), ("CIEND", interval.end)):
                if type(info[name]) is not str or re.fullmatch(r"-?[0-9]+,-?[0-9]+", info[name]) is None:
                    raise ValueError
                points.append(Breakpoint(interval.contig, position,
                    confidence_interval=tuple(int(v) for v in info[name].split(","))))
            if tuple(points) != identity.breakpoints:
                raise ValueError
            envelopes.append(GenomicInterval(interval.contig,
                interval.start + points[0].confidence_interval[0], interval.end + points[1].confidence_interval[1]))
        except (ValueError, TypeError, KeyError, IndexError):
            raise AdmissionInventoryError("invalid supplied structural envelope binding") from None
    # Cover every supported contributor, not a convenient point/intersection.
    return GenomicInterval(identity.interval.contig, min(b.start for b in envelopes), max(b.end for b in envelopes))


def build_admission_inventory(
    inputs: SuppliedInputs, *, alleles: tuple[AlleleEvidence, ...] | None = None,
) -> AdmissionInventory:
    """Build complete accounting, optionally retaining a complete enriched allele set.

    ``None`` builds the initial unannotated projection. An explicitly supplied
    tuple must contain exactly all generated allele IDs and unchanged identities
    and bound CallEvidence. Missing provided alleles are errors, not unsupported
    records. No target/population/annotation filter selects membership here.

    Required keys are ``(derived allele_id, rule_id)`` or, for unresolved native
    conversion, ``(source_allele_id, rule_id)``. A no-target rule outcome still
    owes that key. Pure non-variant records have source-backed non-applicability.
    """
    try:
        if type(inputs) is SuppliedSourceSetInputs:
            return _build_source_set(inputs, alleles)
        inputs = _typed(inputs, SuppliedCallInputs)
        headers, records, definitions = _bound_records(inputs)
        normal = inputs.normalization
        role = ReferenceRole(normal.reference_role)
        lengths = {c.contig: c.length for c in inputs.qc_report.declared_contigs}
        source_format = _source_format(inputs)
        mitochondrial_contigs = set()
        if source_format is not None:
            mitochondrial_contigs = set(source_format["mitochondrial_contigs"])
            if (len(mitochondrial_contigs) != len(source_format["mitochondrial_contigs"])
                    or not mitochondrial_contigs <= set(lengths)):
                raise ValueError
        class_inventory = build_variant_class_inventory(
            inputs.qc_report, allele_balance_screen=inputs.class_inventory.allele_balance_screen.enabled,
        )
        if inputs.class_inventory != class_inventory:
            raise ValueError
        mapped = {m.source.source_allele_id: m for m in normal.source_to_normalized}
        quarantined = {q.source.source_allele_id: q for q in normal.quarantined_source_alleles}
        actual_sources = {**{key: m.source for key, m in mapped.items()},
                          **{key: q.source for key, q in quarantined.items()}}
        expected_sources = {}
        for row in records:
            r = row.qc
            original = SourceSmallVariantRecord(r.record_number, r.contig, r.position, r.reference, r.alternates)
            for n, alt in enumerate(r.alternates, 1):
                key = source_allele_id(r.source_record_id, n)
                expected_sources[key] = SourceAllele(key, r.source_record_id, inputs.input_sha256, original, n, alt)
        if (actual_sources != expected_sources or normal.total_source_records != len(records)
                or normal.total_source_alleles != len(expected_sources)):
            raise ValueError
        screen_candidates = {(c.source_record_id, c.sample_call.sample_name): c
                             for c in class_inventory.allele_balance_screen.calls
                             if c.outcome is ScreenOutcome.CANDIDATE}
        normalized_by_id = {a.normalized_allele_id: a for a in normal.normalized_alleles}
        norm_sha = hashlib.sha256(variant_normalization_as_json(normal).encode("ascii")).hexdigest()
        class_sha = hashlib.sha256(variant_inventory_as_json(class_inventory).encode("ascii")).hexdigest()
        accounting, derived, additional_unresolved = [], [], []
        groups = defaultdict(list)
        for row in records:
            r = row.qc
            if source_format is not None and source_format["repeat_identity"] == longtr.FORMAT:
                sources = tuple(expected_sources[k] for k in row.source_allele_ids)
                blocked = tuple(q.reason.value for s in sources if (q := quarantined.get(s.source_allele_id)) is not None
                                and q.reason in (QuarantineReason.REFERENCE_MISMATCH, QuarantineReason.REFERENCE_UNRESOLVABLE))
                views, links, limits = (), {s.source_allele_id: [] for s in sources}, {s.source_allele_id: [] for s in sources}
                try:
                    if blocked:
                        raise ValueError
                    views, links, limits = longtr.build_views(row, sources, profile=source_format,
                        profile_sha256=inputs.source_format_sha256, role=role, reference_id=normal.reference_identity_id,
                        reference_sha256=normal.reference_manifest_sha256, normalized=normalized_by_id, mappings=mapped)
                except (ValueError, TypeError, KeyError, IndexError, OverflowError):
                    limits = {s.source_allele_id: list(blocked or ("longtr-source-interpretation-unresolved",)) for s in sources}
                derived.extend(views)
                for s in sources:
                    key = s.source_allele_id
                    required = _CLASS_RULES["repeat"]
                    accounting.append(SourceAlleleAccounting(s, (WgsClass.REPEAT,),
                        mapped[key].normalized_allele_id if key in mapped else None, tuple(links[key]),
                        "unresolved" if limits[key] else "derived", required,
                        tuple((rule, "source-class:repeat") for rule in RULE_IDS if rule not in required),
                        tuple(dict.fromkeys(limits[key])), any(
                            qc.technical_flags_for(r, sample, s.allele_number)
                            for sample in r.sample_calls
                        ) or any(
                            c.source.source_fields["read_validation_required"]
                            for a in views for c in a.calls or ()
                        )))
                continue
            represented, unsupported, _ = _record_classes(r, definitions, lengths[r.contig])
            info = row.source.source_fields["info"]
            svtype = info.get("SVTYPE")
            native_repeat_hint = (svtype in {"STR", "VNTR"}
                                  or any(dict(c.class_format_fields).get("REPCN") not in (None, ".") for c in r.sample_calls))
            repeat_hint = "RU" in info and info["RU"] != "." or native_repeat_hint
            cn_hint = any(dict(c.class_format_fields).get("CN") not in (None, ".") for c in r.sample_calls)
            structural_hint = ("SVTYPE" in info and svtype not in (".", "STR", "VNTR") or cn_hint
                               or any(key in info and info[key] != "." for key in ("CIPOS", "CIEND", "IMPRECISE")))
            screen_record = any((r.source_record_id, c.sample_name) in screen_candidates for c in r.sample_calls)
            for key in row.source_allele_ids:
                original = expected_sources[key]
                alt = original.alternate
                kind = r.variant_classes[original.allele_number - 1]
                symbolic_type = alt[1:-1] if alt.startswith("<") and alt.endswith(">") else None
                repeat = repeat_hint or REPEAT_SYMBOL.fullmatch(alt) is not None
                structural = (structural_hint or symbolic_type in SV_TYPES
                              or CN_SYMBOL.fullmatch(alt) is not None or BREAKEND.fullmatch(alt) is not None)
                # Tandem-repeat annotations describe sequence context; RU does
                # not replace a concrete REF/ALT with a native repeat identity.
                if {"MOTIF", "START", "INEXACT_ALLELE"} <= set(info):
                    # Known repeat-dialect markers without its bound profile are
                    # unresolved native evidence, not ordinary exact SNV calls.
                    rule_class, classes = "repeat", {WgsClass.REPEAT}
                elif kind in SEQUENCE_CLASSES and not structural:
                    rule_class, classes = "sequence", {WgsClass.SMALL_VARIANT}
                elif repeat and structural:
                    rule_class, classes = "unknown", represented | unsupported
                elif repeat:
                    rule_class, classes = "repeat", {WgsClass.REPEAT}
                elif structural:
                    rule_class, classes = "structural", {WgsClass.STRUCTURAL}
                    if cn_hint or symbolic_type in {"DEL", "DUP", "CNV"}:
                        classes.add(WgsClass.COPY_NUMBER)
                elif kind is qc.VariantClass.NON_VARIANT:
                    rule_class, classes = "non-variant", set()
                else:
                    rule_class, classes = "unknown", represented | unsupported
                if r.contig in mitochondrial_contigs or r.contig in {"M", "chrM", "MT"}:
                    classes.add(WgsClass.MITOCHONDRIAL)
                if screen_record:
                    classes.add(WgsClass.MOSAIC)
                additional_repeat = rule_class == "sequence" and native_repeat_hint
                if additional_repeat:
                    classes.add(WgsClass.REPEAT)
                reasons = []
                quarantine = quarantined.get(key)
                if quarantine is not None and quarantine.reason in {
                    QuarantineReason.REFERENCE_MISMATCH, QuarantineReason.REFERENCE_UNRESOLVABLE,
                }:
                    reasons.append(quarantine.reason.value)
                normalized_id = mapped[key].normalized_allele_id if key in mapped else None
                derived_ids, native_calls = [], []
                interpretation = {"adapter": SCHEMA, "policy_sha256": inputs.policy_sha256,
                                  "source_format_sha256": inputs.source_format_sha256,
                                  "source_format": source_format}
                if not reasons and rule_class == "sequence":
                    if normalized_id is None:
                        reasons.append("sequence-normalization-unresolved")
                    elif WgsClass.MITOCHONDRIAL in classes and r.contig not in mitochondrial_contigs:
                        reasons.append("mitochondrial-source-format-binding-required")
                    elif r.contig in mitochondrial_contigs or screen_record:
                        n = normalized_by_id[normalized_id]
                        sequence = SequenceAllele(role, normal.reference_identity_id, n.contig, n.position, n.reference, n.alternate)
                        for sample in r.sample_calls:
                            support = ReadSupport(sample.depth, sample.allelic_depths, original.allele_number)
                            screen = screen_candidates.get((r.source_record_id, sample.sample_name))
                            if r.contig in mitochondrial_contigs:
                                identity = MitochondrialAllele(sequence, None, support, sample.genotype, sample.ploidy)
                            elif screen is not None:
                                identity = MosaicAllele(sequence, screen.allele_balance, support, sample.genotype, sample.ploidy)
                            else:
                                identity = sequence
                            native_calls.append((sample, identity, screen))
                    else:
                        derived_ids.append(normalized_id)
                        groups[normalized_id].append(original)
                elif not reasons and rule_class in {"structural", "repeat"}:
                    try:
                        end = info.get("END")
                        if (len(r.alternates) != 1 or definitions.get(("INFO", "END")) != ("1", "Integer")
                                or type(end) is not str or not end.isdigit()
                                or not r.position <= int(end) <= lengths[r.contig]):
                            raise ValueError
                        interval = GenomicInterval(r.contig, r.position, int(end))
                        if rule_class == "structural":
                            from genosapiens.supplied_structural import parse_structural_record, INTERVAL_SEMANTICS
                            interpretation.update(interval_semantics=INTERVAL_SEMANTICS, interval_precision="unknown")
                            native_calls.extend((sample, identity, None) for sample, identity in parse_structural_record(
                                r, definitions, lengths[r.contig], role=role, reference_id=normal.reference_identity_id,
                                vcf_version=row.source.source_release))
                        else:
                            if (source_format is None or source_format["repeat_identity"] != REPEAT_IDENTITY_FORMAT
                                    or REPEAT_SYMBOL.fullmatch(alt) is None or svtype not in (None, ".", "STR", "VNTR")
                                    or "SVTYPE" in info and definitions.get(("INFO", "SVTYPE")) != ("1", "String")
                                    or any(definitions.get((ns, name)) != contract for ns, name, contract in (
                                        ("INFO", "REPID", ("1", "String")), ("INFO", "RU", ("1", "String")),
                                        ("FORMAT", "REPCN", ("1", "String"))))
                                    or type(info.get("REPID")) is not str or info["REPID"] in {"", "."}
                                    or any(c.isspace() for c in info["REPID"])
                                    or type(info.get("RU")) is not str or re.fullmatch(r"[ACGT]+", info["RU"]) is None):
                                raise ValueError
                            interpretation.update(interval_semantics="source-format-closed-repeat-region",
                                                  repeat_identity=REPEAT_IDENTITY_FORMAT, locus_id_source="INFO/REPID",
                                                  repeat_representation="locus_call",
                                                  individual_allele_correspondence="unresolved")
                            for sample in r.sample_calls:
                                sizes = dict(sample.class_format_fields).get("REPCN")
                                counts = None
                                if sizes not in (None, "."):
                                    try:
                                        if (type(sizes) is not str
                                                or re.fullmatch(r"(?:[0-9]+|\.)(?:/(?:[0-9]+|\.))*", sizes) is None):
                                            raise ValueError
                                        counts = tuple(RepeatSize(None if v == "." else int(v), RepeatSizeUnit.REPEAT_UNITS)
                                                       for v in sizes.split("/"))
                                    except (ValueError, OverflowError):
                                        reasons.append("repeat-size-measurements-unresolved")
                                identity = RepeatAllele(role, normal.reference_identity_id, info["REPID"], interval,
                                                        info["RU"], counts, sample.genotype, sample.ploidy)
                                native_calls.append((sample, identity, None))
                    except (ValueError, TypeError, KeyError, OverflowError):
                        native_calls = []
                        reasons.append("unsupported-structural-source-format" if rule_class == "structural"
                                       else "unsupported-or-unbound-repeat-source-format")
                elif not reasons and rule_class == "unknown":
                    reasons.append("unsupported-native-source-semantics")

                for sample_number, (sample, identity, screen) in enumerate(native_calls, 1):
                    allele_id = "allele-" + hashlib.sha256(json.dumps(
                        [SCHEMA, original.source_allele_id, sample_number, identity.kind.value],
                        separators=(",", ":"), ensure_ascii=True).encode("ascii")).hexdigest()[:24]
                    bound_interpretation = {**interpretation, "sample_name": sample.sample_name,
                        "normalization_sha256": norm_sha if normalized_id is not None else None,
                        "normalization_mapping": asdict(mapped[key]) if key in mapped else None,
                        "screen": _plain(screen) if screen is not None else None,
                        "screen_inventory_sha256": class_sha if screen is not None else None}
                    if isinstance(identity, RepeatAllele):
                        raw_sizes = dict(sample.class_format_fields).get("REPCN")
                        bound_interpretation["size_measurement_status"] = (
                            "unknown" if raw_sizes in (None, ".") else "unresolved"
                            if identity.allele_sizes is None else "retained")
                    source = _call_source(row, original, sample, read_required=screen is not None,
                                          interpretation=bound_interpretation)
                    if normalized_id is None:
                        source_identity = identity
                    else:
                        sequence = SequenceAllele(role, normal.reference_identity_id, r.contig, r.position, r.reference, alt)
                        source_identity = replace(identity, sequence=sequence) if type(identity) in (MosaicAllele, MitochondrialAllele) else sequence
                    identities = [IdentityMapping(key, IdentityStage.SOURCE, source_identity,
                                                  normal.reference_manifest_sha256, source)]
                    parents = (key,)
                    if normalized_id is not None:
                        sequence = identity.sequence if type(identity) in (MosaicAllele, MitochondrialAllele) else identity
                        identities.append(IdentityMapping(normalized_id, IdentityStage.NORMALIZED, sequence,
                            normal.reference_manifest_sha256, source, parents, VARIANT_NORMALIZATION_SCHEMA))
                        parents = (normalized_id,)
                    native_id = "native-" + allele_id
                    identities.append(IdentityMapping(native_id, IdentityStage.NATIVE, identity,
                                                      normal.reference_manifest_sha256, source, parents, SCHEMA))
                    call = CallEvidence(source, sample,
                        None if r.filter_status is qc.FilterStatus.NOT_APPLIED else r.filters,
                        row.source.source_fields["quality"])
                    derived.append(AlleleEvidence(allele_id, native_id, tuple(identities), (call,)))
                    derived_ids.append(allele_id)

                status = "unresolved" if reasons else "non_variant" if rule_class == "non-variant" else "derived"
                required = _CLASS_RULES[rule_class]
                accounting.append(SourceAlleleAccounting(
                    original, tuple(sorted(classes, key=lambda c: c.value)), normalized_id,
                    tuple(derived_ids), status, required,
                    tuple((rule, "source-class:" + rule_class) for rule in RULE_IDS
                          if rule not in required and not (additional_repeat and rule == RULE_IDS[3])),
                    tuple(dict.fromkeys(reasons)), screen_record or any(
                        qc.technical_flags_for(r, sample, original.allele_number)
                        for sample in r.sample_calls
                    ),
                ))
                if additional_repeat:
                    additional_unresolved.append(UnresolvedSourceEntry(
                        "unresolved-repeat-" + key, original.source_record_id, key, normalized_id,
                        ("additional-native-repeat-identity-unresolved",), (RULE_IDS[3],),
                    ))

        accounting = tuple(accounting)
        row_by_id = {r.qc.source_record_id: r for r in records}
        # Preserve the VepRequest.from_normalization identity convention, but
        # index its source groups once instead of rescanning the full result for
        # each allele. This projection's own transformation is identified below.
        for normalized in normal.normalized_alleles:
            key = normalized.normalized_allele_id
            if key not in groups:
                continue
            mappings = tuple(mapped[s.source_allele_id] for s in groups[key])
            sequence = SequenceAllele(role, normal.reference_identity_id, normalized.contig,
                                      normalized.position, normalized.reference, normalized.alternate)
            provenance = SourcedRecord("genosapiens.variant-normalization", "current",
                                      key, norm_sha, None, {
                "normalized": {**asdict(normalized), "kind": normalized.kind.value},
                "source_to_normalized": [asdict(m) for m in mappings],
                "input_sha256": normal.input_sha256,
                "reference_concordance_sha256": normal.reference_concordance_sha256,
            })
            identities, calls = [], []
            for mapping in mappings:
                s = mapping.source
                r = s.record
                identities.append(IdentityMapping(
                    s.source_allele_id, IdentityStage.SOURCE,
                    SequenceAllele(role, normal.reference_identity_id, r.contig, r.position, r.reference, s.alternate),
                    normal.reference_manifest_sha256,
                    replace(provenance, source_record_id=s.source_allele_id, source_fields=asdict(s)),
                ))
                row = row_by_id[s.source_record_id]
                for sample in row.qc.sample_calls:
                    call_source = _call_source(row, s, sample)
                    calls.append(CallEvidence(call_source, sample,
                        None if row.qc.filter_status is qc.FilterStatus.NOT_APPLIED else row.qc.filters,
                        row.source.source_fields["quality"]))
            identities.append(IdentityMapping(key, IdentityStage.NORMALIZED, sequence,
                normal.reference_manifest_sha256, provenance, tuple(i.identity_id for i in identities), VARIANT_NORMALIZATION_SCHEMA))
            native_id = "native-" + key
            identities.append(IdentityMapping(native_id, IdentityStage.NATIVE, sequence,
                normal.reference_manifest_sha256, provenance, (key,), SCHEMA))
            derived.append(AlleleEvidence(key, native_id, tuple(identities), tuple(calls)))
        derived = _retain_enrichment(derived, alleles)
        unresolved = tuple(UnresolvedSourceEntry(
            "unresolved-" + i.source.source_allele_id, i.source.source_record_id, i.source.source_allele_id,
            i.normalized_allele_id, i.reasons, i.required_rule_ids,
        ) for i in accounting if i.status == "unresolved") + tuple(additional_unresolved)
        keys = tuple(sorted(
            {(subject, rule) for i in accounting for subject in i.derived_allele_ids for rule in i.required_rule_ids}
            | {(i.source_allele_id, rule) for i in unresolved for rule in i.required_rule_ids}
        ))
        return AdmissionInventory(
            inputs.input_sha256, len(inputs.source_bytes), inputs.policy_sha256,
            role, normal.reference_identity_id, normal.reference_manifest_sha256,
            normal.reference_concordance_sha256,
            hashlib.sha256(qc.vcf_qc_as_json(inputs.qc_report).encode("ascii")).hexdigest(), norm_sha,
            class_sha, inputs.source_format_sha256,
            headers, records, accounting, tuple(derived), unresolved, keys,
        )
    except (ValueError, TypeError, RuntimeError, AttributeError, KeyError, IndexError,
            OverflowError, EOFError, OSError, zlib.error):
        raise AdmissionInventoryError("invalid or incomplete supplied-call inventory inputs") from None


def validate_admission_inventory(
    inventory: AdmissionInventory, *, inputs: SuppliedInputs,
    alleles: tuple[AlleleEvidence, ...] | None = None,
) -> None:
    """Rebuild against independent original inputs, never the emitted subset.

    Supply the independently held complete enriched tuple if annotations were
    attached. A hash on an internally consistent but shortened inventory is not
    evidence of complete source accounting.
    """
    expected = build_admission_inventory(inputs, alleles=alleles)
    try:
        inventory = _typed(inventory, AdmissionInventory)
    except (ValueError, TypeError, RuntimeError, AttributeError, OverflowError):
        raise AdmissionInventoryError("invalid supplied-call inventory structure") from None
    if inventory != expected:
        raise AdmissionInventoryError("supplied-call inventory differs from bound inputs")


def validate_admission_accounting(
    inventory: AdmissionInventory, *, inputs: SuppliedInputs,
    source_record_ids: tuple[str, ...], allele_ids: tuple[str, ...],
    unresolved_entry_ids: tuple[str, ...], evaluated_rule_keys: tuple[RuleEvaluationKey, ...],
    alleles: tuple[AlleleEvidence, ...] | None = None,
) -> tuple[UnresolvedSourceEntry, ...]:
    """Reconcile consumer accounting *in addition to* independently replayed outcomes.

    Project one key per complete allele/rule assessment, including no-target
    outcomes. Per-target decisions may be many-to-one at this boundary; the
    evaluator owns validating their full target/evidence coverage. Project each
    key once. Missing, extra or duplicate subjects/keys are validation errors.
    Return all unresolved source limits; an empty return does NOT establish
    partner-search completion or clear unresolved annotation/rule outcomes.
    """
    validate_admission_inventory(inventory, inputs=inputs, alleles=alleles)
    try:
        pairs = (
            (source_record_ids, tuple(r.qc.source_record_id for r in inventory.records), tuple[str, ...]),
            (allele_ids, tuple(a.allele_id for a in inventory.alleles), tuple[str, ...]),
            (unresolved_entry_ids, tuple(u.entry_id for u in inventory.unresolved_entries), tuple[str, ...]),
            (evaluated_rule_keys, inventory.required_rule_keys, tuple[RuleEvaluationKey, ...]),
        )
        for actual, expected, hint in pairs:
            actual = _typed(actual, hint)
            if len(actual) != len(set(actual)) or set(actual) != set(expected):
                raise ValueError
    except (ValueError, TypeError, RuntimeError, AttributeError):
        raise AdmissionInventoryError("incomplete or invalid admission accounting") from None
    return inventory.unresolved_entries


def admission_inventory_as_json(inventory: AdmissionInventory) -> str:
    """Deterministic controlled serialization, not original-input validation."""
    try:
        value = _typed(inventory, AdmissionInventory)
        return json.dumps({"classification": "controlled", "llm_policy": "no-llm", **_plain(value)},
                          sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False) + "\n"
    except (ValueError, TypeError, RuntimeError, AttributeError, OverflowError):
        raise AdmissionInventoryError("invalid admission inventory serialization") from None
