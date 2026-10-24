"""Source-indexed genotype properties and pair phase, not hypothesis qualification.

The public evaluator replays complete supplied-call accounting. No ALT, annotation,
dosage baseline, population frequency, or inheritance label establishes carriage.
These properties do not certify model completeness, partner-search coverage,
reference acceptance, technical validation, or competitive eligibility.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping
from dataclasses import dataclass
import hashlib
from itertools import product
import re
from typing import Literal

from genosapiens.admission_inventory import (
    AdmissionInventory, SuppliedInputs, SuppliedSourceSetManifest, admission_inventory_as_json,
    build_admission_inventory,
)
from genosapiens.allele_identity import MitochondrialAllele, MosaicAllele, RepeatAllele, SequenceAllele, StructuralAllele
from genosapiens.evidence_model import HypothesisPhase, _Controlled, _typed
from genosapiens.evidence_serialization import _canonical, _loads
from genosapiens.vcf_contracts import POLICY_SHA256 as VCF_POLICY_SHA256
from genosapiens.vcf_qc import GENOTYPE, SampleCallQc, technical_flags_for
from genosapiens import longtr_admission as longtr


SCHEMA = "genosapiens.genotype-properties"
POLICY_BYTES = (_canonical({
    "schema": "genosapiens.genotype-properties-policy",
    "genotype": "designated-sample; original-ALT-index; partial-copy-bounds; no-ALT-carriage-inference",
    "copies": "intersect-supported-source-bounds; never-sum-duplicate-records",
    "ploidy": "source-GT-slots; lone-missing-GT-unknown; no-sex-or-locus-baseline-inference",
    "phase": "diploid-single-copy-members; same-record-distinct-copies-or-shared-source-block",
    "distinct_members": "conservative-sequence/event-anchor-guard; changed-measurements-do-not-establish-distinct-members",
    "blocks": "input/sample/contig/namespace/key; PS-and-PID-are-distinct-namespaces",
    "source_set": "explicit-selected-source-samples; original-labels; no-cross-file-phase-block",
    "exact_event_views": "intersect-source-copy-bounds-without-summing; propagate-conflicts-and-obligations; retain-view-lineage",
    "implicit_PS": "VCFv4.2/v4.3-phased-GT-without-PS; same-contig-only",
    "PGT": "source-bound-producer-contract; complete-GT-multiset-concordance; never-fills-missing-GT",
    "vcf_evidence_policy_sha256": VCF_POLICY_SHA256,
    "conflicts": "retain-opposing-pair-relations-and-source-genotype-conflicts",
    "repeat": "unbound-locus-copies-unknown; source-replayed-LongTR-individual-projections-use-ALT-index",
    "repeat_parent": "distinct-source-ALT-copy-counts-are-not-conflicting-estimates-of-one-member",
    "longtr_phase": "GT-pipe-is-source-local; no-cross-record-PS-or-PGT-interpretation",
    "low_balance": "replayed-existing-screen-only; unchanged-GT; mandatory-validation-retained",
    "scope": "supplied-GT-properties-only; no-parental-origin-or-model-completeness",
}) + "\n").encode("ascii")
POLICY_SHA256 = hashlib.sha256(POLICY_BYTES).hexdigest()


class GenotypePropertiesError(ValueError):
    """Static diagnostics only; supplied values must not reach logs."""


class GenotypeCapacityError(GenotypePropertiesError):
    """Execution capacity failure, never a smaller successful property set."""


class _BoundedPhaseBases(list):
    def __init__(self, maximum, initial=()):
        if type(maximum) is not int or maximum < 0:
            raise GenotypePropertiesError("invalid phase basis capacity")
        self.maximum = maximum
        super().__init__()
        for value in initial:
            self.append(value)

    def append(self, value):
        if len(self) >= self.maximum:
            raise GenotypeCapacityError("phase basis capacity exceeded")
        super().append(value)


@dataclass(frozen=True, slots=True, repr=False)
class PhaseOrientation(_Controlled):
    namespace: Literal["GT/PS", "GT/implicit-PS", "PGT/PID"]
    block: str | None
    haplotype: int
    contract_sha256: str | None = None


@dataclass(frozen=True, slots=True, repr=False)
class SourceGenotypeProperties(_Controlled):
    source_record_id: str
    source_allele_id: str
    allele_number: int
    contig: str
    sample: SampleCallQc
    derived_allele_ids: tuple[str, ...]
    source_status: Literal["derived", "unresolved", "non_variant"]
    reported_ploidy: int | None
    alternate_present: bool | None
    minimum_alt_copies: int | None
    maximum_alt_copies: int | None
    orientations: tuple[PhaseOrientation, ...]
    phase_conflict: bool
    phase_limits: tuple[str, ...]
    read_validation_required: bool
    input_sha256: str | None = None


@dataclass(frozen=True, slots=True, repr=False)
class AlleleGenotypeProperties(_Controlled):
    allele_id: str
    native_identity_id: str
    source_allele_ids: tuple[str, ...]
    alternate_present: bool | None
    minimum_alt_copies: int | None
    maximum_alt_copies: int | None
    reported_ploidies: tuple[int, ...]
    conflicting_source_genotypes: bool
    phase_conflict_bases: tuple[PairPhaseBasis, ...]
    screened_low_balance_candidate: bool
    read_validation_required: bool
    limits: tuple[str, ...]
    event_source_allele_ids: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True, repr=False)
class PairPhaseBasis(_Controlled):
    source_allele_ids: tuple[str, str]
    relation: HypothesisPhase
    rule: Literal["same-record-distinct-copies", "shared-phase-block", "source-genotype-conflict"]
    namespace: str | None
    block: str | None


@dataclass(frozen=True, slots=True, repr=False)
class PairPhaseProperties(_Controlled):
    allele_ids: tuple[str, str]
    phase: HypothesisPhase
    bases: tuple[PairPhaseBasis, ...]
    limits: tuple[str, ...]


@dataclass(frozen=True, slots=True, repr=False)
class SuppliedGenotypeProperties(_Controlled):
    inventory_sha256: str
    input_sha256: str
    reference_manifest_sha256: str
    sample_name: str
    policy_sha256: str
    source_properties: tuple[SourceGenotypeProperties, ...]
    allele_properties: tuple[AlleleGenotypeProperties, ...]
    pair_phases: tuple[PairPhaseProperties, ...]
    source_set: SuppliedSourceSetManifest | None = None
    ledger_source_samples: tuple[tuple[str, str], ...] = ()

    def __post_init__(self):
        _Controlled.__post_init__(self)
        if (
            self.ledger_source_samples != tuple(sorted(self.ledger_source_samples))
            or len({key for key, _sample in self.ledger_source_samples})
            != len(self.ledger_source_samples)
            or any(not key or not sample for key, sample in self.ledger_source_samples)
            or self.source_set is not None and self.ledger_source_samples
        ):
            raise GenotypePropertiesError("invalid ledger genotype source binding")


def source_sample_name(properties, input_sha256):
    """Resolve explicit case/sample assignment, retaining original source labels."""
    if properties.ledger_source_samples:
        if input_sha256 == properties.input_sha256:
            return properties.sample_name
        try:
            return next(
                sample for key, sample in properties.ledger_source_samples
                if key == input_sha256
            )
        except StopIteration:
            raise ValueError("unbound ledger genotype source") from None
    if properties.source_set is None:
        if input_sha256 != properties.input_sha256:
            raise ValueError("unbound genotype input")
        return properties.sample_name
    manifest = properties.source_set
    if manifest.sha256 != properties.input_sha256 or manifest.sample_id != properties.sample_name:
        raise ValueError("invalid genotype source-set binding")
    try:
        return manifest.sample_for(input_sha256)
    except StopIteration:
        raise ValueError("unbound genotype source-set member") from None


def selected_source_calls(allele, properties, *, _sources=None):
    """Check actual per-file call correspondence after upstream full-source replay."""
    calls = tuple(c for c in allele.calls or () if c.sample is not None
                  and c.sample.sample_name == source_sample_name(properties, c.source.source_sha256))
    sources = (_sources if _sources is not None else
               {s.source_allele_id: s for s in properties.source_properties})
    for call in calls:
        source = sources[call.source.source_fields["source_allele_id"]]
        expected_input = source.input_sha256 or properties.input_sha256
        if ((not properties.ledger_source_samples
             and call.source.source_sha256 != expected_input)
                or call.sample != source.sample
                or allele.allele_id not in source.derived_allele_ids):
            raise ValueError("genotype call differs from selected source")
    return calls


def _tokens(value):
    return () if value is None else tuple(None if t == "." else int(t) for t in re.split(r"[/|]", value))


def _phase_event(identity):
    """Conservative distinct-member guard, never a replacement native identity.

    Changed GT, fraction, depth, dosage, or uncertainty alone cannot make two
    descriptions of the same sequence/event anchors into distinct pair members.
    Native evidence and uncertainty are never merged or rewritten by this key.
    """
    if isinstance(identity, (MitochondrialAllele, MosaicAllele)):
        identity = identity.sequence
    if isinstance(identity, SequenceAllele):
        return ("sequence", identity)
    if isinstance(identity, StructuralAllele):
        return ("structural", identity.reference_role, identity.reference_identity_id,
                identity.sv_type, identity.interval,
                tuple((b.contig, b.position, b.orientation) for b in identity.breakpoints), identity.exact_sequence)
    return ("repeat", identity.reference_role, identity.reference_identity_id,
            identity.locus_id, identity.interval, identity.repeat_unit)


def _allele_event(allele):
    native = next(i.identity for i in allele.identities if i.identity_id == allele.native_identity_id)
    if longtr.is_projection(allele):
        normalized = tuple(i.identity for i in allele.identities if i.stage.value == "normalized")
        return _phase_event(normalized[0] if normalized else native.exact_sequence)
    return _phase_event(native)


def assess_pair_phase(left, right, *, left_properties, right_properties,
                      left_sources, right_sources, unprojected_repeats, maximum_phase_bases=None):
    """Evaluate one closed distinct-event pair using the shared phase kernel.

    Domain membership, sample/source authentication and parent exclusion belong
    to the supplying genotype owner. Only the two members' exact source facts
    are resident here; all source bases and their ordering survive unchanged.
    """
    props = {left: left_properties, right: right_properties}
    contributors = {left: left_sources, right: right_sources}
    initial = (*props[left].phase_conflict_bases, *props[right].phase_conflict_bases)
    bases = list(initial) if maximum_phase_bases is None else _BoundedPhaseBases(maximum_phase_bases, initial)
    for a, b in product(contributors[left], contributors[right]):
        source_ids = (a.source_allele_id, b.source_allele_id)
        if a.phase_conflict or b.phase_conflict:
            bases.append(PairPhaseBasis(source_ids, HypothesisPhase.CONFLICTING,
                                        "source-genotype-conflict", None, None))
        if (a.contig != b.contig or a.reported_ploidy != 2 or b.reported_ploidy != 2
                or (a.minimum_alt_copies, a.maximum_alt_copies) != (1, 1)
                or (b.minimum_alt_copies, b.maximum_alt_copies) != (1, 1)
                or left in unprojected_repeats
                or right in unprojected_repeats):
            continue
        if (a.input_sha256 is not None and a.input_sha256 == b.input_sha256
                and a.sample.sample_name == b.sample.sample_name
                and a.source_record_id == b.source_record_id and a.allele_number != b.allele_number):
            bases.append(PairPhaseBasis(source_ids, HypothesisPhase.IN_TRANS,
                                        "same-record-distinct-copies", None, None))
        for x, y in product(a.orientations, b.orientations):
            if (a.input_sha256 is not None and a.input_sha256 == b.input_sha256
                    and a.sample.sample_name == b.sample.sample_name
                    and (x.namespace, x.block, x.contract_sha256) == (y.namespace, y.block, y.contract_sha256)):

                bases.append(PairPhaseBasis(source_ids,
                    HypothesisPhase.IN_CIS if x.haplotype == y.haplotype else HypothesisPhase.IN_TRANS,
                    "shared-phase-block", x.namespace, x.block))
    bases = tuple(sorted(set(bases), key=_canonical))
    relations = {b.relation for b in bases}
    if (any(props[key].conflicting_source_genotypes for key in (left, right))
            or HypothesisPhase.CONFLICTING in relations or len(relations) > 1):
        phase = HypothesisPhase.CONFLICTING
    else:
        phase = next(iter(relations), HypothesisPhase.UNKNOWN)
    limits = ("parental-origin-not-established",)
    if phase is HypothesisPhase.UNKNOWN:
        limits += ("no-consistent-shared-source-phase-basis",)
    return PairPhaseProperties((left, right), phase, bases, limits)


def evaluate_source_genotype(*, original, record, sample, derived, source_status, release, read_required):
    """Interpret one source ALT's actual sample facts through the shared kernel.

    Closure of derived member identities and technical obligations belongs to
    the supplying inventory owner. No other source or context member is needed
    to interpret GT/PGT/PS, copies, ploidy and source-local phase orientation.
    """
    gt = _tokens(sample.genotype)
    ploidy = None if sample.genotype in (None, ".") else len(gt)
    lower = gt.count(original.allele_number) if ploidy is not None else None
    upper = lower + gt.count(None) if lower is not None else None
    present = None if lower is None or lower == 0 and upper > 0 else lower > 0
    if source_status == "non_variant":
        lower = upper = None
        present = None
    orientations, limits = [], list(sample.phase_reasons)
    conflict = sample.phase_state.value == "conflicting"
    if sample.pgt is not None:
        if not sample.pgt_pid_declared or sample.pid is None:
            limits.append("PGT/PID-declaration-or-block-unavailable")
        elif GENOTYPE.fullmatch(sample.pgt) is None:
            limits.append("PGT-unsupported-syntax")
        else:
            pgt = _tokens(sample.pgt)
            if not gt or None in gt or None in pgt:
                limits.append("PGT-cannot-complete-missing-GT")
            elif sorted(gt) != sorted(pgt):
                conflict = True
                limits.append("PGT-GT-copy-conflict")
            elif sample.pid_contract_sha256 is None:
                limits.append("PGT/PID-producer-contract-unavailable")
            elif len(gt) == 2 and lower == 1 and "|" in sample.pgt and "/" not in sample.pgt:
                orientations.append(PhaseOrientation("PGT/PID", sample.pid, pgt.index(original.allele_number), sample.pid_contract_sha256))
    if ploidy == 2 and None not in gt and lower == 1 and "|" in sample.genotype:
        namespace = "GT/implicit-PS" if sample.phase_set is None else "GT/PS"
        block = None if sample.phase_set is None else str(int(sample.phase_set))
        orientations.append(PhaseOrientation(namespace, block, gt.index(original.allele_number)))
    if ploidy != 2 or None in gt or lower != 1:
        limits.append("pair-orientation-requires-complete-diploid-single-ALT-copy")
    if release == "VCFv4.1":  # Closed LongTR profile replay already passed.
        orientations = []
        limits.append("longtr-phase-local-to-source-record; no-cross-record-phase-set")
    return SourceGenotypeProperties(
        original.source_record_id, original.source_allele_id, original.allele_number,
        record.contig, sample, derived, source_status, ploidy, present, lower, upper,
        tuple(orientations), conflict, tuple(dict.fromkeys(limits)), read_required or bool(sample.phase_reasons), original.input_sha256,
    )


def evaluate_allele_genotype(*, allele_id, native_identity_id, own_calls, calls, locus_call, screened,
                            maximum_phase_bases=None):
    """Aggregate exact-event source facts; never sum conflicting source copies.

    Both own calls and the complete exact-event contributor set are supplied by
    the globally closed owner. This is the existing allele-property decision,
    separated from rich inventory materialization, not a new interpretation.
    """
    # A parent references every original ALT slot, whose different copy counts
    # are not competing measurements of one individual allele. Keep the raw
    # slot properties and independent ploidy/phase conflicts, but do not
    # intersect sibling ALT counts into a fabricated genotype contradiction.
    known = () if locus_call else tuple(s for s in calls if s.minimum_alt_copies is not None)
    lower = max((s.minimum_alt_copies for s in known), default=None)
    upper = min((s.maximum_alt_copies for s in known), default=None)
    ploidies = tuple(sorted({s.reported_ploidy for s in calls if s.reported_ploidy is not None}))
    conflict = len(ploidies) > 1 or lower is not None and lower > upper
    present = None if lower is None or lower == 0 and upper > 0 else lower > 0
    limits = []
    if conflict:
        lower = upper = present = None
        limits.append("conflicting-source-genotypes; copies-not-summed")
    phase_conflicts = [] if maximum_phase_bases is None else _BoundedPhaseBases(maximum_phase_bases)
    orientations_by_block = defaultdict(lambda: defaultdict(list))
    for call in calls:
        if call.phase_conflict:
            phase_conflicts.append(PairPhaseBasis((call.source_allele_id, call.source_allele_id),
                HypothesisPhase.CONFLICTING, "source-genotype-conflict", None, None))
        for orientation in call.orientations:
            orientations_by_block[call.input_sha256, call.sample.sample_name, call.contig,
                                  orientation.namespace, orientation.block, orientation.contract_sha256][orientation.haplotype].append(call)
    for (_, _, _, namespace, block, _), orientations in orientations_by_block.items():
        for a, b in product(orientations[0], orientations[1]):
            phase_conflicts.append(PairPhaseBasis(tuple(sorted((a.source_allele_id, b.source_allele_id))),
                HypothesisPhase.CONFLICTING, "source-genotype-conflict", namespace, block))
    if locus_call:
        lower = upper = present = None
        limits.append("source-ALT-not-bound-to-native-repeat-sizes")
    if {s.source_allele_id for s in calls} != {s.source_allele_id for s in own_calls}:
        limits.append("exact-event-source-views-assessed-together; copies-not-summed")
    if any(s.reported_ploidy is None or None in _tokens(s.sample.genotype) for s in calls):
        limits.append("missing-source-GT-properties-retained")
    if screened:
        limits.append("screen-candidate-not-confirmed-carriage; source-GT-unchanged")
    return AlleleGenotypeProperties(
        allele_id, native_identity_id, tuple(s.source_allele_id for s in own_calls),
        present, lower, upper, ploidies, conflict, tuple(sorted(set(phase_conflicts), key=_canonical)), screened,
        screened or any(s.read_validation_required for s in calls), tuple(limits),
        tuple(s.source_allele_id for s in calls),
    )


def _evaluate(
    inventory: AdmissionInventory,
    sample_name: str,
    pairs,
    *,
    inventory_sha256: str | None = None,
    ledger_source_samples: tuple[tuple[str, str], ...] = (),
    call_source_samples: Mapping[str, str] | None = None,
) -> SuppliedGenotypeProperties:
    records = {r.qc.source_record_id: r.qc for r in inventory.records}
    samples = ({inventory.input_sha256: sample_name} if inventory.source_set is None else
               {m.input_sha256: m.sample_name for m in inventory.source_set.members})
    releases = {r.source.source_sha256: r.source.source_release for r in inventory.records}
    native = {a.allele_id: next(i.identity for i in a.identities if i.identity_id == a.native_identity_id)
               for a in inventory.alleles}
    projected = {a.allele_id: longtr.is_projection(a) for a in inventory.alleles}
    parents = {a.allele_id for a in inventory.alleles if longtr.is_parent(a)}
    events = {a.allele_id: _allele_event(a) for a in inventory.alleles}
    selected_samples = samples if call_source_samples is None else call_source_samples
    selected_calls = {
        a.allele_id: tuple(c for c in a.calls or () if c.sample is not None
                          and c.sample.sample_name == selected_samples[c.source.source_sha256])
        for a in inventory.alleles
    }
    obligations = {
        (a.allele_id, c.source.source_fields["source_allele_id"]): c.source.source_fields["read_validation_required"]
        for a in inventory.alleles for c in selected_calls[a.allele_id]
    }
    sources = []
    for entry in inventory.source_alleles:
        original = entry.source
        record = records[original.source_record_id]
        selected_sample = (original.sample_name if hasattr(original, "sample_name")
                           else samples[original.input_sha256])
        sample = next(
            s for s in record.sample_calls if s.sample_name == selected_sample
        )
        # Native projections are sample-specific; normalized sequence identities
        # may instead aggregate source contributors from several samples.
        derived = tuple(key for key in entry.derived_allele_ids if selected_calls[key])
        read_required = any(obligations[key, original.source_allele_id] for key in derived
                            ) if derived else bool(technical_flags_for(
                                record, sample, original.allele_number,
                            ))
        sources.append(evaluate_source_genotype(original=original,record=record,sample=sample,
            derived=derived,source_status=entry.status,release=releases[original.input_sha256],
            read_required=read_required))
    sources = tuple(sorted(sources, key=lambda s: s.source_allele_id))
    by_allele = defaultdict(list)
    for source in sources:
        for key in source.derived_allele_ids:
            by_allele[key].append(source)
    # Only exact sequence events establish equivalence for copy aggregation. The
    # conservative structural/repeat anchor guard alone is not equality proof.
    by_event = defaultdict(dict)
    for key, calls in by_allele.items():
        if events[key][0] == "sequence":
            by_event[events[key]].update((s.source_allele_id, s) for s in calls)
    contributors = {key: tuple(sorted(by_event[events[key]].values(), key=lambda s: s.source_allele_id))
                    if events[key][0] == "sequence" else tuple(calls) for key, calls in by_allele.items()}
    properties = []
    for allele in sorted(inventory.alleles, key=lambda a: a.allele_id):
        own_calls = by_allele[allele.allele_id]
        if not own_calls:
            continue
        calls = contributors[allele.allele_id]
        locus_call = isinstance(native[allele.allele_id], RepeatAllele) and not projected[allele.allele_id]
        properties.append(evaluate_allele_genotype(allele_id=allele.allele_id,
            native_identity_id=allele.native_identity_id,own_calls=own_calls,calls=calls,
            locus_call=locus_call,screened=isinstance(native[allele.allele_id], MosaicAllele)))
    props = {p.allele_id: p for p in properties}
    unprojected_repeats = {key for key in props
        if isinstance(native[key], RepeatAllele) and not projected[key]}
    pair_results = []
    for left, right in sorted(pairs):
        if left not in props or right not in props or left in parents or right in parents or events[left] == events[right]:
            raise ValueError
        pair_results.append(assess_pair_phase(left, right,
            left_properties=props[left], right_properties=props[right],
            left_sources=contributors[left], right_sources=contributors[right],
            unprojected_repeats=unprojected_repeats))
    return SuppliedGenotypeProperties(
        inventory_sha256 or hashlib.sha256(
            admission_inventory_as_json(inventory).encode("ascii")
        ).hexdigest(), inventory.input_sha256,
        inventory.reference_manifest_sha256, sample_name, POLICY_SHA256, sources,
        tuple(properties), tuple(pair_results), inventory.source_set,
        ledger_source_samples,
    )


@dataclass(frozen=True, slots=True)
class _LedgerQcView:
    source_record_id: str
    contig: str
    sample_calls: tuple[SampleCallQc, ...]


@dataclass(frozen=True, slots=True)
class _LedgerRecordView:
    qc: _LedgerQcView
    source: object


@dataclass(frozen=True, slots=True)
class _LedgerSourceView:
    source_record_id: str
    source_allele_id: str
    allele_number: int
    input_sha256: str
    sample_name: str


@dataclass(frozen=True, slots=True)
class _LedgerSourceEntryView:
    source: _LedgerSourceView
    derived_allele_ids: tuple[str, ...]
    status: Literal["derived"] = "derived"


@dataclass(frozen=True, slots=True)
class _LedgerInventoryView:
    records: tuple[_LedgerRecordView, ...]
    source_alleles: tuple[_LedgerSourceEntryView, ...]
    alleles: tuple[object, ...]
    input_sha256: str
    reference_manifest_sha256: str
    source_set: None = None


def _identity_contig(allele, source_allele_id: str) -> str:
    identity = next(
        item.identity for item in allele.identities
        if item.identity_id == source_allele_id
    )
    if isinstance(identity, (MitochondrialAllele, MosaicAllele)):
        identity = identity.sequence
    if isinstance(identity, (StructuralAllele, RepeatAllele)):
        return identity.interval.contig
    return identity.contig


def _canonical_pair_keys(pairs):
    # Validate ordinary callers too, but retain the original canonical bounded
    # sequence passed by complete-context enumeration rather than a second list.
    _typed(pairs, tuple[tuple[str, str], ...])
    canonical_pairs = pairs if all(a < b for a, b in pairs) else tuple(tuple(sorted(p)) for p in pairs)
    if len(set(canonical_pairs)) != len(canonical_pairs) or any(a == b for a, b in canonical_pairs):
        raise ValueError
    return canonical_pairs


def ledger_call_binding(allele, call, *, sample_name, source_manifest_sha256):
    """Validate one unchanged ledger call's original source/sample relation."""
    fields = call.source.source_fields or {}
    source_id = fields.get("source_allele_id")
    record_id = fields.get(
        "source_record_id", call.source.source_record_id
    )
    allele_number = fields.get("allele_number")
    sample_binding = fields.get("merged_intake_sample")
    selected_sample = sample_name
    if sample_binding is not None:
        if (set(sample_binding) != {"manifest_sha256", "sample_id", "source_sample_name"}
                or sample_binding["manifest_sha256"] != source_manifest_sha256
                or sample_binding["sample_id"] != sample_name
                or not isinstance(sample_binding["source_sample_name"], str)
                or not sample_binding["source_sample_name"]):
            raise ValueError
        selected_sample = sample_binding["source_sample_name"]
    if (
        call.sample is None
        or call.sample.sample_name != selected_sample
        or type(source_id) is not str
        or type(record_id) is not str
        or not isinstance(allele_number, int)
        or isinstance(allele_number, bool)
        or allele_number < 1
        or source_id not in {
            item.identity_id for item in allele.identities
        }
    ):
        raise ValueError
    source_input = fields.get("source_input_sha256")
    if type(source_input) is not str or len(source_input) != 64:
        raise ValueError
    return source_id, record_id, allele_number, source_input


def assess_ledger_genotypes(
    members,
    *,
    sample_name: str,
    source_manifest_sha256: str,
    reference_manifest_sha256: str,
    pairs: tuple[tuple[str, str], ...] = (),
    policy_bytes: bytes = POLICY_BYTES,
    policy_sha256: str = POLICY_SHA256,
) -> SuppliedGenotypeProperties:
    """Evaluate the existing genotype/phase kernel over closed ledger members.

    The caller must first authenticate and globally close the context membership.
    Each retained call carries its original sample, source-record/ALT identity,
    genotype, phase fields and technical obligation; this adapter only rebuilds
    the small context-local views consumed by ``_evaluate``.
    """

    try:
        from genosapiens.wgs_admission import LedgerAdmissionInput

        if (
            type(members) is not tuple
            or not members
            or any(type(value) is not LedgerAdmissionInput for value in members)
            or any(value.allele is None for value in members)
            or type(sample_name) is not str
            or not sample_name
            or type(source_manifest_sha256) is not str
            or type(reference_manifest_sha256) is not str
            or policy_bytes != POLICY_BYTES
            or policy_sha256 != hashlib.sha256(policy_bytes).hexdigest()
        ):
            raise ValueError
        alleles = tuple(value.allele for value in members)
        if len({allele.allele_id for allele in alleles}) != len(alleles):
            raise ValueError
        canonical_pairs = _canonical_pair_keys(pairs)

        source_members = defaultdict(list)
        record_members = defaultdict(list)
        source_samples = {}
        for allele in alleles:
            for call in allele.calls or ():
                source_id, record_id, allele_number, source_input = ledger_call_binding(
                    allele,call,sample_name=sample_name,source_manifest_sha256=source_manifest_sha256)
                source_members[source_id].append((allele, call))
                record_members[record_id].append((allele, call))
                previous = source_samples.setdefault(
                    call.source.source_sha256, call.sample.sample_name
                )
                if previous != call.sample.sample_name:
                    raise ValueError
                if (type(source_input) is not str or len(source_input) != 64
                        or source_samples.setdefault(source_input, call.sample.sample_name) != call.sample.sample_name):
                    raise ValueError
        expected_sources = {
            source_id for value in members for source_id in value.source_allele_ids
        }
        if set(source_members) != expected_sources:
            raise ValueError

        records = []
        for record_id, values in sorted(record_members.items()):
            first_allele, first_call = values[0]
            if any(
                call.sample != first_call.sample
                or (call.source.source_fields or {}).get(
                    "source_record_id", call.source.source_record_id
                ) != record_id
                for _allele, call in values
            ):
                raise ValueError
            source = first_call.source
            release = source.source_release
            records.append(_LedgerRecordView(
                _LedgerQcView(
                    record_id,
                    _identity_contig(first_allele, (
                        source.source_fields or {}
                    )["source_allele_id"]),
                    (first_call.sample,),
                ),
                type(source)(
                    source.source_dataset, release, source.source_record_id,
                    source.source_fields["source_input_sha256"], source.public_manifest_sha256,
                    source.source_fields,
                ),
            ))

        entries = []
        for source_id, values in sorted(source_members.items()):
            first_allele, first_call = values[0]
            fields = first_call.source.source_fields or {}
            if any(call.sample != first_call.sample for _allele, call in values):
                raise ValueError
            entries.append(_LedgerSourceEntryView(
                _LedgerSourceView(
                    fields.get(
                        "source_record_id", first_call.source.source_record_id
                    ), source_id,
                    fields["allele_number"], fields["source_input_sha256"],
                    first_call.sample.sample_name,
                ),
                tuple(sorted({allele.allele_id for allele, _call in values})),
            ))
        inventory = _LedgerInventoryView(
            tuple(records), tuple(entries), alleles, source_manifest_sha256,
            reference_manifest_sha256,
        )
        digest = hashlib.sha256(_canonical({
            "schema": "genosapiens.ledger-genotype-inventory",
            "source_manifest_sha256": source_manifest_sha256,
            "reference_manifest_sha256": reference_manifest_sha256,
            "members": tuple(value.as_bytes().decode("ascii") for value in members),
        }).encode("ascii")).hexdigest()
        return _evaluate(
            inventory, sample_name, canonical_pairs,
            inventory_sha256=digest,
            ledger_source_samples=tuple(sorted(source_samples.items())),
            call_source_samples=source_samples,
        )
    except Exception:
        raise GenotypePropertiesError(
            "ledger genotype property replay failed"
        ) from None


def assess_supplied_genotypes(
    inputs: SuppliedInputs, *, sample_name: str,
    pairs: tuple[tuple[str, str], ...] = (),
    policy_bytes: bytes = POLICY_BYTES, policy_sha256: str = POLICY_SHA256,
    expected_result: SuppliedGenotypeProperties | None = None,
) -> SuppliedGenotypeProperties:
    """Replay original VCF/QC/normalization, then evaluate all supplied ALT slots.

    ``pairs`` requests phase assessment only; it is not partner-search coverage or
    a candidate-generation gate. Existing QC/normalization require a nonempty VCF;
    a preparation failure is not a successful empty assessment. The fixed policy
    pin is an engineering contract, not empirical release acceptance.
    Source ALT copy counts, including mitochondrial GT, are not organism/locus
    dosage or evidence of required model count. No sex or parental origin is read.
    """
    try:
        if (type(policy_bytes) is not bytes or policy_bytes != POLICY_BYTES
                or policy_sha256 != hashlib.sha256(policy_bytes).hexdigest()
                or type(sample_name) is not str):
            raise ValueError
        inventory = build_admission_inventory(inputs)
        if (inventory.source_set is not None and sample_name != inventory.source_set.sample_id
                or inventory.source_set is None and any(sample_name not in {s.sample_name for s in r.qc.sample_calls}
                                                        for r in inventory.records)):
            raise ValueError
        canonical_pairs = _canonical_pair_keys(pairs)
        result = _evaluate(inventory, sample_name, canonical_pairs)
        if expected_result is not None and result != _typed(expected_result, SuppliedGenotypeProperties):
            raise ValueError
        return result
    except Exception:
        raise GenotypePropertiesError("supplied genotype property replay failed") from None


def genotype_properties_as_json(result: SuppliedGenotypeProperties) -> bytes:
    try:
        return (_canonical({"schema": SCHEMA, "classification": "controlled", "llm_policy": "no-llm",
                            "result": _typed(result, SuppliedGenotypeProperties)}) + "\n").encode("ascii")
    except Exception:
        raise GenotypePropertiesError("cannot encode supplied genotype properties") from None


def genotype_properties_from_json(content: bytes) -> SuppliedGenotypeProperties:
    """Structural decoding only; use assess_supplied_genotypes for source replay."""
    try:
        data = _loads(content)
        if set(data) != {"schema", "classification", "llm_policy", "result"} or (
            data["schema"], data["classification"], data["llm_policy"]
        ) != (SCHEMA, "controlled", "no-llm"):
            raise ValueError
        return _typed(data["result"], SuppliedGenotypeProperties, decode=True)
    except Exception:
        raise GenotypePropertiesError("invalid supplied genotype property JSON") from None
