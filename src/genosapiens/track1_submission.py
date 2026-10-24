"""Source-replayed exact-identity selection for the Track 1 submission adapter.

This module selects hypotheses only. It does not assign EPCR, write the official
CSV, run preflight, freeze a policy, or claim that public conformance matches the
private truth set. Every candidate and score is rebuilt from independently pinned
inputs before representation and eligibility are assessed.
"""

from __future__ import annotations

from collections.abc import Mapping
from contextlib import ExitStack, contextmanager
from dataclasses import dataclass
import hashlib
import re
from typing import Callable, Iterator, Literal

from genosapiens import admission_links
from genosapiens.admission_links import PinnedPublicSource
from genosapiens.allele_identity import (
    AlleleKind,
    MitochondrialAllele,
    MosaicAllele,
    ReferenceRole,
    SequenceAllele,
)
from genosapiens.evidence_fusion import FusionRanking, rank_hypotheses
from genosapiens.evidence_model import (
    AlleleEvidence,
    Completeness,
    EnumerationStatus,
    EvidenceDocument,
    IdentityStage,
    _Controlled,
    _exact,
    _ids,
    _sha,
    _typed,
)
from genosapiens.evidence_serialization import _canonical
from genosapiens.genotype_hypotheses import (
    ConstructedHypothesis,
    HypothesisGeneration,
    HypothesisReplay,
    replay_genotype_hypotheses,
)
from genosapiens.genotype_properties import selected_source_calls
from genosapiens.hpo_ic import HpoICStore
from genosapiens.reference_concordance import (
    ReferenceIdentityContract,
    ReferenceRoleBinding,
    RescueReferenceAuthorization,
    _OpenIndexedFasta,
    _parse_fai_bytes,
    verify_reference_construction,
    verify_reference_manifest_contract,
)
from genosapiens.submission_rendering import (
    CANONICAL_CONTIGS,
    SubmissionRenderingError,
    load_submission_contig_renderer,
)
from genosapiens.variant_normalization import NormalizedAllele, _variant_kind


SCHEMA = "genosapiens.track1-submission-selection"
PUBLIC_REFERENCE_IDENTITY = "grch38-no-alt-hs38d1-public-native-v1"
_DNA = re.compile(r"[ACGT]+")
_ERRORS = (
    ValueError,
    TypeError,
    RuntimeError,
    AttributeError,
    KeyError,
    IndexError,
    OverflowError,
    RecursionError,
    StopIteration,
    OSError,
)


class Track1SubmissionError(ValueError):
    """Static diagnostics only; candidate identities and paths remain controlled."""


@dataclass(frozen=True, slots=True, repr=False)
class SubmissionRankingReplay:
    """Independent source-bound ranker inputs; a ranking result is never trusted alone."""

    document: EvidenceDocument
    ic_store: HpoICStore
    profile_evidence_id: str | None
    policy_bytes: bytes
    policy_sha256: str
    derivation_bytes: bytes
    expected_derivation_sha256: str
    expected_result: FusionRanking | None = None

    def __repr__(self):
        return "<controlled evidence; no-llm>"


@dataclass(frozen=True, slots=True, repr=False)
class SubmissionAlleleTuple(_Controlled):
    """One literal evaluator allele tuple in the configured canonical name view."""

    contig: str
    position: int
    reference: str
    alternate: str

    def __post_init__(self):
        _Controlled.__post_init__(self)
        canonical = self.contig.removeprefix("chr")
        if (
            self.contig != "chr" + canonical
            or canonical not in CANONICAL_CONTIGS
            or self.position < 1
            or _DNA.fullmatch(self.reference) is None
            or _DNA.fullmatch(self.alternate) is None
            or self.reference == self.alternate
        ):
            raise Track1SubmissionError("invalid exact submission allele tuple")


@dataclass(frozen=True, slots=True, repr=False)
class AuthorizedReadEvidenceTrace(_Controlled):
    """Validated upstream proof bindings, never a read-execution interface."""

    reference_authorization: RescueReferenceAuthorization
    parent_prepared_sha256: str
    parent_early_snapshot_sha256: str
    added_read_inventory_sha256: str
    read_identity_report_sha256: str
    producer_stack_sha256: str
    selected_plan_sha256: str
    stage_output_sha256: str
    read_delta_sha256: str
    merged_intake_sha256: str

    def __post_init__(self):
        _Controlled.__post_init__(self)
        authorization = self.reference_authorization
        _ids((authorization.reference_identity_id,))
        for value in (
            authorization.reference_manifest_sha256,
            authorization.input_sha256,
            authorization.validated_input_sha256,
            authorization.concordance_report_sha256,
            self.parent_prepared_sha256,
            self.parent_early_snapshot_sha256,
            self.added_read_inventory_sha256,
            self.read_identity_report_sha256,
            self.producer_stack_sha256,
            self.selected_plan_sha256,
            self.stage_output_sha256,
            self.read_delta_sha256,
            self.merged_intake_sha256,
        ):
            _sha(value)
        if (
            authorization.authorized is not True
            or authorization.total_records <= 0
            or authorization.matching_records != authorization.total_records
        ):
            raise Track1SubmissionError("invalid authorized read-evidence trace")


@dataclass(frozen=True, slots=True, repr=False)
class SubmissionSourceTrace(_Controlled):
    """One selected call linked through its source identity and authority."""

    source_allele_id: str
    source_identity_record_id: str
    source_identity_sha256: str
    source_record_id: str
    source_input_sha256: str
    sample_name: str
    authority: Literal["supplied_source", "authorized_read_evidence"] = "supplied_source"
    authorization: AuthorizedReadEvidenceTrace | None = None

    def __post_init__(self):
        _Controlled.__post_init__(self)
        for value in (
            self.source_allele_id,
            self.source_identity_record_id,
            self.source_record_id,
            self.sample_name,
        ):
            _ids((value,))
        _sha(self.source_identity_sha256)
        _sha(self.source_input_sha256)
        if (
            self.authority == "supplied_source"
            and self.authorization is not None
            or self.authority == "authorized_read_evidence"
            and (
                self.authorization is None
                or self.source_input_sha256 != self.authorization.stage_output_sha256
            )
        ):
            raise Track1SubmissionError("invalid submission source authority")


@dataclass(frozen=True, slots=True, repr=False)
class SubmissionAlleleTrace(_Controlled):
    """Source-to-normalized-to-rendered mapping plus the additional public REF binding."""

    allele_id: str
    native_identity_id: str
    native_kind: AlleleKind
    normalized_identity_id: str
    rendered: SubmissionAlleleTuple
    source_reference_role: ReferenceRole
    source_reference_identity_id: str
    source_reference_manifest_sha256: str
    source_traces: tuple[SubmissionSourceTrace, ...]
    submission_reference_role: Literal["public_annotation_benchmark"]
    submission_reference_identity_id: str
    submission_reference_manifest_sha256: str
    contig_aliases_sha256: str

    def __post_init__(self):
        _Controlled.__post_init__(self)
        for value in (
            self.allele_id,
            self.native_identity_id,
            self.normalized_identity_id,
            self.source_reference_identity_id,
            self.submission_reference_identity_id,
        ):
            _ids((value,))
        _sha(self.source_reference_manifest_sha256)
        _sha(self.submission_reference_manifest_sha256)
        _sha(self.contig_aliases_sha256)
        if (
            self.native_kind not in (AlleleKind.SEQUENCE, AlleleKind.MITOCHONDRIAL, AlleleKind.MOSAIC)
            or not self.source_traces
            or self.source_traces != tuple(sorted(self.source_traces, key=lambda item: item.source_allele_id))
            or any(
                item.authority == "authorized_read_evidence"
                and item.authorization.reference_authorization.reference_manifest_sha256
                != self.source_reference_manifest_sha256
                for item in self.source_traces
            )
        ):
            raise Track1SubmissionError("invalid exact submission allele trace")


@dataclass(frozen=True, slots=True, repr=False)
class SubmissionRepresentationAssessment(_Controlled):
    allele_id: str
    native_identity_id: str
    status: Literal["eligible", "audit_only"]
    reasons: tuple[str, ...]
    trace: SubmissionAlleleTrace | None

    def __post_init__(self):
        _Controlled.__post_init__(self)
        _ids((self.allele_id, self.native_identity_id))
        _ids(self.reasons)
        if (
            (self.status == "eligible") != (not self.reasons and self.trace is not None)
            or self.trace is not None
            and (self.trace.allele_id, self.trace.native_identity_id)
            != (self.allele_id, self.native_identity_id)
        ):
            raise Track1SubmissionError("invalid submission representation assessment")


@dataclass(frozen=True, slots=True, repr=False)
class SubmissionGroup(_Controlled):
    """One exact singleton or unordered pair represented by its best eligible interpretation."""

    primary_hypothesis_id: str
    source_rank: int
    source_score: str
    alleles: tuple[SubmissionAlleleTrace, ...]
    alternative_hypothesis_ids: tuple[str, ...]
    partial_member: bool = False

    @property
    def identity(self) -> tuple[SubmissionAlleleTuple, ...]:
        return tuple(item.rendered for item in self.alleles)

    def __post_init__(self):
        _Controlled.__post_init__(self)
        _ids((self.primary_hypothesis_id,))
        _ids(self.alternative_hypothesis_ids)
        _exact(self.source_score)
        if (
            self.source_rank < 1
            or len(self.alleles) not in (1, 2)
            or self.alleles != tuple(sorted(self.alleles, key=_allele_trace_key))
            or len(set(self.identity)) != len(self.identity)
            or self.primary_hypothesis_id in self.alternative_hypothesis_ids
            or type(self.partial_member) is not bool
            or self.partial_member and len(self.alleles) != 1
        ):
            raise Track1SubmissionError("invalid submission identity group")


@dataclass(frozen=True, slots=True, repr=False)
class SubmissionExclusion(_Controlled):
    hypothesis_id: str
    source_rank: int
    source_score: str
    reasons: tuple[str, ...]

    def __post_init__(self):
        _Controlled.__post_init__(self)
        _ids((self.hypothesis_id,))
        _ids(self.reasons)
        _exact(self.source_score)
        if self.source_rank < 1 or not self.reasons:
            raise Track1SubmissionError("invalid submission exclusion")


@dataclass(frozen=True, slots=True, repr=False)
class SubmissionReferenceBinding(_Controlled):
    contract_id: str
    contract_sha256: str
    manifest_sha256: str
    role: Literal["public_annotation_benchmark"]
    identity_id: str
    fasta_sha256: str
    fai_sha256: str
    contig_aliases_sha256: str

    def __post_init__(self):
        _Controlled.__post_init__(self)
        _ids((self.contract_id, self.identity_id))
        for value in (
            self.contract_sha256,
            self.manifest_sha256,
            self.fasta_sha256,
            self.fai_sha256,
            self.contig_aliases_sha256,
        ):
            _sha(value)


@dataclass(frozen=True, slots=True, repr=False)
class SubmissionSelection(_Controlled):
    generation_sha256: str
    ranking_sha256: str
    ranking_policy_sha256: str
    ranking_order: tuple[str, ...]
    reference: SubmissionReferenceBinding
    representation_assessments: tuple[SubmissionRepresentationAssessment, ...]
    selected: tuple[SubmissionGroup, ...]
    eligible_tail: tuple[SubmissionGroup, ...]
    excluded: tuple[SubmissionExclusion, ...]
    schema: Literal["genosapiens.track1-submission-selection"] = SCHEMA
    analysis_mode: Literal["primary"] = "primary"
    maximum_groups: Literal[10] = 10
    epcr_assigned: Literal[False] = False
    csv_emitted: Literal[False] = False
    controlled_use_ready: Literal[False] = False

    def __post_init__(self):
        _Controlled.__post_init__(self)
        for value in (self.generation_sha256, self.ranking_sha256, self.ranking_policy_sha256):
            _sha(value)
        _ids(self.ranking_order)
        if (
            len(self.selected) > self.maximum_groups
            or len(self.selected) != min(self.maximum_groups, len(self.selected) + len(self.eligible_tail))
            or self.representation_assessments
            != tuple(sorted(self.representation_assessments, key=lambda item: item.allele_id))
            or len({item.allele_id for item in self.representation_assessments})
            != len(self.representation_assessments)
            or any(group.source_rank > following.source_rank for group, following in zip(
                (*self.selected, *self.eligible_tail),
                (*self.selected, *self.eligible_tail)[1:],
            ))
        ):
            raise Track1SubmissionError("invalid deterministic submission selection")
        groups = (*self.selected, *self.eligible_tail)
        identities = tuple(group.identity for group in groups)
        if len(set(identities)) != len(identities):
            raise Track1SubmissionError("duplicate exact submission identity group")
        ranked = {identifier: rank for rank, identifier in enumerate(self.ranking_order, 1)}
        accounted = (
            *(group.primary_hypothesis_id for group in groups),
            *(identifier for group in groups for identifier in group.alternative_hypothesis_ids),
            *(item.hypothesis_id for item in self.excluded),
        )
        if (
            len(accounted) != len(self.ranking_order)
            or set(accounted) != set(self.ranking_order)
            or any(ranked[group.primary_hypothesis_id] != group.source_rank for group in groups)
            or any(ranked[item.hypothesis_id] != item.source_rank for item in self.excluded)
            or any(
                tuple(sorted(group.alternative_hypothesis_ids, key=ranked.__getitem__))
                != group.alternative_hypothesis_ids
                or any(ranked[identifier] <= group.source_rank for identifier in group.alternative_hypothesis_ids)
                for group in groups
            )
        ):
            raise Track1SubmissionError("submission selection does not account for the ranking")

    def as_bytes(self) -> bytes:
        return (_canonical({
            "schema": self.schema,
            "classification": "controlled",
            "llm_policy": "no-llm",
            "epcr_assigned": False,
            "csv_emitted": False,
            "controlled_use_ready": False,
            "result": self,
        }) + "\n").encode("utf-8")


def _submission_audit_rows(
    selection: SubmissionSelection,
    incidental_hypothesis_ids: tuple[str, ...],
) -> tuple[dict[str, object], ...]:
    selection = _typed(selection, SubmissionSelection)
    _ids(incidental_hypothesis_ids)
    incidental = set(incidental_hypothesis_ids)
    exclusions = {item.hypothesis_id: item for item in selection.excluded}
    if (
        len(incidental) != len(incidental_hypothesis_ids)
        or not incidental <= set(exclusions)
    ):
        raise Track1SubmissionError("incidental hypotheses must be unique exclusions")
    rows: list[dict[str, object]] = []
    for group in selection.selected:
        rows.append({
            "hypothesis_id": group.primary_hypothesis_id,
            "category": "submitted",
            "genotype_status": "incomplete_partial_member" if group.partial_member else "as_assessed",
            "source_rank": group.source_rank,
            "source_score": group.source_score,
            "reasons": (),
        })
        rows.extend({
            "hypothesis_id": identifier,
            "category": "secondary",
            "source_rank": selection.ranking_order.index(identifier) + 1,
            "source_score": None,
            "reasons": ("alternate-interpretation-of-submitted-identity",),
        } for identifier in group.alternative_hypothesis_ids)
    for group in selection.eligible_tail:
        rows.append({
            "hypothesis_id": group.primary_hypothesis_id,
            "category": "secondary",
            "source_rank": group.source_rank,
            "source_score": group.source_score,
            "reasons": ("outside-primary-row-limit",),
        })
        rows.extend({
            "hypothesis_id": identifier,
            "category": "secondary",
            "source_rank": selection.ranking_order.index(identifier) + 1,
            "source_score": None,
            "reasons": ("alternate-interpretation-of-secondary-identity",),
        } for identifier in group.alternative_hypothesis_ids)
    rows.extend({
        "hypothesis_id": item.hypothesis_id,
        "category": "incidental" if item.hypothesis_id in incidental else "excluded",
        "source_rank": item.source_rank,
        "source_score": item.source_score,
        "reasons": item.reasons,
    } for item in selection.excluded)
    rows.sort(key=lambda item: (item["source_rank"], item["hypothesis_id"]))
    if tuple(item["hypothesis_id"] for item in rows) != selection.ranking_order:
        raise Track1SubmissionError("submission audit does not account for the ranking")
    return tuple(rows)


def submission_audit_as_json(
    selection: SubmissionSelection,
    *,
    incidental_hypothesis_ids: tuple[str, ...] = (),
) -> str:
    """Return the complete controlled machine-readable selection audit."""
    rows = _submission_audit_rows(selection, incidental_hypothesis_ids)
    return _canonical({
        "schema": "genosapiens.track1-submission-audit",
        "classification": "controlled",
        "llm_policy": "no-llm",
        "selection_sha256": hashlib.sha256(selection.as_bytes()).hexdigest(),
        "categories": ("submitted", "excluded", "secondary", "incidental"),
        "records": rows,
        "selection": selection,
    }) + "\n"


def submission_audit_as_markdown(
    selection: SubmissionSelection,
    *,
    incidental_hypothesis_ids: tuple[str, ...] = (),
) -> str:
    """Return a human-readable controlled/no-LLM audit of every ranked hypothesis."""
    rows = _submission_audit_rows(selection, incidental_hypothesis_ids)
    lines = [
        "# Track 1 submission audit",
        "",
        "Classification: `controlled`; LLM policy: `no-llm`.",
        "",
        f"Selection SHA-256: `{hashlib.sha256(selection.as_bytes()).hexdigest()}`",
        "",
        "| Rank | Hypothesis | Category | Source score | Reasons |",
        "| ---: | --- | --- | ---: | --- |",
    ]
    for item in rows:
        score = item["source_score"] if item["source_score"] is not None else "see selection"
        reasons = ", ".join(item["reasons"]) or "none"
        if item.get("genotype_status") == "incomplete_partial_member":
            reasons += "; incomplete recessive partial member; partner unresolved"
        lines.append(
            f"| {item['source_rank']} | `{item['hypothesis_id']}` | "
            f"{item['category']} | {score} | {reasons} |"
        )
    return "\n".join(lines) + "\n"


def verify_selected_submission_reference(
    selection: SubmissionSelection,
    *,
    reference: PinnedPublicSource,
    expected_public_identity_id: str = PUBLIC_REFERENCE_IDENTITY,
) -> None:
    """Reopen the pinned public-native reference and recheck every selected REF."""
    try:
        selection = _typed(selection, SubmissionSelection)
        _verify_selected_groups_reference(
            selection.reference,
            selection.selected,
            reference=reference,
            expected_public_identity_id=expected_public_identity_id,
        )
    except _ERRORS:
        raise Track1SubmissionError(
            "selected submission public-reference replay failed"
        ) from None


def _verify_selected_groups_reference(
    binding: SubmissionReferenceBinding,
    groups: tuple[SubmissionGroup, ...],
    *,
    reference: PinnedPublicSource,
    expected_public_identity_id: str,
) -> None:
    """Replay the small final group tuple using the shared representation kernel."""

    try:
        with _open_submission_reference(
            reference,
            expected_public_identity_id=expected_public_identity_id,
        ) as opened:
            if opened.binding != binding:
                raise ValueError
            for group in groups:
                for trace in group.alleles:
                    allele = trace.rendered
                    canonical = allele.contig.removeprefix("chr")
                    public_name = opened.public_names.get(canonical)
                    if (
                        public_name is None
                        or opened.fetch_public(
                            public_name, allele.position, len(allele.reference)
                        )
                        != allele.reference
                    ):
                        raise ValueError
    except _ERRORS:
        raise Track1SubmissionError(
            "selected submission public-reference replay failed"
        ) from None


def verify_sharded_selected_submission_reference(
    selection,
    *,
    reference: PinnedPublicSource,
    expected_public_identity_id: str = PUBLIC_REFERENCE_IDENTITY,
) -> None:
    """Recheck selected REF values from a bounded selection manifest."""

    try:
        from genosapiens.wgs_scientific_products import ScientificSelectionProduct

        if type(selection) is not ScientificSelectionProduct:
            raise ValueError
        _verify_selected_groups_reference(
            selection.reference,
            selection.selected,
            reference=reference,
            expected_public_identity_id=expected_public_identity_id,
        )
    except _ERRORS:
        raise Track1SubmissionError(
            "sharded selected submission public-reference replay failed"
        ) from None


@dataclass(frozen=True, slots=True)
class _OpenSubmissionReference:
    contract: ReferenceIdentityContract
    binding: SubmissionReferenceBinding
    public_names: Mapping[str, str]
    fetch_public: Callable[[str, int, int], str]


@dataclass(frozen=True, slots=True)
class _EligibleCandidate:
    constructed: ConstructedHypothesis
    source_rank: int
    source_score: str
    alleles: tuple[SubmissionAlleleTrace, ...]

    @property
    def identity(self):
        return tuple(item.rendered for item in self.alleles)


def _allele_tuple_key(value: SubmissionAlleleTuple):
    return CANONICAL_CONTIGS.index(value.contig.removeprefix("chr")), value.position, value.reference, value.alternate


def _allele_trace_key(value: SubmissionAlleleTrace):
    return _allele_tuple_key(value.rendered)


@contextmanager
def _open_submission_reference(
    pin: PinnedPublicSource,
    *,
    expected_public_identity_id: str,
) -> Iterator[_OpenSubmissionReference]:
    if type(pin) is not PinnedPublicSource:
        raise ValueError
    _ids((expected_public_identity_id,))
    with ExitStack() as stack:
        files = admission_links._Files(stack)
        data, manifest = files.source(pin)
        contract = ReferenceIdentityContract.from_mapping(data)
        verify_reference_manifest_contract(contract, manifest, contract_sha256=pin.contract_sha256)
        outputs = {item.path: item for item in manifest.outputs}
        report_output = outputs[contract.expected_output("construction_report").path]
        files.open(pin.root / report_output.path, report_output.sha256, report_output.size_bytes)
        verify_reference_construction(
            contract,
            manifest,
            pin.root,
            contract_sha256=pin.contract_sha256,
        )
        public_role = contract.reference_role("public_annotation_benchmark")
        if public_role.identity_id != expected_public_identity_id:
            raise ValueError
        renderer = load_submission_contig_renderer(
            identity_path=pin.contract_path,
            reference_manifest_path=pin.manifest_path,
            reference_root=pin.root,
            reference_binding=ReferenceRoleBinding(public_role.role, public_role.identity_id),
        )
        if renderer.reference_manifest_sha256 != pin.manifest_sha256:
            raise ValueError

        fasta = outputs[contract.expected_output(public_role.fasta_output_role).path]
        fai = outputs[contract.expected_output(public_role.fai_output_role).path]
        fasta_stream = files.open(pin.root / fasta.path, fasta.sha256, fasta.size_bytes)
        fai_stream = files.open(pin.root / fai.path, fai.sha256, fai.size_bytes)
        reference = _OpenIndexedFasta(fasta_stream, _parse_fai_bytes(fai_stream.read()))
        names = {}
        for source_name, (canonical, _length) in renderer._names.items():
            if canonical in names:
                raise ValueError
            names[canonical] = source_name
        yield _OpenSubmissionReference(
            contract,
            SubmissionReferenceBinding(
                contract.contract_id,
                pin.contract_sha256,
                pin.manifest_sha256,
                "public_annotation_benchmark",
                public_role.identity_id,
                fasta.sha256,
                fai.sha256,
                renderer.contig_aliases_sha256,
            ),
            names,
            reference.fetch,
        )


def _native_sequence(identity):
    if type(identity) is SequenceAllele:
        return identity
    if type(identity) in (MitochondrialAllele, MosaicAllele):
        return identity.sequence
    return None


def _source_ancestors(allele: AlleleEvidence, normalized_id: str):
    identities = {item.identity_id: item for item in allele.identities}
    pending = list(identities[normalized_id].derived_from_identity_ids)
    result = set()
    while pending:
        key = pending.pop()
        item = identities[key]
        if item.stage is IdentityStage.SOURCE:
            result.add(key)
        else:
            pending.extend(item.derived_from_identity_ids)
    return result


def _audit_only(allele, reason, trace=None):
    return SubmissionRepresentationAssessment(
        allele.allele_id,
        allele.native_identity_id,
        "audit_only",
        (reason,),
        trace,
    )


def _assess_representation_with_genotypes(
    allele: AlleleEvidence,
    genotypes,
    reference: _OpenSubmissionReference,
    renderer_cache: dict[tuple[str, str], object],
    reference_pin: PinnedPublicSource,
) -> SubmissionRepresentationAssessment:
    def source_facts():
        return (selected_source_calls(allele,genotypes),
                {item.source_allele_id:item for item in genotypes.source_properties},genotypes.input_sha256)
    return assess_representation_source_facts(allele,source_facts,reference,renderer_cache,reference_pin)


def assess_representation_source_facts(allele, source_facts, reference, renderer_cache, reference_pin):
    """Original rendering/trace kernel with a lazy bounded source-call owner.

    Source facts are demanded only after the same native identity guards as
    the expanded adapter. The callback returns selected calls, their indexed
    source properties and the input identity; it cannot certify eligibility.
    """
    native = next(item for item in allele.identities if item.identity_id == allele.native_identity_id)
    sequence = _native_sequence(native.identity)
    if sequence is None:
        return _audit_only(allele, "native-identity-not-faithfully-representable")
    normalized = tuple(item for item in allele.identities if item.stage is IdentityStage.NORMALIZED)
    if len(normalized) != 1 or type(normalized[0].identity) is not SequenceAllele:
        return _audit_only(allele, "normalized-sequence-identity-unavailable")
    normalized_mapping = normalized[0]
    if normalized_mapping.identity != sequence:
        return _audit_only(allele, "normalized-identity-differs-from-native-sequence")
    if (
        native.reference_manifest_sha256 != reference.binding.manifest_sha256
        or normalized_mapping.reference_manifest_sha256 != reference.binding.manifest_sha256
    ):
        return _audit_only(allele, "source-reference-manifest-not-submission-bound")

    identities = {item.identity_id: item for item in allele.identities}
    ancestors = _source_ancestors(allele, normalized_mapping.identity_id)
    try:
        calls, source_properties, input_sha256 = source_facts()
        traces = []
        for call in calls:
            fields = call.source.source_fields
            if not isinstance(fields, Mapping):
                raise ValueError
            source_allele_id = fields["source_allele_id"]
            source_mapping = identities[source_allele_id]
            source_property = source_properties[source_allele_id]
            expected_input = source_property.input_sha256 or input_sha256
            authorization = None
            if call.source.source_dataset == "supplied-vcf":
                trace_input_sha256 = expected_input
                source_binding_valid = call.source.source_sha256 == expected_input
            elif call.source.source_dataset == "genosapiens.wgs-source-ledger":
                trace_input_sha256 = fields.get("source_input_sha256")
                source_binding_valid = (
                    type(trace_input_sha256) is str
                    and len(trace_input_sha256) == 64
                    and all(value in "0123456789abcdef" for value in trace_input_sha256)
                    and fields.get("storage_envelope") == "bounded-ledger-v1"
                    and fields.get("source_record_id") == source_property.source_record_id
                    and fields.get("source_allele_id") == source_property.source_allele_id
                    and call.source.source_sha256 == source_mapping.source.source_sha256
                )
            elif fields.get("storage_envelope") == "bounded-read-delta-v1":
                proof = fields["authorized_read_evidence"]
                authorization = AuthorizedReadEvidenceTrace(
                    _typed(proof["reference_authorizations"][0]["authorization"], RescueReferenceAuthorization, decode=True),
                    *(proof[name] for name in ("parent_prepared_sha256", "parent_early_snapshot_sha256",
                        "added_read_inventory_sha256", "read_identity_report_sha256", "producer_stack_sha256",
                        "selected_plan_sha256", "stage_output_sha256", "read_delta_sha256", "merged_intake_sha256")))
                trace_input_sha256 = fields.get("source_input_sha256")
                source_binding_valid = (
                    trace_input_sha256 == authorization.stage_output_sha256
                    and call.source.source_sha256 == authorization.stage_output_sha256
                    and call.source.source_sha256 == source_mapping.source.source_sha256
                    and fields.get("source_record_id") == source_property.source_record_id)
            else:
                trace_input_sha256 = None
                source_binding_valid = False
            if (
                type(source_allele_id) is not str
                or source_allele_id not in ancestors
                or source_mapping.stage is not IdentityStage.SOURCE
                or source_mapping.reference_manifest_sha256 != reference.binding.manifest_sha256
                or not source_binding_valid
                or source_property.sample != call.sample
            ):
                raise ValueError
            traces.append(SubmissionSourceTrace(
                source_allele_id,
                source_mapping.source.source_record_id,
                source_mapping.source.source_sha256,
                call.source.source_record_id,
                trace_input_sha256,
                call.sample.sample_name,
                "authorized_read_evidence" if authorization is not None else "supplied_source",
                authorization,
            ))
        traces = tuple(sorted(set(traces), key=lambda item: item.source_allele_id))
        if not traces or {item.source_allele_id for item in traces} != ancestors:
            raise ValueError
    except (ValueError, TypeError, KeyError, AttributeError):
        return _audit_only(allele, "source-to-normalized-lineage-unresolved")

    key = (sequence.reference_role.value, sequence.reference_identity_id)
    try:
        renderer = renderer_cache.get(key)
        if renderer is None:
            renderer = load_submission_contig_renderer(
                identity_path=reference_pin.contract_path,
                reference_manifest_path=reference_pin.manifest_path,
                reference_root=reference_pin.root,
                reference_binding=ReferenceRoleBinding(*key),
            )
            renderer_cache[key] = renderer
        normalized_allele = NormalizedAllele(
            normalized_mapping.identity_id,
            sequence.reference_identity_id,
            sequence.reference_role.value,
            sequence.contig,
            sequence.position,
            sequence.reference,
            sequence.alternate,
            _variant_kind(sequence.reference, sequence.alternate),
        )
        view = renderer.render(normalized_allele, submission_contig_prefix="chr")
    except (SubmissionRenderingError, ValueError, TypeError, KeyError, AttributeError):
        return _audit_only(allele, "canonical-contig-or-normalized-representation-unavailable")

    canonical = view.contig.removeprefix("chr")
    public_name = reference.public_names.get(canonical)
    if public_name is None:
        return _audit_only(allele, "canonical-contig-not-present-in-public-reference")
    try:
        public_ref = reference.fetch_public(public_name, sequence.position, len(sequence.reference))
    except (ValueError, KeyError, IndexError, OSError):
        return _audit_only(allele, "public-native-ref-unavailable")
    if public_ref != sequence.reference:
        return _audit_only(allele, "public-native-ref-disagreement")
    if _DNA.fullmatch(sequence.reference) is None or _DNA.fullmatch(sequence.alternate) is None:
        return _audit_only(allele, "ambiguous-sequence-identity")
    rendered = SubmissionAlleleTuple(
        view.contig,
        normalized_allele.position,
        normalized_allele.reference,
        normalized_allele.alternate,
    )
    trace = SubmissionAlleleTrace(
        allele.allele_id,
        allele.native_identity_id,
        native.identity.kind,
        normalized_mapping.identity_id,
        rendered,
        sequence.reference_role,
        sequence.reference_identity_id,
        normalized_mapping.reference_manifest_sha256,
        traces,
        "public_annotation_benchmark",
        reference.binding.identity_id,
        reference.binding.manifest_sha256,
        reference.binding.contig_aliases_sha256,
    )
    return SubmissionRepresentationAssessment(
        allele.allele_id,
        allele.native_identity_id,
        "eligible",
        (),
        trace,
    )


def _assess_representation(
    allele: AlleleEvidence,
    generation: HypothesisGeneration,
    reference: _OpenSubmissionReference,
    renderer_cache: dict[tuple[str, str], object],
    reference_pin: PinnedPublicSource,
) -> SubmissionRepresentationAssessment:
    return _assess_representation_with_genotypes(
        allele, generation.genotypes, reference, renderer_cache, reference_pin,
    )


@dataclass(frozen=True, slots=True, repr=False)
class SubmissionEligibilityInputs:
    locus_call: bool
    allele_count: int
    enumeration_complete: bool
    source_inventory_complete: bool
    members_evaluable: bool
    model_state: Literal["supported", "unresolved", "contradicted"]
    nonqualifying_member: bool
    read_validation_required: bool
    population_assessed: bool
    population_not_excluded: bool
    incomplete: bool
    valid_partial: bool
    representation_reasons: tuple[tuple[str, ...], ...]


def requires_partial_member_checks(hypothesis):
    """Unknown required counts cannot waive an explicitly incomplete model's checks."""
    return hypothesis.completeness is Completeness.INCOMPLETE or hypothesis.context.model == "incomplete"


def submission_eligibility_inputs(item, search, population, representations, *, complete_source_statuses):
    hypothesis = item.hypothesis
    source_complete = (search is not None and search.context_id == item.admission_context_id
                       and search.source_inventory_status in complete_source_statuses)
    population_assessed = population is not None and population.hypothesis == hypothesis
    representation_reasons = tuple(submission_representation_reasons(ref,representations.get(ref.allele_id))
                                   for ref in hypothesis.allele_refs)
    return SubmissionEligibilityInputs(bool(hypothesis.locus_call_refs),len(hypothesis.allele_refs),
        hypothesis.enumeration_status is EnumerationStatus.COMPLETE,source_complete,
        source_complete and set(hypothesis.evidence_subject_ids).issubset(search.evaluable_allele_ids),
        item.model.state,any(member.state == "nonqualifying" for member in item.members),item.read_validation_required,
        population_assessed,population_assessed and population.outcome.value == "not_population_excluded",
        requires_partial_member_checks(hypothesis),
        len(hypothesis.allele_refs) == 1 and hypothesis.context.model == "incomplete" and item.incompleteness_penalty is not None,
        tuple(representation_reasons))


def submission_representation_reasons(allele_ref, assessment):
    if assessment is None or assessment.native_identity_id != allele_ref.identity_id:
        return ("allele-representation-not-assessed",)
    if assessment.status != "eligible":
        return tuple("allele-" + reason for reason in assessment.reasons)
    return ()


def evaluate_submission_eligibility(inputs):
    """Original ordered eligibility rules over complete evaluated dependencies."""
    inputs = _typed(inputs,SubmissionEligibilityInputs)
    if inputs.allele_count < 0 or len(inputs.representation_reasons) != inputs.allele_count:
        raise Track1SubmissionError("incomplete submission eligibility dependencies")
    reasons = []
    if inputs.locus_call:
        reasons.append("locus-call-is-not-an-exact-submission-allele")
    if inputs.allele_count not in (1,2):
        reasons.append("official-allele-cardinality-unsupported")
    if not inputs.enumeration_complete:
        reasons.append("required-enumeration-incomplete")
    if not inputs.source_inventory_complete:
        reasons.append("source-inventory-incomplete")
    elif not inputs.members_evaluable:
        reasons.append("selected-member-outside-declared-domain")
    if inputs.model_state != "supported":
        reasons.append("model-construction-" + inputs.model_state)
    if inputs.nonqualifying_member:
        reasons.append("member-qualification-contradicted")
    if inputs.read_validation_required:
        reasons.append("selected-allele-mandatory-read-validation-unmet")
    if not inputs.population_assessed:
        reasons.append("population-eligibility-not-assessed")
    elif not inputs.population_not_excluded:
        reasons.append("population-excluded")
    if inputs.incomplete and not inputs.valid_partial:
        reasons.append("incomplete-partial-model-invalid")
    for values in inputs.representation_reasons:
        reasons.extend(values)
    return tuple(dict.fromkeys(reasons))


def _eligibility_reasons_from_context(
    item: ConstructedHypothesis,
    search,
    population,
    representations: Mapping[str, SubmissionRepresentationAssessment],
    *,
    complete_source_statuses: tuple[str, ...] = ("complete",),
) -> tuple[str, ...]:
    return evaluate_submission_eligibility(submission_eligibility_inputs(
        item,search,population,representations,complete_source_statuses=complete_source_statuses))


def _eligibility_reasons(
    item: ConstructedHypothesis,
    generation: HypothesisGeneration,
    representations: Mapping[str, SubmissionRepresentationAssessment],
) -> tuple[str, ...]:
    searches = tuple(
        search for search in generation.searches
        if search.context_id == item.admission_context_id
    )
    population = tuple(
        assessment for assessment in generation.population_assessments
        if assessment.hypothesis.hypothesis_id == item.hypothesis.hypothesis_id
    )
    if len(searches) != 1:
        return _eligibility_reasons_from_context(
            item,
            None,
            population[0] if len(population) == 1 else None,
            representations,
        )
    return _eligibility_reasons_from_context(
        item,
        searches[0],
        population[0] if len(population) == 1 else None,
        representations,
    )


def _group_candidates(
    candidates: tuple[_EligibleCandidate, ...],
    exclusions: tuple[SubmissionExclusion, ...],
):
    grouped = {}
    for candidate in candidates:
        grouped.setdefault(candidate.identity, []).append(candidate)
    selected, tail, rejected = [], [], list(exclusions)
    selected_pairs = set()
    selected_singletons = {}
    for members in grouped.values():
        candidate = members[0]
        hypothesis = candidate.constructed.hypothesis
        context = candidate.constructed.admission_context_id
        partial = len(candidate.alleles) == 1 and hypothesis.completeness is Completeness.INCOMPLETE
        if partial and (context, hypothesis.allele_ids[0]) in selected_pairs:
            rejected.extend(SubmissionExclusion(
                item.constructed.hypothesis.hypothesis_id,
                item.source_rank,
                item.source_score,
                ("complete-pair-not-split",),
            ) for item in members)
            continue
        group = SubmissionGroup(
            hypothesis.hypothesis_id,
            candidate.source_rank,
            candidate.source_score,
            candidate.alleles,
            tuple(item.constructed.hypothesis.hypothesis_id for item in members[1:]),
            partial,
        )
        if len(selected) == 10:
            tail.append(group)
            continue
        if len(candidate.alleles) == 2 and hypothesis.completeness is Completeness.COMPLETE:
            for allele_id in hypothesis.allele_ids:
                key = (context, allele_id)
                earlier = selected_singletons.pop(key, None)
                if earlier is not None:
                    selected.remove(earlier[0])
                    rejected.extend(SubmissionExclusion(
                        item.constructed.hypothesis.hypothesis_id,
                        item.source_rank,
                        item.source_score,
                        ("complete-pair-not-split",),
                    ) for item in earlier[1])
                selected_pairs.add(key)
        selected.append(group)
        if partial:
            selected_singletons[(context, hypothesis.allele_ids[0])] = (group, members)
    return tuple(selected), tuple(tail), tuple(sorted(rejected, key=lambda item: item.source_rank))


def _select_replayed(
    generation: HypothesisGeneration,
    ranking: FusionRanking,
    reference: _OpenSubmissionReference,
    *,
    pin: PinnedPublicSource,
) -> SubmissionSelection:
    by_id = {item.hypothesis.hypothesis_id: item for item in generation.hypotheses}
    ledgers = {item.hypothesis_id: item for item in ranking.ledgers}
    if set(by_id) != set(ranking.order) or set(ledgers) != set(ranking.order):
        raise ValueError
    cache = {}
    assessments = tuple(sorted((
        _assess_representation(allele, generation, reference, cache, pin)
        for allele in generation.alleles
    ), key=lambda item: item.allele_id))
    representations = {item.allele_id: item for item in assessments}

    eligible = []
    excluded = []
    for rank, hypothesis_id in enumerate(ranking.order, 1):
        item = by_id[hypothesis_id]
        ledger = ledgers[hypothesis_id]
        reasons = _eligibility_reasons(
            item,
            generation,
            representations,
        )
        if reasons:
            excluded.append(SubmissionExclusion(hypothesis_id, rank, ledger.total, reasons))
            continue
        traces = tuple(representations[allele_ref.allele_id].trace
                       for allele_ref in item.hypothesis.allele_refs)
        if any(trace is None for trace in traces):
            raise ValueError
        traces = tuple(sorted(traces, key=_allele_trace_key))
        if len({trace.rendered for trace in traces}) != len(traces):
            raise ValueError
        eligible.append(_EligibleCandidate(item, rank, ledger.total, traces))
    selected, tail, excluded = _group_candidates(tuple(eligible), tuple(excluded))
    return SubmissionSelection(
        hashlib.sha256(generation.as_bytes()).hexdigest(),
        hashlib.sha256(ranking.as_bytes()).hexdigest(),
        ranking.policy_sha256,
        ranking.order,
        reference.binding,
        assessments,
        selected,
        tail,
        excluded,
    )


def _select_submission_hypotheses_for_synthetic_test(
    *,
    hypothesis_replay: HypothesisReplay,
    ranking_replay: SubmissionRankingReplay,
    reference: PinnedPublicSource,
    expected_public_identity_id: str,
    expected_result: SubmissionSelection | None = None,
) -> SubmissionSelection:
    """Exercise selection with an independently synthetic reference identity."""

    return _select_submission_hypotheses(
        hypothesis_replay=hypothesis_replay,
        ranking_replay=ranking_replay,
        reference=reference,
        expected_public_identity_id=expected_public_identity_id,
        expected_result=expected_result,
    )


def select_submission_hypotheses(
    *,
    hypothesis_replay: HypothesisReplay,
    ranking_replay: SubmissionRankingReplay,
    reference: PinnedPublicSource,
    expected_result: SubmissionSelection | None = None,
) -> SubmissionSelection:
    """Replay and select using the required production submission reference."""

    return _select_submission_hypotheses(
        hypothesis_replay=hypothesis_replay,
        ranking_replay=ranking_replay,
        reference=reference,
        expected_public_identity_id=PUBLIC_REFERENCE_IDENTITY,
        expected_result=expected_result,
    )


def _select_submission_hypotheses(
    *,
    hypothesis_replay: HypothesisReplay,
    ranking_replay: SubmissionRankingReplay,
    reference: PinnedPublicSource,
    expected_public_identity_id: str,
    expected_result: SubmissionSelection | None,
) -> SubmissionSelection:
    """Shared implementation; only the private synthetic-test boundary varies identity.

    No caller-supplied eligibility, order, score, rendered tuple, or standalone
    authorization claim is accepted.
    """

    try:
        if type(hypothesis_replay) is not HypothesisReplay or type(ranking_replay) is not SubmissionRankingReplay:
            raise ValueError
        if ranking_replay.ic_store != hypothesis_replay.ic_store:
            raise ValueError
        generation = replay_genotype_hypotheses(hypothesis_replay)
        ranking = rank_hypotheses(
            ranking_replay.document,
            replay=hypothesis_replay,
            ic_store=ranking_replay.ic_store,
            profile_evidence_id=ranking_replay.profile_evidence_id,
            policy_bytes=ranking_replay.policy_bytes,
            policy_sha256=ranking_replay.policy_sha256,
            derivation_bytes=ranking_replay.derivation_bytes,
            expected_derivation_sha256=ranking_replay.expected_derivation_sha256,
            expected_result=ranking_replay.expected_result,
        )
        with _open_submission_reference(
            reference,
            expected_public_identity_id=expected_public_identity_id,
        ) as opened:
            result = _select_replayed(
                generation,
                ranking,
                opened,
                pin=reference,
            )
        if expected_result is not None and result.as_bytes() != _typed(expected_result, SubmissionSelection).as_bytes():
            raise ValueError
        return result
    except _ERRORS:
        raise Track1SubmissionError("source-replayed submission selection failed") from None
