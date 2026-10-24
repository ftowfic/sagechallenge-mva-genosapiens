"""Canonical supplied-call diagnostics, never rescue execution or authorization.

Every scientific claim is replayed from independently selected source inputs.
Reports retain per-file evidence and selected samples; they never synthesize a
combined VCF, QC report, concordance report, or raw-file identity.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction
import hashlib
from typing import Literal

from genosapiens.admission_inventory import (
    SuppliedCallInputs, SuppliedSourceBinding, SuppliedSourceSetInputs,
    admission_inventory_as_json, build_admission_inventory,
)
from genosapiens.candidate_admission import _bytes
from genosapiens.evidence_fusion import (
    FusionRanking, _PreparedFusion, _is_prepared_fusion, rank_hypotheses,
)
from genosapiens.evidence_model import (
    Completeness, EnumerationStatus, EvidenceDocument, EvidenceReference,
    _Controlled, _exact, _ids, _sha, _typed,
)
from genosapiens.genotype_hypotheses import (
    HypothesisGeneration, HypothesisReplay, replay_genotype_hypotheses,
)
from genosapiens.genotype_properties import selected_source_calls
from genosapiens.hpo_ic import HpoICStore
from genosapiens.reference_concordance import (
    PublicValidationConcordanceReport, ReferenceConcordanceError,
    ReferenceConcordanceReport,
    ReferenceConstructionReport, _is_order_only_report,
    concordance_report_as_json,
)
from genosapiens.variant_inventory import (
    AssessmentStatus, VariantClassInventory, WgsClass,
    build_variant_class_inventory, variant_inventory_as_json,
)
from genosapiens.variant_normalization import variant_normalization_as_json
from genosapiens.vcf_qc import VcfQcError, VcfQcReport, vcf_qc_as_json


SCHEMA = "genosapiens.track1-diagnostics"
State = Literal["present", "absent", "unknown", "not_evaluated"]
ConcordanceStatus = Literal[
    "reported_accepted", "failed", "blocked_order_adjudication", "not_evaluated",
]
QcStatus = Literal["passed", "failed", "unknown", "not_evaluated"]
_INDICATION_ORDER = (
    "vcf-qc-failure",
    "absent-credible-complete-hypotheses",
    "incomplete-recessive-hypotheses",
    "phenotype-supported-unassessed-mechanism",
    "explicit-operator-request",
)


class DiagnosticsError(ValueError):
    """Static errors only; invalid evidence is not an assessed negative."""


@dataclass(frozen=True, slots=True, repr=False)
class MechanismRule(_Controlled):
    """Exact source-requirement token mapping selected in the policy."""

    field: Literal["allelic_requirement", "consequence_requirement", "mechanism_direction"]
    value: str
    variant_class: WgsClass

    def __post_init__(self):
        _Controlled.__post_init__(self)
        _ids((self.value,))


@dataclass(frozen=True, slots=True, repr=False)
class DiagnosticPolicy(_Controlled):
    """Explicit development settings; this module intentionally has no default."""

    policy_id: str
    maximum_uncallable_fraction: str | None
    maximum_filtered_fraction: str | None
    minimum_credible_score: str | None
    mechanism_rules: tuple[MechanismRule, ...]

    def __post_init__(self):
        _Controlled.__post_init__(self)
        _ids((self.policy_id,))
        for value in (self.maximum_uncallable_fraction, self.maximum_filtered_fraction):
            if value is not None and not 0 <= Fraction(value) <= 1:
                raise DiagnosticsError("invalid diagnostic QC threshold")
        if self.minimum_credible_score is not None:
            _exact(self.minimum_credible_score)
        if len(set(self.mechanism_rules)) != len(self.mechanism_rules):
            raise DiagnosticsError("duplicate diagnostic mechanism rule")

    @property
    def sha256(self) -> str:
        return hashlib.sha256(_bytes({"schema": SCHEMA, "policy": self})).hexdigest()


@dataclass(frozen=True, slots=True, repr=False)
class SourceDiagnosticEvidence:
    """One source's actual stage products, keyed by its raw validated-file hash.

    Exceptions record a failed stage without retaining their potentially
    sensitive messages. They are inputs only and never enter report bytes.
    """

    input_sha256: str
    selected_sample_name: str
    concordance: (
        ReferenceConcordanceReport
        | PublicValidationConcordanceReport
        | ReferenceConcordanceError
        | None
    )
    qc: VcfQcReport | VcfQcError | None
    class_inventory: VariantClassInventory | None

    def __post_init__(self):
        try:
            _sha(self.input_sha256)
            _ids((self.selected_sample_name,))
            if not isinstance(
                self.concordance,
                (
                    ReferenceConcordanceReport,
                    PublicValidationConcordanceReport,
                    ReferenceConcordanceError,
                    type(None),
                ),
            ):
                raise ValueError
            if not isinstance(self.qc, (VcfQcReport, VcfQcError, type(None))):
                raise ValueError
            if self.class_inventory is not None and type(self.class_inventory) is not VariantClassInventory:
                raise ValueError
        except Exception:
            raise DiagnosticsError("invalid source diagnostic evidence") from None

    def __repr__(self):
        return "<controlled evidence; no-llm>"


@dataclass(frozen=True, slots=True, repr=False)
class RankingReplay:
    """Explicit source-bound fusion inputs; no persisted score is trusted alone."""

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
class SourceDiagnostic(_Controlled):
    manifest_index: int
    manifest_member_sha256: str | None
    input_sha256: str
    selected_sample_name: str
    reference_role: str | None
    reference_identity_id: str | None
    reference_manifest_sha256: str | None
    source_format_sha256: str | None
    concordance_report: (
        ReferenceConcordanceReport | PublicValidationConcordanceReport | None
    )
    concordance_status: ConcordanceStatus
    qc_report: VcfQcReport | None
    qc_status: QcStatus
    qc_measurements: Mapping[str, int | None]
    class_inventory: VariantClassInventory | None
    blockers: tuple[str, ...]

    def __post_init__(self):
        _Controlled.__post_init__(self)
        if self.manifest_index < 0:
            raise DiagnosticsError("invalid source diagnostic index")
        for value in (
            self.manifest_member_sha256, self.input_sha256,
            self.reference_manifest_sha256, self.source_format_sha256,
        ):
            if value is not None:
                _sha(value)
        _ids((self.selected_sample_name,))
        if self.reference_identity_id is not None:
            _ids((self.reference_identity_id,))


@dataclass(frozen=True, slots=True, repr=False)
class MechanismIndication(_Controlled):
    context_id: str
    rule: MechanismRule
    phenotype_proof_ids: tuple[str, ...]
    requirement_refs: tuple[EvidenceReference, ...]
    source_input_sha256: str
    inventory_status: AssessmentStatus
    inventory_scope: str


@dataclass(frozen=True, slots=True, repr=False)
class TechnicalObligation(_Controlled):
    source_input_sha256: str
    sample_name: str
    allele_id: str
    native_identity_id: str
    source_record_id: str
    source_allele_id: str
    reasons: tuple[str, ...]
    scope: Literal["supplied-exact-event"] = "supplied-exact-event"
    status: Literal["unmet"] = "unmet"


@dataclass(frozen=True, slots=True, repr=False)
class DiagnosticReport(_Controlled):
    policy: DiagnosticPolicy
    policy_sha256: str
    construction_report: ReferenceConstructionReport | None
    construction_status: Literal["reported_accepted", "failed", "not_evaluated"]
    source_manifest_sha256: str | None
    source_roster_status: Literal["verified", "not_evaluated"]
    source_diagnostics: tuple[SourceDiagnostic, ...]
    concordance_status: ConcordanceStatus
    qc_status: QcStatus
    hypothesis_evaluation: Literal["evaluated", "evaluated_empty", "not_evaluated"]
    generation_sha256: str | None
    hypothesis_count: int | None
    ranking_status: Literal["evaluated", "evaluated_empty", "not_evaluated"]
    ranking_sha256: str | None
    ranking_policy_sha256: str | None
    credibility_threshold: str | None
    credible_complete_status: State
    credible_complete_ids: tuple[str, ...]
    incomplete_recessive_ids: tuple[str, ...]
    mechanism_indications: tuple[MechanismIndication, ...]
    operator_request_id: str | None
    indications: tuple[str, ...]
    prerequisite_blockers: tuple[str, ...]
    selected_hypothesis_blockers: tuple[tuple[str, tuple[str, ...]], ...]
    technical_obligations: tuple[TechnicalObligation, ...]
    schema: Literal["genosapiens.track1-diagnostics"] = SCHEMA
    validation_status: Literal["engineering-only-not-controlled-use-frozen"] = "engineering-only-not-controlled-use-frozen"
    rescue_authorization: Literal["not_issued"] = "not_issued"

    def as_bytes(self) -> bytes:
        return _bytes({"classification": "controlled", "llm_policy": "no-llm", "report": self})


def _hash_report(
    report: ReferenceConcordanceReport | PublicValidationConcordanceReport,
) -> str:
    return hashlib.sha256(concordance_report_as_json(report).encode("ascii")).hexdigest()


def _source_roster(
    replay: HypothesisReplay | None,
    source_inputs: SuppliedCallInputs | SuppliedSourceSetInputs | None,
    analysis_sample_name: str | None,
):
    if replay is not None:
        if source_inputs is not None and source_inputs != replay.inputs:
            raise ValueError
        if analysis_sample_name is not None and analysis_sample_name != replay.sample_name:
            raise ValueError
        source_inputs = replay.inputs
        analysis_sample_name = replay.sample_name
    if source_inputs is None:
        if analysis_sample_name is not None:
            raise ValueError
        return None, ()
    inputs = source_inputs
    if analysis_sample_name is None:
        raise ValueError
    if type(inputs) is SuppliedCallInputs:
        return None, ((None, inputs, analysis_sample_name),)
    if type(inputs) is not SuppliedSourceSetInputs:
        raise ValueError
    manifest = inputs.manifest
    if analysis_sample_name != manifest.sample_id:
        raise ValueError
    by_hash = {source.input_sha256: source for source in inputs.sources}
    return manifest.sha256, tuple(
        (member, by_hash[member.input_sha256], member.sample_name)
        for member in manifest.members
    )


def _member_matches(member: SuppliedSourceBinding, source: SuppliedCallInputs) -> bool:
    normal = source.normalization
    inventory = build_admission_inventory(source)
    return (
        member.input_sha256 == source.input_sha256
        and member.input_size_bytes == len(source.source_bytes)
        and member.reference_role.value == normal.reference_role
        and member.reference_identity_id == normal.reference_identity_id
        and member.reference_manifest_sha256 == normal.reference_manifest_sha256
        and member.reference_concordance_sha256 == source.qc_report.reference_concordance_sha256
        and member.vcf_qc_sha256 == hashlib.sha256(vcf_qc_as_json(source.qc_report).encode("ascii")).hexdigest()
        and member.normalization_sha256 == hashlib.sha256(
            variant_normalization_as_json(normal).encode("ascii")
        ).hexdigest()
        and member.class_inventory_sha256 == hashlib.sha256(
            variant_inventory_as_json(source.class_inventory).encode("ascii")
        ).hexdigest()
        and member.source_format_sha256 == source.source_format_sha256
        and member.inventory_sha256 == hashlib.sha256(
            admission_inventory_as_json(inventory).encode("ascii")
        ).hexdigest()
    )


def _source_diagnostic(
    evidence: SourceDiagnosticEvidence, *, index: int,
    expected: tuple[SuppliedSourceBinding | None, SuppliedCallInputs, str] | None,
    policy: DiagnosticPolicy,
) -> SourceDiagnostic:
    member = source = None
    if expected is not None:
        member, source, sample_name = expected
        if ((evidence.input_sha256, evidence.selected_sample_name)
                != (source.input_sha256, sample_name)
                or evidence.qc != source.qc_report
                or evidence.class_inventory != source.class_inventory
                or member is not None and not _member_matches(member, source)):
            raise ValueError

    concordance_failed = isinstance(evidence.concordance, ReferenceConcordanceError)
    report = None if concordance_failed else evidence.concordance
    if not isinstance(
        report,
        (ReferenceConcordanceReport, PublicValidationConcordanceReport, type(None)),
    ):
        raise ValueError
    concordance_status: ConcordanceStatus = (
        "failed" if concordance_failed else
        "not_evaluated" if report is None else
        "reported_accepted" if report.accepted else
        "blocked_order_adjudication"
        if isinstance(report, ReferenceConcordanceReport)
        and _is_order_only_report(report)
        else "failed"
    )
    blockers = []
    if concordance_status != "reported_accepted":
        blockers.append("reference-concordance-" + concordance_status)

    qc_failed = isinstance(evidence.qc, VcfQcError)
    qc_report = None if qc_failed else _typed(evidence.qc, VcfQcReport | None)
    measurements: dict[str, int | None] = {}
    qc_status: QcStatus = "failed" if qc_failed else "not_evaluated"
    if qc_report is not None:
        if qc_report.input_sha256 != evidence.input_sha256:
            raise ValueError
        try:
            summary = next(
                item for item in qc_report.sample_summaries
                if item.sample_name == evidence.selected_sample_name
            )
        except StopIteration:
            raise ValueError from None
        filtered = sum(record.filter_status.value == "filtered" for record in qc_report.records)
        measurements = {
            "file_sample_count": qc_report.sample_count,
            "records": qc_report.total_records,
            "selected_sample_calls": summary.records,
            "selected_sample_uncallable_calls": summary.uncallable_records,
            "selected_sample_unknown_callability_calls": summary.unknown_callability_records,
            "filtered_records": filtered,
        }
        qc_status = "passed"
        if qc_report.sample_count != 1:
            qc_status = "failed"
            blockers.append("source-input-sample-count-not-one")
        for limit, lower, upper, denominator in (
            (policy.maximum_uncallable_fraction, summary.uncallable_records,
             summary.uncallable_records + summary.unknown_callability_records, summary.records),
            (policy.maximum_filtered_fraction, filtered, filtered, qc_report.total_records),
        ):
            if limit is not None:
                if Fraction(lower, denominator) > Fraction(limit):
                    qc_status = "failed"
                elif Fraction(upper, denominator) > Fraction(limit) and qc_status != "failed":
                    qc_status = "unknown"
        if (report is None or report.validated_input_sha256 != evidence.input_sha256
                or qc_report.total_records != report.total_records
                or qc_report.reference_concordance_sha256 != _hash_report(report)):
            blockers.append("qc-concordance-binding-unverified")
        if evidence.class_inventory is None or evidence.class_inventory != build_variant_class_inventory(
            qc_report, allele_balance_screen=evidence.class_inventory.allele_balance_screen.enabled,
        ):
            raise ValueError
    elif evidence.class_inventory is not None:
        raise ValueError
    if qc_status != "passed":
        blockers.append("vcf-qc-" + qc_status)

    normal = None if source is None else source.normalization
    reference_role = normal.reference_role if normal is not None else (None if report is None else report.reference_role)
    reference_identity_id = normal.reference_identity_id if normal is not None else (
        None if report is None else report.reference_identity_id
    )
    reference_manifest_sha256 = normal.reference_manifest_sha256 if normal is not None else (
        None if report is None else report.reference_manifest_sha256
    )
    if report is not None and (
        report.validated_input_sha256 != evidence.input_sha256
        or report.reference_role != reference_role
        or report.reference_identity_id != reference_identity_id
        or report.reference_manifest_sha256 != reference_manifest_sha256
    ):
        blockers.append("source-reference-binding-mismatch")
    member_sha = None if member is None else hashlib.sha256(_bytes(member)).hexdigest()
    return SourceDiagnostic(
        index, member_sha, evidence.input_sha256, evidence.selected_sample_name,
        reference_role, reference_identity_id, reference_manifest_sha256,
        None if source is None else source.source_format_sha256,
        report, concordance_status, qc_report, qc_status, measurements,
        evidence.class_inventory, tuple(dict.fromkeys(blockers)),
    )


def _aggregate_status(values, *, not_evaluated):
    if not values or all(value == "not_evaluated" for value in values):
        return "not_evaluated"
    if "failed" in values:
        return "failed"
    if "blocked_order_adjudication" in values:
        return "blocked_order_adjudication"
    if "unknown" in values:
        return "unknown"
    if "not_evaluated" in values:
        return not_evaluated
    return values[0]


def _technical_obligations(generation: HypothesisGeneration) -> tuple[TechnicalObligation, ...]:
    alleles = {item.allele_id: item for item in generation.alleles}
    properties = {item.allele_id: item for item in generation.genotypes.allele_properties}
    source_properties = {
        item.source_allele_id: item for item in generation.genotypes.source_properties
    }
    calls = {}
    for allele in generation.alleles:
        for call in selected_source_calls(allele, generation.genotypes):
            source_id = call.source.source_fields["source_allele_id"]
            source_calls = calls.setdefault(source_id, {})
            if (
                allele.allele_id in source_calls
                and source_calls[allele.allele_id] != call
            ):
                raise ValueError
            source_calls[allele.allele_id] = call

    obligations = []
    for allele_id, prop in properties.items():
        source_ids = prop.event_source_allele_ids or prop.source_allele_ids
        for source_id in source_ids:
            source = source_properties[source_id]
            if not source.read_validation_required:
                continue
            source_calls = calls.get(source_id, {})
            call = source_calls.get(allele_id)
            if call is None and len(source_calls) == 1:
                call = next(iter(source_calls.values()))
            if call is None or call.source.source_fields["read_validation_required"] is not True:
                raise ValueError
            reasons = tuple(call.source.source_fields["technical_flags"])
            if not reasons:
                reasons = ("source-allele-read-validation-required",)
            obligations.append(TechnicalObligation(
                source.input_sha256 or generation.genotypes.input_sha256,
                source.sample.sample_name, allele_id, alleles[allele_id].native_identity_id,
                source.source_record_id, source_id, reasons,
            ))
        if prop.read_validation_required != any(
            source_properties[source_id].read_validation_required for source_id in source_ids
        ):
            raise ValueError
    obligations = tuple(sorted(set(obligations), key=lambda item: (
        item.allele_id, item.source_input_sha256, item.source_allele_id,
    )))
    for item in generation.hypotheses:
        matching = any(
            obligation.allele_id in item.hypothesis.evidence_subject_ids
            for obligation in obligations
        )
        if item.read_validation_required != matching:
            raise ValueError
    return obligations


def _evaluate_diagnostics(
    *, policy: DiagnosticPolicy, expected_policy_sha256: str,
    construction_report: ReferenceConstructionReport | ReferenceConcordanceError | None,
    source_evidence: tuple[SourceDiagnosticEvidence, ...],
    source_inputs: SuppliedCallInputs | SuppliedSourceSetInputs | None = None,
    analysis_sample_name: str | None = None,
    hypothesis_replay: HypothesisReplay | None = None,
    ranking_replay: RankingReplay | None = None,
    selected_hypothesis_ids: tuple[str, ...] = (),
    operator_request_id: str | None = None,
    prepared_fusion: _PreparedFusion | None = None,
) -> DiagnosticReport:
    """Replay diagnostic claims while keeping indications separate from blockers."""
    try:
        policy = _typed(policy, DiagnosticPolicy)
        if policy.sha256 != expected_policy_sha256:
            raise ValueError
        source_evidence = _typed(source_evidence, tuple[SourceDiagnosticEvidence, ...])
        selected_hypothesis_ids = _typed(selected_hypothesis_ids, tuple[str, ...])
        _ids(selected_hypothesis_ids)
        if operator_request_id is not None:
            operator_request_id = _typed(operator_request_id, str)
            _ids((operator_request_id,))
        if hypothesis_replay is not None and type(hypothesis_replay) is not HypothesisReplay:
            raise ValueError
        if ranking_replay is not None and type(ranking_replay) is not RankingReplay:
            raise ValueError
        if ranking_replay is not None and hypothesis_replay is None:
            raise ValueError
        if prepared_fusion is not None and (
            not _is_prepared_fusion(prepared_fusion)
            or hypothesis_replay != prepared_fusion.replay
            or ranking_replay is not None
        ):
            raise ValueError

        construction_failed = isinstance(construction_report, ReferenceConcordanceError)
        construction_report = None if construction_failed else _typed(
            construction_report, ReferenceConstructionReport | None,
        )
        construction_status = (
            "failed" if construction_failed else
            "not_evaluated" if construction_report is None else
            "reported_accepted" if construction_report.accepted else "failed"
        )
        blockers = []
        if construction_status != "reported_accepted":
            blockers.append("reference-construction-" + construction_status)

        manifest_sha256, expected_sources = _source_roster(
            hypothesis_replay,
            source_inputs,
            analysis_sample_name,
        )
        if expected_sources and tuple(
            (item.input_sha256, item.selected_sample_name) for item in source_evidence
        ) != tuple((source.input_sha256, sample) for _member, source, sample in expected_sources):
            raise ValueError
        diagnostics = tuple(
            _source_diagnostic(
                item, index=index,
                expected=expected_sources[index] if expected_sources else None,
                policy=policy,
            )
            for index, item in enumerate(source_evidence)
        )
        if expected_sources and len(diagnostics) != len(expected_sources):
            raise ValueError
        concordance_status = _aggregate_status(
            tuple(item.concordance_status for item in diagnostics),
            not_evaluated="failed",
        )
        qc_status = _aggregate_status(
            tuple(item.qc_status for item in diagnostics),
            not_evaluated="unknown",
        )
        if not diagnostics:
            blockers.extend(("reference-concordance-not_evaluated", "vcf-qc-not_evaluated"))
        else:
            if concordance_status != "reported_accepted":
                blockers.append("reference-concordance-" + concordance_status)
            if qc_status != "passed":
                blockers.append("vcf-qc-" + qc_status)
            if any("qc-concordance-binding-unverified" in item.blockers for item in diagnostics):
                blockers.append("qc-concordance-binding-unverified")
            if any("source-reference-binding-mismatch" in item.blockers for item in diagnostics):
                blockers.append("source-reference-binding-mismatch")
        for item in diagnostics:
            report = item.concordance_report
            expected_identity = (
                construction_report.public_annotation_benchmark_identity_id
                if construction_report is not None
                and isinstance(report, PublicValidationConcordanceReport)
                else None
                if construction_report is None
                else construction_report.challenge_alignment_identity_id
            )
            reported_reference_contigs = (
                report.available_reference_contigs
                if isinstance(report, PublicValidationConcordanceReport)
                else None if report is None else report.reference_contigs
            )
            if construction_report is not None and report is not None and (
                expected_identity != report.reference_identity_id
                or construction_report.total_records != reported_reference_contigs
            ):
                blockers.append("reference-report-binding-mismatch")
        blockers = list(dict.fromkeys(blockers))

        indications = []
        if qc_status == "failed":
            indications.append("vcf-qc-failure")

        generation = None
        evaluation = "not_evaluated"
        if hypothesis_replay is not None and not blockers:
            generation = (prepared_fusion.generation if prepared_fusion is not None
                          else replay_genotype_hypotheses(hypothesis_replay))
            evaluation = "evaluated" if generation.hypotheses else "evaluated_empty"

        incomplete = []
        mechanisms = []
        obligations = ()
        if generation is not None:
            obligations = _technical_obligations(generation)
            for item in generation.hypotheses:
                if item.hypothesis.completeness is Completeness.INCOMPLETE and any(
                    fact.family in ("recessive", "x_recessive")
                    for fact in item.model.requirement_facts
                ):
                    incomplete.append(item.hypothesis.hypothesis_id)
            if incomplete:
                indications.append("incomplete-recessive-hypotheses")

            inventories = tuple(
                (item.input_sha256, item.class_inventory) for item in diagnostics
                if item.class_inventory is not None
            )
            for seed in generation.admission.contexts:
                proofs = tuple(
                    proof.proof_id for proof in generation.admission.phenotype_proofs
                    if proof.proof_id in seed.phenotype_proof_ids
                    and proof.outcome == "admitted"
                    and seed.context.context_id in proof.context_ids
                    and proof.information_content > 0
                )
                if not proofs:
                    continue
                for rule in policy.mechanism_rules:
                    refs = tuple(sorted({
                        ref for requirement in seed.context.requirements
                        if getattr(requirement, rule.field) == rule.value
                        for ref in seed.context.evidence_refs
                        if ref.record_id in requirement.record_ids
                    }, key=lambda ref: (ref.evidence_id, ref.record_id)))
                    statuses = tuple(
                        (input_sha256, next(entry for entry in inventory.classes
                                            if entry.variant_class is rule.variant_class))
                        for input_sha256, inventory in inventories
                    )
                    if refs and statuses and all(entry.status in (
                        AssessmentStatus.NOT_ASSESSED, AssessmentStatus.UNSUPPORTED,
                    ) for _input_sha256, entry in statuses):
                        mechanisms.extend(MechanismIndication(
                            seed.context.context_id, rule, proofs, refs, input_sha256,
                            entry.status, entry.scope,
                        ) for input_sha256, entry in statuses)
            if mechanisms:
                indications.append("phenotype-supported-unassessed-mechanism")

        ranking = None
        ranking_status = "not_evaluated"
        if generation is not None and prepared_fusion is not None:
            ranking = prepared_fusion.ranking
            ranking_status = "evaluated" if ranking.ledgers else "evaluated_empty"
        elif generation is not None and ranking_replay is not None:
            if ranking_replay.ic_store != hypothesis_replay.ic_store:
                raise ValueError
            ranking = rank_hypotheses(
                ranking_replay.document, replay=hypothesis_replay,
                ic_store=ranking_replay.ic_store,
                profile_evidence_id=ranking_replay.profile_evidence_id,
                policy_bytes=ranking_replay.policy_bytes,
                policy_sha256=ranking_replay.policy_sha256,
                derivation_bytes=ranking_replay.derivation_bytes,
                expected_derivation_sha256=ranking_replay.expected_derivation_sha256,
                expected_result=ranking_replay.expected_result,
            )
            ranking_status = "evaluated" if ranking.ledgers else "evaluated_empty"

        credible_status: State = "not_evaluated"
        credible = []
        if ranking is not None and policy.minimum_credible_score is not None:
            by_id = {item.hypothesis.hypothesis_id: item for item in generation.hypotheses}
            ledgers = {item.hypothesis_id: item for item in ranking.ledgers}
            for hypothesis_id in ranking.order:
                item = by_id[hypothesis_id]
                if (
                    item.hypothesis.completeness is Completeness.COMPLETE
                    and item.model.state == "supported"
                    and all(member.state == "qualifying" for member in item.members)
                    and _exact(ledgers[hypothesis_id].total) >= _exact(policy.minimum_credible_score)
                ):
                    credible.append(hypothesis_id)
            uncertain = (
                any(search.enumeration_status is not EnumerationStatus.COMPLETE
                    or search.unresolved_source_allele_ids
                    or search.unresolved_relevance_source_allele_ids
                    for search in generation.searches)
                or bool(generation.admission.inventory.unresolved_entries)
                or any(
                    item.hypothesis.completeness is Completeness.UNKNOWN
                    or item.hypothesis.completeness is not Completeness.INCOMPLETE
                    and (item.model.state == "unresolved"
                         or any(member.state == "unresolved" for member in item.members))
                    for item in generation.hypotheses
                )
            )
            credible_status = "present" if credible else "unknown" if uncertain else "absent"
            if credible_status == "absent":
                indications.append("absent-credible-complete-hypotheses")

        by_id = {} if generation is None else {
            item.hypothesis.hypothesis_id: item for item in generation.hypotheses
        }
        selected_blockers = []
        for hypothesis_id in sorted(selected_hypothesis_ids):
            reasons = []
            item = by_id.get(hypothesis_id)
            if item is None:
                reasons.append("selected-hypothesis-not-evaluated")
            else:
                hypothesis = item.hypothesis
                if hypothesis.locus_call_refs:
                    reasons.append("locus-call-is-not-an-exact-submission-allele")
                if hypothesis.enumeration_status is not EnumerationStatus.COMPLETE:
                    reasons.append("required-enumeration-incomplete")
                if item.model.state != "supported":
                    reasons.append("model-construction-" + item.model.state)
                if hypothesis.completeness is Completeness.INCOMPLETE and (
                    len(hypothesis.allele_refs) != 1
                    or hypothesis.context.model != "incomplete"
                    or item.incompleteness_penalty is None
                ):
                    reasons.append("incomplete-partial-model-invalid")
                if any(
                    obligation.allele_id in hypothesis.evidence_subject_ids
                    for obligation in obligations
                ):
                    reasons.append("selected-allele-mandatory-read-validation-unmet")
            if reasons:
                selected_blockers.append((hypothesis_id, tuple(reasons)))

        if operator_request_id is not None:
            indications.append("explicit-operator-request")
        indications = tuple(item for item in _INDICATION_ORDER if item in indications)
        return DiagnosticReport(
            policy, policy.sha256, construction_report, construction_status,
            manifest_sha256, "verified" if expected_sources else "not_evaluated",
            diagnostics, concordance_status, qc_status, evaluation,
            None if generation is None else hashlib.sha256(generation.as_bytes()).hexdigest(),
            None if generation is None else len(generation.hypotheses),
            ranking_status,
            None if ranking is None else hashlib.sha256(ranking.as_bytes()).hexdigest(),
            None if ranking is None else ranking.policy_sha256,
            policy.minimum_credible_score,
            credible_status, tuple(credible), tuple(sorted(incomplete)),
            tuple(sorted(mechanisms, key=lambda item: (
                item.context_id, item.rule.field, item.rule.value, item.source_input_sha256,
            ))),
            operator_request_id, indications, tuple(blockers),
            tuple(selected_blockers), obligations,
        )
    except Exception:
        raise DiagnosticsError("diagnostic evidence or policy validation failed") from None


def evaluate_diagnostics(
    *, policy: DiagnosticPolicy, expected_policy_sha256: str,
    construction_report: ReferenceConstructionReport | ReferenceConcordanceError | None,
    source_evidence: tuple[SourceDiagnosticEvidence, ...],
    source_inputs: SuppliedCallInputs | SuppliedSourceSetInputs | None = None,
    analysis_sample_name: str | None = None,
    hypothesis_replay: HypothesisReplay | None = None,
    ranking_replay: RankingReplay | None = None,
    selected_hypothesis_ids: tuple[str, ...] = (),
    operator_request_id: str | None = None,
) -> DiagnosticReport:
    """Standalone boundary: authenticate imported generation and ranking by replay."""
    return _evaluate_diagnostics(
        policy=policy,
        expected_policy_sha256=expected_policy_sha256,
        construction_report=construction_report,
        source_evidence=source_evidence,
        source_inputs=source_inputs,
        analysis_sample_name=analysis_sample_name,
        hypothesis_replay=hypothesis_replay,
        ranking_replay=ranking_replay,
        selected_hypothesis_ids=selected_hypothesis_ids,
        operator_request_id=operator_request_id,
    )


def _evaluate_prepared_diagnostics(
    prepared_fusion: _PreparedFusion, *, policy: DiagnosticPolicy,
    expected_policy_sha256: str,
    construction_report: ReferenceConstructionReport | ReferenceConcordanceError | None,
    source_evidence: tuple[SourceDiagnosticEvidence, ...],
    selected_hypothesis_ids: tuple[str, ...] = (),
    operator_request_id: str | None = None,
) -> DiagnosticReport:
    """Run-local boundary used only after `_prepared_fusion` authenticates sources."""
    return _evaluate_diagnostics(
        policy=policy,
        expected_policy_sha256=expected_policy_sha256,
        construction_report=construction_report,
        source_evidence=source_evidence,
        hypothesis_replay=prepared_fusion.replay,
        selected_hypothesis_ids=selected_hypothesis_ids,
        operator_request_id=operator_request_id,
        prepared_fusion=prepared_fusion,
    )
