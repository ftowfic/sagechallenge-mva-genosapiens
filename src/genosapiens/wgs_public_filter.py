"""Bounded production adapter for the approved public common-allele filter."""

from __future__ import annotations

import copyreg
from dataclasses import asdict, dataclass, replace
from decimal import Decimal, InvalidOperation
from fractions import Fraction
import multiprocessing
import queue
import ctypes
import hashlib
import json
import os
from pathlib import Path
import pickle
import signal
import time
import tempfile
from types import MappingProxyType
from typing import Callable, Iterable, Iterator, Mapping

from genosapiens.common_allele_filter import (
    CLINICAL_ROLES,
    ClinicalProtectionFacts,
    CommonAlleleFilterPolicy,
    PopulationFilterFacts,
    PopulationNonContributionProof,
    PopulationSubsetFacts,
    PrefilterOccurrence,
    RationalEvidence,
    evaluate_common_allele_filter,
    common_allele_filter_policy_as_json,
    prefilter_bypass_reasons,
)
from genosapiens.track1_early_scope import EarlyScopeSession, EarlyScopeDecision
from genosapiens.common_allele_public_indexes import PublicIndexQuery, SequenceAlias
from genosapiens.clinical_annotation import affirmative_source_conflict
from genosapiens.submission_rendering import SubmissionContigRenderer
from genosapiens.evidence_model import (
    AlleleEvidence,
    IdentityMapping,
    IdentityStage,
    PopulationEvidence,
    SourcedRecord,
)
from genosapiens.gnomad_faf95_adapter import (
    GnomadFaf95AdapterError,
    GnomadFaf95PrecisionQualification,
    adapt_joint_faf95,
)
from genosapiens.population_annotation import PopulationAnnotationResult, PopulationRequest
from genosapiens.population_clinical_session import PopulationClinicalResourceSession
from genosapiens.vcf_qc import PhaseEvidence, PhaseStatus, VariantClass
from genosapiens.wgs_admission import (
    ADMISSION_INVENTORY_SCHEMA,
    _call,
    _normalized_sequence,
    _source_record,
    _source_sequence,
    _technical,
)
from genosapiens.wgs_filter import SourceBoundDecisionStream, ledger_filter_decision
from genosapiens.wgs_ledger import (
    Artifact, FilterDecision, MAX_ROW_BYTES, MAX_SHARD_BYTES, canonical,
    load_model, publish_bytes, publish_stream, read_artifact, iter_artifact_lines,
)
from genosapiens.wgs_preparation import PreparationFilterContext
from genosapiens.wgs_failures import WgsFailure, SecondaryFailure, failure_boundary, static_failure
from genosapiens.wgs_source_routing import (
    NativeSourceContract,
    SourceOccurrenceRoute,
    SourceRoutingContext,
    WgsSourceRoutingError,
)
from genosapiens.wgs_runtime import RuntimePhaseRecorder, runtime_phase


class WgsPublicFilterError(RuntimeError):
    """Static failure while adapting authenticated public filter evidence."""


def _restore_immutable_mapping(value):
    return MappingProxyType(value)


def _reduce_immutable_mapping(value):
    return _restore_immutable_mapping, (dict(value),)


copyreg.pickle(MappingProxyType, _reduce_immutable_mapping)


def _ratio(value: str | None) -> RationalEvidence | None:
    if value in {None, "."}:
        return None
    try:
        fraction = Fraction(Decimal(value))
        return RationalEvidence(fraction.numerator, fraction.denominator)
    except (InvalidOperation, ValueError, ArithmeticError):
        return None


def _raw_flag(value: object) -> bool | None:
    if value in {1, "1", True}:
        return True
    if value in {0, "0", False}:
        return False
    return None


def _subset(records: tuple[PopulationEvidence, ...], name: str, proof=None) -> PopulationSubsetFacts | None:
    values = tuple(
        value for value in records
        if value.subset == name and value.stratum == "adjusted"
    )
    if len(values) != 1:
        return None
    value = values[0]
    fields = value.source.source_fields
    raw_context = fields.get("raw_context") if isinstance(fields, Mapping) else None
    if proof is not None:
        from genosapiens.common_allele_filter import _canonical
        if (type(proof) is not PopulationNonContributionProof or name != 'exomes'
                or proof.raw_variant_ac != value.allele_count or proof.raw_variant_an != value.allele_number
                or not isinstance(raw_context, Mapping)
                or hashlib.sha256(_canonical(dict(raw_context))).hexdigest() != proof.source_fields_sha256):
            raise ValueError('noncontribution source binding mismatch')
    count_status = (
        "valid"
        if value.allele_count is not None
        and value.allele_number is not None
        and value.allele_count <= value.allele_number
        else "missing"
        if value.allele_count is None or value.allele_number is None
        else "inconsistent"
    )
    contribution = (
        "observed"
        if count_status == "valid" and value.allele_number > 0
        else "proven_no_callable"
        if count_status == "valid"
        and value.allele_count == 0
        and value.allele_number == 0
        else "unknown"
    )
    if (proof is not None and count_status == 'missing'
            and value.allele_count in (None, 0) and value.allele_number in (None, 0)):
        contribution = 'proven_no_callable'
    defined = fields.get("source_filters_defined") if isinstance(fields, Mapping) else None
    defined = True if defined == 1 else False if defined == 0 else None
    not_called_key = (
        "REGION_NOT_CALLED_IN_EXOMES"
        if name == "exomes" else "REGION_NOT_CALLED_IN_GENOMES"
    )
    not_called = _raw_flag(
        raw_context.get(not_called_key) if isinstance(raw_context, Mapping) else None
    )
    interval = None
    if name == "exomes":
        failed = _raw_flag(
            raw_context.get("REGION_FAIL_INTERVAL_QC")
            if isinstance(raw_context, Mapping) else None
        )
        interval = None if failed is None else not failed
    return PopulationSubsetFacts(
        name,
        value.allele_count,
        value.allele_number,
        count_status,
        contribution,
        defined,
        value.filters if defined is True else None,
        not_called,
        interval,
        proof,
    )


@failure_boundary("population_adaptation", "population_response_adapt")
def population_filter_facts(
    allele: AlleleEvidence,
    *,
    precision: GnomadFaf95PrecisionQualification,
    population_manifest_sha256: str,
    precision_manifest_sha256: str,
    noncontribution_proofs: tuple[PopulationNonContributionProof, ...] = (),
) -> PopulationFilterFacts:
    """Conservatively adapt one authenticated population annotation."""

    return _population_filter_records(allele.population or (), allele.identities,
        precision=precision, population_manifest_sha256=population_manifest_sha256,
        precision_manifest_sha256=precision_manifest_sha256, noncontribution_proofs=noncontribution_proofs)


def _population_filter_records(records, identities, *, precision, population_manifest_sha256,
        precision_manifest_sha256, noncontribution_proofs=()):
    """Shared policy adapter for authenticated full reports and private facts."""
    if type(noncontribution_proofs) is not tuple or len(noncontribution_proofs) > 1:
        raise ValueError('invalid noncontribution proof scope')
    proof = None if not noncontribution_proofs else noncontribution_proofs[0]
    if proof is not None:
        if type(proof) is not PopulationNonContributionProof:
            raise ValueError('invalid noncontribution proof type')
        from genosapiens.allele_identity import SequenceAllele
        mapping = next((value for value in identities if value.identity_id == proof.normalized_allele_id), None)
        if (type(proof) is not PopulationNonContributionProof or mapping is None
                or type(mapping.identity) is not SequenceAllele
                or mapping.reference_manifest_sha256 != proof.reference_manifest_sha256
                or mapping.identity.reference_identity_id != proof.reference_identity_id
                or (mapping.identity.contig, mapping.identity.position, mapping.identity.reference, mapping.identity.alternate)
                != (proof.contig, proof.position, proof.reference, proof.alternate)):
            raise ValueError('noncontribution allele/reference binding mismatch')
    maxima = tuple(
        value for value in records
        if value.subset == "joint" and value.stratum == "fafmax95"
    )
    if not maxima:
        return PopulationFilterFacts(
            "no_match", None, None, None, None, "missing",
            None, None, None, None, None, None, None, None, "unknown", (),
            population_manifest_sha256, precision_manifest_sha256,
        )
    try:
        joint = adapt_joint_faf95(records, precision)
        interval = joint.precision_interval
        adjusted = tuple(
            value for value in records
            if value.subset == "joint" and value.stratum == "adjusted"
        )
        joint_ac = adjusted[0].allele_count if len(adjusted) == 1 else None
        joint_an = adjusted[0].allele_number if len(adjusted) == 1 else None
        fields = joint.maximum_record.source.source_fields
        raw_context = fields.get("raw_context") if isinstance(fields, Mapping) else None
        discordance = joint.maximum_record.discordant
        discordance_p = _ratio(
            raw_context.get("DISC_UNION_P")
            if isinstance(raw_context, Mapping) else None
        )
        exomes = _subset(records, "exomes", proof)
        genomes = _subset(records, "genomes")
        if exomes is None or genomes is None:
            comparison = "unknown"
        elif exomes.contribution == genomes.contribution == "observed":
            comparison = "applicable"
        elif "proven_no_callable" in {exomes.contribution, genomes.contribution}:
            comparison = "proved_not_applicable"
        else:
            comparison = "unknown"
        return PopulationFilterFacts(
            "exact_match",
            joint.maximum_raw_token,
            joint.maximum_ancestry,
            None if interval is None else RationalEvidence(
                interval.lower_bound.numerator, interval.lower_bound.denominator,
            ),
            None if interval is None else RationalEvidence(
                interval.upper_bound.numerator, interval.upper_bound.denominator,
            ),
            "qualified" if interval is not None else (
                "missing" if joint.precision_status == "missing_qualification"
                else "unassessable"
            ),
            joint.matching_ancestry_ac,
            joint.matching_ancestry_an,
            joint_ac,
            joint_an,
            exomes,
            genomes,
            discordance,
            discordance_p,
            comparison,
            (),
            population_manifest_sha256,
            precision_manifest_sha256,
        )
    except GnomadFaf95AdapterError:
        raise WgsFailure(
            "population_adaptation", "malformed_evidence",
            operation="population_response_adapt", exception_kind="project_failure",
        ) from None


def _tokens(value: object) -> tuple[str, ...]:
    if type(value) is str:
        values = value.replace("|", ",").split(",")
    elif type(value) in {tuple, list}:
        values = value
    else:
        return ()
    normalized = []
    for item in values:
        if type(item) is not str:
            continue
        token = "_".join(item.strip().lower().replace("-", "_").split())
        if token:
            normalized.append(token)
    return tuple(sorted(set(normalized)))


def _complete_significance(value: object) -> bool:
    """Missing members cannot disappear when assertions are unioned."""
    if type(value) is str:
        return all(part.strip() not in {"", "."} for part in value.replace("|", ",").split(","))
    if type(value) in {tuple, list}:
        return bool(value) and all(_complete_significance(part) for part in value)
    return False


@failure_boundary("clinical_adaptation", "clinical_response_adapt")
def clinical_protection_facts(
    status: str,
    records: tuple[Mapping[str, object], ...],
    *,
    clinical_manifest_sha256: str,
    exception_manifest_sha256: str,
) -> ClinicalProtectionFacts:
    significance = set()
    conflicts = set()
    scopes = set()
    exact = interval = False
    incomplete = False
    for record in records:
        source = record.get("source")
        evidence = record.get("evidence")
        record_type = record.get("record_type")
        if not isinstance(evidence, Mapping):
            incomplete = True
            continue
        if source == "clinvar":
            incomplete |= not _complete_significance(evidence.get("clinical_significance"))
            significance.update(_tokens(evidence.get("clinical_significance")))
            conflicts.update(_tokens(evidence.get("conflicting_significance")))
            if affirmative_source_conflict("ClinVar VCF", {
                "CLNSIG": evidence.get("clinical_significance"),
                "CLNSIGCONF": evidence.get("conflicting_significance"),
                "CLNREVSTAT": evidence.get("review_status"),
            }):
                conflicts.add("conflicting_classifications_of_pathogenicity")
            scopes.add("aggregate")
        elif source == "open_targets":
            incomplete |= not _complete_significance(evidence.get("clinical_significances"))
            significance.update(_tokens(evidence.get("clinical_significances")))
            if affirmative_source_conflict("Open Targets evidence_eva", {
                "confidence": evidence.get("confidence"),
            }):
                conflicts.add("conflicting_classifications_of_pathogenicity")
            scopes.add("individual")
        elif source == "common_allele_exceptions":
            exact = exact or record_type == "exact"
            interval = interval or record_type == "interval"
        for token in tuple(significance):
            if "conflict" in token:
                conflicts.add(token)
    query_status = {
        "matches": "matches",
        "no_match": "no_match",
        "retain_unknown": "unassessable",
        "unavailable": "unavailable",
    }.get(status, "unassessable")
    if incomplete and query_status in {"matches", "no_match"}:
        query_status = "unassessable"
    return ClinicalProtectionFacts(
        query_status,
        tuple(sorted(significance)),
        tuple(sorted(conflicts)),
        exact,
        interval,
        "complete" if query_status in {"matches", "no_match"} else "uncertain",
        CLINICAL_ROLES,
        clinical_manifest_sha256,
        exception_manifest_sha256,
        tuple(sorted(scopes)),
    )


@failure_boundary(
    "source_classification", "source_disposition_classify",
    known_errors=(WgsSourceRoutingError,),
)
def _prefilter_occurrence(
    value,
    renderer: SubmissionContigRenderer,
    contract: NativeSourceContract,
    *, route: SourceOccurrenceRoute,
) -> PrefilterOccurrence:
    if type(route) is not SourceOccurrenceRoute:
        raise WgsSourceRoutingError("source occurrence route is invalid")
    normalized = value.normalization.normalized_allele
    variant_class = value.record.variant_classes[value.source.allele_number - 1]
    phase = set()
    from genosapiens.vcf_qc import phase_reasons_for
    if (value.sample.phase_status is PhaseStatus.PHASED
            and value.sample.phase_evidence is PhaseEvidence.GT):
        phase.add("phased_gt")
    if value.sample.phase_set is not None:
        phase.add("phase_set")
    if value.sample.phase_evidence is PhaseEvidence.PGT_PID:
        phase.add("pgt_pid")
    if phase_reasons_for(value.sample, value.source.allele_number):
        phase.add("phase_" + value.sample.phase_state.value)
        if value.sample.pgt is not None or value.sample.pid is not None:
            phase.add("pgt_pid")
    representation = set()
    if variant_class is VariantClass.MNV:
        representation.add("mnv")
    if variant_class is VariantClass.DELINS:
        representation.add("complex_delins")
    if len(value.record.alternates) > 1:
        representation.add("ambiguous_composite")
    obligation = _technical(value)
    obligations = set(obligation.technical_flags)
    if obligation.source_filtered:
        obligations.add("source_filtered")
    if obligation.read_validation_required:
        obligations.add("read_validation_required")
    return PrefilterOccurrence(
        value.record.source_record_id,
        value.source.source_allele_id,
        value.record.record_number,
        value.source.allele_number,
        route.normalized_allele_id,
        route.preparation_disposition,
        route.contig_class,
        variant_class.value,
        normalized.kind.value
        if route.preparation_disposition == "normalized" else None,
        "exact"
        if route.preparation_disposition == "normalized" else "unavailable",
        tuple(sorted(phase)),
        tuple(sorted(representation)),
        route.native_class_tags,
        tuple(sorted(obligations)),
    )


@failure_boundary("population_request", "population_request_build")
def _population_request(value, *, reference_manifest_sha256: str) -> PopulationRequest:
    normalized = _normalized_sequence(value)
    source = _source_record(value, _technical(value))
    source_mapping = IdentityMapping(
        value.source.source_allele_id,
        IdentityStage.SOURCE,
        _source_sequence(value),
        reference_manifest_sha256,
        source,
    )
    normalized_source = SourcedRecord(
        "genosapiens.wgs-normalization-ledger",
        "1",
        value.normalization.normalized_allele.normalized_allele_id,
        value.payload_sha256,
        None,
        {"source_manifest_sha256": value.parent_source_manifest_sha256},
    )
    normalized_id = value.normalization.normalized_allele.normalized_allele_id
    normalized_mapping = IdentityMapping(
        normalized_id,
        IdentityStage.NORMALIZED,
        normalized,
        reference_manifest_sha256,
        normalized_source,
        (source_mapping.identity_id,),
        "genosapiens.variant-normalization",
    )
    native_id = "native-" + normalized_id
    native_mapping = IdentityMapping(
        native_id,
        IdentityStage.NATIVE,
        normalized,
        reference_manifest_sha256,
        normalized_source,
        (normalized_id,),
        ADMISSION_INVENTORY_SCHEMA,
    )
    allele = AlleleEvidence(
        value.source.source_allele_id,
        native_id,
        (source_mapping, normalized_mapping, native_mapping),
        (_call(value, source),),
    )
    return PopulationRequest(allele, normalized_id)


@failure_boundary("clinical_request", "clinical_request_build")
def _clinical_requests(normalized, renderer):
    return tuple(PublicIndexQuery(
        renderer.accepted_names[item[0].record.contig][0],
        item[0].record.position,
        item[0].record.reference,
        item[0].source.alternate,
        source_span=(item[0].record.position,
                     item[0].record.position + len(item[0].record.reference) - 1),
        normalized_span=(item[0].normalization.normalized_allele.position,
                         item[0].normalization.normalized_allele.position
                         + len(item[0].normalization.normalized_allele.reference) - 1),
        normalized_alias=SequenceAlias(
            renderer.accepted_names[
                item[0].normalization.normalized_allele.contig
            ][0],
            item[0].normalization.normalized_allele.position,
            item[0].normalization.normalized_allele.reference,
            item[0].normalization.normalized_allele.alternate,
        ),
    ) for item in normalized)


@dataclass(frozen=True, slots=True)
class _PreparedQueryBatch:
    ordinal: int
    normalized: tuple[str, ...]
    population_requests: tuple[PopulationRequest, ...]
    clinical_queries: tuple[PublicIndexQuery, ...]


@dataclass(frozen=True, slots=True)
class _CompletedQueryBatch:
    ordinal: int
    population: PopulationAnnotationResult | None
    clinical: tuple
    failure: WgsFailure | None


def _adapt_query_batch(
    batch: _PreparedQueryBatch,
    result: _CompletedQueryBatch,
    pending: tuple,
    policy: CommonAlleleFilterPolicy,
    precision: GnomadFaf95PrecisionQualification,
    source_manifest_sha256: str,
    population_manifest_sha256: str,
    precision_manifest_sha256: str,
    clinical_manifest_sha256: str,
    exception_manifest_sha256: str,
):
    """Pure spawned-worker adaptation; return only static failure details."""
    try:
        unavailable_population = PopulationFilterFacts(
            "unavailable", None, None, None, None, "missing",
            None, None, None, None, None, None, None, None, "unknown", (),
            population_manifest_sha256, precision_manifest_sha256,
        )
        unavailable_clinical = ClinicalProtectionFacts(
            "unavailable", (), (), False, False, "unavailable", CLINICAL_ROLES,
            clinical_manifest_sha256, exception_manifest_sha256,
        )
        population_by_source = {}
        clinical_by_source = {}
        if batch.normalized:
            expected = {
                request.allele.allele_id: request.allele
                for request in batch.population_requests
            }
            responses = {}
            for annotated in result.population.alleles:
                key = annotated.allele_id
                if key in responses or key not in expected or (
                    replace(annotated, population=expected[key].population)
                    != expected[key]
                ):
                    raise WgsFailure(
                        "population_lookup", "malformed_evidence",
                        operation="population_response_bind",
                        exception_kind="project_failure",
                    )
                responses[key] = annotated
            if set(responses) != set(expected):
                raise WgsFailure(
                    "population_lookup", "missing_response",
                    operation="population_response_bind",
                    exception_kind="project_failure",
                )
            for key, protected in zip(
                batch.normalized, result.clinical, strict=True,
            ):
                population_by_source[key] = population_filter_facts(
                    responses[key], precision=precision,
                    population_manifest_sha256=population_manifest_sha256,
                    precision_manifest_sha256=precision_manifest_sha256,
                )
                clinical_by_source[key] = clinical_protection_facts(
                    protected.status, protected.records,
                    clinical_manifest_sha256=clinical_manifest_sha256,
                    exception_manifest_sha256=exception_manifest_sha256,
                )
        decisions = []
        for occurrence, scope_decision in pending:
            with failure_boundary("decision_conversion", "decision_convert"):
                deferred = scope_decision.status in ("outside_scope", "unresolved_scope")
                pure = None if deferred else evaluate_common_allele_filter(
                    occurrence,
                    population_by_source.get(
                        occurrence.source_allele_id, unavailable_population,
                    ),
                    clinical_by_source.get(
                        occurrence.source_allele_id, unavailable_clinical,
                    ),
                    policy=policy,
                )
                decisions.append(ledger_filter_decision(
                    occurrence, pure,
                    parent_source_manifest_sha256=source_manifest_sha256,
                    scope_decision=scope_decision, policy_sha256=hashlib.sha256(common_allele_filter_policy_as_json(policy)).hexdigest(),
                    resource_sha256=(population_manifest_sha256, clinical_manifest_sha256,
                        exception_manifest_sha256, precision_manifest_sha256),
                    population_request_keys=int(occurrence.source_allele_id in population_by_source),
                    clinical_query_keys=int(occurrence.source_allele_id in clinical_by_source),
                ))
        return tuple(decisions), None
    except (Exception, KeyboardInterrupt) as error:
        failure = static_failure(
            error, "resource_session", operation="resource_session_work",
        )
        return (), _failure_fields(failure)


@dataclass(frozen=True, slots=True)
class PublicFilterWorkerConfig:
    """Spawn-safe pinned resource descriptors; no session, factory or open handle."""

    joint: object
    all_sites: object
    reference: object
    runtime_image: str
    limits: object
    public_query_limits: object
    clinical_indexes: object
    staging_root: Path

    def __post_init__(self):
        from genosapiens.population_annotation import (
            GnomadJointReference, GnomadAllSitesReference,
            PopulationReferenceCompatibility, PopulationQueryLimits,
        )
        from genosapiens.bounded_wgs_config import PublicQueryLimits
        from genosapiens.common_allele_public_indexes import VerifiedPublicIndexes
        for value, cls in ((self.joint, GnomadJointReference),
                           (self.all_sites, GnomadAllSitesReference),
                           (self.reference, PopulationReferenceCompatibility),
                           (self.limits, PopulationQueryLimits),
                           (self.public_query_limits, PublicQueryLimits),
                           (self.clinical_indexes, VerifiedPublicIndexes)):
            if type(value) is not cls:
                raise WgsPublicFilterError("public worker descriptor is invalid")
        resource_paths = (self.joint.projection_contract_path, self.joint.publication_contract_path,
            self.joint.validation_contract_path, self.joint.manifest_path, self.joint.release_directory,
            self.all_sites.contract_path, self.all_sites.reference_root, self.reference.identity_path,
            self.reference.reference_root, self.reference.manifest_path,
            self.clinical_indexes.projection_path, self.clinical_indexes.sqlite_path)
        if (any(not isinstance(path, Path) or not path.is_absolute() or ".." in path.parts for path in resource_paths)
                or not isinstance(self.staging_root, Path)
                or not self.staging_root.is_absolute() or ".." in self.staging_root.parts
                or self.public_query_limits.concurrent_batches != 10
                or not 0 < self.public_query_limits.batch_rows <= 1024
                or self.limits.max_alleles != self.public_query_limits.batch_rows
                or self.limits.max_output_bytes != self.public_query_limits.max_response_bytes
                or self.limits.timeout_seconds != self.public_query_limits.timeout_seconds):
            raise WgsPublicFilterError("public worker limits are invalid")

    def session(self):
        return PopulationClinicalResourceSession(joint=self.joint, all_sites=self.all_sites,
            reference=self.reference, runtime_image=self.runtime_image, limits=self.limits,
            public_query_limits=self.public_query_limits, clinical_indexes=self.clinical_indexes)


@dataclass(frozen=True, slots=True)
class _StagedQueryResult:
    ordinal: int
    count: int
    blocks: tuple[Artifact, ...]
    query_nanoseconds: int
    adaptation_nanoseconds: int
    failure_fields: tuple[str, ...] | None
    population_request_keys: int = 0
    clinical_query_keys: int = 0


def _failure_fields(error, operation="resource_session_work"):
    failure = static_failure(error, "resource_session", operation=operation)
    return (failure.origin, failure.category, failure.operation, failure.exception_kind,
        tuple((row.failure_code, row.operation, row.exception_kind) for row in failure.secondary_failures),
        failure.secondary_failures_truncated)


def _worker_failure(fields):
    if type(fields) is not tuple or len(fields) != 6 or type(fields[4]) is not tuple or len(fields[4]) > 4:
        raise WgsPublicFilterError("public worker failure is invalid")
    return WgsFailure(fields[0], fields[1], operation=fields[2], exception_kind=fields[3],
        secondary_failures=tuple(SecondaryFailure(*row) for row in fields[4]),
        secondary_failures_truncated=fields[5])


def _run_query_worker(index, config, science, jobs, results):
    """Own one resource session and one process group until explicit shutdown."""
    os.setsid()

    def interrupt(_signum, _frame):
        raise KeyboardInterrupt

    signal.signal(signal.SIGTERM, interrupt)
    manager = None
    opened = False
    try:
        try:
            manager = config.session()
            session = manager.__enter__()
            opened = True
            results.put(("ready", index, None))
        except BaseException as error:
            results.put(("ready", index, _failure_fields(error, "resource_session_open")))
            return
        while True:
            job = jobs.get()
            if job is None:
                break
            ordinal, count, request_artifact = job
            query_ns = adaptation_ns = 0
            population_count = clinical_count = 0
            blocks = ()
            failure = None
            cancelled = False
            try:
                from genosapiens.wgs_ledger import _directory, _read
                if type(request_artifact) is not Artifact or request_artifact.size_bytes > 64 * 1024**2:
                    raise WgsPublicFilterError("public worker request exceeds byte bound")
                with _directory(config.staging_root) as fd, _read(fd, request_artifact.name) as (stream, metadata):
                    if metadata.st_size != request_artifact.size_bytes:
                        raise WgsPublicFilterError("public worker request size changed")
                    content = stream.read(request_artifact.size_bytes + 1)
                if len(content) != request_artifact.size_bytes or hashlib.sha256(content).hexdigest() != request_artifact.sha256:
                    raise WgsPublicFilterError("public worker request digest changed")
                records = []
                offset = 0
                while offset < len(content):
                    if offset + 4 > len(content):
                        raise WgsPublicFilterError("public worker request is truncated")
                    size = int.from_bytes(content[offset:offset + 4], "big")
                    offset += 4
                    if size <= 0 or offset + size > len(content) or len(records) >= config.public_query_limits.batch_rows:
                        raise WgsPublicFilterError("public worker request exceeds bounds")
                    records.append(pickle.loads(content[offset:offset + size]))
                    offset += size
                if len(records) != count:
                    raise WgsPublicFilterError("public worker request count changed")
                del content
                if any(type(row) is not tuple or len(row) != 4
                       or type(row[3]) is not EarlyScopeDecision
                       or (row[1] is None) != (row[2] is None)
                       or row[3].status != "in_scope" and row[1] is not None for row in records):
                    raise WgsPublicFilterError("public worker scope transport is invalid")
                pending = tuple((row[0], row[3]) for row in records)
                normalized = tuple(row for row in records if row[1] is not None)
                batch = _PreparedQueryBatch(ordinal, tuple(row[0].source_allele_id for row in normalized),
                    tuple(row[1] for row in normalized), tuple(row[2] for row in normalized))
                started = time.monotonic_ns()
                population, clinical = None, ()
                if normalized:
                    with failure_boundary("population_lookup", "population_query"):
                        population = session.query_population(batch.population_requests)
                        population_count = len(batch.population_requests)
                    with failure_boundary("clinical_index", "clinical_index_query"):
                        clinical = session.query_clinical(batch.clinical_queries)
                        clinical_count = len(batch.clinical_queries)
                query_ns = time.monotonic_ns() - started
                started = time.monotonic_ns()
                decisions, failure = _adapt_query_batch(batch,
                    _CompletedQueryBatch(ordinal, population, clinical, None), pending, *science)
                adaptation_ns = time.monotonic_ns() - started
                if failure is None:
                    artifacts = []
                    block, size = [], 0

                    def flush():
                        name = f"result-{ordinal:09d}-{len(artifacts):04d}.jsonl"
                        artifacts.append(publish_stream(config.staging_root, name, iter(block)))

                    for decision in decisions:
                        encoded = decision.as_bytes()
                        if len(encoded) > MAX_ROW_BYTES:
                            raise WgsPublicFilterError("public worker decision exceeds row bound")
                        if block and size + len(encoded) > MAX_SHARD_BYTES:
                            flush()
                            block, size = [], 0
                        block.append(encoded)
                        size += len(encoded)
                    if block:
                        flush()
                    blocks = tuple(artifacts)
                del records, pending, normalized, batch, population, clinical, decisions
            except BaseException as error:
                failure = _failure_fields(error)
                cancelled = isinstance(error, KeyboardInterrupt)
            results.put(("result", index, _StagedQueryResult(
                ordinal, count, blocks, query_ns, adaptation_ns, failure,
                population_count, clinical_count)))
            if cancelled:
                break
    except KeyboardInterrupt:
        pass
    finally:
        cleanup = None
        if opened:
            try:
                manager.__exit__(None, None, None)
            except BaseException as error:
                cleanup = _failure_fields(error, "resource_session_close")
        results.put(("closed", index, cleanup))


class _ProcessQueryPool:
    """Two batches per worker; all waits check child health and deadlines."""

    def __init__(self, config, science):
        self.config = config
        self.science = science
        context = multiprocessing.get_context("spawn")
        self.worker_count = config.public_query_limits.concurrent_batches
        self.max_outstanding = 2 * self.worker_count
        self.results = context.Queue(maxsize=self.max_outstanding)
        self.jobs = tuple(context.Queue(maxsize=2) for _ in range(self.worker_count))
        self.workers = tuple(context.Process(target=_run_query_worker,
            args=(index, config, science, self.jobs[index], self.results),
            name=f"wgs-public-query-{index}") for index in range(self.worker_count))
        self.outstanding = {}
        self.completed = {}
        self.openings = {}
        self.closings = {}
        self.closed = False
        self.subreaper_previous = None

    def _receive(self, deadline):
        while True:
            try:
                message = self.results.get(timeout=min(.1, max(.001, deadline - time.monotonic())))
                if (type(message) is not tuple or len(message) != 3
                        or type(message[1]) is not int or not 0 <= message[1] < self.worker_count):
                    raise WgsPublicFilterError("public worker result is malformed")
                kind, index, value = message
                if kind == "ready" and index not in self.openings:
                    self.openings[index] = value
                elif kind == "closed" and index not in self.closings:
                    self.closings[index] = value
                elif kind == "result" and type(value) is _StagedQueryResult:
                    if (value.ordinal not in self.outstanding or value.ordinal in self.completed
                            or value.ordinal % self.worker_count != index
                            or value.count != self.outstanding[value.ordinal][0]
                            or type(value.blocks) is not tuple or len(value.blocks) > value.count
                            or any(type(n) is not int or n < 0 for n in (
                                value.query_nanoseconds, value.adaptation_nanoseconds,
                                value.population_request_keys, value.clinical_query_keys))):
                        raise WgsPublicFilterError("public worker result is malformed")
                    self.completed[value.ordinal] = value
                else:
                    raise WgsPublicFilterError("public worker result is unexpected")
                return
            except queue.Empty:
                if any(worker.exitcode is not None for worker in self.workers):
                    raise WgsFailure("resource_session", "missing_response",
                        operation="resource_session_work", exception_kind="project_failure")
                if time.monotonic() >= deadline:
                    raise WgsFailure("resource_session", "timeout",
                        operation="resource_session_work", exception_kind="project_failure")

    def open(self):
        # Adopt only descendants orphaned by these workers, and reap by group.
        libc = ctypes.CDLL(None, use_errno=True)
        previous = ctypes.c_int()
        if libc.prctl(37, ctypes.byref(previous), 0, 0, 0) or libc.prctl(36, 1, 0, 0, 0):
            raise WgsPublicFilterError("public worker descendant ownership unavailable")
        self.subreaper_previous = previous.value
        for worker in self.workers:
            worker.start()
        deadline = time.monotonic() + self.config.public_query_limits.timeout_seconds
        while len(self.openings) != self.worker_count:
            self._receive(deadline)
        failures = [_worker_failure(value) for _, value in sorted(self.openings.items()) if value is not None]
        if failures:
            for secondary in failures[1:]:
                failures[0].add_secondary(secondary)
            raise failures[0]

    def submit(self, ordinal, expected, content):
        if (self.closed or ordinal in self.outstanding or len(self.outstanding) >= self.max_outstanding
                or not 0 < len(expected) <= self.config.public_query_limits.batch_rows
                or not 0 < len(content) <= 64 * 1024**2):
            raise WgsPublicFilterError("public worker submission exceeds bounds")
        artifact = publish_stream(self.config.staging_root, f"request-{ordinal:09d}.bin", (content,), maximum_bytes=64 * 1024**2)
        self.outstanding[ordinal] = (len(expected), artifact,
            time.monotonic() + 3 * self.config.public_query_limits.timeout_seconds)
        self.jobs[ordinal % self.worker_count].put((ordinal, len(expected), artifact), timeout=1)

    def take(self, ordinal):
        while ordinal not in self.completed:
            self._receive(min(item[2] for item in self.outstanding.values()))
        result = self.completed[ordinal]
        if result.failure_fields is not None:
            failures = [_worker_failure(value.failure_fields)
                        for _, value in sorted(self.completed.items()) if value.failure_fields is not None]
            for secondary in failures[1:]:
                failures[0].add_secondary(secondary)
            raise failures[0]
        return result

    def release(self, ordinal):
        result = self.completed.pop(ordinal)
        _, request, _ = self.outstanding.pop(ordinal)
        for artifact in (request, *result.blocks):
            (self.config.staging_root / artifact.name).unlink()

    def close(self):
        if self.closed:
            return ()
        self.closed = True
        interrupted = bool(self.outstanding)
        for index, worker in enumerate(self.workers):
            if worker.pid is None:
                continue
            if interrupted:
                try:
                    os.killpg(worker.pid, signal.SIGTERM)
                except ProcessLookupError:
                    if worker.is_alive():
                        worker.terminate()
            else:
                try:
                    self.jobs[index].put(None, timeout=.1)
                except queue.Full:
                    worker.terminate()
        deadline = time.monotonic() + 30
        while any(worker.is_alive() for worker in self.workers) and time.monotonic() < deadline:
            try:
                kind, index, value = self.results.get(timeout=.05)
                if kind == "closed":
                    self.closings[index] = value
            except queue.Empty:
                pass
        for worker in self.workers:
            if worker.pid is None:
                continue
            # Include query descendants even if the worker itself died first.
            try:
                os.killpg(worker.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            if worker.is_alive():
                worker.kill()
            worker.join(timeout=1)
        for worker in self.workers:
            if worker.pid is None:
                continue
            while True:
                try:
                    child, _ = os.waitpid(-worker.pid, os.WNOHANG)
                except ChildProcessError:
                    break
                if child == 0:
                    if time.monotonic() >= deadline:
                        raise WgsPublicFilterError("public worker descendants did not exit")
                    time.sleep(.01)
        if self.subreaper_previous is not None:
            ctypes.CDLL(None).prctl(36, self.subreaper_previous, 0, 0, 0)
        while True:
            try:
                kind, index, value = self.results.get_nowait()
                if kind == "closed":
                    self.closings[index] = value
            except queue.Empty:
                break
        for channel in (*self.jobs, self.results):
            channel.cancel_join_thread()
            channel.close()
        return tuple(_worker_failure(value) for _, value in sorted(self.closings.items()) if value is not None)


@dataclass(frozen=True, slots=True)
class PublicEarlyFilterDecisionProvider:
    policy: CommonAlleleFilterPolicy
    precision: GnomadFaf95PrecisionQualification
    source_routing_context: SourceRoutingContext
    worker_config: PublicFilterWorkerConfig
    scope_session: EarlyScopeSession
    phase_recorder: RuntimePhaseRecorder | None = None

    def __post_init__(self):
        if (type(self.policy) is not CommonAlleleFilterPolicy
                or type(self.precision) is not GnomadFaf95PrecisionQualification
                or type(self.source_routing_context) is not SourceRoutingContext
                or type(self.worker_config) is not PublicFilterWorkerConfig
                or type(self.scope_session) is not EarlyScopeSession
                or self.phase_recorder is not None and type(self.phase_recorder) is not RuntimePhaseRecorder):
            raise WgsPublicFilterError("public early-filter provider is invalid")

    @property
    def native_source_contract(self):
        return self.source_routing_context.contract

    @property
    def contig_renderer(self):
        return self.source_routing_context.renderer

    def __call__(self, context):
        if type(context) is not PreparationFilterContext or context.source_cursor is None:
            raise WgsPublicFilterError("public early-filter source cursor is required")
        return SourceBoundDecisionStream(context.source_cursor, self._decisions(context), completed_batches=True)

    def _decisions(self, context):
        if (context.source_cursor.routing != self.source_routing_context
                or self.native_source_contract.parent_source_manifest_sha256 != context.source_manifest.sha256
                or self.scope_session.policy_sha256 != context.scope_policy_sha256
                or self.scope_session.manifest.sha256 != context.scope_manifest_sha256
                or self.scope_session.reference_sha256 != context.source_manifest.source.reference_manifest_sha256
                or self.scope_session.renderer != self.contig_renderer):
            raise WgsPublicFilterError("public early-filter context is invalid")
        from genosapiens.wgs_admission import _source_header
        header = _source_header(context.source_root, context.source_manifest)
        native_context = dict(
            definitions={(ns, name): (number, kind) for ns, name, number, kind in header.class_field_definitions},
            vcf_version=next((line.split("=", 1)[1].strip() for line in header.raw_lines
                if line.startswith("##fileformat=")), None),
            repeat_identity_format=self.native_source_contract.repeat_identity_format,
            header_lines=header.raw_lines)
        self.worker_config.staging_root.mkdir(mode=0o700, parents=True, exist_ok=True)
        from genosapiens.wgs_ledger import _directory
        with _directory(self.worker_config.staging_root):
            pass
        with tempfile.TemporaryDirectory(prefix="query-", dir=self.worker_config.staging_root) as temporary:
            config = replace(self.worker_config, staging_root=Path(temporary))
            pool = _ProcessQueryPool(config, (self.policy, self.precision, context.source_manifest.sha256,
                context.population_manifest_sha256, context.precision_manifest_sha256,
                context.clinical_manifest_sha256, context.exception_manifest_sha256))
            held, pending, encoded = {}, [], bytearray()
            submitted = consumed = 0
            primary = None
            generator_closed = False

            def submit():
                nonlocal submitted, pending, encoded
                held[submitted] = tuple(pending)
                pool.submit(submitted, held[submitted], bytes(encoded))
                submitted += 1
                pending, encoded = [], bytearray()

            def consume():
                nonlocal consumed
                result = pool.take(consumed)
                expected = held[consumed]
                seen = 0
                completed_decisions = []
                counts = dict(population_request_keys=0, clinical_query_keys=0)
                for index, artifact in enumerate(result.blocks):
                    if artifact.name != f"result-{consumed:09d}-{index:04d}.jsonl":
                        raise WgsPublicFilterError("public worker result path is invalid")
                    for line in iter_artifact_lines(config.staging_root, artifact, maximum_line_bytes=MAX_ROW_BYTES):
                        decision = load_model(line, FilterDecision)
                        if seen >= len(expected):
                            raise WgsPublicFilterError("public worker returned extra decisions")
                        occurrence, expected_scope = expected[seen]
                        evidence = json.loads(decision.evidence_json)
                        if evidence["scope_decision"] != json.loads(json.dumps(asdict(expected_scope))):
                            raise WgsPublicFilterError("public worker scope decision changed")
                        for key in counts:
                            counts[key] += evidence["request_counts"][key]
                        if (decision.record_number, decision.allele_number, decision.source_allele_id,
                            decision.preparation_disposition, decision.normalized_allele_id) != (
                                occurrence.record_number, occurrence.allele_number, occurrence.source_allele_id,
                                occurrence.preparation_disposition, occurrence.normalized_allele_id):
                            raise WgsFailure("filter_publication", "malformed_evidence",
                                operation="filter_source_tuple_compare", exception_kind="project_failure")
                        if (decision.parent_source_manifest_sha256 != context.source_manifest.sha256
                                or decision.policy_sha256 != context.policy_sha256):
                            raise WgsPublicFilterError("public worker decision parents changed")
                        seen += 1
                        completed_decisions.append(decision)
                if counts != dict(population_request_keys=result.population_request_keys,
                                  clinical_query_keys=result.clinical_query_keys):
                    raise WgsPublicFilterError("public worker request counts changed")
                if seen != len(expected):
                    raise WgsPublicFilterError("public worker omitted decisions")
                if self.phase_recorder is not None:
                    self.phase_recorder.record_worker_duration("prefilter_adaptation", result.adaptation_nanoseconds)
                    self.phase_recorder.record_worker_duration("prefilter_worker_query", result.query_nanoseconds)
                from genosapiens.wgs_filter import CompletedFilterBatch
                yield CompletedFilterBatch(tuple(completed_decisions),
                    result.population_request_keys, result.clinical_query_keys)
                pool.release(consumed)
                del held[consumed]
                consumed += 1

            try:
                with runtime_phase(self.phase_recorder, "prefilter_session_authentication"):
                    pool.open()
                with runtime_phase(self.phase_recorder, "prefilter_query"):
                    for value, route in context.source_cursor:
                        occurrence = _prefilter_occurrence(value, self.contig_renderer, self.native_source_contract, route=route)
                        scope_decision = self.scope_session.classify(value, route, structural_context=native_context)
                        request = query = None
                        if scope_decision.status == "in_scope" and not prefilter_bypass_reasons(occurrence):
                            request = _population_request(value, reference_manifest_sha256=context.source_manifest.source.reference_manifest_sha256)
                            query, = _clinical_requests(((value, occurrence),), self.contig_renderer)
                        item = pickle.dumps((occurrence, request, query, scope_decision), protocol=5)
                        if len(item) + 4 > 64 * 1024**2:
                            raise WgsPublicFilterError("public worker request exceeds byte bound")
                        if pending and len(encoded) + len(item) + 4 > 64 * 1024**2:
                            submit()
                            if len(pool.outstanding) == pool.max_outstanding:
                                yield from consume()
                        pending.append((occurrence, scope_decision))
                        encoded.extend(len(item).to_bytes(4, "big"))
                        encoded.extend(item)
                        if len(pending) == config.public_query_limits.batch_rows:
                            submit()
                            if len(pool.outstanding) == pool.max_outstanding:
                                yield from consume()
                    if pending:
                        submit()
                    while pool.outstanding:
                        yield from consume()
            except GeneratorExit:
                generator_closed = True
                raise
            except BaseException as error:
                primary = static_failure(error, "resource_session", operation="resource_session_work")
                raise primary from None
            finally:
                try:
                    failures = pool.close()
                except BaseException as error:
                    failures = (static_failure(error, "resource_session", operation="resource_session_close"),)
                if failures:
                    if primary is not None:
                        for failure in failures:
                            primary.add_secondary(failure)
                    elif not generator_closed:
                        for failure in failures[1:]:
                            failures[0].add_secondary(failure)
                        raise failures[0]
