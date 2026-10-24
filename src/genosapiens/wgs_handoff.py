"""Current Track 1 handoff over the complete normalized scientific roots."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from contextlib import ExitStack, closing, nullcontext
import hashlib
import json
import os
from pathlib import Path
import stat
from typing import Literal

from genosapiens.evidence_model import Observation, ScoreLedger, _sha, _exact, _typed
from genosapiens.phenotype_scoring import PhenotypeScore
from genosapiens.track1_preflight import SelectedTechnicalAssessment
from genosapiens.track1_preflight import CompetitiveRunEvidence
from genosapiens.admission_links import PinnedPublicSource
from genosapiens.wgs_ledger import LedgerManifest, FilterManifest, WorkingSetManifest
from genosapiens.wgs_provenance import execution_sha256
from genosapiens.wgs_scientific import _bytes, _load_manifest
from genosapiens.wgs_scientific_artifact import (
    ScientificArtifact, ScientificArtifactReader, scientific_artifact_from_bytes,
)
from genosapiens.wgs_scientific_products import (
    ScientificOrderProduct, ScientificPrimaryProduct,
    ScientificSelectionProduct, scientific_product_from_bytes, product_descriptors,
    open_scientific_products, validate_scientific_products, scientific_report_roots,
)
from genosapiens.wgs_evidence_store import EvidenceReference as RootReference
from genosapiens.wgs_case_resources import bound_reference
from genosapiens.wgs_logical_evidence import ScientificHypothesisDetail, hypothesis_detail, iter_logical_provenance
from genosapiens.wgs_failures import WgsFailure


HANDOFF_SCHEMA = "genosapiens.track1-sharded-mechanism-ready-handoff"
MODES = ("diagnosis_blind_primary",)


class WgsHandoffError(WgsFailure):
    """Static failure of the shared current handoff."""

    def __init__(self, message, *, category="malformed_evidence"):
        super().__init__(
            "final_publication", category,
            operation="final_publish", exception_kind="project_failure",
        )
        RuntimeError.__init__(self, message)


@dataclass(frozen=True, slots=True, repr=False)
class HandoffHypothesis:
    analysis_mode: Literal["diagnosis_blind_primary"]
    scientific: ScientificHypothesisDetail
    mode_artifact_sha256: str
    rank: int
    selected_for_submission: bool
    ledger: ScoreLedger
    phenotype_score: PhenotypeScore
    ordering_total: str
    prior_contribution: str
    prior_reason: str | None
    prior_proofs: tuple[Mapping[str, object], ...]
    technical_assessments: tuple[SelectedTechnicalAssessment, ...]

    def __post_init__(self):
        _sha(self.mode_artifact_sha256)
        if (self.analysis_mode not in MODES or type(self.rank) is not int or self.rank < 1
                or self.ledger.hypothesis_id != self.scientific.ledger.hypothesis_id
                or self.selected_for_submission and self.analysis_mode != MODES[0]
                or _exact(self.ordering_total) != _exact(self.ledger.total) + _exact(self.prior_contribution)
                or (
                    self.prior_proofs or self.prior_reason is not None
                    or _exact(self.prior_contribution) != 0)):
            raise WgsHandoffError("handoff scientific score binding is invalid")

    @property
    def hypothesis_id(self):
        return self.analysis_mode + ":" + self.scientific.ledger.hypothesis_id


@dataclass(frozen=True, slots=True, repr=False)
class HandoffRecord:
    kind: Literal["observation", "hypothesis"]
    record_id: str
    observation: Observation | None
    hypothesis: HandoffHypothesis | None

    def __post_init__(self):
        if ((self.kind, self.observation is not None, self.hypothesis is not None)
                not in (("observation", True, False), ("hypothesis", False, True))
                or self.record_id != (self.observation.scoring.observation_id
                    if self.kind == "observation" else self.hypothesis.hypothesis_id)):
            raise WgsHandoffError("handoff record identity is invalid")

    @property
    def key(self):
        return self.kind + ":" + self.record_id


@dataclass(frozen=True, slots=True, repr=False)
class HandoffArtifact:
    role: str
    relative_path: str
    sha256: str
    size_bytes: int

    def __post_init__(self):
        path = Path(self.relative_path)
        _sha(self.sha256)
        if (not self.role or self.relative_path != path.as_posix() or path.is_absolute()
                or not path.parts or any(part in ("", ".", "..") for part in path.parts)
                or type(self.size_bytes) is not int or self.size_bytes < 1):
            raise WgsHandoffError("handoff artifact binding is invalid")


ARTIFACT_ROLES = frozenset((
    "prepared", "early", "source", "source_accounting", "source_audit",
    "preparation_run_options", "prepared_index", "completed_source_input", "prefilter", "working",
    "evidence", "scientific", "order",
    "primary", "selection", "execution", "release",
    "run_evidence", "preflight", "case", "concordance", "output", "epcr", "csv", "analysis_profile",
))


@dataclass(frozen=True, slots=True, repr=False)
class HandoffScientificProducts:
    scientific: ScientificArtifact
    scientific_directory: Path
    order: ScientificOrderProduct
    order_directory: Path
    primary: ScientificPrimaryProduct
    selection: ScientificSelectionProduct
    selection_directory: Path
    eligibility_root: Path | None = None
    profile_root: Path | None = None


@dataclass(frozen=True, slots=True, repr=False)
class HandoffPublicationAuthority:
    source: LedgerManifest
    prefilter: FilterManifest
    working: WorkingSetManifest
    case_manifest_path: Path
    concordance_content: bytes
    run_evidence: CompetitiveRunEvidence
    public_reference: PinnedPublicSource
    schema_path: Path
    contract_path: Path
    application_version: str
    evaluator_source_path: Path | None = None

    def __post_init__(self):
        if (type(self.source) is not LedgerManifest or type(self.prefilter) is not FilterManifest
                or type(self.working) is not WorkingSetManifest
                or type(self.run_evidence) is not CompetitiveRunEvidence
                or type(self.public_reference) is not PinnedPublicSource
                or not self.application_version or not self.concordance_content
                or any(not path.is_absolute() for path in (
                    self.case_manifest_path, self.schema_path, self.contract_path))):
            raise WgsHandoffError("handoff publication authority is invalid")


def read_handoff_artifact(root, artifact, *, maximum_bytes=8 * 1024 * 1024):
    root = Path(root)
    path = root / artifact.relative_path
    if (not root.is_absolute() or artifact.size_bytes > maximum_bytes
            or any(part.is_symlink() for part in (path, *path.parents))
            or not path.resolve(strict=True).is_relative_to(root.resolve(strict=True))):
        raise WgsHandoffError("handoff artifact path or capacity is invalid")
    descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    with os.fdopen(descriptor, "rb") as handle:
        before = os.fstat(handle.fileno())
        if not stat.S_ISREG(before.st_mode) or before.st_size != artifact.size_bytes:
            raise WgsHandoffError("handoff artifact file or size is invalid")
        content = handle.read(artifact.size_bytes + 1)
        after = os.fstat(handle.fileno())
        if (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns, before.st_ctime_ns) != (
                after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns, after.st_ctime_ns):
            raise WgsHandoffError("handoff artifact changed during read")
    if len(content) != artifact.size_bytes or hashlib.sha256(content).hexdigest() != artifact.sha256:
        raise WgsHandoffError("handoff artifact binding changed")
    return content


def load_handoff_science(root, manifest):
    artifacts = {row.role: row for row in manifest.artifacts}
    scientific = scientific_artifact_from_bytes(read_handoff_artifact(root, artifacts["scientific"]),
        expected_sha256=artifacts["scientific"].sha256)
    values = {}
    for role, cls in (("order", ScientificOrderProduct), ("primary", ScientificPrimaryProduct),
            ("selection", ScientificSelectionProduct)):
        values[role] = scientific_product_from_bytes(read_handoff_artifact(root, artifacts[role]),
            expected_sha256=artifacts[role].sha256, product_type=cls)
        if values[role].scientific_manifest_sha256 != scientific.sha256:
            raise WgsHandoffError("handoff scientific parent changed")
    directory = lambda role: (Path(root) / artifacts[role].relative_path).parent
    return HandoffScientificProducts(scientific, directory("scientific"), values["order"], directory("order"),
        values["primary"],
        values["selection"], directory("selection"), profile_root=directory("analysis_profile"))


def check_scientific_shard_limits(products, config):
    for shard in (*products.scientific.evidence.node_shards, *products.scientific.evidence.root_shards):
        if (shard.size_bytes > config.source_shards.max_uncompressed_bytes
                or shard.row_count > config.source_shards.max_uncompressed_rows):
            raise WgsHandoffError("scientific shard exceeds bound execution limits", category="response_capacity")


ROOT_ROLES = ("case", "primary", "primary_order", "reference", "resources", "selected", "selection")


@dataclass(frozen=True, slots=True, repr=False)
class ShardedHandoffManifest:
    producer_version: str
    schema_sha256: str
    contract_sha256: str
    canonical_artifact_path: str
    run: Mapping[str, object]
    artifacts: tuple[HandoffArtifact, ...]
    scientific_roots: tuple[tuple[str, RootReference], ...]
    observation_count: int
    hypothesis_count_by_mode: tuple[tuple[str, int], ...]
    schema: Literal["genosapiens.track1-sharded-mechanism-ready-handoff"] = HANDOFF_SCHEMA
    version: Literal[6] = 6
    frozen: Literal[True] = True

    def __post_init__(self):
        for value in (self.schema_sha256, self.contract_sha256):
            _sha(value)
        evidence = next((row for row in self.artifacts if row.role == "evidence"), None)
        profile = next((row for row in self.artifacts if row.role == "analysis_profile"), None)
        if (not self.producer_version or self.schema != HANDOFF_SCHEMA or self.version != 6
                or self.frozen is not True
                or self.canonical_artifact_path != "controlled-work/track1/handoff/mechanism-ready-handoff.json"
                or self.artifacts != tuple(sorted(self.artifacts, key=lambda row: row.role))
                or {row.role for row in self.artifacts} not in (ARTIFACT_ROLES, ARTIFACT_ROLES | {"read_intake"})
                or profile is None or profile.relative_path != "analysis/profile/manifest.json"
                or len(self.artifacts) != len({row.role for row in self.artifacts})
                or len({row.relative_path for row in self.artifacts}) != len(self.artifacts)
                or tuple(key for key, _ in self.scientific_roots) != ROOT_ROLES
                or any(type(ref) is not RootReference or evidence is None
                       or ref.manifest_sha256 != evidence.sha256 for _, ref in self.scientific_roots)
                or type(self.observation_count) is not int or self.observation_count < 0
                or tuple(mode for mode, _ in self.hypothesis_count_by_mode) != MODES
                or any(type(count) is not int or count < 1 for _, count in self.hypothesis_count_by_mode)
                or len({count for _, count in self.hypothesis_count_by_mode}) != 1):
            raise WgsHandoffError("normalized handoff manifest is invalid")

    @property
    def handoff_id(self):
        return "track1-handoff:" + self.sha256

    def as_bytes(self):
        content = _bytes({"classification": "controlled", "llm_policy": "no-llm", "manifest": self})
        if len(content) > 8 * 1024 * 1024:
            raise WgsHandoffError("handoff manifest exceeds capacity", category="response_capacity")
        return content

    @property
    def sha256(self):
        return hashlib.sha256(self.as_bytes()).hexdigest()


def handoff_manifest_from_json(content, *, expected_sha256):
    return _load_manifest(content, expected_sha256=expected_sha256,
        result_type=ShardedHandoffManifest, result_key="manifest")


def validate_handoff_roots(manifest, products, reader):
    if type(reader) is not ScientificArtifactReader or reader.artifact != products.scientific:
        raise WgsHandoffError("handoff requires its authenticated scientific owner")
    with reader.reader.read_scope():
        descriptors = product_descriptors(products)
        validate_scientific_products(reader, descriptors,
            expected_prefilter_manifest_sha256=products.primary.prefilter_manifest_sha256)
        if (manifest.run.get('retention_scope') != retention_scope(reader.resources.analysis_profile)
                or manifest.scientific_roots != scientific_report_roots(reader, descriptors)
                or manifest.observation_count != len(reader.resources.phenotypes.observations)
                or manifest.hypothesis_count_by_mode != tuple((mode, reader.artifact.hypothesis_count) for mode in MODES)):
            raise WgsHandoffError("handoff complete scientific closure changed")


def retention_scope(profile):
    """Describe the authenticated selection universe without a partner inference."""
    from genosapiens.variant_retention import POLICY_SHA256
    if profile is None:
        raise WgsHandoffError('handoff requires the retained universe')
    return dict(policy_sha256=POLICY_SHA256, profile_manifest_sha256=profile['manifest_sha256'],
        evaluated_universe='retained_alleles_only',
        identities=profile['counts']['retention_dispositions'],
        source_occurrences=profile['counts']['retention_occurrences'],
        deferred_inventory='ordered_profile_decision_shards',
        deferred_downstream_work=False, outside_partner_absence_established=False)


def write_handoff_manifest(products, *, reader, directory, execution_config, producer_version,
                           schema_sha256, contract_sha256, canonical_artifact_path, run, artifacts):
    from genosapiens.wgs_ledger import publish_bytes
    if type(reader) is not ScientificArtifactReader or reader.artifact != products.scientific:
        raise WgsHandoffError("handoff requires its authenticated scientific owner")
    if dict(reader.parents).get("execution_config") != execution_sha256(execution_config):
        raise WgsHandoffError("handoff execution limits changed")
    check_scientific_shard_limits(products, execution_config)
    with reader.reader.read_scope():
        manifest = ShardedHandoffManifest(producer_version, schema_sha256, contract_sha256,
            canonical_artifact_path, {**run, "retention_scope": retention_scope(reader.resources.analysis_profile),
                **({"scientific_scope": reader.scientific_scope.header} if hasattr(reader, "scientific_scope") else {})},
            tuple(sorted(artifacts, key=lambda row: row.role)),
            scientific_report_roots(reader, product_descriptors(products)),
            len(reader.resources.phenotypes.observations),
            tuple((mode, reader.artifact.hypothesis_count) for mode in MODES))
        validate_handoff_roots(manifest, products, reader)
        Path(directory).mkdir(mode=0o700, parents=True, exist_ok=True)
        publish_bytes(Path(directory), "mechanism-ready-handoff.json", manifest.as_bytes())
        return manifest


def project_handoff_records(products, *, reader, run_evidence):
    """Explicit export of complete primary records, restoring each context once."""
    if type(reader) is not ScientificArtifactReader or reader.artifact != products.scientific:
        raise WgsHandoffError("handoff requires its authenticated scientific owner")
    selected = {row.primary_hypothesis_id: row for row in products.selection.selected}
    primary_sha256 = products.primary.sha256
    with reader.reader.read_scope():
        validate_scientific_products(reader, product_descriptors(products),
            expected_prefilter_manifest_sha256=products.primary.prefilter_manifest_sha256)
        for observation in reader.resources.phenotypes.observations:
            yield HandoffRecord("observation", observation.scoring.observation_id, observation, None)
        count = 0
        case = reader.case
        for row in case.contexts():
            with case._restore(reader.reader, row, reader.reference, reader.reference_pin, reader.renderers,
                    retain_identity_index=True) as context:
                positions = None if context.universe.pairs is None else context.universe.pairs.domain.positions()
                with closing(context.order.identity_locators()) as locators:
                    for hypothesis_id, locator in locators:
                        primary = case.logical_record(context, locator, context.scores, positions=positions)
                        scientific = hypothesis_detail(reader, primary)
                        hypothesis = scientific.constructed.hypothesis
                        if hypothesis.hypothesis_id != hypothesis_id:
                            raise WgsHandoffError("handoff scientific identity changed")
                        group = selected.get(hypothesis_id)
                        if group is not None and group.partial_member != (hypothesis.completeness.value == "incomplete"):
                            raise WgsHandoffError("handoff selected genotype completeness changed")
                        technical = tuple(row for row in run_evidence.technical_assessments
                                          if row.allele_id in hypothesis.evidence_subject_ids)
                        ledger = scientific.ledger
                        ties = tuple(item.value for item in ledger.tie_keys[:-1])
                        rank = case.rank_of((context.ordinal, *locator), total=ledger.total, ties=ties)
                        projected = HandoffHypothesis(MODES[0], scientific, primary_sha256, rank, group is not None,
                            ledger, scientific.evidence.phenotype_score, ledger.total, "0", None, (), technical)
                        record = HandoffRecord("hypothesis", projected.hypothesis_id, None, projected)
                        if len(_bytes(record)) > 256 * 1024**2:
                            raise WgsHandoffError("logical handoff record exceeds capacity", category="response_capacity")
                        yield record
                        count += 1
        if count != products.scientific.hypothesis_count:
            raise WgsHandoffError("handoff omits sealed scientific hypotheses")


@dataclass(frozen=True, slots=True, repr=False)
class VerifiedHandoffBundle:
    manifest: ShardedHandoffManifest
    root: Path
    directory: Path
    execution_config: object
    products: HandoffScientificProducts
    evidence: object
    prepared: object
    source: object
    release: object
    run_evidence: object
    case: object
    scratch_directory: Path
    read_intake: object | None = None

    @property
    def source_inputs(self):
        from genosapiens.case_manifest import CaseFileClass
        original = {row.sha256 for row in self.case.inputs if row.file_class is CaseFileClass.VCF}
        if self.read_intake is not None:
            original.update(stage.output_artifact.sha256 for delta in self.read_intake.manifest.deltas
                            for stage in delta.stage_results)
        return tuple(sorted(original))

    def open_science(self):
        return open_scientific_products(self.products.scientific, product_descriptors(self.products),
            directory=self.products.scientific_directory,
            expected_parent_bindings=dict(self.products.scientific.parent_bindings),
            scratch_directory=self.scratch_directory,eligibility_root=self.products.eligibility_root,
            profile_root=self.products.profile_root)

    def records(self):
        with self.open_science() as reader:
            validate_handoff_roots(self.manifest, self.products, reader)
            yield from project_handoff_records(self.products, reader=reader, run_evidence=self.run_evidence)

    def provenance(self):
        with self.open_science() as reader:
            validate_handoff_roots(self.manifest, self.products, reader)
            yield from iter_logical_provenance(reader)


def validate_handoff_allele(allele, *, source_inputs, reference_identity_id,
                            reference_manifest_sha256):
    """Retain source/native lineage checks for every persisted allele."""
    from genosapiens.evidence_model import IdentityStage
    from genosapiens.allele_identity import SequenceAllele
    mappings = {row.identity_id: row for row in allele.identities}
    sources = {key: row for key, row in mappings.items() if row.stage is IdentityStage.SOURCE}
    native = mappings[allele.native_identity_id]
    if (not sources or not allele.calls
            or sum(row.stage is IdentityStage.NATIVE for row in mappings.values()) != 1
            or any(row.derived_from_identity_ids for row in sources.values())):
        raise WgsHandoffError("handoff original identity lineage is incomplete")
    ancestors = set()
    pending = list(native.derived_from_identity_ids)
    while pending:
        key = pending.pop()
        if key not in ancestors:
            ancestors.add(key)
            pending.extend(mappings[key].derived_from_identity_ids)
    if (set(sources) - ancestors or ancestors != set(mappings) - {native.identity_id}):
        raise WgsHandoffError("handoff native lineage has missing or extra contributors")
    normalized = [row for row in mappings.values() if row.stage is IdentityStage.NORMALIZED]
    for mapping in mappings.values():
        identity = getattr(mapping.identity, "sequence", mapping.identity)
        if (identity.reference_identity_id != reference_identity_id
                or mapping.reference_manifest_sha256 != reference_manifest_sha256
                or identity.reference_role != getattr(native.identity, "sequence", native.identity).reference_role):
            raise WgsHandoffError("handoff identity reference binding changed")
    if isinstance(native.identity, SequenceAllele) and normalized and (
            len(normalized) != 1 or normalized[0].identity != native.identity):
        raise WgsHandoffError("handoff normalized and native geometry disagree")
    interpretation = (native.source.source_fields or {}).get("native_interpretation") or {}
    if interpretation.get("repeat_representation") == "allele_projection":
        if interpretation.get("parent_identity_id") not in sources:
            raise WgsHandoffError("conditional source projection parent changed")
        _sha(interpretation.get("source_format_sha256"))
        _sha(interpretation.get("protocol_sha256"))
    declarations = set()
    for call in allele.calls:
        fields = call.source.source_fields or {}
        source_id = fields.get("source_allele_id")
        source = sources.get(source_id)
        if (source is None or call.sample is None
                or fields.get("source_input_sha256", call.source.source_sha256) not in source_inputs
                or call.source.source_sha256 != source.source.source_sha256
                or call.source.source_record_id != source.source.source_record_id):
            raise WgsHandoffError("handoff source call is undeclared or misbound")
        declaration = (source_id, call.sample.sample_name)
        if declaration in declarations:
            raise WgsHandoffError("handoff source call is duplicated")
        declarations.add(declaration)
    if {source_id for source_id, _ in declarations} != set(sources):
        raise WgsHandoffError("handoff omits original source calls")


def validate_handoff_bundle(manifest, *, root, directory, scratch_path,
                            validated_reader=None):
    """Rehash metadata and authenticate the complete normalized scientific closure.

    Source VCFs and public datasets remain upstream-attested. All persisted
    scientific and provenance records used in the handoff are authenticated.
    """
    from genosapiens.bounded_wgs_config import load_bounded_wgs_execution_config
    from genosapiens.case_manifest import CaseFileClass, load_case_manifest, case_manifest_as_bytes
    from genosapiens.track1_preflight import CompetitiveRunEvidence, ReleasePreflight
    from genosapiens.track1_installation import InstalledVersions, load_run_versions as load_base_release
    from genosapiens.track1_publication import load_final_output_manifest, FinalRole
    from genosapiens.wgs_ledger import (
        PreparedCaseManifest, LedgerManifest, FilterManifest, WorkingSetManifest, load_model,
    )
    from genosapiens.wgs_prepared_case import (
        LedgerBoundEarlySnapshot, EARLY_SNAPSHOT_SCHEMA, _decode_envelope,
        load_prepared_index,
    )
    from genosapiens.wgs_preparation import (
        load_source_accounting, load_source_audit, load_preparation_run_options,
        load_completed_source_input, resolve_completed_source,
    )
    root = Path(root)
    artifacts = {row.role: row for row in manifest.artifacts}
    contents = {role: read_handoff_artifact(root, artifact) for role, artifact in artifacts.items()}
    config = load_bounded_wgs_execution_config(contents["execution"], artifacts["execution"].sha256)
    prepared = load_model(contents["prepared"], PreparedCaseManifest)
    source = load_model(contents["source"], LedgerManifest)
    completed_input = load_completed_source_input(contents["completed_source_input"])
    completed_source = resolve_completed_source(completed_input,
        expected_source=source.source, audit_requested=completed_input.audit_requested)
    source_reports = {row.name: row.sha256 for row in source.summary.reports}
    accounting = load_source_accounting(contents["source_accounting"])
    audit = load_source_audit(contents["source_audit"])
    options = load_preparation_run_options(contents["preparation_run_options"])
    prepared_index = load_prepared_index(contents["prepared_index"])
    filtered = load_model(contents["prefilter"], FilterManifest)
    working = load_model(contents["working"], WorkingSetManifest)
    early = _decode_envelope(contents["early"], LedgerBoundEarlySnapshot,
                            schema=EARLY_SNAPSHOT_SCHEMA, kind="ledger_bound_early_snapshot")
    release = load_base_release(contents["release"], expected_sha256=artifacts["release"].sha256)
    case = load_case_manifest(root / artifacts["case"].relative_path)
    run = _typed(json.loads(contents["run_evidence"]), CompetitiveRunEvidence, decode=True)
    preflight = _typed(json.loads(contents["preflight"]), ReleasePreflight, decode=True)
    output = load_final_output_manifest(contents["output"])
    products = load_handoff_science(root, manifest)
    check_scientific_shard_limits(products, config)
    evidence = products.scientific.evidence
    if evidence.as_bytes() != contents["evidence"]:
        raise WgsHandoffError("handoff normalized evidence declaration changed")
    expected_parents = {
        "prepared_case": prepared.sha256, "source_manifest": source.sha256,
        "prefilter_manifest": filtered.sha256, "working_set": working.sha256,
        "execution_config": artifacts["execution"].sha256,
        "policy_bindings": release.policy_bindings.sha256,
    }
    from genosapiens.wgs_analysis_eligibility import scientific_eligibility_component,eligibility_component_parents
    eligibility=scientific_eligibility_component(products.scientific,products.scientific_directory,eligibility_root=getattr(products,"eligibility_root",None))
    if eligibility is not None:
        expected_parents.update(eligibility_component_parents(eligibility))
    from genosapiens.wgs_profile_projection import scientific_profile_component, profile_component_parents
    profile = scientific_profile_component(products.scientific, products.scientific_directory,
        profile_root=products.profile_root)
    if profile is None or profile['manifest_sha256'] != artifacts['analysis_profile'].sha256:
        raise WgsHandoffError("handoff analysis profile differs from scientific parents")
    if validated_reader is None:
        # This independent consumer requests complete handoff validation. Shared
        # publication owners already authenticate their copied profile bytes.
        from genosapiens.wgs_profile_projection import verify_analysis_profile_component
        verify_analysis_profile_component(products.scientific_directory.parent, profile,
            full=True, scratch_directory=Path(scratch_path).parent, profile_root=products.profile_root)
    expected_parents.update(profile_component_parents(profile))
    from genosapiens.wgs_association_projection import verify_assessment_projections, assessment_parents
    projections = verify_assessment_projections(products.scientific_directory.parent,
        dict(products.scientific.parent_bindings), full=False, scratch_directory=Path(scratch_path).parent)
    if projections is not None:
        expected_parents.update(assessment_parents(*projections))
    source_inputs = {item.sha256 for item in case.inputs if item.file_class is CaseFileClass.VCF}
    read_intake = None
    if "read_intake" in artifacts:
        from genosapiens.wgs_read_rescue import (
            load_merged_intake, MergedReadIntake, MergedIntakeIndex, validate_read_intake_run,
        )
        intake_manifest = load_merged_intake(contents["read_intake"], expected_sha256=artifacts["read_intake"].sha256)
        read_intake = MergedReadIntake(intake_manifest, (root / artifacts["read_intake"].relative_path).parent)
        validate_read_intake_run(intake_manifest, run_evidence=run, release=release)
        primary_binding = next((row for row in intake_manifest.source_bindings if row.prepared_sha256 == prepared.sha256), None)
        if (intake_manifest.case_manifest_sha256 != artifacts["case"].sha256
                or {row.input_sha256 for row in intake_manifest.source_bindings} != source_inputs
                or primary_binding is None
                or (primary_binding.early_snapshot_sha256, primary_binding.source_manifest_sha256,
                    primary_binding.prefilter_manifest_sha256, primary_binding.working_set_sha256,
                    primary_binding.input_sha256, primary_binding.original_occurrences)
                != (early.sha256, source.sha256, filtered.sha256, working.sha256,
                    source.source.input_sha256, source.summary.total_physical_alternate_slots)
                or any((delta.case_id, delta.sample_id) != (case.case_id, case.sample_id)
                       for delta in intake_manifest.deltas)):
            raise WgsHandoffError("handoff read intake source authority changed")
        expected_parents["read_intake"] = intake_manifest.sha256
        source_inputs.update(stage.output_artifact.sha256 for delta in intake_manifest.deltas for stage in delta.stage_results)
    elif run.rescue_execution_status != "not_performed":
        raise WgsHandoffError("handoff rescue run has no persisted intake")
    if (prepared.completed_source_input.sha256 != completed_input.sha256
            or prepared_index.completed_source_input_sha256 != completed_input.sha256
            or early.completed_source_input_sha256 != completed_input.sha256
            or prepared_index.completed_source_state_root != completed_input.source_state_root
            or prepared.source_manifest.sha256 != source.sha256
            or prepared.source_accounting.sha256 != accounting.sha256
            or prepared.source_audit.sha256 != audit.sha256
            or prepared.preparation_run_options.sha256 != options.sha256
            or prepared.prepared_index.sha256 != prepared_index.sha256
            or early.source_accounting_sha256 != accounting.sha256
            or early.source_audit_sha256 != audit.sha256
            or early.preparation_run_options_sha256 != options.sha256
            or early.prepared_index_sha256 != prepared_index.sha256
            or prepared_index.source_manifest_sha256 != source.sha256
            or prepared_index.filter_manifest_sha256 != filtered.sha256
            or prepared_index.working_manifest_sha256 != working.sha256
            or prepared_index.source_accounting_sha256 != accounting.sha256
            or prepared_index.source_audit_sha256 != audit.sha256
            or prepared_index.preparation_run_options_sha256 != options.sha256
            or accounting.source_manifest_sha256 != source.sha256
            or accounting.input_sha256 != source.source.input_sha256
            or accounting.header_sha256 != source.header.sha256
            or accounting.qc_summary_sha256 != source_reports.get("qc.json")
            or accounting.normalization_summary_sha256 != source_reports.get("normalization.json")
            or accounting.class_summary_sha256 != source_reports.get("classes.json")
            or accounting.total_records != source.summary.total_records
            or accounting.total_physical_alternate_slots != source.summary.total_physical_alternate_slots
            or accounting.total_biological_alternates != source.summary.total_biological_alternates
            or audit.accounting_sha256 != accounting.sha256
            or audit.source_manifest_sha256 != source.sha256
            or audit.run_options_sha256 != options.sha256
            or audit.requested != options.verify_source_ledger
            or audit.status != ("passed" if options.verify_source_ledger else "not_requested")
            or audit.input_sha256 != source.source.input_sha256
            or audit.reference_manifest_sha256 != source.source.reference_manifest_sha256
            or audit.executable_stack_sha256 != completed_input.producer.source.executable_stack_sha256
            or audit.runtime_image_digest != completed_input.producer.runtime_image_digest
            or run.as_bytes() != contents["run_evidence"]
            or preflight.as_bytes() != contents["preflight"]
            or early.as_bytes() != contents["early"]
            or case_manifest_as_bytes(case) != contents["case"]
            or output.competitive_status != "ready"
            or output.input_manifest_sha256 != prepared.sha256
            or output.early_snapshot_sha256 != early.sha256
            or output.run_id != run.run_id
            or prepared.prefilter_manifest.sha256 != filtered.sha256
            or prepared.working_sets[0].sha256 != working.sha256
            or prepared.case_manifest_sha256 != artifacts["case"].sha256
            or source.source.case_manifest_sha256 != artifacts["case"].sha256
            or prepared.original_input_sha256 != source.source.input_sha256
            or source.source.input_sha256 not in source_inputs
            or type(release) is not InstalledVersions and prepared.base_release_sha256 != release.sha256
            or prepared.runtime_image_digest != release.runtime_image_digest
            or prepared.execution_limits_sha256 != artifacts["execution"].sha256
            or release.policy_bindings.bounded_wgs_execution_sha256 != artifacts["execution"].sha256
            or run.input_manifest_sha256 != prepared.sha256
            or run.early_snapshot_sha256 != artifacts["early"].sha256
            or early.prepared_case_sha256 != prepared.sha256
            or run.primary_snapshot_sha256 != products.primary.sha256
            or run.selection_sha256 != products.selection.sha256
            or run.release_manifest_sha256 != release.sha256
            or run.accepted_reference_status != "accepted" or run.concordance_status != "complete"
            or run.accepted_reference_manifest_sha256 != source.source.reference_manifest_sha256
            or run.accepted_reference_manifest_sha256 != release.reference_manifest_sha256
            or run.accepted_reference_manifest_sha256 != case.reference_bundle.manifest_sha256
            or run.concordance_report_sha256 != artifacts["concordance"].sha256
            or prepared.concordance_sha256 != artifacts["concordance"].sha256
            or preflight.run_evidence_sha256 != artifacts["run_evidence"].sha256
            or preflight.release_manifest_sha256 != release.sha256
            or preflight.selection_sha256 != products.selection.sha256
            or preflight.epcr_assignment_sha256 != artifacts["epcr"].sha256
            or run.epcr_assignment_sha256 != artifacts["epcr"].sha256
            or not preflight.ready or preflight.competitive_csv_sha256 != artifacts["csv"].sha256
            or evidence.sha256 != products.scientific.evidence.sha256
            or products.scientific.policy_sha256 != release.policy_bindings.fusion_policy_sha256):
        raise WgsHandoffError("handoff prepared authority or readiness binding changed")
    outputs = {row.role: row for row in output.outputs}
    if read_intake is not None:
        bound = outputs.get(FinalRole.READ_INTAKE.value)
        artifact = artifacts["read_intake"]
        if bound is None or (bound.sha256, bound.size_bytes, bound.relative_path) != (
                artifact.sha256, artifact.size_bytes, artifact.relative_path):
            raise WgsHandoffError("handoff read intake output declaration changed")
    for role, final_role in (("evidence", FinalRole.PROVENANCE), ("scientific", FinalRole.HYPOTHESES),
            ("order", FinalRole.RANKING), ("selection", FinalRole.SELECTION),
            ("preflight", FinalRole.PREFLIGHT), ("epcr", FinalRole.EPCR), ("csv", FinalRole.CSV)):
        bound = outputs[final_role.value]
        artifact = artifacts[role]
        if (bound.sha256, bound.size_bytes, bound.relative_path) != (
                artifact.sha256, artifact.size_bytes, artifact.relative_path):
            raise WgsHandoffError("handoff final output declaration changed")
    expected_run = {
        "run_id": run.run_id, "prepared_manifest_sha256": prepared.sha256,
        "source_input_sha256": source.source.input_sha256,
        "reference_identity_id": run.accepted_reference_identity_id,
        "reference_manifest_sha256": run.accepted_reference_manifest_sha256,
        "concordance_report_sha256": run.concordance_report_sha256,
        "configuration_sha256": case.configuration.sha256,
        "image_digest": release.runtime_image_digest,
        "release_manifest_sha256": release.sha256,
        "output_manifest_sha256": artifacts["output"].sha256,
    }
    expected_run['retention_scope'] = retention_scope(profile)
    if dict(manifest.run) != expected_run:
        raise WgsHandoffError("handoff run declarations changed")
    owner = (open_scientific_products(products.scientific, product_descriptors(products),
        directory=products.scientific_directory, expected_parent_bindings=expected_parents,
        scratch_directory=Path(scratch_path).parent,eligibility_root=products.eligibility_root,
        profile_root=products.profile_root) if validated_reader is None else nullcontext(validated_reader))
    with owner as reader, reader.reader.read_scope():
        if type(reader) is not ScientificArtifactReader or reader.parents != expected_parents:
            raise WgsHandoffError("handoff verified scientific parents changed")
        validate_handoff_roots(manifest, products, reader)
        if validated_reader is not None:
            # Rehearsal links the actual owner's authenticated leaves. The
            # publication path separately copies and hashes immutable bytes.
            for shard in (*evidence.node_shards, *evidence.root_shards):
                original = Path(reader.directory) / shard.relative_path
                staged = products.scientific_directory / shard.relative_path
                a, b = original.lstat(), staged.lstat()
                if (not stat.S_ISREG(a.st_mode) or not stat.S_ISREG(b.st_mode)
                        or (a.st_dev, a.st_ino, a.st_size) != (b.st_dev, b.st_ino, b.st_size)):
                    raise WgsHandoffError("handoff staged scientific leaf changed")
        with ExitStack() as stack:
            intake_index = (None if read_intake is None else stack.enter_context(MergedIntakeIndex(
                read_intake, execution_config=config,
                scratch_path=Path(scratch_path).with_name(Path(scratch_path).name + ".intake"))))
            with closing(reader.iter_member_records()) as members:
                for member in members:
                    allele = member["member"].allele
                    validate_handoff_allele(allele, source_inputs=source_inputs,
                        reference_identity_id=run.accepted_reference_identity_id,
                        reference_manifest_sha256=run.accepted_reference_manifest_sha256)
                    if intake_index is not None:
                        intake_index.validate_allele(allele)
    return VerifiedHandoffBundle(manifest, root, Path(directory), config, products,
        evidence, prepared, source, release, run, case, Path(scratch_path).parent, read_intake)
