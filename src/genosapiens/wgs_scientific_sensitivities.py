"""Explicit, immutable optional analyses of one sealed primary universe.

Copy the normalized DAG without expanding rich pairs, then use the existing
source, score-view and ordering kernels. The derivative has its own writer,
completion descriptor and strict replay; no primary file is opened for writing.
"""
from contextlib import ExitStack, nullcontext
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path

from genosapiens.evidence_model import _exact, _sha, _typed
from genosapiens.evidence_serialization import _canonical, _loads
from genosapiens.wgs_case_evidence import CompleteScientificCase
from genosapiens.wgs_case_resources import bound_reference, component_reference
from genosapiens.wgs_context_order import CompleteContextOrder
from genosapiens.wgs_evidence_replay import EvidenceReplay
from genosapiens.wgs_evidence_store import EvidenceManifest, EvidenceReader, EvidenceWriter, EvidenceStoreError
from genosapiens.wgs_scientific_artifact import ScientificArtifactReader
from genosapiens.wgs_scientific_products import (
    ScientificPrimaryProduct, ScientificSensitivityProduct, ScientificSensitivityProducts, scientific_products,
)
from genosapiens.wgs_score_views import CompleteContextScoreView, VerifiedContextScoreView
from genosapiens.wgs_sensitivity_evidence import (
    CaseAblationEvidence, CaseDiagnosisEvidence, ContextAblationEvidence, ContextDiagnosisEvidence,
)

MODES = ("dominant_term_sensitivity", "diagnosis_informed_sensitivity")
SCHEMA = "genosapiens.wgs-scientific-sensitivities"
KIND = "scientific-sensitivities-v1"
DESCRIPTOR = "scientific-sensitivities.json"
MAX_DESCRIPTOR_BYTES = 8 * 1024**2


@dataclass(frozen=True, slots=True, repr=False)
class ScientificSensitivities:
    evidence: EvidenceManifest
    root_sha256: str
    scientific_manifest_sha256: str
    primary_manifest_sha256: str
    candidate_universe_sha256: str
    context_count: int
    hypothesis_count: int
    schema: str = SCHEMA
    version: int = 1

    def __post_init__(self):
        for digest in (self.root_sha256, self.scientific_manifest_sha256,
                       self.primary_manifest_sha256, self.candidate_universe_sha256):
            _sha(digest)
        if (type(self.evidence) is not EvidenceManifest or self.schema != SCHEMA
                or type(self.version) is not int or self.version != 1
                or any(type(n) is not int or not 0 <= n < 2**63 for n in (self.context_count, self.hypothesis_count))):
            raise EvidenceStoreError("invalid optional scientific publication")
        self.as_bytes()

    def as_bytes(self):
        content = (_canonical({"classification": "controlled", "llm_policy": "no-llm", "sensitivities": self}) + "\n").encode("ascii")
        if len(content) > MAX_DESCRIPTOR_BYTES:
            raise EvidenceStoreError("optional scientific descriptor exceeds capacity", capacity=True)
        return content

    @property
    def sha256(self):
        return hashlib.sha256(self.as_bytes()).hexdigest()


def sensitivities_from_bytes(content, *, expected_sha256):
    _sha(expected_sha256)
    if type(content) is not bytes or len(content) > MAX_DESCRIPTOR_BYTES or hashlib.sha256(content).hexdigest() != expected_sha256:
        raise EvidenceStoreError("optional scientific descriptor binding differs")
    value = _loads(content)
    if (type(value) is not dict or set(value) != {"classification", "llm_policy", "sensitivities"}
            or value["classification"] != "controlled" or value["llm_policy"] != "no-llm"):
        raise EvidenceStoreError("invalid optional scientific envelope")
    result = _typed(value["sensitivities"], ScientificSensitivities, decode=True)
    if result.as_bytes() != content:
        raise EvidenceStoreError("noncanonical optional scientific descriptor")
    return result


def _validate_primary(reader, product):
    if type(reader) is not ScientificArtifactReader or type(product) is not ScientificPrimaryProduct:
        raise EvidenceStoreError("verified immutable primary reader required")
    reader._require_open()
    if scientific_products(reader, expected_prefilter_manifest_sha256=product.prefilter_manifest_sha256).primary != product:
        raise EvidenceStoreError("optional scientific primary product changed")


class CompleteScientificSensitivities(CompleteScientificCase):
    """Share the primary context/rank accessors, owning only optional indexes."""
    def __init__(self, writer, primary, product, scientific_sha256, reader, reference, reference_pin, renderers):
        self.writer, self.connection, self.resources = writer, writer.connection, primary.resources
        self.name, self.stack, self.closed, self.created = "sensitivities", ExitStack(), False, False
        self.primary, self.primary_product = primary, product
        self.table = "scientific_sensitivity_views"
        self.published = False
        self.orders = {}
        for field in ("context_count", "hypothesis_count", "identities"):
            setattr(self, field, getattr(primary, field))
        try:
            self.resources.require(writer)
            self.connection.execute(f"CREATE TABLE {self.table} (ordinal INTEGER PRIMARY KEY,views TEXT)")
            self.created = True
            order = primary.orders["primary"]
            locator = next(order.tape.rows(start=0, count=min(1, self.hypothesis_count)), None)
            leader = (nullcontext(None) if locator is None else
                self._restore(reader, primary.context(locator[0]), reference, reference_pin, renderers))
            with leader as context:
                self.ablation = CaseAblationEvidence(writer, self.resources.phenotypes, order,
                    None if context is None else context.scores, prior_policy=self.resources.prior_policy, name=self.name + ".ablation")
            self.diagnosis = CaseDiagnosisEvidence(writer, self.resources.genes,
                diagnosis=self.resources.diagnosis, prior_policy=self.resources.prior_policy, name=self.name + ".diagnosis")
            self.diagnosis_matched_hypotheses = 0
            for mode in MODES:
                self.orders[mode] = self._order(mode)
            for row in primary.contexts():
                with self._restore(reader, row, reference, reference_pin, renderers) as context:
                    evidence = (ContextAblationEvidence(writer, self.ablation, context.scores, name=context.name + ".ablation"),
                        ContextDiagnosisEvidence(writer, self.diagnosis, context.scores, name=context.name + ".diagnosis"))
                    if evidence[1].proofs:
                        self.diagnosis_matched_hypotheses += context.hypothesis_count
                    views = {"ordinal": context.ordinal, "primary": context.reference, "modes": {}}
                    for mode, source in zip(MODES, evidence, strict=True):
                        with CompleteContextScoreView(writer, context.scores, name=context.name + ".view." + mode,
                                **source.view_arguments(writer, context.scores)) as scores:
                            with CompleteContextOrder(writer, scores, name=context.name + ".order." + mode) as context_order:
                                self.orders[mode].add_context(context.ordinal, context_order)
                                views["modes"][mode] = {"source": source.reference,
                                    "scores": component_reference(scores, "complete-context-score-view-v1"),
                                    "order": component_reference(context_order, "complete-context-order-v1")}
                    self.connection.execute(f"INSERT INTO {self.table} VALUES (?,?)", (context.ordinal, _canonical(views)))
            for mode in MODES:
                self.orders[mode].publish()
                if (self.orders[mode].context_count != self.context_count or self.orders[mode].hypothesis_count != self.hypothesis_count):
                    raise EvidenceStoreError("optional scientific universe coverage changed")
            digest = writer.add_sequence(self.name + ".views", "scientific-case-views-v1",
                (json.loads(content) for content, in self.connection.execute(f"SELECT views FROM {self.table} ORDER BY ordinal")))
            self.header = {"schema": SCHEMA, "version": 1, "primary": primary.reference,
                "scientific_manifest_sha256": scientific_sha256, "primary_manifest_sha256": product.sha256,
                "candidate_universe_sha256": product.candidate_universe_sha256,
                "assessment_domain_sha256": self.resources.assessment_domain_sha256,
                "views": {"name": self.name + ".views", "kind": "scientific-case-views-v1", "sha256": digest},
                "ablation": self.ablation.reference, "diagnosis": self.diagnosis.reference,
                "orders": {mode: self.orders[mode].reference for mode in MODES},
                "context_count": self.context_count, "hypothesis_count": self.hypothesis_count}
            self.sha256 = writer.add(self.name, KIND, self.header)
            self.published = True
        except BaseException:
            writer.failed = True
            self.close()
            raise

    @property
    def reference(self):
        self._require_open()
        if not self.published:
            raise EvidenceStoreError("optional scientific case is incomplete")
        return component_reference(self, KIND)

    def context(self, ordinal):
        return self.primary.context(ordinal)

    def context_scores(self, reader, context, mode, stack):
        self._require_open()
        if mode not in MODES or context.resources is not self.resources or reader.connection is not self.connection:
            raise EvidenceStoreError("foreign optional scientific score view")
        content, = self.connection.execute(f"SELECT views FROM {self.table} WHERE ordinal=?", (context.ordinal,)).fetchone()
        view = json.loads(content)["modes"][mode]
        source = (ContextAblationEvidence(EvidenceReplay(reader), self.ablation, context.scores, name=context.name + ".ablation")
            if mode == MODES[0] else
            ContextDiagnosisEvidence(EvidenceReplay(reader), self.diagnosis, context.scores, name=context.name + ".diagnosis"))
        if source.reference != view["source"]:
            raise EvidenceStoreError("optional scientific source changed")
        return stack.enter_context(VerifiedContextScoreView(reader,
            bound_reference(reader, view["scores"], "complete-context-score-view-v1"), context.scores,
            **source.view_arguments(EvidenceReplay(reader), context.scores)))


def _copied_primary(original, reader, stack):
    """Replay the exact primary roots within the derivative's own connection."""
    copy = ScientificArtifactReader(original.artifact, directory=reader.directory,
        expected_parent_bindings=original.parents, expected_reference_binding=original.expected_reference_binding,
        scratch_directory=reader.directory)
    stack.callback(copy.close)
    copy.reader, copy.entered = reader, True
    from genosapiens.wgs_association_projection import verify_assessment_projections
    copy.assessment_projections=verify_assessment_projections(reader.directory.parent,
        original.parents,full=False)
    copy._restore_primary()
    return copy


def _retain_components(primary, root):
    from genosapiens.wgs_analysis_eligibility import eligibility_retained_files
    from genosapiens.wgs_profile_projection import profile_retained_files
    from genosapiens.wgs_association_projection import assessment_retained_files
    from genosapiens.wgs_ledger import publish_stream
    from itertools import chain
    component = primary.resources.analysis_eligibility
    eligibility = () if component is None else eligibility_retained_files(primary.directory.parent, component,eligibility_root=primary.eligibility_root)
    component = primary.resources.analysis_profile
    profile = () if component is None else profile_retained_files(primary.directory.parent,component,profile_root=primary.profile_root)
    for source, relative, digest, size in chain(eligibility,profile,assessment_retained_files(primary.directory.parent, primary.parents)):
        target = root / relative
        target.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
        with source.open("rb") as handle:
            observed = hashlib.sha256()
            def chunks():
                while content := handle.read(1024**2):
                    observed.update(content)
                    yield content
            publish_stream(target.parent, target.name, chunks(), maximum_bytes=size)
        if observed.hexdigest() != digest or target.stat().st_size != size:
            raise EvidenceStoreError("optional retained scientific component changed")


def write_scientific_sensitivities(primary, product, *, output_directory):
    _validate_primary(primary, product)
    root = Path(output_directory)
    if (not root.is_absolute() or any(p.is_symlink() for p in (root, *root.parents))
            or any(root.is_relative_to(source) for source in (primary.directory,
                primary.directory.parent / "eligibility", primary.directory.parent / "associations",
                primary.profile_root))):
        raise EvidenceStoreError("invalid optional scientific output directory")
    root.mkdir(mode=0o700, parents=True, exist_ok=True)
    if any(root.iterdir()):
        raise EvidenceStoreError("optional scientific output must be empty")
    _retain_components(primary, root)
    with EvidenceWriter(root / "evidence", scratch_directory=root) as writer:
        writer.copy_from(primary.reader)
        with writer.read_view() as reader, ExitStack() as stack:
            copy = _copied_primary(primary, reader, stack)
            case = stack.enter_context(CompleteScientificSensitivities(writer, copy.case, product,
                primary.artifact.sha256, reader, copy.reference, copy.reference_pin, copy.renderers))
            values = dict(root_sha256=case.sha256, scientific_manifest_sha256=primary.artifact.sha256,
                primary_manifest_sha256=product.sha256, candidate_universe_sha256=product.candidate_universe_sha256,
                context_count=case.context_count, hypothesis_count=case.hypothesis_count)
        evidence = writer.finish()
    artifact = ScientificSensitivities(evidence=evidence, **values)
    from genosapiens.wgs_ledger import publish_bytes
    publish_bytes(root, DESCRIPTOR, artifact.as_bytes())
    return artifact


class ScientificSensitivitiesReader:
    def __init__(self, artifact, primary, product, *, directory, scratch_directory):
        _validate_primary(primary, product)
        if (type(artifact) is not ScientificSensitivities or artifact.scientific_manifest_sha256 != primary.artifact.sha256
                or artifact.primary_manifest_sha256 != product.sha256 or artifact.candidate_universe_sha256 != product.candidate_universe_sha256):
            raise EvidenceStoreError("foreign optional scientific primary")
        self.artifact, self.primary, self.product = artifact, primary, product
        self.directory, self.scratch = Path(directory), Path(scratch_directory)
        self.stack, self.closed, self.entered = ExitStack(), False, False

    def __enter__(self):
        if self.closed or self.entered:
            raise EvidenceStoreError("optional scientific reader is unavailable")
        self.entered = True
        try:
            from genosapiens.wgs_association_projection import verify_assessment_projections
            verify_assessment_projections(self.directory, self.primary.parents, full=False, scratch_directory=self.scratch)
            self.reader = self.stack.enter_context(EvidenceReader(self.artifact.evidence, self.directory / "evidence", scratch_directory=self.scratch))
            with self.reader.read_scope(), self.reader.complete_root_scope():
                self.copy = _copied_primary(self.primary, self.reader, self.stack)
                root = self.reader.reference("sensitivities", kind=KIND)
                if root.sha256 != self.artifact.root_sha256:
                    raise EvidenceStoreError("optional scientific root changed")
                self.case = self.stack.enter_context(CompleteScientificSensitivities(EvidenceReplay(self.reader), self.copy.case,
                    self.product, self.primary.artifact.sha256, self.reader, self.copy.reference, self.copy.reference_pin, self.copy.renderers))
                if any(getattr(self.case, key) != getattr(self.artifact, key) for key in ("context_count", "hypothesis_count")):
                    raise EvidenceStoreError("optional scientific counts changed")
            return self
        except BaseException:
            self.close()
            raise

    def _require_open(self):
        if self.closed or not self.entered:
            raise EvidenceStoreError("optional scientific reader is closed")
        self.reader._check_files()
        self.primary._require_open()

    def ordered(self, mode):
        self._require_open()
        with self.reader.read_scope():
            yield from self.case.ordered(mode)

    def record_at_rank(self, rank, *, mode):
        self._require_open()
        with self.reader.read_scope():
            return self.case.record_at_rank(self.reader, rank, self.copy.reference, self.copy.reference_pin, self.copy.renderers, mode=mode)

    def products(self):
        self._require_open()
        with self.reader.read_scope():
            case, resources, primary = self.case, self.copy.resources, self.product
            def root(value):
                return bound_reference(self.reader, value, value["kind"])
            ablation = case.ablation
            common = dict(scientific_manifest_sha256=self.artifact.sha256, root=root(case.header["views"]),
                primary_manifest_sha256=primary.sha256, prefilter_manifest_sha256=primary.prefilter_manifest_sha256,
                candidate_universe_sha256=primary.candidate_universe_sha256, ranking_policy_sha256=primary.ranking_policy_sha256,
                hypothesis_count=case.hypothesis_count, prior_policy=resources.prior_policy)
            dominant = ScientificSensitivityProduct(**common, mode=MODES[0], order_root=root(case.orders[MODES[0]].reference),
                source_root=root(ablation.reference), status=ablation.status, reason=ablation.reason,
                dominant_term=ablation.selection, removed_observation_ids=ablation.removed, diagnosis=None)
            reason = ("diagnosis_not_supplied" if resources.diagnosis is None else "empty_primary_universe"
                if not case.hypothesis_count else "no_prior_match" if not case.diagnosis_matched_hypotheses
                else "zero_prior_strength" if _exact(resources.prior_policy.strength) == 0 else None)
            diagnosis = ScientificSensitivityProduct(**common, mode=MODES[1], order_root=root(case.orders[MODES[1]].reference),
                source_root=root(case.diagnosis.reference), status="not_supplied" if resources.diagnosis is None else
                "not_applicable" if reason else "applied", reason=reason, dominant_term=None, removed_observation_ids=(), diagnosis=resources.diagnosis)
            return ScientificSensitivityProducts(dominant, diagnosis)

    def close(self):
        self.stack.close()
        self.closed = True

    def __exit__(self, *_):
        self.close()
