"""Task 7.5: pure, fixed-universe sensitivity core over task 9.1 evidence.

This is not candidate generation, a full primary ranker, or a controlled-use
release. Callers supply completed seven-component ledgers, explicit executable
phenotype/ordering contracts and a prior strength before disclosing a diagnosis.
Task 9.6 must derive and freeze numerical policy and validate release readiness.
No function accepts an acquisition, ranking, rescue, or I/O callback.
"""

from __future__ import annotations

from collections import OrderedDict, deque
from collections.abc import Mapping
from dataclasses import dataclass, replace
from fractions import Fraction
import hashlib
import json
import re
from types import MappingProxyType
from typing import Literal

from genosapiens.disease_gene_evidence import (
    assertion_native_disease_id, is_disease_identifier,
)
from genosapiens.disease_profiles import (
    DiseaseProfileStore,
)
from genosapiens.evidence_model import (
    ComponentName, DiseaseGeneEvidence, EvidenceModelError, GenotypeHypothesis,
    ScoreComponent, ScoreLedger, ScoringView, TieKey, _Controlled, _exact, _ids,
    _sha, _typed, _json_value,
)
from genosapiens.evidence_serialization import _canonical, _plain
from genosapiens.hpo_ic import HpoGraph, HpoICStore, build_hpo_graph, build_hpo_ic, hpo_ic_as_bytes
from genosapiens.ontology_bundle import ParsedOntology
from genosapiens.phenotype_scoring import (
    POLICY, POLICY_SHA256, PhenotypeScore, _score_phenotypes_batch,
)


DOMINANT_TIE = "positive-post-cap-sum-desc-canonical-hpo-id-asc"
ORDERING = "exact-total-desc-declared-ties-null-last-id-asc"
TRANSFORM = "phenotype-raw-score-affine-clamp"
_ERRORS = (ValueError, TypeError, RuntimeError, AttributeError, KeyError, OverflowError, RecursionError)
_AUTHENTICATED_RESOURCES_SEAL = object()


class SensitivityError(EvidenceModelError):
    """Static diagnostics; explicit fields and serialized bytes remain sensitive."""


def _bytes(value):
    return (_canonical(value) + "\n").encode("utf-8")


def _hash(value):
    return hashlib.sha256(_bytes(value)).hexdigest()


@dataclass(frozen=True, slots=True, repr=False)
class PhenotypeTransform(_Controlled):
    """Declarative exact affine/clamp, not fitted defaults or executable code.

    Raw is raw_evidence_units / the scorer's pinned scale. Missing evidence stays
    null; fusion uses neutral, without weight redistribution. All parameters are
    required. Unsupported transforms must be implemented/validated explicitly,
    never approximated by subtracting an old phenotype contribution.
    """

    slope: str
    intercept: str
    lower: str
    upper: str
    neutral: str
    weight: str

    def __post_init__(self):
        _Controlled.__post_init__(self)
        slope, _, lower, upper, neutral, weight = map(_exact, (
            self.slope, self.intercept, self.lower, self.upper, self.neutral, self.weight,
        ))
        if slope <= 0 or lower >= upper or not lower <= neutral <= upper or weight < 0:
            raise SensitivityError("invalid phenotype transform policy")

    @property
    def parameters(self):
        return MappingProxyType({
            **_plain(self), "operation": TRANSFORM, "scorer_sha256": POLICY_SHA256,
            "raw_scale": POLICY["scale"], "missing": "null-with-declared-neutral",
        })

    @property
    def sha256(self):
        return _hash(self.parameters)

    def component(self, score: PhenotypeScore) -> ScoreComponent:
        """Reusable primary/ablation transform; dependency audits stay in the scorer."""
        try:
            if type(score) is not PhenotypeScore or score.ledger["policy_sha256"] != POLICY_SHA256:
                raise SensitivityError("incompatible phenotype score")
            units = score.ledger["raw_evidence_units"]
            raw, normalized = self.numeric_values(units)
            return ScoreComponent(
                "phenotype", str(raw) if raw is not None else None,
                str(normalized) if normalized is not None else None, self.neutral, self.weight,
                TRANSFORM, self.sha256, self.parameters,
            )
        except _ERRORS:
            raise SensitivityError("cannot transform phenotype score") from None

    def numeric_values(self, units):
        """Shared exact transform after the owner binds the phenotype scorer."""
        if units is not None and type(units) is not int:
            raise SensitivityError("invalid phenotype raw units")
        raw = Fraction(units, POLICY["scale"]) if units is not None else None
        normalized = None if raw is None else min(_exact(self.upper), max(
            _exact(self.lower), raw * _exact(self.slope) + _exact(self.intercept),
        ))
        return raw, normalized


@dataclass(frozen=True, slots=True, repr=False)
class RankTie(_Controlled):
    """Ordered tie contract: fixed primary metadata or a recomputed component.

    Component ties read its normalized-or-neutral value. Fixed keys must be
    phenotype-independent. Numeric strings are never silently sorted as text.
    """

    name: str
    source: Literal[ComponentName, "fixed"]
    kind: Literal["exact", "text"]
    direction: Literal["ascending", "descending"]

    def __post_init__(self):
        _Controlled.__post_init__(self)
        _ids((self.name,))
        if self.source != "fixed" and self.kind != "exact":
            raise SensitivityError("component tie must be exact numeric")


def order_ledgers(
    ledgers: tuple[ScoreLedger, ...], ties: tuple[RankTie, ...], *,
    totals: Mapping[str, str] | None = None,
) -> tuple[str, ...]:
    """Bounded reusable ordering of existing scores, not candidate selection."""
    try:
        ledgers, ties = _typed(ledgers, tuple[ScoreLedger, ...]), _typed(ties, tuple[RankTie, ...])
        _ids(tuple(item.hypothesis_id for item in ledgers))
        _ids(tuple(item.name for item in ties))
        if totals is None:
            totals = {item.hypothesis_id: item.total for item in ledgers}
        elif not isinstance(totals, Mapping) or set(totals) != {item.hypothesis_id for item in ledgers}:
            raise SensitivityError("secondary total membership mismatch")
        totals = {key: _exact(value) for key, value in totals.items()}
        values = {}
        for ledger in ledgers:
            if tuple((key.name, key.direction) for key in ledger.tie_keys) != tuple(
                (key.name, key.direction) for key in ties
            ):
                raise SensitivityError("primary tie contract mismatch")
            components = {item.name: item for item in ledger.components}
            for rule, key in zip(ties, ledger.tie_keys, strict=True):
                value = _exact(key.value) if key.value is not None and rule.kind == "exact" else key.value
                if rule.source != "fixed":
                    component = components[rule.source]
                    expected = _exact(component.normalized if component.normalized is not None else component.neutral_value)
                    if value != expected:
                        raise SensitivityError("stale component tie value")
                values[ledger.hypothesis_id, rule.name] = value
        ordered = sorted(ledgers, key=lambda item: item.hypothesis_id)
        for rule in reversed(ties):
            present = [item for item in ordered if values[item.hypothesis_id, rule.name] is not None]
            absent = [item for item in ordered if values[item.hypothesis_id, rule.name] is None]
            ordered = sorted(present, key=lambda item: values[item.hypothesis_id, rule.name],
                             reverse=rule.direction == "descending") + absent
        return tuple(item.hypothesis_id for item in sorted(ordered, key=lambda item: -totals[item.hypothesis_id]))
    except _ERRORS:
        raise SensitivityError("invalid primary ordering contract") from None


@dataclass(frozen=True, slots=True, repr=False)
class PriorPolicy(_Controlled):
    """Caller-supplied validation policy pin, NOT a task 9.6 release attestation."""

    policy_id: str
    strength: str
    derivation_sha256: str
    derivation_scope: Literal["public", "independently_synthetic"]

    def __post_init__(self):
        _Controlled.__post_init__(self)
        _ids((self.policy_id,))
        _sha(self.derivation_sha256)
        if _exact(self.strength) < 0:
            raise SensitivityError("invalid diagnosis prior policy")


@dataclass(frozen=True, slots=True, repr=False)
class DisclosedDiagnosis(_Controlled):
    """Structured identifier only; syntax does not establish source eligibility.

    Declared catalog Symbol namespaces use canonical percent encoding, not raw
    narrative, URLs, arbitrary namespaces, or case/whitespace normalization.
    Exact native/equivalence evidence is still required by the prior evaluator.
    """
    disease_id: str

    def __post_init__(self):
        _Controlled.__post_init__(self)
        if not is_disease_identifier(self.disease_id):
            raise SensitivityError("diagnosis must be a structured disease identifier")


def score_primary_phenotypes(
    view: ScoringView, hypotheses: tuple[GenotypeHypothesis, ...], *,
    ic_store: HpoICStore, profile_store: DiseaseProfileStore | None,
    prepared_scores=None,
) -> Mapping[str, PhenotypeScore]:
    """Diagnosis firewall: closed observations and existing typed hypotheses only.

    Narrative, manual gene lists, diagnosis objects and evidence documents are
    not scoring inputs. No claim about how the caller generated hypotheses.
    """
    try:
        view = _typed(view, ScoringView)
        hypotheses = _typed(hypotheses, tuple[GenotypeHypothesis, ...])
        _ids(tuple(item.hypothesis_id for item in hypotheses))
        if (type(ic_store) is not HpoICStore
                or profile_store is not None and type(profile_store) is not DiseaseProfileStore):
            raise SensitivityError("typed resources required")
        scores = _score_phenotypes_batch(
            view.observations, (item.context.disease_id for item in hypotheses),
            ic_store=ic_store, profile_store=profile_store,
            prepared=prepared_scores,
        )
        return MappingProxyType({item.hypothesis_id: score
                                 for item, score in zip(hypotheses, scores, strict=True)})
    except _ERRORS:
        raise SensitivityError("invalid diagnosis-blind scoring input") from None


def replace_phenotype_component(ledger: ScoreLedger, score: PhenotypeScore, *,
                                transform: PhenotypeTransform, ties: tuple[RankTie, ...]) -> ScoreLedger:
    """Replace one component and its dynamic ties without altering fixed evidence."""
    old = next(value for value in ledger.components if value.name == "phenotype")
    phenotype = replace(transform.component(score), evidence_refs=old.evidence_refs, source_values=old.source_values)
    components = tuple(phenotype if value.name == "phenotype" else value for value in ledger.components)
    total = sum(((-1 if value.name == "technical_and_conflict_penalties" else 1) * _exact(value.weight)
        * _exact(value.normalized if value.normalized is not None else value.neutral_value)
        for value in components), Fraction(0))
    by_name = {value.name: value for value in components}
    keys = tuple(key if rule.source == "fixed" else TieKey(key.name,
        by_name[rule.source].normalized if by_name[rule.source].normalized is not None else by_name[rule.source].neutral_value,
        key.direction) for rule, key in zip(ties, ledger.tie_keys, strict=True))
    return replace(ledger, components=components, total=str(total), tie_keys=keys)


def _ic_fingerprint(store):
    if (type(store) is not HpoICStore or not store.ready or type(store.graph) is not HpoGraph
            or type(store.graph.ontology) is not ParsedOntology):
        raise SensitivityError("ready pinned IC resource required")
    content = hpo_ic_as_bytes(store)
    return _hash({"store": json.loads(content),
                  "active_aliases": store.graph.ontology.active_aliases})


def _ic_hash(store):
    resource_hash = _ic_fingerprint(store)
    content = hpo_ic_as_bytes(store)
    rebuilt = build_hpo_ic(build_hpo_graph(store.graph.ontology), (dict(item.fields) for item in store.annotations),
                           source_binding=store.source_binding)
    if hpo_ic_as_bytes(rebuilt) != content:
        raise SensitivityError("IC graph or values differ from source evidence")
    # The fingerprint includes executable alias lookup as well as serialized
    # terms, IC values, and closures.
    return resource_hash


@dataclass(frozen=True, slots=True, repr=False, init=False)
class _AuthenticatedSensitivityResources:
    """Coordinator-local proof for an already checksum-loaded public IC store."""

    ic_store: HpoICStore
    ic_sha256: str
    store_sha256: str
    _seal: object
    _binding: tuple[object, ...]

    def __repr__(self):
        return "<authenticated public sensitivity resources>"


def _authenticate_loaded_ic_store(
    store: HpoICStore, *, expected_store_sha256: str,
) -> _AuthenticatedSensitivityResources:
    """Bind the exact immutable store loaded from producer-validated pinned bytes."""
    try:
        _sha(expected_store_sha256)
        content = hpo_ic_as_bytes(store)
        if hashlib.sha256(content).hexdigest() != expected_store_sha256:
            raise SensitivityError("loaded IC store checksum mismatch")
        resource_hash = _hash({
            "store": json.loads(content),
            "active_aliases": store.graph.ontology.active_aliases,
        })
        result = object.__new__(_AuthenticatedSensitivityResources)
        graph = store.graph
        ontology = graph.ontology
        binding = (
            store,
            graph,
            ontology,
            ontology.terms,
            ontology.active_aliases,
            ontology._alias_targets,
            graph.parents,
            graph.ancestors,
            graph.resolutions,
            store.annotations,
            store.diseases,
            store.disease_counts,
            store.ic,
            store.unresolved_required,
            store.source_binding,
            resource_hash,
            expected_store_sha256,
        )
        for name, value in (
            ("ic_store", store), ("ic_sha256", resource_hash),
            ("store_sha256", expected_store_sha256),
            ("_seal", _AUTHENTICATED_RESOURCES_SEAL), ("_binding", binding),
        ):
            object.__setattr__(result, name, value)
        return result
    except _ERRORS:
        raise SensitivityError("cannot authenticate loaded IC store") from None


def _authenticated_ic_hash(resources, store: HpoICStore) -> str:
    binding = getattr(resources, "_binding", ())
    if type(resources) is _AuthenticatedSensitivityResources:
        graph = resources.ic_store.graph
        ontology = graph.ontology
        expected = (
            resources.ic_store,
            graph,
            ontology,
            ontology.terms,
            ontology.active_aliases,
            ontology._alias_targets,
            graph.parents,
            graph.ancestors,
            graph.resolutions,
            resources.ic_store.annotations,
            resources.ic_store.diseases,
            resources.ic_store.disease_counts,
            resources.ic_store.ic,
            resources.ic_store.unresolved_required,
            resources.ic_store.source_binding,
            resources.ic_sha256,
            resources.store_sha256,
        )
    else:
        expected = ()
    if (
        type(resources) is not _AuthenticatedSensitivityResources
        or getattr(resources, "_seal", None) is not _AUTHENTICATED_RESOURCES_SEAL
        or len(binding) != len(expected)
        or any(actual is not original for actual, original in zip(binding, expected, strict=True))
        or resources.ic_store is not store
    ):
        raise SensitivityError("authenticated IC resource binding mismatch")
    return resources.ic_sha256


def select_dominant_term(score: PhenotypeScore) -> Mapping[str, object]:
    """Pure positive-unit reduction; not a certificate of primary completion.

    Each contribution is counted once, never expanded through observation_ids.
    The bounded primary producer authenticates complete score membership.
    """
    try:
        if type(score) is not PhenotypeScore:
            raise SensitivityError("typed phenotype score required")
        totals, seen = {}, set()
        for item in score.ledger["contributions"]:
            term, units, identifier = item["observation_term"], item["units"], item["contribution_id"]
            if (type(term) is not str or re.fullmatch(r"HP:[0-9]{7}", term) is None
                    or type(units) is not int or units < 0 or identifier in seen
                    or item["kind"] not in ("support", "penalty")):
                raise SensitivityError("invalid phenotype contribution")
            seen.add(identifier)
            if item["kind"] == "support" and units > 0:
                totals[term] = totals.get(term, 0) + units
        denominator = sum(totals.values())
        term = min(totals, key=lambda key: (-totals[key], key)) if totals else None
        return MappingProxyType({
            "canonical_hpo_id": term, "positive_units": totals.get(term, 0),
            "all_positive_units": denominator,
            "percentage_of_positive_support": str(Fraction(100 * totals[term], denominator)) if term else None,
            "tie_policy": DOMINANT_TIE,
        })
    except _ERRORS:
        raise SensitivityError("cannot select dominant phenotype term") from None


_DIAGNOSIS_MATCH_SEAL = object()


class PreparedDiagnosisMatchContext:
    """Exact gene-evidence/diagnosis owner, with bounded immutable match facts."""
    def __init__(self, *, seal=None):
        if seal is not _DIAGNOSIS_MATCH_SEAL:
            raise SensitivityError("diagnosis matching requires preparation")
        self.seal = seal
        self.cache = OrderedDict()
        self.cache_bytes = 0

    def require(self, evidence, diagnosis):
        if (self.seal is not _DIAGNOSIS_MATCH_SEAL or evidence is not self.evidence
                or diagnosis is not self.diagnosis):
            raise SensitivityError("prepared diagnosis source owner changed")

    def remember(self, key, value):
        value = _json_value(_plain(value))
        size = len(_bytes(value)) + len(_bytes(key))
        if size <= 64 * 1024 * 1024:
            while self.cache and (len(self.cache) >= 1024 or self.cache_bytes + size > 64 * 1024 * 1024):
                _, (_, old_size) = self.cache.popitem(last=False)
                self.cache_bytes -= old_size
            self.cache[key] = (value, size)
            self.cache_bytes += size
        return value

    def cached(self, key):
        value = self.cache.get(key)
        if value is not None:
            self.cache.move_to_end(key)
            return value[0]
        return None

    def paths(self, start):
        key = ('paths', start)
        cached = self.cached(key)
        if cached is not None:
            return cached
        reached, pending = {start: ()}, deque([start])
        while pending:
            current = pending.popleft()
            for step in self.edges.get(current, ()):
                target = step['to']
                if target not in reached:
                    reached[target] = (*reached[current], step)
                    pending.append(target)
        return self.remember(key, reached)


def prepare_diagnosis_match_context(*, disease_gene_evidence, diagnosis):
    _typed(disease_gene_evidence, tuple[DiseaseGeneEvidence, ...])
    _typed(diagnosis, DisclosedDiagnosis)
    _ids(tuple(public.evidence_id for public in disease_gene_evidence))
    owner = PreparedDiagnosisMatchContext(seal=_DIAGNOSIS_MATCH_SEAL)
    owner.evidence, owner.diagnosis = disease_gene_evidence, diagnosis
    owner.evidence_sha256, owner.diagnosis_sha256 = _hash(disease_gene_evidence), _hash(diagnosis)
    edges, assertions = {}, {}
    for public in disease_gene_evidence:
        for assertion in public.store.assertions:
            assertions.setdefault((public.evidence_id, assertion.record_id), []).append(assertion)
        for relation in public.store.relations:
            if relation.kind == 'equivalent' and relation.subject_id and relation.object_id and relation.record_ids:
                for a, b in ((relation.subject_id, relation.object_id), (relation.object_id, relation.subject_id)):
                    edges.setdefault(a, []).append({
                        'from': a, 'to': b, 'record_ids': relation.record_ids,
                        'evidence_id': public.evidence_id, 'public_manifest_sha256': public.public_manifest_sha256,
                        'evidence_sha256': public.store.evidence_sha256,
                    })
    owner.edges = _json_value({key: tuple(sorted(rows, key=lambda step: (step['to'], _bytes(step))))
                              for key, rows in edges.items()})
    owner.assertions = MappingProxyType({key: tuple(rows) for key, rows in assertions.items()})
    owner.diagnosis_paths = owner.paths(diagnosis.disease_id)
    return owner


def _prior_matches_sources(hypothesis, disease_gene_evidence, diagnosis, *, prepared=None):
    _typed(hypothesis,GenotypeHypothesis)
    return prior_matches_context(hypothesis.context,hypothesis.evidence_refs,disease_gene_evidence,diagnosis,prepared=prepared)


def prior_matches_context(context, evidence_refs, disease_gene_evidence, diagnosis, *, prepared=None):
    """Original source-match kernel; allele identity/model outcomes are not inputs."""
    from genosapiens.evidence_model import HypothesisContext, EvidenceReference
    _typed(context,HypothesisContext)
    _typed(evidence_refs,tuple[EvidenceReference,...])
    if prepared is None:
        prepared = prepare_diagnosis_match_context(disease_gene_evidence=disease_gene_evidence,diagnosis=diagnosis)
    prepared.require(disease_gene_evidence,diagnosis)
    proofs = []
    key = ('matches', _hash((context,evidence_refs)))
    cached = prepared.cached(key)
    if cached is not None:
        return cached
    anchor = context.native_disease_id or context.disease_id
    anchor_paths, diagnosis_paths = prepared.paths(anchor), prepared.diagnosis_paths
    context_anchor_path = prepared.paths(context.disease_id).get(anchor)
    for public in disease_gene_evidence:
        store = public.store
        refs = {ref.record_id for ref in (*context.evidence_refs, *evidence_refs)
                if ref.evidence_id == public.evidence_id and ref.public_manifest_sha256 == public.public_manifest_sha256}
        target_assertions = tuple(assertion for record_id in refs
            for assertion in prepared.assertions.get((public.evidence_id,record_id),())
            if (context.gene_id is None or assertion.gene_id == context.gene_id)
            and (assertion.locus_id is None or assertion.locus_id == context.locus_id)
            and (context.locus_id is None or any(
                binding.source_ref.evidence_id == public.evidence_id
                and binding.source_ref.public_manifest_sha256 == public.public_manifest_sha256
                and binding.source_ref.record_id == (assertion.locus_record_id or assertion.record_id)
                for binding in context.locus_bindings
            )))
        witnesses = tuple(sorted(item.record_id for item in target_assertions
                                if assertion_native_disease_id(store, item) == anchor
                                and item.record_id in store.views.get(context.disease_id, ())))
        # A native anchor needs a target-bound mapped witness or explicit equivalence.
        if (context.native_disease_id and anchor != context.disease_id
                and not witnesses and context_anchor_path is None):
            continue
        for assertion in target_assertions:
            native = assertion_native_disease_id(store, assertion)
            if not native or native not in anchor_paths or native not in diagnosis_paths:
                continue
            proofs.append({
                "evidence_id": public.evidence_id, "public_manifest_sha256": public.public_manifest_sha256,
                "evidence_sha256": store.evidence_sha256, "assertion_record_id": assertion.record_id,
                "native_disease_id": native, "target_gene_id": assertion.gene_id, "target_locus_id": context.locus_id,
                "match": "exact_native" if native == diagnosis.disease_id else "source_equivalent",
                "context_equivalence_path": anchor_paths[native], "diagnosis_equivalence_path": diagnosis_paths[native],
                "context_anchor_equivalence_path": context_anchor_path,
                "context_native_witness_record_ids": witnesses,
                "locus_bindings": context.locus_bindings,
            })
    return prepared.remember(key, tuple(sorted(proofs, key=_bytes)))
