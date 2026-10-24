"""Task 7.3: deterministic, observation-only phenotype compatibility.

This development policy is synthetically validated, NOT the task 9.6 public
ranking freeze or a clinical probability model. No files, networks, notes,
diagnoses, genes, or inheritance hypotheses are accepted by the scorer.

One assignment node represents one resolved supplied term across ALL scopes
and polarities. Distinct supplied context variants share that node's budget;
duplicate records retain provenance without increasing it. Assignment maximizes
nonnegative evidence magnitude, not signed agreement (which would evade negative
evidence). Selected positive and negative allocations remain separately visible,
including unresolved contradictions. Correlation caps are applied AFTER exact
assignment, not claimed to be a joint constrained-assignment optimum.
"""

from __future__ import annotations

from collections import OrderedDict
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass, field, fields, is_dataclass, replace
from fractions import Fraction
import hashlib
import json
import math
import re
from types import MappingProxyType

from genosapiens.disease_gene_evidence import is_disease_identifier
from genosapiens.disease_profiles import DiseaseProfileStore, disease_profiles_as_bytes
from genosapiens.hpo_ic import HpoICStore, hpo_ic_as_bytes
from genosapiens.ontology_bundle import OBSOLETE_RESOLUTION_POLICY


SCHEMA = "genosapiens.phenotype-scoring"
_TOKEN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.:-]{0,127}$")
_HPO = re.compile(r"^HP:[0-9]{7}$")
_SHA = re.compile(r"^[0-9a-f]{64}$")


class PhenotypeScoringError(ValueError):
    """Invalid structured evidence; errors never interpolate case content."""


def _plain(value):
    if isinstance(value, Mapping):
        return {key: _plain(child) for key, child in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(child) for child in value]
    return value


def _bytes(value):
    try:
        return (json.dumps(_plain(value), sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode()
    except (ValueError, TypeError, RecursionError):
        raise PhenotypeScoringError("evidence must be finite JSON") from None


def _freeze(value):
    if isinstance(value, Mapping):
        return MappingProxyType({key: _freeze(child) for key, child in value.items()})
    if isinstance(value, (list, tuple)):
        return tuple(_freeze(child) for child in value)
    return value


def _digest(value):
    return hashlib.sha256(_bytes(value)).hexdigest()


# Category bounds are definitions, not estimated patient probabilities. Onset
# windows use years since birth, 365.25 days/year and half-open postnatal ranges.
# Prenatal lower bounds are intentionally open, not invented gestational ages.
POLICY = _freeze({
    "id": "phenotype-compatibility-development",
    "validation": "synthetic-only; task-9.6-public-derivation-and-freeze-pending",
    "scale": 1_000_000,
    "semantic": "resnik_ic/(1+resnik_ic); null -> neutral zero",
    "assignment": "maximum-sum-integer-magnitude; one canonical term per side across all scopes/statuses",
    "tie": "sorted term IDs; Hungarian ascending columns; strict predecessor improvement; omit zero edges",
    "rounding": "uncapped: original-value floors then largest original fractional remainders to floor(total); explicit cap: proportional scaling first; ties by allocation key",
    "observation_duplicates": "same canonical term and scoring context share provenance, not weight",
    "observation_numeric_identity": "same ID requires identical canonical raw JSON; scoring normalizes int/float/signed zero; smallest observation ID represents equivalent variants; all distinct raw inputs retained",
    "observation_variants": "equal budget across distinct supplied scope/subject/status/context variants",
    "source_variants": "minimum factor over all applicable/unknown-applicability annotations in all P groups; exclude only conclusively sex-inapplicable variants; no eligible variants -> 0",
    "source_conflicts": "among aggregation-eligible variants, positive versus negative/unknown/unresolved polarity: no positive coercion",
    "missing_frequency_positive": 1,
    "missing_frequency_negative": 1,
    "frequency_positive": "(1+lower_bound)/2; zero frequency -> zero",
    "frequency_negative": "lower_bound; missing/unsupported -> neutral one, never imputed source frequency",
    "negative_weight": 0.5,
    "negative_entailment": "absent term must equal or be ancestor of profile term; siblings are not contradictions",
    "onset_positive": "evaluate every known constraint even alongside unknown alternatives: overlap -> 1; known disjoint or present before minimum onset age -> 0.5; unknown -> neutral 1",
    "onset_negative": "known age before lower -> 0; within/open window -> 0.5; at/after finite upper -> 1",
    "context_negative": "all sourced onset alternatives, sex, modifiers must be evaluable; unknown requirement -> 0",
    "unspecified_context": "no known requirement -> neutral 1, not an assertion of age/sex independence",
    "sex_positive": "known incompatible -> zero effect retained outside aggregation; unknown -> neutral 1, retained inside aggregation",
    "quality": "supplied quality in [0,1], missing -> neutral 1; propagated source -> 0.5, missing -> neutral 1",
    "contradictions": "retain both polarities and provenance unresolved; shared term budget, no adjudication",
    "scope_weights": {"proband": 1, "perinatal": 1, "family_history": 1},
    "family_history": "contextual phenotype only; no parental disease, segregation, inherited/de-novo or mechanism inference",
    "caps": "positive/negative separately; source reuse then same-domain dependency then growth then perinatal; proportional integer allocation",
    "source_reuse_cap": "largest selected single-term allocation for a repeated native source record",
    "dependency_cap_units": 1_000_000,
    "growth_cap_units": 1_000_000,
    "perinatal_cap_units": 1_000_000,
    "growth_roots": ["HP:0001507"],
    "growth_terms": ["HP:0001518", "HP:0001510", "HP:0004322"],
    "growth_mapping": "HPO growth abnormality ancestry plus small-for-gestational-age, growth-delay and short-stature terms; supplied growth domain also eligible",
    "unlisted_features": "unknown; no profile-completeness or all-domain checklist penalty",
    "interpretation": "uncalibrated compatibility only; no isolated finding is diagnostic; no candidate-set normalization",
    "mapping_provenance": {
        "source": "https://hpo.jax.org/browse/term/",
        "verified_on": "2026-09-06",
        "method": "ToolUniverse HPO_get_term live JAX definitions, not an immutable HPO release pin; engineering transforms are separate",
        "age_convention": "years since birth; 365.25 days/year; juvenile 5-15 completed years -> [5,16); pediatric excludes first 28 days",
        "prenatal_coarsening": "embryonal/fetal both map to before birth; no gestational duration inferred; distinct source terms retained",
    },
    "frequency_bounds": {
        "HP:0040280": [1, 1], "HP:0040281": [0.80, 0.99],
        "HP:0040282": [0.30, 0.79], "HP:0040283": [0.05, 0.29],
        "HP:0040284": [0.01, 0.04], "HP:0040285": [0, 0],
    },
    "frequency_names": {
        "Obligate": "HP:0040280", "Very frequent": "HP:0040281",
        "Frequent": "HP:0040282", "Occasional": "HP:0040283",
        "Very rare": "HP:0040284", "Excluded": "HP:0040285",
    },
    "onset_windows": {
        "HP:0030674": [None, 0], "HP:0011460": [None, 0], "HP:0011461": [None, 0],
        "HP:0003577": [0, 0], "HP:0003623": [0, 28 / 365.25],
        "HP:0003593": [28 / 365.25, 1], "HP:0011463": [1, 5],
        "HP:0003621": [5, 16], "HP:0410280": [28 / 365.25, 16],
        "HP:0003581": [16, None], "HP:0011462": [16, 40],
        "HP:0003596": [40, 60], "HP:0003584": [60, None],
    },
    "sex_mapping": {"F": "female", "FEMALE": "female", "M": "male", "MALE": "male"},
})
POLICY_SHA256 = _digest(POLICY)


@dataclass(frozen=True, slots=True, repr=False)
class ObservationProvenance:
    input_sha256: str
    source_locator: str
    input_type: str

    def __post_init__(self):
        if (not isinstance(self.input_sha256, str) or not _SHA.fullmatch(self.input_sha256)
                or not isinstance(self.source_locator, str)
                or len(self.source_locator) > 256
                or not re.fullmatch(r"(?:line:[1-9][0-9]*|record:[A-Za-z0-9_.:-]+)", self.source_locator)
                or self.input_type not in ("hpo-list", "structured", "synthetic")):
            raise PhenotypeScoringError("observation provenance is invalid")


@dataclass(frozen=True, slots=True, repr=False)
class ScoringObservation:
    """Closed scoring view, not an ingestion or handoff schema.

    Provenance binds all populated fields. Age is this observation subject's age
    at assessment, NEVER the proband's age implicitly copied to family history.
    Onset accepts only an HPO identifier, not unrestricted narrative. Optional
    subject_id is an opaque within-case identity; absent family identity stays
    unknown. Domain/dependency labels are supplied structured context, not genes.
    """

    observation_id: str
    source_hpo_id: str
    scope: str
    status: str
    provenance: ObservationProvenance
    subject_id: str | None = None
    onset: str | None = None
    age_years: float | None = None
    sex: str | None = None
    domain: str | None = None
    dependency_cluster: str | None = None
    quality: float | None = None

    def __post_init__(self):
        tokens = (self.observation_id, self.subject_id, self.domain, self.dependency_cluster)
        if (any(value is not None and (not isinstance(value, str) or not _TOKEN.fullmatch(value)) for value in tokens)
                or self.observation_id is None
                or not isinstance(self.source_hpo_id, str) or not _HPO.fullmatch(self.source_hpo_id)
                or not isinstance(self.scope, str) or self.scope not in POLICY["scope_weights"]
                or self.status not in ("present", "explicitly_absent")
                or not isinstance(self.provenance, ObservationProvenance)
                or self.onset is not None and (not isinstance(self.onset, str) or not _HPO.fullmatch(self.onset))
                or self.sex not in (None, "female", "male")
                or self.dependency_cluster is not None and self.domain is None):
            raise PhenotypeScoringError("structured scoring observation is invalid")
        for value, low, high in ((self.age_years, -1, 200), (self.quality, 0, 1)):
            if value is not None and (type(value) not in (int, float) or not math.isfinite(value) or not low <= value <= high):
                raise PhenotypeScoringError("observation numeric context is invalid")


@dataclass(frozen=True, slots=True, repr=False)
class PhenotypeScore:
    ledger: Mapping[str, object] = field(repr=False)

    @property
    def total(self) -> float:
        return self.ledger["total_units"] / POLICY["scale"]

    def as_bytes(self) -> bytes:
        """Controlled/no-LLM JSON bytes; the caller owns storage and dispatch."""
        return _bytes(self.ledger)


def maximum_weight_matching(weights: Mapping[tuple[str, str], int]) -> tuple[tuple[str, str], ...]:
    """Exact rectangular Hungarian assignment with zero-weight unmatched dummies.

    Integer arithmetic has no epsilon or floating-point tie ambiguity. Input is
    a sparse edge map; absent edges have weight zero. O(L^2 * (R+L)) time.
    """
    if not isinstance(weights, Mapping) or any(
        not isinstance(pair, tuple) or len(pair) != 2
        or any(not isinstance(key, str) or not key for key in pair)
        or type(weight) is not int or weight < 0 for pair, weight in weights.items()
    ):
        raise PhenotypeScoringError("matching requires nonnegative integer edge weights")
    left = sorted({a for a, _ in weights})
    right = sorted({b for _, b in weights})
    n, real = len(left), len(right)
    m = real + n
    u, v, p, way = [0] * (n + 1), [0] * (m + 1), [0] * (m + 1), [0] * (m + 1)
    for i in range(1, n + 1):
        p[0], j0 = i, 0
        minimum, used = [None] * (m + 1), [False] * (m + 1)
        while True:
            used[j0] = True
            i0, delta, j1 = p[j0], None, 0
            for j in range(1, m + 1):
                if used[j]:
                    continue
                cost = -weights.get((left[i0 - 1], right[j - 1]), 0) if j <= real else 0
                reduced = cost - u[i0] - v[j]
                if minimum[j] is None or reduced < minimum[j]:
                    minimum[j], way[j] = reduced, j0
                if delta is None or minimum[j] < delta:
                    delta, j1 = minimum[j], j
            for j in range(m + 1):
                if used[j]:
                    u[p[j]] += delta
                    v[j] -= delta
                elif minimum[j] is not None:
                    minimum[j] -= delta
            j0 = j1
            if p[j0] == 0:
                break
        while j0:
            j1 = way[j0]
            p[j0], j0 = p[j1], j1
    return tuple(sorted((left[p[j] - 1], right[j - 1]) for j in range(1, real + 1)
                        if p[j] and weights.get((left[p[j] - 1], right[j - 1]), 0) > 0))


def _allocate(values: Mapping[str, Fraction], target: int | None = None) -> dict[str, int]:
    total = sum(values.values(), Fraction(0))
    if target is None:
        target, scaled = int(total), values
    else:
        scaled = {key: value * target / total if total else Fraction(0) for key, value in values.items()}
    result = {key: int(value) for key, value in scaled.items()}
    order = sorted(scaled, key=lambda key: (-(scaled[key] - result[key]), key))
    for key in order[:target - sum(result.values())]:
        result[key] += 1
    return result


def _frequency(fields):
    raw_id, raw_name = fields["frequency_id"], fields["frequency_name"]
    identifier = raw_id.strip() if raw_id else None
    bounds, fraction, exact_lower, reason = None, None, None, "missing_neutral"
    if identifier in POLICY["frequency_bounds"]:
        bounds, reason = POLICY["frequency_bounds"][identifier], "hpo_category_bounds"
        exact_lower = Fraction(str(bounds[0]))
    elif identifier and re.fullmatch(r"[0-9]+/[0-9]+", identifier):
        numerator, denominator = map(int, identifier.split("/"))
        if denominator == 0 or numerator > denominator:
            raise PhenotypeScoringError("profile frequency fraction is invalid")
        fraction = {"numerator": numerator, "denominator": denominator}
        value = Fraction(numerator, denominator)
        exact_lower = value
        bounds, reason = (float(value), float(value)), "source_fraction"
    elif identifier and re.fullmatch(r"(?:[0-9]+(?:\.[0-9]+)?|\.[0-9]+)%", identifier):
        value = Fraction(identifier[:-1]) / 100
        if not 0 <= value <= 1:
            raise PhenotypeScoringError("profile frequency percentage is invalid")
        exact_lower = value
        bounds, reason = (float(value), float(value)), "source_percentage"
    elif identifier:
        reason = "unsupported_source_value_neutral"
        if "/" in identifier or "%" in identifier:
            raise PhenotypeScoringError("profile frequency syntax is invalid")
    name_id = POLICY["frequency_names"].get(raw_name)
    # Names are source evidence too: retain a conflicting recognized name as a
    # bound alternative, rather than silently preferring the numeric ID.
    alternatives = [bounds] if bounds is not None else []
    if name_id:
        alternatives.append(POLICY["frequency_bounds"][name_id])
    lower_values = [exact_lower] if exact_lower is not None else []
    if name_id:
        lower_values.append(Fraction(str(POLICY["frequency_bounds"][name_id][0])))
    lower = min(lower_values) if lower_values else None
    if lower is not None:
        positive, negative = ((1 + lower) / 2 if lower else Fraction(0)), lower
    else:
        positive = negative = Fraction(1)
    return {
        "source_id": raw_id, "source_name": raw_name, "bounds": bounds,
        "source_fraction": fraction, "name_bounds": POLICY["frequency_bounds"].get(name_id),
        "reason": reason, "bound_conflict": len({tuple(item) for item in alternatives}) > 1,
        "lower_bound_used": {"numerator": lower.numerator, "denominator": lower.denominator} if lower is not None else None,
        "positive_factor": float(positive), "negative_factor": float(negative),
        "policy_sha256": POLICY_SHA256,
    }, positive, negative


def _context(observation, fields):
    source_onset, source_sex, modifiers = fields["onset"], fields["sex"], fields["modifiers"]
    windows = [(term, POLICY["onset_windows"].get(term)) for term in sorted(set(source_onset or ()), key=lambda item: item or "")]
    observed = POLICY["onset_windows"].get(observation.onset)
    positive, negative = Fraction(1), Fraction(1)
    reasons = []
    if windows:
        known_windows = [(term, window) for term, window in windows if window is not None]
        if len(known_windows) != len(windows) or observation.age_years is None:
            negative = Fraction(0)
            reasons.append("negative_onset_context_unknown")
        else:
            age = observation.age_years
            factors = [0 if low is not None and age < low else 1 if high is not None and age >= high else 0.5
                       for _, (low, high) in windows]
            negative = Fraction(str(min(factors)))
            reasons.append("negative_onset_all_alternatives")
        if known_windows and observed is not None:
            a, b = observed
            compatible = []
            for _, (low, high) in known_windows:
                # Congenital onset is the point zero; other windows are half-open.
                if a == b:
                    overlap = (low is None or a >= low) and (high is None or a < high or low == high == a)
                elif low == high:
                    overlap = (a is None or low >= a) and (b is None or low < b)
                else:
                    overlap = (b is None or low is None or low < b) and (high is None or a is None or a < high)
                compatible.append(overlap)
            positive = Fraction(1) if all(compatible) else Fraction(1, 2)
            reasons.append("positive_onset_known_constraints")
        else:
            reasons.append("positive_onset_unknown_neutral")
        if observation.age_years is not None and any(
            low is not None and observation.age_years < low for _, (low, _) in known_windows
        ):
            positive = Fraction(1, 2)
            reasons.append("present_before_sourced_minimum_onset_age")
    else:
        reasons.append("no_sourced_onset_requirement_neutral")
    normalized_sex = POLICY["sex_mapping"].get(source_sex.upper()) if source_sex else None
    sex_applicability = "unknown"
    if source_sex:
        if normalized_sex is None or observation.sex is None:
            negative = Fraction(0)
            reasons.append("negative_sex_context_unknown")
        elif normalized_sex != observation.sex:
            sex_applicability = "inapplicable"
            positive = negative = Fraction(0)
            reasons.append("sex_inapplicable")
        else:
            sex_applicability = "applicable"
            reasons.append("sex_compatible")
    else:
        reasons.append("no_sourced_sex_requirement_neutral")
    if modifiers:
        negative = Fraction(0)
        reasons.append("uninterpreted_modifier_blocks_context_penalty")
    return {
        "source_onset": source_onset, "mapped_windows": windows,
        "source_sex": source_sex, "normalized_sex": normalized_sex, "source_modifiers": modifiers,
        "sex_applicability": sex_applicability,
        "observed_onset_window": observed, "reasons": reasons,
        "positive_factor": float(positive), "negative_factor": float(negative),
        "policy_sha256": POLICY_SHA256,
    }, positive, negative


def score_phenotypes(
    observations: Sequence[ScoringObservation], disease_id: str, *,
    ic_store: HpoICStore, profile_store: DiseaseProfileStore | None,
) -> PhenotypeScore:
    """Score one existing candidate without candidate generation or gene priors.

    Uses task 7.1/7.2 stores, checks their graph binding, and retains all relevant
    source groups/annotations. Both HPO source roles require matching complete
    pins unless BOTH supplied stores explicitly declare synthetic=True. An
    explicitly unavailable profile store or a missing disease profile yields
    null raw evidence and a declared neutral total, not exclusion. Returned
    evidence is always controlled and no-LLM, even when callers use independently
    synthetic inputs.
    """
    return _score_phenotypes(observations, disease_id, ic_store=ic_store, profile_store=profile_store)


_SCORING_RESOURCE_SEAL = object()


@dataclass(frozen=True, slots=True, init=False, repr=False)
class _ScoringResourceBinding:
    ic_store: HpoICStore
    profile_store: DiseaseProfileStore | None
    graph_sha256: str
    ic_store_sha256: str
    annotation_index: Mapping
    _seal: object
    _owners: tuple[object, ...]


def _deeply_immutable(value):
    if value is None or type(value) in (str, int, float, bool):
        return True
    if type(value) in (tuple, frozenset):
        return all(_deeply_immutable(child) for child in value)
    if type(value) is MappingProxyType:
        return all(_deeply_immutable(key) and _deeply_immutable(child) for key, child in value.items())
    if is_dataclass(value) and value.__dataclass_params__.frozen:
        return all(_deeply_immutable(getattr(value, item.name)) for item in fields(value))
    return False


def _resource_owners(ic_store, profile_store):
    owners = (ic_store, ic_store.graph, ic_store.graph.ontology)
    if profile_store is not None:
        owners += (profile_store,)
    return (*owners, *(getattr(owner, item.name) for owner in owners for item in fields(owner)))


def _scoring_graph_sha256(ic_store):
    return _digest({"policy": OBSOLETE_RESOLUTION_POLICY,
                    "terms": [asdict(term) for _, term in sorted(ic_store.graph.ontology.terms.items())]})


_PREPARED_SCORES_SEAL = object()


@dataclass(frozen=True, slots=True, repr=False)
class PreparedPhenotypeScores:
    """Caller-owned, bounded score cache for one observation/resource/prior view.

    The producer computes complete metadata and observation preparation once.
    Consumers use a 256-disease LRU; no process-global scores or parent reuse.
    """
    observations: tuple[ScoringObservation, ...]
    ic_store: HpoICStore
    profile_store: DiseaseProfileStore | None
    prior_identity: str
    key: tuple[str, ...]
    resource_binding: _ScoringResourceBinding
    observation_owners: tuple[object, ...]
    observation_context: tuple
    scores: OrderedDict
    seal: object

    def score(self, disease_id: str) -> PhenotypeScore:
        owners = tuple(owner for item in self.observations
                       for owner in (item, *(getattr(item, part.name) for part in fields(item))))
        resource_owners = (self.resource_binding.graph_sha256,
                           self.resource_binding.ic_store_sha256, self.resource_binding.annotation_index,
                           *_resource_owners(self.ic_store, self.profile_store))
        if (self.seal is not _PREPARED_SCORES_SEAL
                or len(owners) != len(self.observation_owners)
                or any(a is not b for a, b in zip(owners, self.observation_owners, strict=True))
                or len(resource_owners) != len(self.resource_binding._owners)
                or any(a is not b for a, b in zip(resource_owners, self.resource_binding._owners, strict=True))):
            raise PhenotypeScoringError("prepared score cache parent mismatch")
        if disease_id in self.scores:
            value = self.scores.pop(disease_id)
            self.scores[disease_id] = value
            return value
        value = _score_phenotypes(
            self.observations, disease_id, ic_store=self.ic_store,
            profile_store=self.profile_store,
            _resource_binding=self.resource_binding,
            _observation_context=self.observation_context,
        )
        self.scores[disease_id] = value
        if len(self.scores) > 256:
            self.scores.popitem(last=False)
        return value


def prepare_phenotype_scores(
    observations, *, ic_store: HpoICStore,
    profile_store: DiseaseProfileStore | None, prior_identity: str,
    _skip_verification=False, _ic_store_sha256=None,
) -> PreparedPhenotypeScores:
    """Create a case/view-scoped cache; never accept mutable resource owners."""
    observations = tuple(observations)
    if (not prior_identity or type(prior_identity) is not str
            or any(type(item) is not ScoringObservation for item in observations)
            or not isinstance(ic_store, HpoICStore) or not ic_store.ready
            or not _skip_verification and (
                not _deeply_immutable(ic_store) or not _deeply_immutable(profile_store))):
        raise PhenotypeScoringError("prepared scoring requires immutable bound resources")
    binding = object.__new__(_ScoringResourceBinding)
    graph_sha256 = (profile_store.graph_sha256 if _skip_verification and profile_store is not None
        else _scoring_graph_sha256(ic_store))
    ic_store_sha256 = (_ic_store_sha256 if _skip_verification and _ic_store_sha256 is not None
        else hashlib.sha256(hpo_ic_as_bytes(ic_store)).hexdigest())
    annotation_index = MappingProxyType(_profile_annotations(profile_store, None))
    for name, value in (
        ("ic_store", ic_store), ("profile_store", profile_store),
        ("graph_sha256", graph_sha256), ("ic_store_sha256", ic_store_sha256),
        ("annotation_index", annotation_index),
        ("_seal", _SCORING_RESOURCE_SEAL),
        ("_owners", (graph_sha256, ic_store_sha256, annotation_index, *_resource_owners(ic_store, profile_store))),
    ):
        object.__setattr__(binding, name, value)
    key = (
        _digest([asdict(item) for item in observations]), POLICY_SHA256,
        graph_sha256, ic_store_sha256,
        "unavailable" if profile_store is None else "live-profile-verification-omitted" if _skip_verification
        else hashlib.sha256(disease_profiles_as_bytes(profile_store)).hexdigest(),
        prior_identity,
    )
    return PreparedPhenotypeScores(
        observations, ic_store, profile_store, prior_identity, key, binding,
        tuple(owner for item in observations
              for owner in (item, *(getattr(item, part.name) for part in fields(item)))),
        _prepare_observation_context(observations, ic_store.graph),
        OrderedDict(), _PREPARED_SCORES_SEAL,
    )


def _score_phenotypes_batch(observations, disease_ids, *, ic_store, profile_store,
                            prepared: PreparedPhenotypeScores | None = None):
    """Reuse complete preparation only under an exact case/view binding."""
    disease_ids = tuple(disease_ids)
    if not disease_ids:
        return ()
    if prepared is not None:
        if (type(prepared) is not PreparedPhenotypeScores
                or prepared.ic_store is not ic_store
                or prepared.profile_store is not profile_store
                or prepared.observations != tuple(observations)):
            raise PhenotypeScoringError("prepared score cache parent mismatch")
        return tuple(prepared.score(disease) for disease in disease_ids)
    if (not isinstance(ic_store, HpoICStore) or not ic_store.ready
            or not _deeply_immutable(ic_store) or not _deeply_immutable(profile_store)):
        return tuple(score_phenotypes(observations, disease, ic_store=ic_store, profile_store=profile_store)
                     for disease in disease_ids)
    batch = prepare_phenotype_scores(
        observations, ic_store=ic_store, profile_store=profile_store,
        prior_identity="diagnosis-blind-primary",
    )
    return tuple(batch.score(disease) for disease in disease_ids)


def _prepare_observation_context(observations, graph):
    """Resolve one immutable observation view for all diseases in its cache."""
    observation_rows, units, variants = {}, {}, {}
    growth = set(POLICY["growth_terms"])
    for item in sorted(observations, key=lambda item: item.observation_id):
        row = asdict(item)
        if item.observation_id in observation_rows and _bytes(observation_rows[item.observation_id]["input"]) != _bytes(row):
            raise PhenotypeScoringError("observation ID is reused with conflicting content")
        resolution = graph.resolve(item.source_hpo_id)
        term = resolution.canonical_id or item.source_hpo_id
        ancestors = graph.ancestors.get(term, ())
        is_growth = item.domain == "growth" or term in growth or bool(set(POLICY["growth_roots"]) & set(ancestors))
        domain = item.domain or ("growth" if is_growth else None)
        semantic = {key: value for key, value in row.items() if key not in ("observation_id", "provenance", "source_hpo_id")}
        semantic.update(canonical_term=term, domain=domain)
        for numeric in ("age_years", "quality"):
            if semantic[numeric] is not None:
                semantic[numeric] = float(semantic[numeric]) or 0.0
        key = _digest(semantic)
        if key not in variants:
            # Sorted IDs choose stable raw provenance; scoring uses only the
            # canonical numeric context, never the first raw numeric spelling.
            normalized = replace(item, age_years=semantic["age_years"], quality=semantic["quality"])
            variants[key] = {"observation": normalized, "ids": set(), "domain": domain, "growth": is_growth, "semantic": semantic}
        variants[key]["ids"].add(item.observation_id)
        units.setdefault(term, set()).add(key)
        observation_rows[item.observation_id] = {"input": row, "resolution": asdict(resolution), "variant_id": key}

    conflicts = []
    ordered = sorted(observation_rows)
    for index, left in enumerate(ordered):
        a = observation_rows[left]
        for right in ordered[index + 1:]:
            b = observation_rows[right]
            x, y = a["input"], b["input"]
            if x["status"] == y["status"]:
                continue
            positive, absent = (a, b) if x["status"] == "present" else (b, a)
            pterm, nterm = positive["resolution"]["canonical_id"], absent["resolution"]["canonical_id"]
            related = nterm in graph.ancestors.get(pterm, ()) if pterm and nterm else x["source_hpo_id"] == y["source_hpo_id"]
            if not related:
                continue
            same_subject = (x["subject_id"] == y["subject_id"] if x["subject_id"] and y["subject_id"]
                            else x["scope"] != "family_history" and y["scope"] != "family_history")
            family_unknown = x["scope"] == y["scope"] == "family_history" and (not x["subject_id"] or not y["subject_id"])
            if same_subject or family_unknown:
                conflicts.append({"observation_ids": [left, right], "status": "unresolved",
                                  "subject_relation": "same_subject" if same_subject else "unknown_family_subject",
                                  "treatment": "retain_both_signed_allocations_no_contextual_adjudication"})

    return (
        MappingProxyType({key: _freeze(value) for key, value in observation_rows.items()}),
        MappingProxyType({key: frozenset(value) for key, value in units.items()}),
        MappingProxyType({key: MappingProxyType({**value, "ids": frozenset(value["ids"])})
                          for key, value in variants.items()}),
        tuple(_freeze(value) for value in conflicts), tuple(ordered),
    )


def _profile_annotations(profile_store, resource_binding):
    if resource_binding is not None:
        return resource_binding.annotation_index
    rows = () if profile_store is None else profile_store.annotations
    annotations = {item.annotation_id: item for item in rows}
    if len(annotations) != len(rows):
        raise PhenotypeScoringError("profile annotation IDs are duplicated")
    return annotations


def _score_phenotypes(
    observations, disease_id, *, ic_store, profile_store, _resource_binding=None,
    _observation_context=None,
):
    if (not isinstance(observations, Sequence) or isinstance(observations, (str, bytes))
            or any(not isinstance(item, ScoringObservation) for item in observations)
            or not is_disease_identifier(disease_id, allow_native_token=True)
            or not isinstance(ic_store, HpoICStore) or not ic_store.ready
            or profile_store is not None and not isinstance(profile_store, DiseaseProfileStore)):
        raise PhenotypeScoringError("scoring requires structured observations and ready source stores")
    graph = ic_store.graph
    if _resource_binding is None:
        graph_hash = _scoring_graph_sha256(ic_store)
        ic_hash = None
    else:
        if (type(_resource_binding) is not _ScoringResourceBinding
                or getattr(_resource_binding, "_seal", None) is not _SCORING_RESOURCE_SEAL
                or _resource_binding.ic_store is not ic_store
                or _resource_binding.profile_store is not profile_store):
            raise PhenotypeScoringError("prepared scoring resource binding mismatch")
        expected = (_resource_binding.graph_sha256, _resource_binding.ic_store_sha256, _resource_binding.annotation_index,
                    *_resource_owners(ic_store, profile_store))
        if (len(_resource_binding._owners) != len(expected)
                or any(actual is not original for actual, original in
                       zip(_resource_binding._owners, expected, strict=True))):
            raise PhenotypeScoringError("prepared scoring resource binding mismatch")
        graph_hash = _resource_binding.graph_sha256
        ic_hash = _resource_binding.ic_store_sha256
    if profile_store is not None and graph_hash != profile_store.graph_sha256:
        raise PhenotypeScoringError("IC and profile graph provenance differ")
    if (not isinstance(ic_store.source_binding, Mapping) or not ic_store.source_binding
            or profile_store is not None and (
                not isinstance(profile_store.source_binding, Mapping) or not profile_store.source_binding)):
        raise PhenotypeScoringError("source bindings are required")
    joint_synthetic = (ic_store.source_binding.get("synthetic") is True
                       and (profile_store is None or profile_store.source_binding.get("synthetic") is True))
    if not joint_synthetic:
        required_roles = {"hpo_ontology": "hpo_ontology", "hpo_disease_annotations": "hpo_annotations"}
        sources = ic_store.source_binding.get("sources")
        if not isinstance(sources, (list, tuple)) or len(sources) != len(required_roles):
            raise PhenotypeScoringError("both IC HPO source roles are required exactly once")
        seen_roles = set()
        for source in sources:
            if not isinstance(source, Mapping):
                raise PhenotypeScoringError("IC source binding is malformed")
            role = source.get("role")
            if not isinstance(role, str) or role not in required_roles or role in seen_roles:
                raise PhenotypeScoringError("IC HPO source role is missing, unknown, or duplicated")
            seen_roles.add(role)
            bound = None if profile_store is None else profile_store.source_binding.get(required_roles[role])
            for pin in (source,) if profile_store is None else (source, bound):
                if (not isinstance(pin, Mapping) or pin.get("role") != role
                        or not isinstance(pin.get("sha256"), str) or not _SHA.fullmatch(pin["sha256"])
                        or not isinstance(pin.get("release"), str) or not pin["release"].strip()):
                    raise PhenotypeScoringError("HPO source pin requires its role, SHA-256, and nonempty release")
            if profile_store is not None and any(source[key] != bound[key] for key in ("sha256", "release")):
                raise PhenotypeScoringError("IC and profile HPO source pins differ")
    if _observation_context is None:
        _observation_context = _prepare_observation_context(observations, graph)
    observation_rows, units, variants, conflicts, ordered = _observation_context

    profile = None if profile_store is None else profile_store.profiles.get(disease_id)
    annotations = _profile_annotations(profile_store, _resource_binding)
    terms, groups, evidence = {}, [], {}
    for group in sorted(profile.groups if profile else (), key=lambda item: item.group_id):
        groups.append(asdict(group))
        if not group.annotation_ids or any(key not in annotations for key in group.annotation_ids):
            raise PhenotypeScoringError("profile group has unresolved evidence references")
        for key in group.annotation_ids:
            annotation = annotations[key]
            resolved = graph.resolve(annotation.normalized_source_hpo_id) if annotation.normalized_source_hpo_id else None
            if resolved != annotation.resolution or (resolved.canonical_id if resolved else None) != group.canonical_hpo_id:
                raise PhenotypeScoringError("profile term resolution does not match its graph")
            if annotation.source_fields["aspect"] != group.aspect:
                raise PhenotypeScoringError("profile group aspect differs from source")
            evidence[key] = {
                "annotation_id": key, "source_fields": annotation.source_fields,
                "resolution": asdict(resolved) if resolved else None,
                "polarity": annotation.polarity, "diagnostics": annotation.diagnostics,
            }
        if group.aspect == "P" and group.canonical_hpo_id is not None:
            terms.setdefault(group.canonical_hpo_id, set()).update(group.annotation_ids)

    edges, weights = [], {}
    for term, keys in sorted(units.items()):
        for target, members in sorted(terms.items()):
            resnik = ic_store.resnik(term, target)
            ic = resnik.information_content
            if ic is not None and (not math.isfinite(ic) or ic < 0):
                raise PhenotypeScoringError("IC evidence is invalid")
            semantic = Fraction(str(ic)) / (1 + Fraction(str(ic))) if ic is not None else Fraction(0)
            allocations, exact_values = [], {}
            for key in sorted(keys):
                variant = variants[key]
                item = variant["observation"]
                negative = item.status == "explicitly_absent"
                entailed = term in graph.ancestors[target]
                effects, factors = [], []
                for member in sorted(members):
                    annotation = annotations[member]
                    fields = annotation.source_fields
                    frequency, fp, fn = _frequency(fields)
                    context, cp, cn = _context(item, fields)
                    included = context["sex_applicability"] != "inapplicable"
                    quality = Fraction(1, 2) if fields["is_ontology_propagated"] is True else Fraction(1)
                    factor = (fn * cn if negative else fp * cp) * quality
                    reasons = []
                    if annotation.polarity != "positive":
                        factor = Fraction(0)
                        reasons.append("source_not_unambiguously_positive")
                    if negative and not entailed:
                        factor = Fraction(0)
                        reasons.append("negative_not_entailed_by_profile_term")
                    effects.append({"annotation_id": member, "frequency": frequency, "context": context,
                                    "applicability": context["sex_applicability"], "included_in_aggregation": included,
                                    "source_quality_raw": fields["is_ontology_propagated"], "source_quality_factor": float(quality),
                                    "factor": float(factor), "reasons": reasons})
                    if included:
                        factors.append(factor)
                factor = min(factors, default=Fraction(0))
                quality = Fraction(str(item.quality)) if item.quality is not None else Fraction(1)
                magnitude = semantic * factor * quality * (Fraction(1, 2) if negative else 1) / len(keys)
                exact_values[key] = magnitude * POLICY["scale"]
                allocations.append({
                    "variant_id": key, "observation_ids": sorted(variant["ids"]),
                    "scope": item.scope, "domain": variant["domain"], "status": item.status,
                    "kind": "penalty" if negative else "support", "contextual": item.scope != "proband",
                    "source_effects": effects, "source_factor": float(factor),
                    "observation_quality_raw": observation_rows[item.observation_id]["input"]["quality"],
                    "observation_quality_factor": float(quality),
                    "shared_term_budget_divisor": len(keys), "negative_weight": 0.5 if negative else None,
                    "exact_units": {"numerator": exact_values[key].numerator, "denominator": exact_values[key].denominator},
                })
            rounded = _allocate(exact_values)
            for allocation in allocations:
                allocation["units"] = rounded[allocation["variant_id"]]
            weight = sum(rounded.values())
            weights[term, target] = weight
            edges.append({"observation_term": term, "profile_term": target, "resnik": asdict(resnik),
                          "semantic_factor": float(semantic), "allocations": allocations, "weight_units": weight})
    matching = maximum_weight_matching(weights)
    selected = set(matching)
    contributions = []
    for edge in edges:
        pair = edge["observation_term"], edge["profile_term"]
        edge["selected"] = pair in selected
        edge["decision"] = "assigned" if pair in selected else "zero_weight" if edge["weight_units"] == 0 else "one_to_one_competition"
        if pair not in selected:
            continue
        for allocation in edge["allocations"]:
            key = allocation["variant_id"]
            contributions.append({"contribution_id": _digest((pair, key)), "observation_term": pair[0],
                                  "profile_term": pair[1], "variant_id": key, "scope": allocation["scope"],
                                  "domain": allocation["domain"], "kind": allocation["kind"],
                                  "before_caps_units": allocation["units"], "units": allocation["units"],
                                  "cap_ids": []})

    # Each cap acts only on its own members. Overlap never merges growth,
    # perinatal, or independent domains into one transitive constellation cap.
    cap_sets = {}
    source_members = {}
    for contribution in contributions:
        key, kind = contribution["variant_id"], contribution["kind"]
        variant = variants[key]
        item = variant["observation"]
        native_records = {tuple(annotations[member].source_fields[name] for name in
                                ("source_family", "source_dataset", "source_release", "source_record_id"))
                          for member in terms[contribution["profile_term"]]}
        for record in native_records:
            source_members.setdefault((record, kind), []).append(contribution)
        if item.dependency_cluster:
            cap_sets.setdefault(("1_dependency", variant["domain"], item.dependency_cluster, kind), []).append(contribution)
        if variant["growth"]:
            cap_sets.setdefault(("2_growth", kind), []).append(contribution)
        if item.scope == "perinatal":
            cap_sets.setdefault(("3_perinatal", kind), []).append(contribution)
    for (record, kind), members in source_members.items():
        if len({item["profile_term"] for item in members}) > 1:
            cap_sets[("0_source_reuse", *record, kind)] = members
    caps = []
    for cap_key, members in sorted(cap_sets.items()):
        limit = POLICY["dependency_cap_units"] if cap_key[0] == "1_dependency" else POLICY["growth_cap_units"] if cap_key[0] == "2_growth" else POLICY["perinatal_cap_units"]
        if cap_key[0] == "0_source_reuse":
            limit = max(sum(item["before_caps_units"] for item in members if item["profile_term"] == term)
                        for term in {item["profile_term"] for item in members})
        before = {item["contribution_id"]: item["units"] for item in members}
        after = _allocate({key: Fraction(value) for key, value in before.items()}, min(sum(before.values()), limit))
        cap_id = _digest(cap_key)
        caps.append({"cap_id": cap_id, "key": cap_key, "limit_units": limit, "before_units": before,
                     "after_units": after, "removed_units": sum(before.values()) - sum(after.values())})
        for item in members:
            item["units"] = after[item["contribution_id"]]
            item["cap_ids"].append(cap_id)
    support = sum(item["units"] for item in contributions if item["kind"] == "support")
    penalty = sum(item["units"] for item in contributions if item["kind"] == "penalty")
    subtotals = []
    for scope, domain in sorted({(item.scope, variants[observation_rows[item.observation_id]["variant_id"]]["domain"])
                                 for item in observations}, key=lambda pair: (pair[0], pair[1] or "")):
        members = [item for item in contributions if (item["scope"], item["domain"]) == (scope, domain)]
        positive = sum(item["units"] for item in members if item["kind"] == "support")
        negative = sum(item["units"] for item in members if item["kind"] == "penalty")
        subtotals.append({"scope": scope, "domain": domain, "support_units": positive,
                          "penalty_units": negative, "total_units": positive - negative})
    for item in contributions:
        item["percentage_of_positive_support"] = 100 * item["units"] / support if support and item["kind"] == "support" else None
    result = {
        "schema": SCHEMA, "classification": "controlled", "llm_policy": "no-llm",
        "disease_id": disease_id, "policy": POLICY, "policy_sha256": POLICY_SHA256,
        "source_binding": {"ic": ic_store.source_binding, "ic_store_sha256": ic_hash if ic_hash is not None else hashlib.sha256(hpo_ic_as_bytes(ic_store)).hexdigest(),
                           "profiles": None if profile_store is None else profile_store.source_binding,
                           "graph_sha256": graph_hash,
                           "candidate_profile_sha256": _digest({"groups": groups, "annotations": evidence})},
        "observations": [observation_rows[key] for key in ordered],
        "observation_binding_sha256": _digest([observation_rows[key] for key in ordered]),
        "observation_variants": [{"variant_id": key, "context": variants[key]["semantic"],
                                  "representative_observation_id": variants[key]["observation"].observation_id,
                                  "observation_ids": sorted(variants[key]["ids"])} for key in sorted(variants)],
        "contradictions": conflicts, "profile_groups": groups,
        "profile_annotations": [evidence[key] for key in sorted(evidence)],
        "profile_status": "unavailable" if profile is None else "no_resolved_phenotype_terms" if not terms else "available",
        "raw_evidence_units": support - penalty if any(edge["resnik"]["information_content"] is not None for edge in edges) else None,
        "missing_evidence_neutral_units": 0, "unlisted_features": "unknown",
        "edges": edges, "matching": matching, "matching_objective_units": sum(weights[pair] for pair in matching),
        "contributions": sorted(contributions, key=lambda item: item["contribution_id"]),
        "caps": caps, "subtotals": subtotals,
        "support_units": support, "penalty_units": penalty, "total_units": support - penalty,
        "inferred_inheritance": None, "inferred_parental_disease": None, "inferred_genomic_mechanism": None,
        "diagnostic_claim": None,
    }
    # Validate the complete result with the ordinary canonical JSON rules.
    # The C encoder traverses native dict/list/tuple/scalar values directly.
    # Only immutable or custom Mapping instances need a shallow conversion.
    # Keep the upstream observation serializer and its contracts unchanged.
    def mapping_default(value):
        if isinstance(value, Mapping):
            return dict(value.items())
        raise TypeError('unsupported phenotype result value')
    try:
        json.dumps(result, sort_keys=True, separators=(",", ":"), allow_nan=False,
                   default=mapping_default)
    except (ValueError, TypeError, RecursionError):
        raise PhenotypeScoringError("evidence must be finite JSON") from None
    def freeze_result(value):
        # Exact native leaves cannot implement Mapping or sequence behavior.
        # Preserve ordinary handling for subclasses and custom containers.
        kind = type(value)
        if kind in (str, int, float, bool, type(None)):
            return value
        if kind is dict or isinstance(value, Mapping):
            return MappingProxyType({key: freeze_result(child) for key, child in value.items()})
        if kind in (list, tuple) or isinstance(value, (list, tuple)):
            return tuple(freeze_result(child) for child in value)
        return value
    return PhenotypeScore(freeze_result(result))
