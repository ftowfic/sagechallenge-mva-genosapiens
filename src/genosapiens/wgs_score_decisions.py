"""Exact scoring dependencies and computed outcomes, without source duplication.

The source adapter must close atomic ownership before projecting these inputs.
All numeric operations call the ordinary fusion/phenotype kernels. Source audit,
hypothesis identity and explicit tie values remain in their authoritative owners.
"""
from collections import OrderedDict
from dataclasses import dataclass
import hashlib
import json

from genosapiens.evidence_fusion import (
    ChargeKind, ComponentName, FusionPolicy, _OWNER, _COMPONENTS, _canonical_evidence,
    _raw_component, component_numeric_value, generated_penalty_values, technical_penalty_values,
    total_score_value,
    policy_as_json,
)
from genosapiens.evidence_model import Completeness, HypothesisPhase, _exact, _sha, _typed
from genosapiens.evidence_serialization import _canonical
from genosapiens.phenotype_scoring import POLICY_SHA256 as PHENOTYPE_POLICY_SHA256
from genosapiens.wgs_evidence_store import EvidenceStoreError
from genosapiens.wgs_evidence_replay import EvidenceReplay

MAX_SCORE_INPUT_BYTES = 64 * 1024
MAX_SCORE_CLASSES = 1_000_000
HOT_SCORE_CLASSES = 4096


@dataclass(frozen=True, slots=True, repr=False)
class ScoreChargeValue:
    kind: ChargeKind
    value: str | None
    member_id: str | None

    @property
    def component(self):
        return _OWNER[self.kind]


@dataclass(frozen=True, slots=True, repr=False)
class ScoreDecisionInputs:
    member_count: int
    paired: bool
    completeness: Completeness
    phase: HypothesisPhase
    phenotype_units: int | None
    charges: tuple[ScoreChargeValue, ...]


@dataclass(frozen=True, slots=True, repr=False)
class ComponentValues:
    name: ComponentName
    raw: str | None
    normalized: str | None


@dataclass(frozen=True, slots=True, repr=False)
class ScoreDecision:
    components: tuple[ComponentValues, ...]
    total: str
    generated_penalties: tuple[tuple[str, str, str | None, str | None], ...]


def score_decision_inputs(value):
    """Project a fully canonical source input, preserving ordered member slots.

    Lexical member order retains the weakest-member identity tie. No evidence
    ownership or proof conflict is lost: the ordinary canonicalizer runs first.
    """
    value = _canonical_evidence(value)
    if value.phenotype_score.ledger["policy_sha256"] != PHENOTYPE_POLICY_SHA256:
        raise EvidenceStoreError("score phenotype policy binding changed")
    members = tuple(sorted(value.hypothesis.evidence_subject_ids))
    slots = {key: str(i) for i, key in enumerate(members)}
    return ScoreDecisionInputs(len(members), len(value.hypothesis.allele_refs) == 2,
        value.hypothesis.completeness, value.hypothesis.phase, value.phenotype_score.ledger["raw_evidence_units"],
        tuple(ScoreChargeValue(c.kind, c.value, None if c.member_id is None else slots[c.member_id]) for c in value.charges))


def evaluate_score_decision(inputs, policy):
    inputs = _typed(inputs, ScoreDecisionInputs)
    if inputs.member_count not in (1, 2):
        raise EvidenceStoreError("bounded score requires one or two evidence subjects")
    members = tuple(str(i) for i in range(inputs.member_count))
    for charge in inputs.charges:
        if (charge.member_id is not None and charge.member_id not in members
                or charge.kind in ("technical", "conflict") and charge.value is not None and _exact(charge.value) < 0):
            raise EvidenceStoreError("invalid score charge scope or value")
    phenotype_raw, phenotype_normalized = policy.phenotype.numeric_values(inputs.phenotype_units)
    values = {"phenotype": ComponentValues("phenotype",
        None if phenotype_raw is None else str(phenotype_raw),
        None if phenotype_normalized is None else str(phenotype_normalized))}
    member_values = {}
    for settings in policy.components:
        if settings.name == "technical_and_conflict_penalties":
            continue
        raw, member_raw, _ = _raw_component(settings, inputs.charges, members)
        normalized = component_numeric_value(settings, raw, member_raw=member_raw)
        values[settings.name] = ComponentValues(settings.name, None if raw is None else str(raw),
                                                None if normalized is None else str(normalized))
        member_values[settings.name] = member_raw
    molecular = next(c for c in policy.components if c.name == "molecular")
    technical = next(c for c in policy.components if c.name == "technical_and_conflict_penalties")
    base_raw, _, charges = _raw_component(technical, inputs.charges, members)
    generated = generated_penalty_values(completeness=inputs.completeness, phase=inputs.phase, paired=inputs.paired,
        policy=policy, molecular=molecular, member_raw=member_values["molecular"])
    penalty_raw, _ = technical_penalty_values(technical, members, base_raw, charges,
                                             tuple(amount for _, amount, *_ in generated))
    normalized = component_numeric_value(technical, penalty_raw)
    values[technical.name] = ComponentValues(technical.name, None if penalty_raw is None else str(penalty_raw),
                                            None if normalized is None else str(normalized))
    settings = {c.name: (c.weight, c.neutral_value) for c in policy.components}
    settings["phenotype"] = policy.phenotype.weight, policy.phenotype.neutral
    total = total_score_value((name, settings[name][0],
        values[name].normalized if values[name].normalized is not None else settings[name][1]) for name in _COMPONENTS)
    return ScoreDecision(tuple(values[name] for name in _COMPONENTS), str(total),
        tuple((kind, str(amount), member, None if weakest is None else str(weakest))
              for kind, amount, member, weakest in generated))


class ScoreDecisionCatalog:
    def __init__(self, writer, policy, *, policy_sha256, name):
        self.writer, self.connection, self.name = writer, writer.connection, name
        self.policy = _typed(policy, FusionPolicy)
        _sha(policy_sha256)
        if hashlib.sha256(policy_as_json(self.policy)).hexdigest() != policy_sha256:
            raise EvidenceStoreError("score decision policy digest changed")
        self.policy_sha256 = policy_sha256
        self.table = "score_decisions_" + hashlib.sha256(name.encode()).hexdigest()
        self.connection.execute(f"CREATE TABLE {self.table} (ordinal INTEGER PRIMARY KEY,digest TEXT UNIQUE,inputs TEXT,result TEXT)")
        self.count = self.evaluations = 0
        self.hot, self.closed, self.published = OrderedDict(), False, False

    def evaluate(self, inputs):
        if self.closed or self.published:
            raise EvidenceStoreError("score decision catalog is closed")
        inputs = _typed(inputs, ScoreDecisionInputs)
        content = _canonical(inputs)
        if len(content) > MAX_SCORE_INPUT_BYTES:
            raise EvidenceStoreError("score decision inputs exceed capacity", capacity=True)
        digest = hashlib.sha256(content.encode("ascii")).hexdigest()
        if digest in self.hot:
            return self.hot[digest]
        row = self.connection.execute(f"SELECT ordinal,inputs FROM {self.table} WHERE digest=?", (digest,)).fetchone()
        if row is None:
            if self.count >= MAX_SCORE_CLASSES:
                raise EvidenceStoreError("score decision class count exceeds capacity", capacity=True)
            result = evaluate_score_decision(inputs, self.policy)
            ordinal = self.count
            self.connection.execute(f"INSERT INTO {self.table} VALUES (?,?,?,?)", (ordinal,digest,content,_canonical(result)))
            self.count += 1
            self.evaluations += 1
        else:
            ordinal, previous = row
            if previous != content:
                raise EvidenceStoreError("score decision input digest collision")
        self.hot[digest] = ordinal
        if len(self.hot) > HOT_SCORE_CLASSES:
            self.hot.popitem(last=False)
        return ordinal

    def result(self, ordinal):
        if self.closed or type(ordinal) is not int or not 0 <= ordinal < self.count:
            raise EvidenceStoreError("score decision ordinal is unavailable")
        content, = self.connection.execute(f"SELECT result FROM {self.table} WHERE ordinal=?", (ordinal,)).fetchone()
        return _typed(json.loads(content), ScoreDecision, decode=True)

    def records(self):
        if self.closed:
            raise EvidenceStoreError("score decision catalog is closed")
        for ordinal,inputs,result in self.connection.execute(
                f"SELECT ordinal,inputs,result FROM {self.table} ORDER BY ordinal"):
            yield ordinal,_typed(json.loads(inputs),ScoreDecisionInputs,decode=True),_typed(json.loads(result),ScoreDecision,decode=True)

    def ordinal_for_inputs(self, inputs):
        if self.closed:
            raise EvidenceStoreError("score decision catalog is closed")
        content = _canonical(_typed(inputs, ScoreDecisionInputs))
        if len(content) > MAX_SCORE_INPUT_BYTES:
            raise EvidenceStoreError("score decision inputs exceed capacity", capacity=True)
        digest = hashlib.sha256(content.encode("ascii")).hexdigest()
        row = self.connection.execute(f"SELECT ordinal,inputs FROM {self.table} WHERE digest=?", (digest,)).fetchone()
        if row is None or row[1] != content:
            raise EvidenceStoreError("missing or foreign computed score inputs")
        return row[0]

    def publish(self):
        if self.closed or self.published:
            raise EvidenceStoreError("score decision catalog is closed")
        rows_sha = self.writer.add_sequence(self.name + ".outcomes", "score-computed-outcomes-v1", (
            {"ordinal": i, "inputs_sha256": digest, "inputs": json.loads(inputs), "result": json.loads(result)}
            for i, digest, inputs, result in self.connection.execute(f"SELECT * FROM {self.table} ORDER BY ordinal")))
        self.sha256 = self.writer.add(self.name, "score-decision-catalog-v1", {
            "schema": "genosapiens.wgs-score-decision-catalog", "version": 1,
            "policy": self.policy, "policy_sha256": self.policy_sha256, "outcome_count": self.count,
            "outcomes": {"name": self.name + ".outcomes", "sha256": rows_sha},
        })
        self.published = True
        return self.sha256

    def close(self):
        if not self.closed:
            self.connection.execute("DROP TABLE " + self.table)
        self.hot.clear()
        self.closed = True

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()


class VerifiedScoreDecisions(ScoreDecisionCatalog):
    def __init__(self, reader, reference, *, policy, policy_sha256):
        if reference.kind != "score-decision-catalog-v1":
            raise EvidenceStoreError("foreign score decision root")
        header = reader.reconstruct(reference)
        if (set(header) != {"schema", "version", "policy", "policy_sha256", "outcome_count", "outcomes"}
                or header["schema"] != "genosapiens.wgs-score-decision-catalog"
                or type(header["version"]) is not int or header["version"] != 1
                or header["policy_sha256"] != policy_sha256
                or _typed(header["policy"], FusionPolicy, decode=True) != policy
                or type(header["outcome_count"]) is not int or not 0 <= header["outcome_count"] <= MAX_SCORE_CLASSES
                or set(header["outcomes"]) != {"name", "sha256"}
                or header["outcomes"]["name"] != reference.name + ".outcomes"):
            raise EvidenceStoreError("invalid score decision header")
        child = reader.reference(header["outcomes"]["name"], kind="score-computed-outcomes-v1")
        if child.sha256 != header["outcomes"]["sha256"]:
            raise EvidenceStoreError("score decision child binding changed")
        super().__init__(EvidenceReplay(reader), policy, policy_sha256=policy_sha256, name=reference.name)
        try:
            count = 0
            for row in reader.iter_array(child, max_item_bytes=2 * MAX_SCORE_INPUT_BYTES):
                if (set(row) != {"ordinal", "inputs_sha256", "inputs", "result"}
                        or type(row["ordinal"]) is not int or row["ordinal"] != count):
                    raise EvidenceStoreError("score outcome coverage changed")
                inputs = _typed(row["inputs"], ScoreDecisionInputs, decode=True)
                if hashlib.sha256(_canonical(inputs).encode("ascii")).hexdigest() != row["inputs_sha256"]:
                    raise EvidenceStoreError("score input digest changed")
                ordinal = self.evaluate(inputs)
                if ordinal != count or self.result(ordinal) != _typed(row["result"], ScoreDecision, decode=True):
                    raise EvidenceStoreError("computed score outcome does not replay")
                count += 1
            if count != header["outcome_count"]:
                raise EvidenceStoreError("score outcome count changed")
            self.publish()
        except BaseException:
            self.close()
            raise
