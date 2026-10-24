"""Disk-backed equivalence of explicit inputs to the existing pair kernels.

Every new input class is evaluated. Cache hits have identical complete decision
dependencies, not merely matching outputs, coordinates or consequence labels.
Identity/source/phase-basis records stay in the context owner and must be joined
back for reconstruction; they are not replaced by another member's evidence.
"""
from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass, replace
import hashlib
import json

from genosapiens.evidence_model import Completeness, HypothesisContext, HypothesisPhase, _typed
from genosapiens.evidence_serialization import _canonical
from genosapiens.genotype_model import (
    ModelDecisionMember, ModelRequirementDecision, evaluate_model_requirement, summarize_model_requirements,
)
from genosapiens.partner_qualification import (
    QualificationDecisionInputs, evaluate_qualification_decision,
)
from genosapiens.genotype_hypotheses import qualify_whole_model, construction_carriage_supported
from genosapiens.wgs_evidence_store import EvidenceReader, EvidenceStoreError, MAX_RECONSTRUCT_BYTES
from genosapiens.wgs_pair_outcomes import MAX_CATALOG_ENTRIES

MAX_DECISION_BYTES = 64 * 1024
# The header retains the complete context, not one decision input/result.
# Use the codec's existing whole-value bound on both publication and replay.
MAX_CATALOG_HEADER_BYTES = MAX_RECONSTRUCT_BYTES
HOT_OUTCOMES = 4096
_PHASES = (HypothesisPhase.UNKNOWN, HypothesisPhase.IN_CIS,
           HypothesisPhase.IN_TRANS, HypothesisPhase.CONFLICTING)


@dataclass(frozen=True, slots=True, repr=False)
class MemberDecisionInputs:
    model: ModelDecisionMember
    qualification: QualificationDecisionInputs


@dataclass(frozen=True, slots=True, repr=False)
class PairDecision:
    requirement_decisions: tuple[ModelRequirementDecision, ...]
    model_state: str
    completeness: Completeness
    model_reasons: tuple[str, ...]
    qualification_states: tuple[str, str]
    qualification_reasons: tuple[tuple[str, ...], tuple[str, ...]]
    constructible: bool
    read_validation_required: bool


@dataclass(frozen=True, slots=True, repr=False)
class _MemberClass:
    owner: object
    ordinal: int


def evaluate_pair_decision(context, left, right, phase):
    """Assemble evaluated decisions exclusively through the existing owners."""
    requirements = tuple(evaluate_model_requirement(model=context.model,
        locus_id=context.locus_id, has_locus_binding=bool(context.locus_bindings),
        requirement=req, members=(left.model, right.model), phase=phase) for req in context.requirements)
    state, complete, reasons = summarize_model_requirements(requirements, context.requirements)
    qualified = tuple(qualify_whole_model(*evaluate_qualification_decision(replace(
        member.qualification, phase=phase)), state) for member in (left, right))
    constructible = all(value[0] != "nonqualifying" for value in qualified) \
        and construction_carriage_supported((left.model.genotype, right.model.genotype))
    return PairDecision(requirements, state, complete, reasons,
        tuple(value[0] for value in qualified), tuple(value[1] for value in qualified), constructible,
        any(member.qualification.read_validation_required for member in (left, right)))


class PairDecisionCatalog:
    """One context's private projection inside its sole normalized evidence owner."""
    def __init__(self, writer, context):
        self.writer = writer
        self.context = _typed(context, HypothesisContext)
        if self.context.model != "compound_heterozygous" or writer.connection is None:
            raise EvidenceStoreError("pair decision catalog requires its closed pair context")
        self.prefix = "pair_" + hashlib.sha256(_canonical(self.context).encode("ascii")).hexdigest()
        self.connection = writer.connection
        self.connection.execute(f"CREATE TABLE {self.prefix}_members (id INTEGER PRIMARY KEY, digest TEXT UNIQUE, content TEXT)")
        self.connection.execute(f"CREATE TABLE {self.prefix}_outcomes (id INTEGER PRIMARY KEY, left_id INTEGER, right_id INTEGER, "
            "phase INTEGER, content TEXT, UNIQUE(left_id,right_id,phase))")
        self.members = self.outcomes = self.evaluations = 0
        self.hot = OrderedDict()
        self.closed = False
        self.published = False

    def member(self, model, qualification):
        if self.closed or self.published:
            raise EvidenceStoreError("pair decision catalog is closed")
        value = _typed(MemberDecisionInputs(model, qualification), MemberDecisionInputs)
        if (qualification.model != self.context.model
                or qualification.context_requirements != self.context.requirements
                or qualification.genotype != model.genotype):
            raise EvidenceStoreError("pair decision member dependencies disagree")
        # Phase is supplied separately for each pair. No caller-provided phase
        # can split or alias the invariant member class.
        value = replace(value, qualification=replace(value.qualification, phase=None))
        content = _canonical(value)
        if len(content) > MAX_DECISION_BYTES:
            raise EvidenceStoreError("pair decision input exceeds capacity", capacity=True)
        digest = hashlib.sha256(content.encode("ascii")).hexdigest()
        found = self.connection.execute(f"SELECT id,content FROM {self.prefix}_members WHERE digest=?", (digest,)).fetchone()
        if found is None:
            if self.members == MAX_CATALOG_ENTRIES:
                raise EvidenceStoreError("pair member class capacity exceeded", capacity=True)
            ordinal = self.members
            self.connection.execute(f"INSERT INTO {self.prefix}_members VALUES (?,?,?)", (ordinal, digest, content))
            self.members += 1
        else:
            ordinal, previous = found
            if previous != content:
                raise EvidenceStoreError("pair member class digest collision")
        return _MemberClass(self, ordinal)

    def evaluate(self, left, right, phase):
        if (self.closed or self.published or any(type(value) is not _MemberClass or value.owner is not self for value in (left, right))
                or type(phase) is not HypothesisPhase or phase not in _PHASES):
            raise EvidenceStoreError("foreign pair decision dependencies")
        key = (left.ordinal, right.ordinal, _PHASES.index(phase))
        if key in self.hot:
            self.hot.move_to_end(key)
            return self.hot[key]
        row = self.connection.execute(f"SELECT id FROM {self.prefix}_outcomes WHERE left_id=? AND right_id=? AND phase=?", key).fetchone()
        if row is None:
            if self.outcomes == MAX_CATALOG_ENTRIES:
                raise EvidenceStoreError("pair outcome class capacity exceeded", capacity=True)
            values = []
            for member in (left, right):
                row = self.connection.execute(f"SELECT content FROM {self.prefix}_members WHERE id=?", (member.ordinal,)).fetchone()
                if row is None:
                    raise EvidenceStoreError("missing pair member class")
                values.append(_typed(json.loads(row[0]), MemberDecisionInputs, decode=True))
            result = evaluate_pair_decision(self.context, *values, phase)
            content = _canonical(result)
            if len(content) > MAX_DECISION_BYTES:
                raise EvidenceStoreError("pair decision result exceeds capacity", capacity=True)
            ordinal = self.outcomes
            self.connection.execute(f"INSERT INTO {self.prefix}_outcomes VALUES (?,?,?,?,?)", (ordinal, *key, content))
            self.outcomes += 1
            self.evaluations += 1
        else:
            ordinal, = row
        self.hot[key] = ordinal
        if len(self.hot) > HOT_OUTCOMES:
            self.hot.popitem(last=False)
        return ordinal

    def result(self, ordinal):
        if self.closed or type(ordinal) is not int or not 0 <= ordinal < self.outcomes:
            raise EvidenceStoreError("pair decision result is unavailable")
        row = self.connection.execute(f"SELECT content FROM {self.prefix}_outcomes WHERE id=?", (ordinal,)).fetchone()
        return _typed(json.loads(row[0]), PairDecision, decode=True)

    def member_inputs(self, ordinal):
        if self.closed or type(ordinal) is not int or not 0 <= ordinal < self.members:
            raise EvidenceStoreError("pair member input is unavailable")
        row = self.connection.execute(f"SELECT content FROM {self.prefix}_members WHERE id=?",(ordinal,)).fetchone()
        return _typed(json.loads(row[0]),MemberDecisionInputs,decode=True)

    def outcome_inputs(self, ordinal):
        if self.closed or type(ordinal) is not int or not 0 <= ordinal < self.outcomes:
            raise EvidenceStoreError("pair outcome input is unavailable")
        left,right,phase = self.connection.execute(f"SELECT left_id,right_id,phase FROM {self.prefix}_outcomes WHERE id=?",(ordinal,)).fetchone()
        return left,right,_PHASES[phase]

    def publish(self, name, *, retain_projection=False):
        if self.closed or self.published:
            raise EvidenceStoreError("pair decision catalog is closed")
        header = {
            "schema": "genosapiens.wgs-pair-decision-catalog", "version": 1,
            "context": self.context, "member_count": self.members, "outcome_count": self.outcomes,
            "members_name": name + ".members", "members_sha256": "0" * 64,
            "outcomes_name": name + ".outcomes", "outcomes_sha256": "0" * 64}
        # Real digests replace equal-length hexadecimal placeholders below.
        # Include the newline charged by EvidenceReader.reconstruct.
        if len(_canonical(header)) + 1 > MAX_CATALOG_HEADER_BYTES:
            self.writer.failed = True
            raise EvidenceStoreError("pair decision catalog header exceeds capacity", capacity=True,
                                     operation="pair_decision_header")
        def members():
            for ordinal, digest, content in self.connection.execute(f"SELECT id,digest,content FROM {self.prefix}_members ORDER BY id"):
                yield {"ordinal": ordinal, "sha256": digest, "inputs": json.loads(content)}
        def outcomes():
            for ordinal, left, right, phase, content in self.connection.execute(
                    f"SELECT id,left_id,right_id,phase,content FROM {self.prefix}_outcomes ORDER BY id"):
                yield {"ordinal": ordinal, "left": left, "right": right, "phase": _PHASES[phase].value,
                       "result": json.loads(content)}
        member_digest = self.writer.add_sequence(name + ".members", "pair-decision-members-v1", members())
        outcome_digest = self.writer.add_sequence(name + ".outcomes", "pair-decision-results-v1", outcomes())
        header.update(members_sha256=member_digest, outcomes_sha256=outcome_digest)
        digest = self.writer.add(name, "pair-decision-catalog-v1", header)
        self.published = True
        self.hot.clear()
        if not retain_projection:
            self.close()
        return digest

    def close(self):
        if not self.closed:
            self.connection.execute(f"DROP TABLE {self.prefix}_outcomes")
            self.connection.execute(f"DROP TABLE {self.prefix}_members")
            self.hot.clear()
            self.closed = True


class VerifiedPairDecisions:
    """Replay every distinct decision dependency class before exposing results."""
    def __init__(self, reader, reference, *, expected_context):
        if type(reader) is not EvidenceReader or reference.kind != "pair-decision-catalog-v1":
            raise EvidenceStoreError("foreign pair decision catalog")
        self.reader, self.reference = reader, reference
        self.connection = reader.connection
        self.prefix = "verified_pair_" + reference.sha256
        self.tables = []
        self.verified = False
        try:
            try:
                header = reader.reconstruct(reference, max_bytes=MAX_CATALOG_HEADER_BYTES)
            except EvidenceStoreError as error:
                if error.category == "response_capacity":
                    raise EvidenceStoreError("pair decision catalog header exceeds capacity", capacity=True,
                                             operation="pair_decision_header") from None
                raise
            if (set(header) != {"schema", "version", "context", "member_count", "outcome_count", "members_name",
                    "members_sha256", "outcomes_name", "outcomes_sha256"}
                    or header["schema"] != "genosapiens.wgs-pair-decision-catalog"
                    or type(header["version"]) is not int or header["version"] != 1
                    or any(type(header[key]) is not int or not 0 <= header[key] <= MAX_CATALOG_ENTRIES
                           for key in ("member_count", "outcome_count"))):
                raise EvidenceStoreError("invalid pair decision catalog")
            context = _typed(header["context"], HypothesisContext, decode=True)
            if context != expected_context or context.model != "compound_heterozygous":
                raise EvidenceStoreError("foreign pair decision context")
            self.connection.execute(f"CREATE TABLE {self.prefix}_members (id INTEGER PRIMARY KEY, digest TEXT UNIQUE, content TEXT)")
            self.tables.append(self.prefix + "_members")
            self.connection.execute(f"CREATE TABLE {self.prefix}_outcomes (id INTEGER PRIMARY KEY, left_id INTEGER, right_id INTEGER, "
                "phase INTEGER, content TEXT, UNIQUE(left_id,right_id,phase))")
            self.tables.append(self.prefix + "_outcomes")
            member_ref = reader.reference(header["members_name"], kind="pair-decision-members-v1")
            outcome_ref = reader.reference(header["outcomes_name"], kind="pair-decision-results-v1")
            if member_ref.sha256 != header["members_sha256"] or outcome_ref.sha256 != header["outcomes_sha256"]:
                raise EvidenceStoreError("pair decision catalog references changed")
            count = 0
            for row in reader.iter_array(member_ref, max_item_bytes=MAX_DECISION_BYTES + 1024):
                if set(row) != {"ordinal", "sha256", "inputs"} or type(row["ordinal"]) is not int or row["ordinal"] != count:
                    raise EvidenceStoreError("pair decision member coverage changed")
                inputs = _typed(row["inputs"], MemberDecisionInputs, decode=True)
                content = _canonical(inputs)
                if (hashlib.sha256(content.encode("ascii")).hexdigest() != row["sha256"]
                        or inputs.qualification.phase is not None
                        or inputs.qualification.model != context.model
                        or inputs.qualification.context_requirements != context.requirements
                        or inputs.qualification.genotype != inputs.model.genotype):
                    raise EvidenceStoreError("pair decision member dependencies changed")
                self.connection.execute(f"INSERT INTO {self.prefix}_members VALUES (?,?,?)", (count, row["sha256"], content))
                count += 1
            if count != header["member_count"]:
                raise EvidenceStoreError("pair decision member count changed")
            count = 0
            for row in reader.iter_array(outcome_ref, max_item_bytes=MAX_DECISION_BYTES + 1024):
                if (set(row) != {"ordinal", "left", "right", "phase", "result"}
                        or type(row["ordinal"]) is not int or row["ordinal"] != count
                        or any(type(row[key]) is not int or not 0 <= row[key] < header["member_count"] for key in ("left", "right"))):
                    raise EvidenceStoreError("pair decision result coverage changed")
                phase = HypothesisPhase(row["phase"])
                if phase not in _PHASES:
                    raise EvidenceStoreError("invalid pair decision phase")
                members = tuple(_typed(json.loads(self.connection.execute(
                    f"SELECT content FROM {self.prefix}_members WHERE id=?", (row[key],)).fetchone()[0]),
                    MemberDecisionInputs, decode=True) for key in ("left", "right"))
                result = _typed(row["result"], PairDecision, decode=True)
                if result != evaluate_pair_decision(context, *members, phase):
                    raise EvidenceStoreError("pair decision result does not replay")
                self.connection.execute(f"INSERT INTO {self.prefix}_outcomes VALUES (?,?,?,?,?)",
                    (count, row["left"], row["right"], _PHASES.index(phase), _canonical(result)))
                count += 1
            if count != header["outcome_count"]:
                raise EvidenceStoreError("pair decision result count changed")
            self.member_count, self.outcome_count = header["member_count"], count
            self.verified = True
        except BaseException:
            self.close()
            raise

    def result(self, ordinal):
        if not self.verified or type(ordinal) is not int or not 0 <= ordinal < self.outcome_count:
            raise EvidenceStoreError("pair decision result is unavailable")
        row = self.connection.execute(f"SELECT content FROM {self.prefix}_outcomes WHERE id=?", (ordinal,)).fetchone()
        return _typed(json.loads(row[0]), PairDecision, decode=True)

    def member_inputs(self, ordinal):
        if not self.verified or type(ordinal) is not int or not 0 <= ordinal < self.member_count:
            raise EvidenceStoreError("pair member input is unavailable")
        row = self.connection.execute(f"SELECT content FROM {self.prefix}_members WHERE id=?",(ordinal,)).fetchone()
        return _typed(json.loads(row[0]),MemberDecisionInputs,decode=True)

    def outcome_inputs(self, ordinal):
        if not self.verified or type(ordinal) is not int or not 0 <= ordinal < self.outcome_count:
            raise EvidenceStoreError("pair outcome input is unavailable")
        left,right,phase = self.connection.execute(f"SELECT left_id,right_id,phase FROM {self.prefix}_outcomes WHERE id=?",(ordinal,)).fetchone()
        return left,right,_PHASES[phase]

    def close(self):
        self.verified = False
        for table in self.tables:
            self.connection.execute("DROP TABLE " + table)
        self.tables.clear()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()
