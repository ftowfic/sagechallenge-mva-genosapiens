"""Computed member population evidence for a complete hypothesis universe.

Each required member/context assessment is evaluated once. Pair assessment is
the original maximum/any aggregation of those computed member results; phase
and partner identity do not change the member population kernel's inputs.
"""
from dataclasses import dataclass
import hashlib
import json
import sqlite3

from genosapiens.evidence_model import GenotypeHypothesis, _Controlled, _typed
from genosapiens.evidence_serialization import _canonical
from genosapiens.population_policy import (
    ContextKey, MemberAssessment, PopulationOutcome, summarize_population_members,
    evaluate_population_member_inputs,
)
from genosapiens.wgs_evidence_replay import EvidenceReplay
from genosapiens.wgs_evidence_store import EvidenceStoreError

MAX_ASSESSMENT_BYTES = 64 * 1024**2
MAX_SETTINGS_CACHE_BYTES = 1024**2


@dataclass(frozen=True, slots=True, repr=False)
class PopulationAssessmentValues(_Controlled):
    selected_context: ContextKey | None
    members: tuple[MemberAssessment, ...]
    population_penalty: str | None
    population_contribution: str | None
    outcome: PopulationOutcome
    reason: str


class CompleteContextPopulation:
    def __init__(self, writer, universe, policy, *, name):
        self.connection, self.writer, self.universe, self.policy = writer.connection, writer, universe, policy
        self.name, self.prefix = name, "context_population_" + hashlib.sha256(name.encode()).hexdigest()
        self.tables, self.closed, self.assessment_count = [], False, 0
        self.settings_cache = None
        self.settings_cache_hits = self.settings_decodes = 0
        if (universe.connection is not writer.connection or policy.connection is not writer.connection
                or universe.closed or policy.closed):
            raise EvidenceStoreError("foreign context population owners")
        try:
            for suffix, definition in (
                ("members", "context TEXT,id TEXT,selected INTEGER,content TEXT,PRIMARY KEY(context,id)"),
                ("contexts", "context TEXT PRIMARY KEY,settings TEXT,exceptions TEXT"),
            ):
                table = self.prefix + "_" + suffix
                self.connection.execute(f"CREATE TABLE {table} ({definition})")
                self.tables.append(table)
            # Singleton contexts can differ by model. Retain exactly the same
            # source member and global exception ordinals for each such key.
            for value in universe.singletons():
                key = ContextKey.from_context(value.hypothesis.context)
                for ref in value.hypothesis.evidence_subject_refs:
                    self._evaluate(key, ref.allele_id, expected_native=ref.identity_id)
            if universe.pairs is not None:
                key = ContextKey.from_context(universe.pairs.context)
                for _, allele_id, native, _, parent in universe.pairs.domain.members():
                    if not parent:
                        self._evaluate(key, allele_id, expected_native=native)
            self.members_sha256 = writer.add_sequence(name + ".members", "context-population-members-v1", self.rows())
            self.sha256 = writer.add(name, "context-population-evidence-v1", {
                "schema": "genosapiens.wgs-context-population-evidence", "version": 1,
                "hypotheses": {"name": universe.name, "sha256": universe.sha256},
                "context_evidence": {"name": universe.members.name, "sha256": universe.members.sha256},
                "policy": policy.reference, "assessment_count": self.assessment_count,
                "hypothesis_count": universe.hypothesis_count,
                "members": {"name": name + ".members", "sha256": self.members_sha256},
            })
        except BaseException as error:
            writer.failed = True
            self.close()
            if isinstance(error, (ValueError, TypeError, KeyError, sqlite3.Error)):
                raise EvidenceStoreError("invalid complete context population evidence") from None
            raise

    def _evaluate(self, key, allele_id, *, expected_native):
        encoded_key = _canonical(key)
        if self.connection.execute(f"SELECT 1 FROM {self.prefix}_members WHERE context=? AND id=?",
                                   (encoded_key, allele_id)).fetchone():
            return
        allele = self.universe.members.member(allele_id).allele
        if allele.native_identity_id != expected_native:
            raise EvidenceStoreError("population hypothesis member identity changed")
        self.policy._require_open()
        binding = (self.policy.sha256, self.policy.policy_sha256, encoded_key)
        cached = self.settings_cache
        if cached is not None and cached[0] is self.policy and cached[1] == binding:
            settings, exceptions = cached[2:4]
            self.settings_cache_hits += 1
        else:
            row = self.connection.execute(f"SELECT settings,exceptions FROM {self.prefix}_contexts WHERE context=?",
                                          (encoded_key,)).fetchone()
            if row is None:
                settings, exceptions = self.policy.inputs(key)
                row = (_canonical(settings), _canonical(exceptions))
                self.connection.execute(f"INSERT INTO {self.prefix}_contexts VALUES (?,?,?)",
                    (encoded_key, *row))
            else:
                from genosapiens.population_policy import ContextPolicy, CuratedException
                settings = _typed(json.loads(row[0]), ContextPolicy | None, decode=True)
                exceptions = _typed(json.loads(row[1]), tuple[tuple[int, CuratedException], ...], decode=True)
                self.settings_decodes += 1
            size = len(encoded_key)+sum(len(value) for value in row)+len(binding[0])+len(binding[1])
            self.settings_cache = ((self.policy,binding,settings,exceptions,size)
                if size <= MAX_SETTINGS_CACHE_BYTES else None)
        assessment = self.universe.members.population_work_assessment(key, allele, self.policy)
        if assessment is None:
            assessment, = evaluate_population_member_inputs(key, (allele,), settings=settings, indexed_exceptions=exceptions)
        content = _canonical(assessment)
        if len(content) > MAX_ASSESSMENT_BYTES:
            raise EvidenceStoreError("population member evidence exceeds capacity", capacity=True)
        self.connection.execute(f"INSERT INTO {self.prefix}_members VALUES (?,?,?,?)",
            (encoded_key, allele_id, int(settings is not None), content))
        self.assessment_count += 1

    def _require_open(self):
        if self.closed:
            raise EvidenceStoreError("context population owner is closed")

    def rows(self):
        self._require_open()
        for context, _, selected, content in self.connection.execute(
                f"SELECT * FROM {self.prefix}_members ORDER BY context,id"):
            yield {"context": json.loads(context), "selected": bool(selected), "assessment": json.loads(content)}

    def member(self, key, allele_id):
        self._require_open()
        row = self.connection.execute(f"SELECT selected,content FROM {self.prefix}_members WHERE context=? AND id=?",
            (_canonical(key), allele_id)).fetchone()
        if row is None:
            raise EvidenceStoreError("missing computed population member")
        return bool(row[0]), _typed(json.loads(row[1]), MemberAssessment, decode=True)

    def assessment(self, hypothesis):
        self._require_open()
        hypothesis = _typed(hypothesis, GenotypeHypothesis)
        key = ContextKey.from_context(hypothesis.context)
        members, selected_values = [], set()
        for ref in sorted(hypothesis.evidence_subject_refs, key=lambda ref: ref.allele_id):
            selected, member = self.member(key, ref.allele_id)
            if member.native_identity_id != ref.identity_id:
                raise EvidenceStoreError("population assessment member identity changed")
            members.append(member)
            selected_values.add(selected)
        if len(selected_values) != 1 or len({m.allele_id for m in members}) != len(members):
            raise EvidenceStoreError("population assessment member scope changed")
        selected, = selected_values
        penalty, contribution, outcome = summarize_population_members(tuple(members))
        return PopulationAssessmentValues(key if selected else None, tuple(members), penalty, contribution, outcome,
            "population_only;not_all_checks_or_scientific_acceptance" if selected else "no_exact_context_policy")

    def close(self):
        self.settings_cache = None
        for table in reversed(self.tables):
            self.connection.execute("DROP TABLE " + table)
        self.tables.clear()
        self.closed = True

    def __enter__(self):
        self._require_open()
        return self

    def __exit__(self, *_):
        self.close()


class VerifiedContextPopulation(CompleteContextPopulation):
    def __init__(self, reader, reference, universe, policy):
        if reference.kind != "context-population-evidence-v1":
            raise EvidenceStoreError("foreign context population evidence")
        reader._require_reference(reference)
        super().__init__(EvidenceReplay(reader), universe, policy, name=reference.name)
