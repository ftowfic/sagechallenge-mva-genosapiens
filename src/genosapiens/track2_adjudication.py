"""Append-only adjudication and non-overridable blocking-fact enforcement."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
import hashlib
import json
import re
from typing import Any, Mapping


DOMAINS = {"endpoint", "mechanism", "subtype", "target"}
ROLES = {"curator", "release_approver", "safety_reviewer", "scientific_reviewer"}
INDEPENDENCE = {"independent", "not_reviewed", "self_review_disclosed"}
ACTIONS = {"initial", "override"}
DECISIONS = {"accepted", "rejected", "revised"}


class AdjudicationError(ValueError):
    """Raised when a review decision or override is absent or invalid."""


@dataclass(frozen=True, slots=True)
class AdjudicationRecord:
    adjudication_id: str
    subject_id: str
    domain: str
    action: str
    reviewer_role: str
    reviewer_count: int
    independence: str
    prior_state: str
    decision: str
    reason: str
    timestamp: str
    evidence_ids: tuple[str, ...]
    override_of: str | None
    blocking_fact: bool
    reviewed_content_sha256: str | None = None

    def __post_init__(self) -> None:
        try:
            reviewed_at = datetime.fromisoformat(self.timestamp.replace("Z", "+00:00"))
        except (TypeError, ValueError) as error:
            raise AdjudicationError("adjudication timestamp must use ISO format") from error
        if (
            not self.adjudication_id
            or not self.subject_id
            or self.domain not in DOMAINS
            or self.action not in ACTIONS
            or self.reviewer_role not in ROLES
            or not isinstance(self.reviewer_count, int)
            or isinstance(self.reviewer_count, bool)
            or self.reviewer_count < 1
            or self.independence not in INDEPENDENCE
            or not self.prior_state
            or self.decision not in DECISIONS
            or not self.reason
            or (self.reviewed_content_sha256 is not None
                and (len(self.reviewed_content_sha256) != 64
                     or any(char not in "0123456789abcdef" for char in self.reviewed_content_sha256)))
            or reviewed_at.tzinfo is None
            or not self.evidence_ids
            or len(self.evidence_ids) != len(set(self.evidence_ids))
        ):
            raise AdjudicationError("adjudication record is incomplete")
        if self.reviewer_count == 1 and self.independence != "self_review_disclosed":
            raise AdjudicationError(
                "solo review must be disclosed and cannot be represented as independent"
            )
        if self.reviewer_count > 1 and self.independence == "not_reviewed":
            raise AdjudicationError("completed adjudication cannot be not_reviewed")
        if (self.action == "override") != (self.override_of is not None):
            raise AdjudicationError("override action requires explicit prior adjudication")

    def as_track2_payload(self) -> dict:
        return {
            "decision": self.decision,
            "independence": self.independence,
            "prior_state": self.prior_state,
            "reason": self.reason,
            "reviewer_role": self.reviewer_role,
            "reviewed_content_sha256": self.reviewed_content_sha256,
            "timestamp": self.timestamp,
        }


@dataclass(frozen=True, slots=True)
class AdjudicationLedger:
    records: tuple[AdjudicationRecord, ...] = ()

    def latest(self, subject_id: str, domain: str) -> AdjudicationRecord | None:
        matches = tuple(
            record
            for record in self.records
            if record.subject_id == subject_id and record.domain == domain
        )
        return matches[-1] if matches else None


def make_adjudication_id(
    subject_id: str, domain: str, action: str, timestamp: str, prior_id: str | None
) -> str:
    content = "\0".join((subject_id, domain, action, timestamp, prior_id or ""))
    return "adjudication:" + hashlib.sha256(content.encode()).hexdigest()[:24]


def append_adjudication(
    ledger: AdjudicationLedger, record: AdjudicationRecord
) -> AdjudicationLedger:
    if record.adjudication_id in {item.adjudication_id for item in ledger.records}:
        raise AdjudicationError("adjudication ID must be unique")
    previous = ledger.latest(record.subject_id, record.domain)
    if previous is None:
        if (
            record.action != "initial"
            or record.override_of is not None
            or record.prior_state != "unreviewed"
        ):
            raise AdjudicationError("initial adjudication cannot be an override")
    else:
        if record.action != "override" or record.override_of != previous.adjudication_id:
            raise AdjudicationError("changed adjudication requires a recorded override")
        if previous.blocking_fact:
            raise AdjudicationError("non-overridable blocking fact cannot be overridden")
        if record.prior_state != previous.decision:
            raise AdjudicationError("override prior state does not match the ledger")
    return AdjudicationLedger((*ledger.records, record))


def assert_adjudicated_state(
    ledger: AdjudicationLedger,
    subject_id: str,
    domain: str,
    declared_state: str,
) -> None:
    latest = ledger.latest(subject_id, domain)
    if latest is None or latest.decision != declared_state:
        raise AdjudicationError("declared state has an unrecorded override")


@dataclass(frozen=True, slots=True)
class ReviewedInputBinding:
    """Bind a reviewed scientific input to the latest subject/domain decision."""

    subject_id: str
    domain: str
    adjudication_id: str
    declared_state: str
    content_sha256: str


def reviewed_input_sha256(value: Mapping[str, Any]) -> str:
    from genosapiens.track2_scientific_storage import json_identity
    return json_identity(value)['sha256']


def assert_reviewed_input(
    ledger: AdjudicationLedger,
    binding: ReviewedInputBinding,
    value: Mapping[str, Any],
) -> None:
    assert_adjudicated_state(ledger, binding.subject_id, binding.domain, binding.declared_state)
    latest = ledger.latest(binding.subject_id, binding.domain)
    if (latest is None or latest.adjudication_id != binding.adjudication_id
            or binding.content_sha256 != reviewed_input_sha256(value)
            or latest.reviewed_content_sha256 != binding.content_sha256
            or latest.blocking_fact or latest.decision not in {"accepted", "revised"}):
        raise AdjudicationError("scientific input lacks a matching nonblocking adjudication")


# Machine decisions have their own authority. They never enter the human ledger.
AUTOMATIC_POLICY = "diagnosis-automatic-acceptance/1.0.0"
AUTOMATIC_CHECKS = {
    "input": ("integrity", "source_binding", "scope", "selection"),
    "mechanism": ("genotype", "quality", "inheritance", "transcripts", "direction", "functional_sources"),
    "plan": ("integrity", "mechanisms", "applicability", "prospective_scope"),
    "discovery_plan": ("integrity", "mechanisms", "hypotheses", "prospective_scope"),
}


def automatic_decision(*, domain, subject, evidence, checks, execution, timestamp):
    """Record a reproducible policy result without asserting human acceptance.

    Scientific producers supply the checked facts. Consumers must replay those
    producers as well as this decision. A hash alone does not prove a fact.
    """
    if domain not in AUTOMATIC_CHECKS or set(checks) != set(AUTOMATIC_CHECKS[domain]):
        raise AdjudicationError("Automatic acceptance checks are incomplete.")
    if not isinstance(subject, Mapping) or not subject or not isinstance(evidence, Mapping) or not evidence:
        raise AdjudicationError("Automatic acceptance lacks its subject or evidence.")
    if (set(execution) != {"code_sha256", "image"}
            or not isinstance(execution['image'],str)
            or re.fullmatch(r'.+@sha256:[0-9a-f]{64}',execution['image']) is None):
        raise AdjudicationError("Automatic acceptance lacks its execution identity.")
    digest = execution["code_sha256"]
    if not isinstance(digest, str) or len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest):
        raise AdjudicationError("Automatic acceptance has an invalid code identity.")
    try:
        parsed = datetime.fromisoformat(timestamp.replace("Z", "+00:00"))
    except (AttributeError, TypeError, ValueError) as error:
        raise AdjudicationError("Automatic acceptance requires an actual timestamp.") from error
    if parsed.tzinfo is None:
        raise AdjudicationError("Automatic acceptance requires a timezone.")
    for row in checks.values():
        if (not isinstance(row, Mapping) or set(row) != {"state", "reasons"}
                or row["state"] not in {"accepted", "conflicting", "unresolved"}
                or not isinstance(row["reasons"], list) or not row["reasons"]
                or any(not isinstance(r, str) or not r for r in row["reasons"])):
            raise AdjudicationError("Automatic acceptance has an invalid check.")
    states = {row["state"] for row in checks.values()}
    outcome = "conflicting" if "conflicting" in states else "unresolved" if "unresolved" in states else "accepted"
    result={"schema": "genosapiens.track2-automatic-decision/v2" if domain=='discovery_plan' else "genosapiens.track2-automatic-decision/v1", "authority": "machine",
        "policy": 'diagnosis-exploratory-planning/1.0.0' if domain=='discovery_plan' else AUTOMATIC_POLICY,
        "domain": domain, "subject_sha256": reviewed_input_sha256(subject),
        "evidence_sha256": reviewed_input_sha256(evidence), "checks": dict(checks),
        "execution": dict(execution), "timestamp": timestamp, "outcome": outcome,
        "reasons": sorted({r for row in checks.values() for r in row["reasons"]}),
        "evidence_strength": "policy_prediction" if domain == "mechanism" else "technical_acceptance"}
    from genosapiens.json_schema import validate_json_schema
    validate_json_schema(result,json.loads((Path(__file__).resolve().parents[2]/
        'schemas/track2/automatic-decision.schema.json').read_bytes()))
    return result


def validate_automatic_decision(record, *, domain, subject, evidence, checks, execution):
    """Replay the decision against independently derived current facts."""
    if not isinstance(record, Mapping) or "timestamp" not in record:
        raise AdjudicationError("Automatic acceptance record is missing.")
    expected = automatic_decision(domain=domain, subject=subject, evidence=evidence,
        checks=checks, execution=execution, timestamp=record["timestamp"])
    if record != expected:
        raise AdjudicationError("Automatic acceptance differs from its subject, evidence, or policy.")
    return expected
