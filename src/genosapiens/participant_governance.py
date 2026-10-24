"""Private, offline participant attestations, not verification of real-world status.

The release authority attests the complete roster, private registered identities,
applicable terms/eligibility and review declarations. Never publish this manifest.
Use the release policy's approval_public_key_sha256 (hash of exact public PEM
bytes), not a key/hash nominated by the manifest. The caller supplies the current
public rules contract SHA, trusted clock, and current replacement manifest before
EVERY authorization. Offline signatures cannot discover revocation or omitted
people. Replacement distribution and host-account identity are operator duties.
"""

from __future__ import annotations

import base64
from dataclasses import asdict, dataclass, field, fields
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import subprocess
import tempfile
from typing import Any

from genosapiens.bounded_json import load_bounded_json, read_bounded_bytes


PARTICIPANT_GOVERNANCE_SCHEMA = "genosapiens.participant-governance/v1"
PARTICIPANT_GOVERNANCE_ATTESTATION_SCHEMA = (
    "genosapiens.participant-governance-attestation/v1"
)
REQUIRED_REVIEW_DOMAINS = frozenset({
    "subtype", "mechanism", "endpoint", "pharmacology", "exposure", "safety",
    "claims", "citations", "privacy", "package",
})
PARTICIPANT_ROLES = frozenset({"team_member", "submitter", "controlled_reviewer"})
_PARTICIPANT_ID = re.compile(r"participant:[0-9a-f]{32}")
_SHA256 = re.compile(r"[0-9a-f]{64}")
_REQUIRED_ACCEPTANCES = (
    "individual_registration", "age_eligible", "huggingface_terms_accepted",
    "hackathon_terms_accepted", "no_recontact_accepted", "no_resharing_accepted",
    "applicable_law_accepted", "professional_obligations_accepted",
)
_SUBMITTER_ACCEPTANCES = (
    "submission_material_cc_by_accepted", "organizer_rerun_accepted",
    "public_name_attribution_accepted",
)


class ParticipantGovernanceError(ValueError):
    """A non-sensitive, fail-closed governance rejection."""


@dataclass(frozen=True, slots=True, repr=False)
class Participant:
    participant_id: str
    registered_identity: str
    host_uid: int
    roles: tuple[str, ...]
    individual_registration: bool
    age_eligible: bool
    huggingface_terms_accepted: bool
    hackathon_terms_accepted: bool
    gated_data_authorized: bool
    no_recontact_accepted: bool
    no_resharing_accepted: bool
    applicable_law_accepted: bool
    professional_obligations_accepted: bool
    submission_material_cc_by_accepted: bool
    organizer_rerun_accepted: bool
    public_name_attribution_accepted: bool


@dataclass(frozen=True, slots=True, repr=False)
class StructuredPackageReview:
    package_sha256: str
    reviewer_participant_id: str
    author_participant_ids: tuple[str, ...]
    scope_domains: tuple[str, ...]
    independence_status: str
    reviewed_at: str
    conflicts: tuple[str, ...]  # An explicit empty tuple attests no other conflicts.


@dataclass(frozen=True, slots=True, repr=False)
class ParticipantGovernanceManifest:
    rules_sha256: str
    issued_at: str
    expires_at: str
    roster_complete: bool
    participants: tuple[Participant, ...]
    team_member_ids: tuple[str, ...]
    controlled_review_roster: tuple[str, ...]
    package_review: StructuredPackageReview | None
    schema: str = PARTICIPANT_GOVERNANCE_SCHEMA
    _authenticated_sha256: str | None = field(default=None, init=False, compare=False)

    def as_dict(self) -> dict[str, Any]:
        """Return PRIVATE payload fields, excluding local authentication state."""
        return {
            **{f.name: getattr(self, f.name) for f in fields(self)
               if f.init and f.name not in {"participants", "package_review"}},
            "participants": [asdict(p) for p in self.participants],
            "package_review": asdict(self.package_review) if self.package_review else None,
        }


def _canonical_bytes(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":"),
                       ensure_ascii=True, allow_nan=False) + "\n").encode("ascii")


def participant_governance_as_json(manifest: ParticipantGovernanceManifest) -> str:
    """Serialize private payload only; compact sorted ASCII JSON plus one LF.

    Lists representing sets must already be sorted and unique. Payload timestamps
    use UTC datetime.isoformat(), including +00:00, not Z. Serialization does not
    validate or authorize an unsigned manifest.
    """
    return _canonical_bytes(manifest.as_dict()).decode("ascii")


def _require(condition: bool) -> None:
    if not condition:
        raise ParticipantGovernanceError("participant governance manifest is invalid")


def _matches(pattern: re.Pattern[str], value: Any) -> bool:
    return type(value) is str and pattern.fullmatch(value) is not None


def _strings(value: Any) -> tuple[str, ...]:
    _require(type(value) is list and all(type(item) is str for item in value))
    _require(value == sorted(set(value)))
    return tuple(value)


def _record(value: Any, record_type: type) -> None:
    _require(type(value) is dict and set(value) == {
        f.name for f in fields(record_type) if f.init
    })


def _timestamp(value: Any) -> datetime:
    _require(type(value) is str)
    parsed = datetime.fromisoformat(value)
    _require(parsed.tzinfo is not None and parsed.utcoffset() is not None)
    _require(value == parsed.astimezone(timezone.utc).isoformat())
    return parsed


def _uid(value: Any) -> bool:
    return type(value) is int and 0 <= value < 2**32 - 1


def _parse_manifest(data: Any, checked_at: datetime,
                    expected_rules_sha256: str) -> ParticipantGovernanceManifest:
    _record(data, ParticipantGovernanceManifest)
    _require(data["schema"] == PARTICIPANT_GOVERNANCE_SCHEMA)
    _require(_matches(_SHA256, expected_rules_sha256)
             and data["rules_sha256"] == expected_rules_sha256)
    _require(isinstance(checked_at, datetime) and checked_at.tzinfo is not None
             and checked_at.utcoffset() is not None)
    issued, expires = _timestamp(data["issued_at"]), _timestamp(data["expires_at"])
    if not issued <= checked_at < expires:
        raise ParticipantGovernanceError("participant governance is not current")
    _require(data["roster_complete"] is True)
    team = _strings(data["team_member_ids"])
    roster = _strings(data["controlled_review_roster"])
    _require(bool(team) and type(data["participants"]) is list)
    participants = []
    for row in data["participants"]:
        _record(row, Participant)
        roles = _strings(row["roles"])
        _require(bool(roles) and set(roles) <= PARTICIPANT_ROLES)
        _require(_matches(_PARTICIPANT_ID, row["participant_id"])
                 and row["participant_id"] != "participant:" + "0" * 32)
        identity = row["registered_identity"]
        _require(type(identity) is str and 0 < len(identity) <= 1024
                 and identity == identity.strip() and identity.isprintable())
        _require(_uid(row["host_uid"]))
        for name in (*_REQUIRED_ACCEPTANCES, *_SUBMITTER_ACCEPTANCES,
                     "gated_data_authorized"):
            _require(type(row[name]) is bool)
        _require(all(row[name] is True for name in _REQUIRED_ACCEPTANCES))
        _require("team_member" in roles or "controlled_reviewer" in roles)
        if "submitter" in roles:
            _require("team_member" in roles
                     and all(row[name] is True for name in _SUBMITTER_ACCEPTANCES))
        if "controlled_reviewer" in roles:
            _require(row["gated_data_authorized"] is True)
        participants.append(Participant(**{**row, "roles": roles}))
    ids = [p.participant_id for p in participants]
    _require(ids == sorted(set(ids)))
    _require(len({p.host_uid for p in participants}) == len(participants))
    _require(len({p.registered_identity.casefold() for p in participants})
             == len(participants))
    _require(team == tuple(p.participant_id for p in participants
                           if "team_member" in p.roles))
    _require(roster == tuple(p.participant_id for p in participants
                             if "controlled_reviewer" in p.roles))
    review = None
    if data["package_review"] is not None:
        row = data["package_review"]
        _record(row, StructuredPackageReview)
        authors = _strings(row["author_participant_ids"])
        domains = _strings(row["scope_domains"])
        conflicts = _strings(row["conflicts"])
        _require(all(0 < len(c) <= 1024 and c == c.strip() and c.isprintable()
                     for c in conflicts))
        submitters = {p.participant_id for p in participants if "submitter" in p.roles}
        _require(bool(authors) and set(authors) <= submitters)
        _require(_matches(_SHA256, row["package_sha256"])
                 and row["reviewer_participant_id"] in roster)
        _require(bool(domains) and set(domains) <= REQUIRED_REVIEW_DOMAINS)
        status = row["independence_status"]
        _require(status in ("independent", "self_review_disclosed"))
        is_author = row["reviewer_participant_id"] in authors
        _require((status == "self_review_disclosed") == is_author)
        _require(status != "independent" or not conflicts)
        _require(_timestamp(row["reviewed_at"]) <= issued)
        review = StructuredPackageReview(**{
            **row, "author_participant_ids": authors, "scope_domains": domains,
            "conflicts": conflicts,
        })
    return ParticipantGovernanceManifest(**{
        **data, "participants": tuple(participants), "team_member_ids": team,
        "controlled_review_roster": roster, "package_review": review,
    })


def _run_openssl(arguments: list[str]) -> bytes:
    try:
        result = subprocess.run(["openssl", *arguments], capture_output=True,
                                check=False, timeout=30)
    except (OSError, subprocess.TimeoutExpired):
        raise ParticipantGovernanceError("participant governance attestation failed") from None
    if result.returncode != 0:
        raise ParticipantGovernanceError("participant governance attestation failed")
    return result.stdout


def load_signed_participant_governance(
    path: str | Path,
    verification_key_path: str | Path,
    expected_verification_key_sha256: str,
    checked_at: datetime,
    expected_rules_sha256: str,
) -> tuple[ParticipantGovernanceManifest, str]:
    """Authenticate a complete private roster at the caller's trusted check time.

    The envelope has exactly schema, algorithm="Ed25519", manifest and signature
    (canonical base64) fields. The signature covers canonical manifest bytes.
    Returns (immutable manifest, SHA256(canonical payload)), NOT envelope SHA.
    Review may be null for mounts; partial review never authorizes a package.
    """
    try:
        for source in (path, verification_key_path):
            artifact = Path(source)
            _require(not artifact.is_symlink() and artifact.is_file())
        key = read_bounded_bytes(verification_key_path, maximum_bytes=64 * 1024)
        if (not _matches(_SHA256, expected_verification_key_sha256)
                or hashlib.sha256(key).hexdigest() != expected_verification_key_sha256):
            raise ParticipantGovernanceError("participant governance authority is invalid")
        envelope, encoded = load_bounded_json(path)
        _require(type(envelope) is dict and set(envelope) == {
            "schema", "algorithm", "manifest", "signature",
        })
        _require(envelope["schema"] == PARTICIPANT_GOVERNANCE_ATTESTATION_SCHEMA
                 and envelope["algorithm"] == "Ed25519")
        _require(_canonical_bytes(envelope) == encoded)
        _require(type(envelope["signature"]) is str)
        signature = base64.b64decode(envelope["signature"], validate=True)
        _require(len(signature) == 64
                 and base64.b64encode(signature).decode("ascii") == envelope["signature"])
        content = _canonical_bytes(envelope["manifest"])
        # Verify the exact pinned bytes, not a key path that could change after hashing.
        with (tempfile.NamedTemporaryFile() as key_file,
              tempfile.NamedTemporaryFile() as message_file,
              tempfile.NamedTemporaryFile() as signature_file):
            for handle, value in ((key_file, key), (message_file, content),
                                  (signature_file, signature)):
                handle.write(value)
                handle.flush()
            der = _run_openssl(["pkey", "-pubin", "-in", key_file.name, "-outform", "DER"])
            _require(len(der) == 44 and der.startswith(bytes.fromhex("302a300506032b6570032100")))
            _run_openssl(["pkeyutl", "-verify", "-pubin", "-inkey", key_file.name,
                          "-rawin", "-in", message_file.name, "-sigfile", signature_file.name])
        manifest = _parse_manifest(envelope["manifest"], checked_at, expected_rules_sha256)
        digest = hashlib.sha256(content).hexdigest()
        object.__setattr__(manifest, "_authenticated_sha256", digest)
        return manifest, digest
    except ParticipantGovernanceError:
        raise
    except (ValueError, TypeError, KeyError, AttributeError, OSError,
            RecursionError, OverflowError):
        raise ParticipantGovernanceError("participant governance manifest is invalid") from None


def _require_authenticated(manifest: ParticipantGovernanceManifest) -> None:
    try:
        if (type(manifest) is ParticipantGovernanceManifest
                and manifest._authenticated_sha256 is not None
                and hashlib.sha256(_canonical_bytes(manifest.as_dict())).hexdigest()
                == manifest._authenticated_sha256):
            return
    except (ValueError, TypeError, AttributeError, RecursionError):
        pass
    raise ParticipantGovernanceError("participant governance authentication required")


def authorize_controlled_access(
    manifest: ParticipantGovernanceManifest, participant_id: str, host_uid: int,
) -> None:
    """Check a freshly loaded manifest against the launcher's ACTUAL host UID.

    Never obtain host_uid from the manifest or an untrusted container environment.
    UID bindings describe one host account namespace; this does not grant a mount.
    """
    _require_authenticated(manifest)
    if _matches(_PARTICIPANT_ID, participant_id) and _uid(host_uid):
        for participant in manifest.participants:
            if (participant.participant_id == participant_id
                    and participant_id in manifest.controlled_review_roster
                    and participant.host_uid == host_uid):
                return
    raise ParticipantGovernanceError("controlled access denied")


def authorize_package_review(
    manifest: ParticipantGovernanceManifest, reviewer_authorization_id: str,
    expected_package_sha256: str, reviewed_at: str | datetime,
) -> str:
    """Bind ReleaseApproval.authorization_id, package and review time to the roster.

    All submitters' reuse acceptances were required during loading, not merely the
    reviewer's. Return the attested independence status, never inferred expertise.
    The approval timestamp may use any aware ISO offset but must be the same instant.
    """
    _require_authenticated(manifest)
    review = manifest.package_review
    try:
        timestamp = datetime.fromisoformat(reviewed_at) if type(reviewed_at) is str else reviewed_at
        if (review is not None and _matches(_PARTICIPANT_ID, reviewer_authorization_id)
                and reviewer_authorization_id == review.reviewer_participant_id
                and reviewer_authorization_id in manifest.controlled_review_roster
                and _matches(_SHA256, expected_package_sha256)
                and review.package_sha256 == expected_package_sha256
                and set(review.scope_domains) == REQUIRED_REVIEW_DOMAINS
                and isinstance(timestamp, datetime) and timestamp.tzinfo is not None
                and timestamp.utcoffset() is not None
                and timestamp == _timestamp(review.reviewed_at)):
            return review.independence_status
    except (ValueError, TypeError, AttributeError, OverflowError):
        pass
    raise ParticipantGovernanceError("package review authorization denied")
