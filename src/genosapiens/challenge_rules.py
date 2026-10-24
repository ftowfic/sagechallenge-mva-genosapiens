"""Pinned, offline-verifiable challenge rule contracts."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import date, datetime, timezone
from pathlib import Path
import re
from typing import Any, Mapping

from genosapiens.bounded_json import BoundedJSONError, load_bounded_json


CHALLENGE_RULES_SCHEMA = "genosapiens.challenge-rules/v1"
SHA256 = re.compile(r"^[0-9a-f]{64}$")
GIT_ID = re.compile(r"^[0-9a-f]{40}$")
EVENT_ID = re.compile(r"^[0-9a-f]{24}$")
SOURCE_ROLES = {
    "announcement_10",
    "methods_template",
    "official_rules",
    "space_config",
    "track2_overview",
    "track2_submission_ui",
}
SOURCE_PATHS = {
    "announcement_10": "discussions/10",
    "methods_template": "static/templates/methods_description_form.xlsx",
    "official_rules": "tabs/rules.py",
    "space_config": "config.py",
    "track2_overview": "tabs/about.py",
    "track2_submission_ui": "tabs/submit_track2.py",
}
DISCLOSURE_FIELDS = {
    "plan_or_tier",
    "provider",
    "relevant_data_handling_setting",
}
REQUIRED_SUBMISSION_FIELDS = {
    "github_repository_url",
    "pitch_video_url",
    "report",
    "team_or_submitter_identity",
}
ACKNOWLEDGEMENT = (
    "This work was made possible through the Hackathon, organized by Sage "
    "Bionetworks in partnership with the MVA Society, Hugging Face, and BEACON "
    "(The Benchmarking, Evaluation, and Assessment Consortium for Science), with "
    "prize sponsorship from AWS and Anthropic. We are deeply grateful to the child "
    "and their family who generously contributed their data and their story to "
    "advance research into this rare disease. We acknowledge their trust in making "
    "this Hackathon possible."
)
JUDGING_CRITERIA = {
    "scientific_rigor": (
        "Scientific Rigor",
        35,
        "Is the variant mechanism characterization sound, and is the proposed drug "
        "repurposing candidate well-supported by that mechanism?",
    ),
    "potential_impact": (
        "Potential Impact",
        25,
        "If independently validated, how much could this contribute to understanding "
        "or diagnosing MVA, for this child or for others living with it?",
    ),
    "innovation": (
        "Innovation",
        25,
        "Did you bring a genuinely creative angle, method, or tool to the problem?",
    ),
    "scalability": (
        "Scalability",
        15,
        "Could your approach realistically be applied to other individuals, other "
        "diseases, or larger datasets - not just this one case?",
    ),
}


class ChallengeRulesError(ValueError):
    """Raised when a challenge rules contract is incomplete or inconsistent."""


def _reject_unknown(data: Mapping[str, Any], allowed: set[str], context: str) -> None:
    unknown = sorted(set(data) - allowed)
    if unknown:
        raise ChallengeRulesError(f"unknown {context} fields: {', '.join(unknown)}")


def _required_string(data: Mapping[str, Any], key: str, context: str) -> str:
    value = data.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ChallengeRulesError(f"{context}.{key} must be a non-empty string")
    return value


@dataclass(frozen=True, slots=True)
class RuleSource:
    role: str
    path: str
    url: str
    source_id: str
    sha256: str
    size_bytes: int
    retrieval_url: str
    extraction: str

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "RuleSource":
        _reject_unknown(
            data,
            {
                "extraction",
                "path",
                "retrieval_url",
                "role",
                "sha256",
                "size_bytes",
                "source_id",
                "url",
            },
            "rule source",
        )
        source = cls(
            role=_required_string(data, "role", "source"),
            path=_required_string(data, "path", "source"),
            url=_required_string(data, "url", "source"),
            source_id=_required_string(data, "source_id", "source"),
            sha256=_required_string(data, "sha256", "source"),
            size_bytes=data.get("size_bytes"),
            retrieval_url=_required_string(data, "retrieval_url", "source"),
            extraction=_required_string(data, "extraction", "source"),
        )
        if source.role not in SOURCE_ROLES:
            raise ChallengeRulesError(f"unsupported rule source role: {source.role}")
        if not source.url.startswith("https://") or not source.retrieval_url.startswith(
            "https://"
        ):
            raise ChallengeRulesError(f"source {source.role} must use HTTPS")
        expected_extraction = (
            "hf_discussion_first_comment_raw"
            if source.role == "announcement_10"
            else "response_bytes"
        )
        if source.extraction != expected_extraction:
            raise ChallengeRulesError(
                f"source {source.role} extraction does not match its pinned source type"
            )
        if not SHA256.fullmatch(source.sha256):
            raise ChallengeRulesError(f"source {source.role} SHA-256 is invalid")
        if (
            not isinstance(source.size_bytes, int)
            or isinstance(source.size_bytes, bool)
            or source.size_bytes <= 0
        ):
            raise ChallengeRulesError(
                f"source {source.role} size_bytes must be a positive integer"
            )
        return source


@dataclass(frozen=True, slots=True)
class JudgingCriterion:
    key: str
    label: str
    weight_percent: int
    organizer_question: str

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "JudgingCriterion":
        _reject_unknown(
            data,
            {"key", "label", "organizer_question", "weight_percent"},
            "judging criterion",
        )
        criterion = cls(
            key=_required_string(data, "key", "judging criterion"),
            label=_required_string(data, "label", "judging criterion"),
            weight_percent=data.get("weight_percent"),
            organizer_question=_required_string(
                data, "organizer_question", "judging criterion"
            ),
        )
        expected = JUDGING_CRITERIA.get(criterion.key)
        if expected != (
            criterion.label,
            criterion.weight_percent,
            criterion.organizer_question,
        ):
            raise ChallengeRulesError("judging criterion differs from the pinned source")
        return criterion


@dataclass(frozen=True, slots=True)
class Track2RuleFacts:
    maximum_submission_attempts: int
    judge_active_submission: str
    designated_team_submitter: bool
    report_formats: tuple[str, ...]
    repository_url_prefix: str
    repository_public_by: str
    pitch_maximum_seconds: int
    pitch_providers: tuple[str, ...]
    methods_template_source: str
    judging_criteria_source: str
    judging_criteria: tuple[JudgingCriterion, ...]
    commercial_llm_disclosure_required_when_used: bool
    commercial_llm_disclosure_fields: tuple[str, ...]
    deadline_display: str
    submission_closes_at_utc_exclusive: str
    required_submission_fields: tuple[str, ...]
    submission_material_license: str
    submission_material_license_source: str
    organizer_rerun_allowed: bool
    public_name_attribution_required: bool
    embargo_begins: str
    embargo_ends: str
    peer_review_submission_during_embargo: str
    conference_preliminary_findings: str
    required_acknowledgement: str
    dataset_citation: str
    controlled_retention_days_after_close: int
    deletion_confirmation_email: str
    requirements_snapshot_freshness_hours: int
    requirements_attestation_public_key_sha256: str

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "Track2RuleFacts":
        _reject_unknown(
            data,
            {
                "commercial_llm_disclosure_fields",
                "commercial_llm_disclosure_required_when_used",
                "conference_preliminary_findings",
                "controlled_retention_days_after_close",
                "dataset_citation",
                "deadline_display",
                "deletion_confirmation_email",
                "designated_team_submitter",
                "embargo_begins",
                "embargo_ends",
                "judge_active_submission",
                "judging_criteria",
                "judging_criteria_source",
                "maximum_submission_attempts",
                "methods_template_source",
                "pitch_maximum_seconds",
                "pitch_providers",
                "peer_review_submission_during_embargo",
                "public_name_attribution_required",
                "required_acknowledgement",
                "required_submission_fields",
                "report_formats",
                "requirements_snapshot_freshness_hours",
                "requirements_attestation_public_key_sha256",
                "repository_public_by",
                "repository_url_prefix",
                "submission_closes_at_utc_exclusive",
                "submission_material_license",
                "submission_material_license_source",
                "organizer_rerun_allowed",
            },
            "Track 2 rule",
        )
        sequence_fields = (
            "commercial_llm_disclosure_fields",
            "pitch_providers",
            "report_formats",
            "required_submission_fields",
        )
        values: dict[str, tuple[str, ...]] = {}
        for key in sequence_fields:
            value = data.get(key)
            if (
                not isinstance(value, list)
                or not value
                or any(not isinstance(item, str) or not item for item in value)
            ):
                raise ChallengeRulesError(f"track2.{key} must be a non-empty string list")
            if len(value) != len(set(value)):
                raise ChallengeRulesError(f"track2.{key} values must be unique")
            values[key] = tuple(value)
        criteria_data = data.get("judging_criteria")
        if not isinstance(criteria_data, list):
            raise ChallengeRulesError("track2.judging_criteria must be a list")
        criteria = tuple(JudgingCriterion.from_mapping(item) for item in criteria_data)
        facts = cls(
            maximum_submission_attempts=data.get("maximum_submission_attempts"),
            judge_active_submission=_required_string(
                data, "judge_active_submission", "track2"
            ),
            designated_team_submitter=data.get("designated_team_submitter"),
            report_formats=values["report_formats"],
            repository_url_prefix=_required_string(
                data, "repository_url_prefix", "track2"
            ),
            repository_public_by=_required_string(
                data, "repository_public_by", "track2"
            ),
            pitch_maximum_seconds=data.get("pitch_maximum_seconds"),
            pitch_providers=values["pitch_providers"],
            methods_template_source=_required_string(
                data, "methods_template_source", "track2"
            ),
            judging_criteria_source=_required_string(
                data, "judging_criteria_source", "track2"
            ),
            judging_criteria=criteria,
            commercial_llm_disclosure_required_when_used=data.get(
                "commercial_llm_disclosure_required_when_used"
            ),
            commercial_llm_disclosure_fields=values[
                "commercial_llm_disclosure_fields"
            ],
            deadline_display=_required_string(data, "deadline_display", "track2"),
            submission_closes_at_utc_exclusive=_required_string(
                data, "submission_closes_at_utc_exclusive", "track2"
            ),
            required_submission_fields=values["required_submission_fields"],
            submission_material_license=_required_string(
                data, "submission_material_license", "track2"
            ),
            submission_material_license_source=_required_string(
                data, "submission_material_license_source", "track2"
            ),
            organizer_rerun_allowed=data.get("organizer_rerun_allowed"),
            public_name_attribution_required=data.get(
                "public_name_attribution_required"
            ),
            embargo_begins=_required_string(data, "embargo_begins", "track2"),
            embargo_ends=_required_string(data, "embargo_ends", "track2"),
            peer_review_submission_during_embargo=_required_string(
                data, "peer_review_submission_during_embargo", "track2"
            ),
            conference_preliminary_findings=_required_string(
                data, "conference_preliminary_findings", "track2"
            ),
            required_acknowledgement=_required_string(
                data, "required_acknowledgement", "track2"
            ),
            dataset_citation=_required_string(data, "dataset_citation", "track2"),
            controlled_retention_days_after_close=data.get(
                "controlled_retention_days_after_close"
            ),
            deletion_confirmation_email=_required_string(
                data, "deletion_confirmation_email", "track2"
            ),
            requirements_snapshot_freshness_hours=data.get(
                "requirements_snapshot_freshness_hours"
            ),
            requirements_attestation_public_key_sha256=_required_string(
                data,
                "requirements_attestation_public_key_sha256",
                "track2",
            ),
        )
        if facts.maximum_submission_attempts != 3:
            raise ChallengeRulesError("Track 2 permits exactly three submission attempts")
        if facts.judge_active_submission != "latest_submitted":
            raise ChallengeRulesError("only the latest Track 2 submission may be judge-active")
        if facts.designated_team_submitter is not True:
            raise ChallengeRulesError("Track 2 requires a designated team submitter")
        if set(facts.report_formats) != {"markdown", "pdf"}:
            raise ChallengeRulesError("Track 2 report formats must be Markdown and PDF")
        if facts.repository_url_prefix != "https://github.com/":
            raise ChallengeRulesError("Track 2 requires an HTTPS GitHub repository URL")
        if facts.repository_public_by != "final_evaluation_begins":
            raise ChallengeRulesError(
                "the repository must be public when final evaluation begins"
            )
        if facts.pitch_maximum_seconds != 180:
            raise ChallengeRulesError("the Track 2 pitch limit must be three minutes")
        if set(facts.pitch_providers) != {"vimeo", "youtube"}:
            raise ChallengeRulesError("the accepted pitch providers must be YouTube and Vimeo")
        if facts.methods_template_source != "methods_template":
            raise ChallengeRulesError("Track 2 must reference the pinned methods template")
        if facts.judging_criteria_source != "track2_overview":
            raise ChallengeRulesError("Track 2 must reference the pinned overview criteria")
        if (
            tuple(item.key for item in facts.judging_criteria)
            != tuple(JUDGING_CRITERIA)
            or sum(item.weight_percent for item in facts.judging_criteria) != 100
        ):
            raise ChallengeRulesError("Track 2 judging criteria are incomplete or reordered")
        if facts.commercial_llm_disclosure_required_when_used is not True:
            raise ChallengeRulesError(
                "commercial LLM disclosure must be required when such a system was used"
            )
        if set(facts.commercial_llm_disclosure_fields) != DISCLOSURE_FIELDS:
            raise ChallengeRulesError(
                "commercial LLM disclosure must include provider, plan/tier, and data handling"
            )
        if facts.deadline_display != "October 24, 2026 at 23:59 UTC":
            raise ChallengeRulesError("Track 2 deadline display differs from the current interface")
        try:
            deadline = datetime.fromisoformat(
                facts.submission_closes_at_utc_exclusive.replace("Z", "+00:00")
            )
        except ValueError as error:
            raise ChallengeRulesError("Track 2 deadline must be an ISO UTC timestamp") from error
        if (
            deadline != datetime(2026, 10, 25, tzinfo=timezone.utc)
            or set(facts.required_submission_fields) != REQUIRED_SUBMISSION_FIELDS
        ):
            raise ChallengeRulesError("Track 2 deadline or required fields have drifted")
        expected_terms = (
            facts.submission_material_license == "CC-BY-4.0"
            and facts.submission_material_license_source == "track2_overview"
            and facts.organizer_rerun_allowed is True
            and facts.public_name_attribution_required is True
            and facts.embargo_begins == "hackathon_close"
            and facts.embargo_ends == "organizer_summary_report_or_preprint_publication"
            and facts.peer_review_submission_during_embargo == "prohibited"
            and facts.conference_preliminary_findings == "prior_written_approval_required"
            and facts.required_acknowledgement == ACKNOWLEDGEMENT
            and facts.dataset_citation
            == "hackathon_synapse_page_reference_at_time_of_publication"
            and facts.controlled_retention_days_after_close == 30
            and facts.deletion_confirmation_email
            == "RarediseaserealkidMVAhackathon2026@synapse.org"
            and facts.requirements_snapshot_freshness_hours == 24
            and SHA256.fullmatch(
                facts.requirements_attestation_public_key_sha256
            )
            is not None
        )
        if not expected_terms:
            raise ChallengeRulesError("Track 2 publication or retention terms have drifted")
        return facts


@dataclass(frozen=True, slots=True)
class RulePrecedence:
    authority: str
    supersedes: str
    subject: str
    effective_value: Any

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "RulePrecedence":
        _reject_unknown(
            data,
            {"authority", "effective_value", "subject", "supersedes"},
            "precedence",
        )
        precedence = cls(
            authority=_required_string(data, "authority", "precedence"),
            supersedes=_required_string(data, "supersedes", "precedence"),
            subject=_required_string(data, "subject", "precedence"),
            effective_value=data.get("effective_value"),
        )
        expected = {
            "track2.maximum_submission_attempts": ("methods_template", 3),
            "track2.repository_public_by": (
                "track2_submission_ui",
                "final_evaluation_begins",
            ),
        }.get(precedence.subject)
        if (
            precedence.authority != "announcement_10"
            or expected is None
            or (precedence.supersedes, precedence.effective_value) != expected
        ):
            raise ChallengeRulesError("announcement 10 precedence is invalid")
        return precedence


@dataclass(frozen=True, slots=True)
class ChallengeRulesContract:
    challenge: str
    space_repository: str
    space_revision: str
    verified_on: str
    sources: tuple[RuleSource, ...]
    track2: Track2RuleFacts
    precedence: tuple[RulePrecedence, ...]
    schema: str = CHALLENGE_RULES_SCHEMA

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "ChallengeRulesContract":
        _reject_unknown(
            data,
            {
                "challenge",
                "precedence",
                "schema",
                "sources",
                "space_repository",
                "space_revision",
                "track2",
                "verified_on",
            },
            "challenge rules",
        )
        schema = _required_string(data, "schema", "rules")
        if schema != CHALLENGE_RULES_SCHEMA:
            raise ChallengeRulesError(
                f"rules schema must be {CHALLENGE_RULES_SCHEMA!r}"
            )
        revision = _required_string(data, "space_revision", "rules")
        if not GIT_ID.fullmatch(revision):
            raise ChallengeRulesError("space_revision must be a full Git commit ID")
        space_repository = _required_string(data, "space_repository", "rules")
        if space_repository != "SageBio/rare-disease-real-kid-mva-hackathon-2026":
            raise ChallengeRulesError("challenge Space repository is invalid")
        verified_on = _required_string(data, "verified_on", "rules")
        try:
            date.fromisoformat(verified_on)
        except ValueError as error:
            raise ChallengeRulesError("verified_on must be an ISO date") from error
        source_data = data.get("sources")
        if not isinstance(source_data, list):
            raise ChallengeRulesError("sources must be a list")
        sources = tuple(RuleSource.from_mapping(item) for item in source_data)
        roles = [source.role for source in sources]
        if set(roles) != SOURCE_ROLES or len(roles) != len(SOURCE_ROLES):
            raise ChallengeRulesError("rules contract must pin each required source exactly once")
        for source in sources:
            if source.path != SOURCE_PATHS[source.role]:
                raise ChallengeRulesError(f"source {source.role} path has drifted")
            base = f"https://huggingface.co/spaces/{space_repository}"
            if source.role == "announcement_10":
                expected_url = f"{base}/discussions/10"
                expected_retrieval = (
                    f"https://huggingface.co/api/spaces/{space_repository}"
                    "/discussions/10"
                )
                if (
                    source.url != expected_url
                    or source.retrieval_url != expected_retrieval
                    or not EVENT_ID.fullmatch(source.source_id)
                ):
                    raise ChallengeRulesError("announcement source identity is invalid")
            else:
                action = "resolve" if source.role == "methods_template" else "raw"
                expected_url = f"{base}/{action}/{revision}/{source.path}"
                if (
                    source.url != expected_url
                    or source.retrieval_url != expected_url
                    or not GIT_ID.fullmatch(source.source_id)
                ):
                    raise ChallengeRulesError(
                        f"source {source.role} must pin its exact Git artifact"
                    )
        track2_data = data.get("track2")
        if not isinstance(track2_data, Mapping):
            raise ChallengeRulesError("track2 must be an object")
        precedence_data = data.get("precedence")
        if not isinstance(precedence_data, list) or len(precedence_data) != 2:
            raise ChallengeRulesError("two explicit precedence rules are required")
        precedence = tuple(RulePrecedence.from_mapping(item) for item in precedence_data)
        if {item.subject for item in precedence} != {
            "track2.maximum_submission_attempts",
            "track2.repository_public_by",
        }:
            raise ChallengeRulesError("required announcement precedence subjects are missing")
        return cls(
            schema=schema,
            challenge=_required_string(data, "challenge", "rules"),
            space_repository=space_repository,
            space_revision=revision,
            verified_on=verified_on,
            sources=sources,
            track2=Track2RuleFacts.from_mapping(track2_data),
            precedence=precedence,
        )


def load_challenge_rules(path: str | Path) -> ChallengeRulesContract:
    rules_path = Path(path)
    try:
        data, _content = load_bounded_json(rules_path)
    except BoundedJSONError as error:
        if "cannot be read" not in str(error):
            raise ChallengeRulesError(f"invalid JSON in {rules_path}: {error}") from error
        raise ChallengeRulesError(f"cannot read challenge rules {rules_path}: {error}") from error
    if not isinstance(data, Mapping):
        raise ChallengeRulesError("challenge rules root must be an object")
    return ChallengeRulesContract.from_mapping(data)
