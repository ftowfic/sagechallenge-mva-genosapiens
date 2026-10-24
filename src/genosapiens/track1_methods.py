"""Deterministic controlled Track 1 methods report and operator checklist."""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
import re
from typing import Literal


SHA256 = re.compile(r"^[0-9a-f]{64}$")
REVISION = re.compile(r"^[0-9a-f]{40,64}$")
SAFE_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,79}$")


class Track1MethodsError(ValueError):
    """Raised when report inputs are malformed or bindings are inconsistent."""


def _text(value: str | None) -> str:
    return "unconfirmed" if value is None else value


def _clean_optional(value: str | None, label: str) -> None:
    if value is not None and (
        not value.strip() or len(value) > 4000 or any(ord(char) < 32 and char not in "\n\t" for char in value)
    ):
        raise Track1MethodsError(f"invalid {label}")


@dataclass(frozen=True, slots=True)
class BoundedWgsMethodsDisclosure:
    """Required policy and qualification disclosure for the bounded WGS path."""

    early_policy_id: str
    early_policy_version: str
    early_policy_sha256: str
    early_population_source: str
    early_population_release: str
    early_population_metric: str
    early_threshold: str
    exception_scope: str
    wgs_coverage: str
    unknown_retention: str
    phased_bypass: str
    native_bypass: str
    downstream_population_limitation: str
    storage_qualification: str
    resource_qualification: str
    predecessor_effect: str
    scope_policy_sha256: str = field(kw_only=True)
    scope_manifest_sha256: str = field(kw_only=True)

    def __post_init__(self) -> None:
        if any(SHA256.fullmatch(value) is None for value in
               (self.early_policy_sha256, self.scope_policy_sha256, self.scope_manifest_sha256)):
            raise Track1MethodsError("invalid bounded WGS policy identity")
        for label, value in (
            ("early policy ID", self.early_policy_id),
            ("early policy version", self.early_policy_version),
            ("early population source", self.early_population_source),
            ("early population release", self.early_population_release),
            ("early population metric", self.early_population_metric),
            ("early threshold", self.early_threshold),
            ("exception scope", self.exception_scope),
            ("WGS coverage", self.wgs_coverage),
            ("unknown retention", self.unknown_retention),
            ("phased bypass", self.phased_bypass),
            ("native bypass", self.native_bypass),
            ("downstream population limitation", self.downstream_population_limitation),
            ("storage qualification", self.storage_qualification),
            ("resource qualification", self.resource_qualification),
            ("predecessor effect", self.predecessor_effect),
        ):
            _clean_optional(value, label)
            if not value.strip():
                raise Track1MethodsError(f"invalid {label}")


def _packet_name_component(value: str, label: str) -> str:
    """Keep submission display names intact and normalize spaces for filenames."""
    component = re.sub(r" +", "-", value.strip(" "))
    if SAFE_NAME.fullmatch(component) is None:
        raise Track1MethodsError(f"{label} is unsafe for a packet filename")
    return component


@dataclass(frozen=True, slots=True)
class Track1MethodsFacts:
    team_name: str | None
    approach_name: str | None
    submitter_registration_confirmed: bool | None
    repository_url: str | None
    code_revision: str
    automation_statement: str | None
    manual_review_statement: str | None
    proprietary_source_statement: str | None
    ai_provider_plan_data_use: str | None
    runtime_and_cost: str | None
    required_acknowledgement: str | None
    dataset_citation: str | None
    retention_coordinator: str | None
    video_status: Literal["operator_confirmed", "unconfirmed"]

    def __post_init__(self) -> None:
        for label, value in (
            ("team name", self.team_name),
            ("approach name", self.approach_name),
            ("repository URL", self.repository_url),
            ("automation statement", self.automation_statement),
            ("manual review statement", self.manual_review_statement),
            ("proprietary source statement", self.proprietary_source_statement),
            ("AI disclosure", self.ai_provider_plan_data_use),
            ("runtime and cost", self.runtime_and_cost),
            ("acknowledgement", self.required_acknowledgement),
            ("dataset citation", self.dataset_citation),
            ("retention coordinator", self.retention_coordinator),
        ):
            _clean_optional(value, label)
        if not REVISION.fullmatch(self.code_revision):
            raise Track1MethodsError("invalid reproducible code revision")
        if self.team_name is not None:
            _packet_name_component(self.team_name, "team name")
        if self.approach_name is not None:
            _packet_name_component(self.approach_name, "approach name")
        if self.repository_url is not None and not self.repository_url.startswith("https://github.com/"):
            raise Track1MethodsError("repository URL must identify GitHub")
        if self.submitter_registration_confirmed not in {True, False, None}:
            raise Track1MethodsError("invalid registration status")


@dataclass(frozen=True, slots=True)
class Track1MethodsBindings:
    run_id: str
    csv_sha256: str | None
    audit_json_sha256: str
    audit_markdown_sha256: str
    preflight_sha256: str
    selected_row_count: int
    software_ready: bool
    public_sources: tuple[str, ...]
    license_notices: tuple[str, ...]

    def __post_init__(self) -> None:
        hashes = (
            self.audit_json_sha256,
            self.audit_markdown_sha256,
            self.preflight_sha256,
        )
        if (
            not self.run_id.startswith("track1-run")
            or any(SHA256.fullmatch(value) is None for value in hashes)
            or self.csv_sha256 is not None and SHA256.fullmatch(self.csv_sha256) is None
            or not 0 <= self.selected_row_count <= 10
            or (self.csv_sha256 is None) != (self.selected_row_count == 0)
            or not self.public_sources
            or self.public_sources != tuple(sorted(set(self.public_sources)))
            or self.license_notices != tuple(sorted(set(self.license_notices)))
        ):
            raise Track1MethodsError("invalid run-bound methods inputs")


@dataclass(frozen=True, slots=True)
class Track1OperatorPacket:
    methods_markdown: str
    methods_sha256: str
    csv_filename: str | None
    report_filename: str | None
    software_ready: bool
    video_status: Literal["operator_confirmed", "unconfirmed"]
    ready_for_operator_upload: bool
    blockers: tuple[str, ...]


def build_track1_operator_packet(
    facts: Track1MethodsFacts,
    bindings: Track1MethodsBindings,
    *,
    bounded_wgs: BoundedWgsMethodsDisclosure | None = None,
    analysis_profile: Literal["focused", "exhaustive"] = "focused",
    profile_source_scope: dict | None = None,
) -> Track1OperatorPacket:
    """Render methods and a checklist; this function never uploads a submission."""
    if (
        type(facts) is not Track1MethodsFacts
        or type(bindings) is not Track1MethodsBindings
        or bounded_wgs is not None
        and type(bounded_wgs) is not BoundedWgsMethodsDisclosure
        or type(analysis_profile) is not str or analysis_profile not in ("focused", "exhaustive")
    ):
        raise Track1MethodsError("typed methods inputs are required")
    source_scope_text = 'unconfirmed. No authenticated source-wide accounting was supplied.'
    if profile_source_scope is not None:
        from genosapiens.wgs_reporting import validate_profile_source_accounting, WgsReportingError
        try:
            scope = validate_profile_source_accounting(profile_source_scope)
            source_scope_text = json.dumps(scope, sort_keys=True, separators=(',', ':'))
        except WgsReportingError:
            raise Track1MethodsError('invalid profile source accounting') from None
    missing = []
    required = {
        "team-name-unconfirmed": facts.team_name,
        "approach-name-unconfirmed": facts.approach_name,
        "registration-unconfirmed": facts.submitter_registration_confirmed is True,
        "repository-url-unconfirmed": facts.repository_url,
        "automation-unconfirmed": facts.automation_statement,
        "manual-review-unconfirmed": facts.manual_review_statement,
        "proprietary-sources-unconfirmed": facts.proprietary_source_statement,
        "ai-provider-plan-data-use-unconfirmed": facts.ai_provider_plan_data_use,
        "runtime-cost-unconfirmed": facts.runtime_and_cost,
        "acknowledgement-unconfirmed": facts.required_acknowledgement,
        "dataset-citation-unconfirmed": facts.dataset_citation,
        "retention-coordinator-unconfirmed": facts.retention_coordinator,
    }
    missing.extend(key for key, value in required.items() if not value)
    if not bindings.software_ready:
        missing.append("software-preflight-not-ready")
    if bindings.selected_row_count == 0:
        missing.append("zero-eligible-hypotheses")
    if facts.video_status != "operator_confirmed":
        missing.append("required-video-unconfirmed")
    blockers = tuple(sorted(set(missing)))

    sources = "\n".join(f"- {item}" for item in bindings.public_sources)
    licenses = (
        "\n".join(f"- {item}" for item in bindings.license_notices)
        if bindings.license_notices
        else "- unconfirmed"
    )
    from genosapiens.track1_association_scope import POLICY_SHA256 as association_scope_sha256
    bounded_section = ""
    if bounded_wgs is not None:
        bounded_section = f"""
## Bounded whole-genome eligibility policy

The global early scope uses the complete Track 1 disease catalog.
It does not select genes from the case diagnosis or HPO terms.
It includes full Open Targets gene spans without padding and all explicit stable regulatory links to catalog genes.
It also includes exact alleles accepted by the shared known-disease route, including gene-free assertions.

- Early scope policy SHA-256: `{bounded_wgs.scope_policy_sha256}`.
- Early scope manifest SHA-256: `{bounded_wgs.scope_manifest_sha256}`.

Source checks, normalization, and accounting cover the complete input.
Outside-scope calls and calls with unresolved scope receive no preparation population or clinical queries.
They receive no annotation or admission unless duplicate closure advances their normalized identity.
These calls retain their original source and scope evidence.
If an in-scope occurrence advances a normalized identity, that identity preserves all original occurrence mappings.

Open Targets boundaries can omit transcript regions recognized by VEP.
Missing gene coordinates and unresolved native geometry reduce coverage.
Scope deferral does not establish benignity, a missing partner, a diagnosis, or a treatment recommendation.
Existing phase and native population protections apply only after scope admission.
Unknown population evidence retains its existing treatment after scope admission.

Diagnostics separate provisional scope counts from effective counts after duplicate closure.
Preparation request-key counts cover successfully committed batches.
They do not measure physical database reads or all queries during the run.
Public scope construction remains separate from case execution and supports reuse.
Track 1 uses a global disease scope that can exceed Track 2's diagnosis-specific scope.
No runtime ratio or genome benchmark establishes acceptance.

- Whole-genome coverage: {bounded_wgs.wgs_coverage}
- Early policy: `{bounded_wgs.early_policy_id}` version `{bounded_wgs.early_policy_version}`.
- Early policy SHA-256: `{bounded_wgs.early_policy_sha256}`.
- Population evidence: {bounded_wgs.early_population_source}, release {bounded_wgs.early_population_release}, metric {bounded_wgs.early_population_metric}.
- Early threshold: {bounded_wgs.early_threshold}.
- Exception scope: {bounded_wgs.exception_scope}.
- Unknown evidence: {bounded_wgs.unknown_retention}.
- Phased/composite bypass: {bounded_wgs.phased_bypass}.
- Native-class bypass: {bounded_wgs.native_bypass}.
- Downstream population limitation: {bounded_wgs.downstream_population_limitation}.
- Storage evidence: {bounded_wgs.storage_qualification}.
- Resource versions: {bounded_wgs.resource_qualification}.
- Effect of early filtering: {bounded_wgs.predecessor_effect}.
- Focused association scope SHA-256: `{association_scope_sha256}`.

Source inventory complete means every original input and source occurrence has an authenticated disposition.
Screened domain completely evaluated means every relevant pair and generated hypothesis in the retained domain receives scientific assessment.
The workflow applies no ranking cutoff within that domain.
Broader rescue not evaluated means deferred or excluded inputs cannot establish benignity or the absence of a partner.

Scientific records preserve source conflicts, negative evidence, phase constraints, and incomplete-genotype penalties.
Context membership closes across all batches before pair enumeration.
Selected members retain their technical checks and original allele representations.
Regulatory evidence preserves the highest finite score per element, target, and method, with stable source-location ties and omitted-row counts.

Full-input validation audits remain separate requested work.
Ordinary execution requires no qualification campaign or previous attempt history.
Completed stages bind their direct code, policy, source, and parent dependencies in the current run manifest.
Changed selection evidence reuses compatible annotation.
Changed report code reuses compatible scientific results.
Accepted reports remain immutable.

Measured work separates annotation, selection, discovery, context formation, scientific assessment, publication, and verification.
Counts distinguish original inputs, retained inputs, deferred inputs, source occurrences, contexts, and pairs.
Missing measurements remain unconfirmed.
Finite measurements do not establish whole-genome runtime or clinical validation.
Stopped checks remain incomplete.
No hard total-runtime or speedup threshold applies.

The early population policy and the retention policy define different decisions.
Early exclusions and later deferrals can change ranks by changing the candidate universe.
Equal-universe scientific parity does not establish universal output equality.
"""
    methods = f"""# Track 1 methods description

Classification: `controlled`.
LLM policy: `no-llm`.

## Submission identity

- Team: {_text(facts.team_name)}
- Approach: {_text(facts.approach_name)}
- Run: `{bindings.run_id}`
- Repository: {_text(facts.repository_url)}
- Reproducible code revision: `{facts.code_revision}`

## Method abstract

The workflow checks supplied whole-genome calls against the reconstructed reference.
It applies quality and population rules before initial annotation and exact public-evidence lookup.
It retains alleles through the four routes below.
It constructs complete or explicitly incomplete genotype hypotheses from retained alleles and ranks them with audit records.
Primary ranking is diagnosis blind.
Optional diagnosis-informed and dominant-term sensitivities require separate requests and publications.
The competitive CSV contains at most ten primary exact singleton or unordered-pair representations.

Source-bound identities, explicit uncertainty, deterministic replay, and public-reference checks support reproducibility.
EPCR values provide an uncalibrated relative ordering.
Public knowledge remains incomplete.
Some native variant classes cannot use the evaluator CSV format.
Missing family or onset information remains unknown.
Parser conformance does not establish agreement with private truth.
{bounded_section}

## Analysis profile and source limitations

Profile metadata: `{analysis_profile}`.
The focused profile uses one shared versioned retention policy.
It retains an allele when one or more routes apply:

1. Exact pathogenic or likely pathogenic clinical evidence, or exact curated disease-relevant damaging functional evidence.
2. Applicable deleterious SIFT or PolyPhen evidence in an established disease gene.
3. A supplied structural call within the declared scope, with its native identity and support limits.
4. A declared coding or splice consequence on a relevant transcript in an established disease gene.

The policy evaluates all applicable transcripts.
A LOW impact label does not defeat a qualifying splice consequence.
An uncertain prediction does not defeat a qualifying missense consequence.
An exact known disease allele can qualify outside a gene interval or without clinical gene metadata.
Conflicting evidence keeps its source meaning.
Retention does not establish pathogenicity, a diagnosis, a mechanism, or competitive eligibility.
Diagnosis labels and phenotype similarity do not select Track 1 alleles.

No qualifying route means deferred_no_retention_evidence.
Unknown effects, ordinary synonymous effects, and generic regulatory overlap do not qualify alone.
Repeat, mitochondrial, and mosaic calls also enter the deferred inventory when no specified route applies.
Deferred alleles enter no contexts, genotype pairs, scores, mechanisms, or drug inputs.
The compact inventory preserves source occurrences, identity, annotation summary, reasons, and policy and source versions.
Reports count identities and occurrences separately.
The profile manifest binds the retained universe and aggregate deferred scope.
A deferred partner can exist in the source genome.
Selection can omit novel noncoding causes and unknown effects without qualifying evidence.

Authenticated source-wide accounting: `{source_scope_text}`.
These source counts describe the pinned resource, independently of case targets and retained contexts.
Complete parsing does not establish target resolution or causality.
Regulatory source records preserve stable, unmapped, and missing targets without invented nearest-gene mappings.
Missing required resources or failed queries stop selection.
They cannot become complete no-match results.
Optional unavailable evidence remains explicitly unavailable.

Population work exclusions require complete source, frequency, protection, and consuming-model checks under the strict >5% policy.
Unknown evidence and positive rescue evidence retain their existing scientific treatment.
Mixed pairs preserve original partners when omission would change scientific values.
Skipped work has no computed hypothesis count or score.
The workflow reconstructs selected records from saved values.
Optimized scientific execution may omit validation and verification checks under user authorization.
The scientific scope report identifies omitted evidence verification.
It does not establish whole-genome runtime.

## Variant pairs and other findings

The workflow emits exact compound-heterozygous pairs when both alleles pass representation and eligibility checks.
Complete pairs stay in one row.
An eligible incomplete recessive hypothesis can contribute one explicitly labeled partial-member prediction.
That prediction retains its incompleteness penalty, unresolved-partner disclosure, and selected-member technical and eligibility checks.
It does not require a no-partner certificate or establish a complete genotype.
Other assessed candidates remain in the complete evidence and ordering records.
Their exhaustive narrative interpretation remains optional.

## Challenge checklist for the human operator

Reviewed organizer guidance dated 2026-10-10 remains separate from software readiness and human completion.
The operator owns prerequisite checks, header-scope confirmation, controlled preparation, analysis, review, submission, and retention closeout.
O1–O5 completion requires operator confirmation.

- O4 submission: At most ten primary-only GRCh38 exact singleton or full-pair rows, with EPCR in (0,1].
- O4 attachments: PDF or Markdown methods and reproducible GitHub code and configuration.
- O4 filenames: Use consistent participant or team names and approach names.
- O4 attempt limit: Six submissions per participant. The best score counts.
- O4 evaluator: The public evaluator includes secondary rows despite contrary submission prose. Competitive output remains primary-only.
- O4 disclosure: Name actual AI tools, models, provider, product tier, settings, and data handling.
- O4 attribution: Include required acknowledgments, dataset citation, and applicable CC BY notices.
- O4 privacy: The controlled workflow remains offline and no-LLM despite conditional organizer permission for external processing.
- O4 public code: Publish reproducible code and configuration by closure or final evaluation.
- O4 restricted data: Do not publish genome-bearing inputs or full controlled audit records.
- O4 dates: The reviewed submission close is October 24, 2026 at 23:59 UTC.
- O4 later dates: Judging continues through December 17. The organizer announces winners on December 18. These dates do not extend submission or retention.
- O4 methods updates: Methods-only repository updates to an existing submission do not consume another scoring attempt under the reviewed organizer reply.
- O4 video: General rules mention a pitch video. Track 1 instructions omit it. Applicability remains unresolved, without a waiver or completion claim.
- O4 incidental findings: Return and public disclosure remain unresolved. Direct family contact and blanket case publication remain unauthorized.
- O4 unanswered procedures: An unanswered procedural question alone does not block software readiness.
- O5 deletion: Delete genome-bearing originals, per-variant derivatives, reconstruction-capable evidence, caches, and genomic-data logs within 30 days of the organizer-defined Hackathon close.
- O5 confirmation: Each participant must provide the required notification and confirmation.
- O5 deadline: Resolve the close-date basis before assigning a calendar deadline. Later judging does not grant an extension.
- O5 retained reports: Reviewed findings, HPO terms, and reports without prohibited genome content may remain after content classification.
- O5 disclosure: Retention permission does not authorize public disclosure. This packet asserts no deletion, notification, or closeout.

Official sources: [Track 1 instructions](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/submit_track1.py), [methods/AI clarification #10](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/discussions/10), [repository clarification #4](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/discussions/4), [timeline #26](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/discussions/26), [retention clarification #2](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/discussions/2) and [incidental findings #20](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/discussions/20).

## Automation and manual review

- Automated output: {_text(facts.automation_statement)}
- Downstream manual review: {_text(facts.manual_review_statement)}

## Data sources

Public sources:

{sources}

Proprietary sources: {_text(facts.proprietary_source_statement)}

License notices:

{licenses}

## Runtime, cost, and AI disclosure

- Runtime and cost: {_text(facts.runtime_and_cost)}
- Commercial generative-AI provider, plan, and data-use setting: {_text(facts.ai_provider_plan_data_use)}

The workflow sends no case-derived phenotype, genomic content, or patient output to an LLM.
The report and packet remain controlled/no-LLM artifacts for authorized human review and official submission.

## Required submission statements

- Registration confirmation: `{"confirmed" if facts.submitter_registration_confirmed is True else "unconfirmed"}`
- Required acknowledgement: {_text(facts.required_acknowledgement)}
- Dataset citation: {_text(facts.dataset_citation)}
- Retention coordinator: {_text(facts.retention_coordinator)}
- Required video: {facts.video_status}

## Artifact bindings

- CSV SHA-256: `{bindings.csv_sha256 if bindings.csv_sha256 is not None else "unavailable-zero-rows"}`
- JSON audit SHA-256: `{bindings.audit_json_sha256}`
- Markdown audit SHA-256: `{bindings.audit_markdown_sha256}`
- Preflight SHA-256: `{bindings.preflight_sha256}`

## Glossary

- Early scope: the public regions and exact known alleles eligible for preparation queries.
- Scope deferral: exclusion from expensive processing before a population decision.
- Source occurrence: one alternate slot in an original source record.
- Normalized identity: the existing sequence allele identity after reference normalization.
- Duplicate closure: retention of all mappings when one occurrence advances a normalized identity.
- Preparation request key: one allele key submitted to a population or clinical provider during preparation.

- Allele: A specific sequence change or native structural event.
- Retained universe: The alleles that satisfy the shared retention policy after upstream checks.
- Deferred inventory: Compact records for alleles without a qualifying retention route.
- EPCR: An uncalibrated relative value used to order eligible predictions.
- Partial-member prediction: An eligible member of an incomplete recessive hypothesis, with its original incompleteness penalty.
- HPO: Human Phenotype Ontology.
- LLM: Large language model.
- VCF: Variant Call Format.
- CSV: Comma-separated values.
- SHA-256: A digest algorithm that identifies exact bytes.

Identifiers, structured source accounting, supplied statements, and official source titles retain their original form for traceability.
This language review does not establish full ASD-STE100 compliance against the official approved-word dictionary.
"""
    content = methods.encode("utf-8")
    names_available = facts.team_name is not None and facts.approach_name is not None
    stem = None if not names_available else (
        _packet_name_component(facts.team_name, "team name")
        + "-" + _packet_name_component(facts.approach_name, "approach name")
    )
    return Track1OperatorPacket(
        methods,
        hashlib.sha256(content).hexdigest(),
        None if stem is None or bindings.csv_sha256 is None else stem + ".csv",
        None if stem is None else stem + "-methods.md",
        bindings.software_ready,
        facts.video_status,
        not blockers,
        blockers,
    )
