"""One terminal outcome for operator phases and public command entry points."""
from dataclasses import asdict, dataclass
import os
from pathlib import Path


@dataclass(frozen=True, slots=True)
class TerminalOutcome:
    status: str
    exit_code: int
    reason_codes: tuple[str, ...]

    def as_dict(self):
        return asdict(self)


def assessment_outcome(readiness, *, require_review=False):
    status=readiness['assessment_status']
    if require_review and not readiness['input_report_ready']:
        return TerminalOutcome('awaiting_review',2,('input_review_required',))
    if status=='not_performed':
        return TerminalOutcome('assessment_not_performed',2,('no_supported_objective_reached_search',))
    if status=='incomplete':
        return TerminalOutcome('assessment_incomplete',2,('required_assessment_incomplete',))
    if status!='complete':
        raise ValueError('Assessment status is undefined.')
    if not readiness['research_report_ready']:
        return TerminalOutcome('awaiting_review',2,('required_review_incomplete',))
    return TerminalOutcome('completed_report',0,('assessment_and_report_checks_complete',))


def publish_outcome(outcome):
    """Write only to the isolated attempt destination supplied by the launcher."""
    from genosapiens.track2_operator_inputs import write_new
    destination=os.environ.get('GENOSAPIENS_TRACK2_OUTCOME')
    if destination:
        write_new(Path(destination),outcome.as_dict())
    return outcome.exit_code


def pending_review(reason):
    return publish_outcome(TerminalOutcome('awaiting_review',2,(reason,)))
