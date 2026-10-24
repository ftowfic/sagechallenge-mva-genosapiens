"""Policy-bound relative EPCR mechanics, not calibration or controlled-use freeze.

The selected primary scores are replay products from task 11.2. The selected
numerical presentation policy remains distinct from a production candidate:
task 9.6 must bind it with the scientific policies before held-out execution and
task 13.1 must bind accepted results before controlled use.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, InvalidOperation, localcontext
from fractions import Fraction
import hashlib
from typing import Literal

from genosapiens.evidence_model import _Controlled, _exact, _ids, _sha, _typed
from genosapiens.evidence_serialization import _canonical, _loads
from genosapiens.track1_submission import SubmissionSelection


SCHEMA = "genosapiens.track1-epcr-policy"
RESULT_SCHEMA = "genosapiens.track1-epcr-assignment"
EXECUTABLE_VERSION = "leader-relative-decimal-exp-rank-decay"
DECIMAL_PRECISION = 80
TIE_RULE = "fixed-multiplicative-rank-decay"
CONFIDENCE_LABEL = "uncalibrated-relative-confidence-not-clinical-probability"
_ERRORS = (
    ValueError,
    TypeError,
    RuntimeError,
    AttributeError,
    KeyError,
    IndexError,
    InvalidOperation,
    OverflowError,
    RecursionError,
)


class Track1EpcrError(ValueError):
    """Static diagnostics only; scores and assignments remain controlled."""


def _bytes(value: object) -> bytes:
    return (_canonical(value) + "\n").encode("utf-8")


def _fraction_decimal(value: Fraction) -> Decimal:
    return Decimal(value.numerator) / Decimal(value.denominator)


def _fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def _decimal_text(value: Decimal) -> str:
    if not value.is_finite() or value <= 0 or value > 1:
        raise Track1EpcrError("EPCR is outside the supported numerical range")
    rendered = format(value, "f")
    if "." in rendered:
        rendered = rendered.rstrip("0").rstrip(".")
    _exact(rendered)
    return rendered


@dataclass(frozen=True, slots=True, repr=False)
class EpcrPolicy(_Controlled):
    """Caller-pinned mechanics; this object is not evidence of a public freeze."""

    policy_id: str
    score_gap_scale: str
    rank_decay: str
    minimum_serialized_epcr: str
    maximum_supported_score_gap: str
    maximum_rows: Literal[10]
    derivation_sha256: str
    derivation_scope: Literal["public", "independently_synthetic"]

    def __post_init__(self):
        _Controlled.__post_init__(self)
        _ids((self.policy_id,))
        _sha(self.derivation_sha256)
        scale = _exact(self.score_gap_scale)
        decay = _exact(self.rank_decay)
        minimum = _exact(self.minimum_serialized_epcr)
        maximum_gap = _exact(self.maximum_supported_score_gap)
        if (
            scale <= 0
            or not 0 < decay < 1
            or not 0 < minimum < 1
            or maximum_gap <= 0
            or self.maximum_rows != 10
        ):
            raise Track1EpcrError("invalid EPCR mechanics policy")
        with localcontext() as context:
            context.prec = DECIMAL_PRECISION
            endpoint_base = (-_fraction_decimal(scale * maximum_gap)).exp()
            for base in (Decimal(1), endpoint_base):
                values = tuple(
                    base * _fraction_decimal(decay) ** offset
                    for offset in range(self.maximum_rows)
                )
                if (
                    values[-1] < _fraction_decimal(minimum)
                    or any(
                        float(left) <= float(right)
                        for left, right in zip(values, values[1:])
                    )
                ):
                    raise Track1EpcrError(
                        "EPCR policy is not float-distinct at its boundaries"
                    )


@dataclass(frozen=True, slots=True, repr=False)
class EpcrValue(_Controlled):
    hypothesis_id: str
    serialized_rank: int
    source_rank: int
    source_score: str
    leader_score_gap: str
    base_epcr: str
    epcr: str

    def __post_init__(self):
        _Controlled.__post_init__(self)
        _ids((self.hypothesis_id,))
        score = _exact(self.source_score)
        gap = _exact(self.leader_score_gap)
        base = _exact(self.base_epcr)
        epcr = _exact(self.epcr)
        if (
            self.serialized_rank < 1
            or self.source_rank < 1
            or gap < 0
            or not 0 < base <= 1
            or not 0 < epcr <= base
            or self.serialized_rank == 1
            and (gap != 0 or base != 1 or epcr != 1)
        ):
            raise Track1EpcrError("invalid EPCR value")


@dataclass(frozen=True, slots=True, repr=False)
class EpcrAssignment(_Controlled):
    selection_sha256: str
    policy_id: str
    policy_sha256: str
    derivation_sha256: str
    derivation_scope: Literal["public", "independently_synthetic"]
    tie_rule: Literal["fixed-multiplicative-rank-decay"]
    decimal_precision: Literal[80]
    values: tuple[EpcrValue, ...]
    confidence_label: Literal[
        "uncalibrated-relative-confidence-not-clinical-probability"
    ] = CONFIDENCE_LABEL
    schema: Literal["genosapiens.track1-epcr-assignment"] = RESULT_SCHEMA
    controlled_use_ready: Literal[False] = False
    csv_emitted: Literal[False] = False

    def __post_init__(self):
        _Controlled.__post_init__(self)
        _sha(self.selection_sha256)
        _ids((self.policy_id,))
        _sha(self.policy_sha256)
        _sha(self.derivation_sha256)
        identifiers = tuple(item.hypothesis_id for item in self.values)
        _ids(identifiers)
        if (
            tuple(item.serialized_rank for item in self.values)
            != tuple(range(1, len(self.values) + 1))
            or any(left.source_rank >= right.source_rank for left, right in zip(self.values, self.values[1:]))
            or any(_exact(left.epcr) <= _exact(right.epcr) for left, right in zip(self.values, self.values[1:]))
        ):
            raise Track1EpcrError("invalid deterministic EPCR assignment")

    def as_bytes(self) -> bytes:
        return _bytes({
            "schema": self.schema,
            "classification": "controlled",
            "llm_policy": "no-llm",
            "confidence_label": self.confidence_label,
            "controlled_use_ready": False,
            "csv_emitted": False,
            "result": self,
        })


def policy_as_json(policy: EpcrPolicy) -> bytes:
    """Return canonical mechanics bytes; there is intentionally no default policy."""

    try:
        return _bytes({
            "schema": SCHEMA,
            "executable_version": EXECUTABLE_VERSION,
            "decimal_precision": DECIMAL_PRECISION,
            "tie_rule": TIE_RULE,
            "policy": _typed(policy, EpcrPolicy),
        })
    except _ERRORS:
        raise Track1EpcrError("cannot encode EPCR mechanics policy") from None


def policy_from_json(content: bytes, *, expected_sha256: str) -> EpcrPolicy:
    """Require independently supplied canonical bytes and their exact hash."""

    try:
        _sha(expected_sha256)
        if type(content) is not bytes or hashlib.sha256(content).hexdigest() != expected_sha256:
            raise ValueError
        data = _loads(content)
        if (
            type(data) is not dict
            or set(data) != {
                "schema", "executable_version", "decimal_precision", "tie_rule", "policy",
            }
            or data["schema"] != SCHEMA
            or data["executable_version"] != EXECUTABLE_VERSION
            or data["decimal_precision"] != DECIMAL_PRECISION
            or data["tie_rule"] != TIE_RULE
        ):
            raise ValueError
        policy = _typed(data["policy"], EpcrPolicy, decode=True)
        if policy_as_json(policy) != content:
            raise ValueError
        return policy
    except _ERRORS:
        raise Track1EpcrError("invalid EPCR mechanics policy") from None


def assign_epcr(
    selection: SubmissionSelection,
    *,
    policy_bytes: bytes,
    policy_sha256: str,
    derivation_bytes: bytes,
    expected_derivation_sha256: str,
    expected_result: EpcrAssignment | None = None,
) -> EpcrAssignment:
    """Assign deterministic relative values to the retained task 11.2 groups."""

    try:
        selection = _typed(selection, SubmissionSelection)
        return _assign_epcr_groups(
            hashlib.sha256(selection.as_bytes()).hexdigest(),
            selection.selected,
            policy_bytes=policy_bytes,
            policy_sha256=policy_sha256,
            derivation_bytes=derivation_bytes,
            expected_derivation_sha256=expected_derivation_sha256,
            expected_result=expected_result,
        )
    except _ERRORS:
        raise Track1EpcrError("source-bound EPCR assignment failed") from None


def _assign_epcr_groups(
    selection_sha256: str,
    groups,
    *,
    policy_bytes: bytes,
    policy_sha256: str,
    derivation_bytes: bytes,
    expected_derivation_sha256: str,
    expected_result: EpcrAssignment | None,
) -> EpcrAssignment:
    """Shared bounded kernel over the already limited selected group tuple."""

    try:
        _sha(selection_sha256)
        policy = policy_from_json(policy_bytes, expected_sha256=policy_sha256)
        _sha(expected_derivation_sha256)
        if (
            type(derivation_bytes) is not bytes
            or hashlib.sha256(derivation_bytes).hexdigest() != expected_derivation_sha256
            or policy.derivation_sha256 != expected_derivation_sha256
        ):
            raise ValueError
        groups = tuple(groups)
        scores = tuple(_exact(group.source_score) for group in groups)
        if (
            len(scores) > policy.maximum_rows
            or any(left < right for left, right in zip(scores, scores[1:]))
        ):
            raise ValueError

        values = []
        if scores:
            leader = scores[0]
            scale = _exact(policy.score_gap_scale)
            maximum_gap = _exact(policy.maximum_supported_score_gap)
            with localcontext() as context:
                context.prec = DECIMAL_PRECISION
                decay = _fraction_decimal(_exact(policy.rank_decay))
                minimum = _fraction_decimal(_exact(policy.minimum_serialized_epcr))
                for offset, (group, score) in enumerate(zip(groups, scores, strict=True)):
                    gap = leader - score
                    if gap > maximum_gap:
                        raise ValueError
                    base = (-_fraction_decimal(scale * gap)).exp()
                    epcr = base * decay**offset
                    if not base.is_finite() or not epcr.is_finite() or epcr < minimum:
                        raise ValueError
                    values.append(EpcrValue(
                        group.primary_hypothesis_id,
                        offset + 1,
                        group.source_rank,
                        group.source_score,
                        _fraction_text(gap),
                        _decimal_text(base),
                        _decimal_text(epcr),
                    ))
        result = EpcrAssignment(
            selection_sha256,
            policy.policy_id,
            policy_sha256,
            expected_derivation_sha256,
            policy.derivation_scope,
            TIE_RULE,
            DECIMAL_PRECISION,
            tuple(values),
        )
        if expected_result is not None and result.as_bytes() != _typed(expected_result, EpcrAssignment).as_bytes():
            raise ValueError
        return result
    except _ERRORS:
        raise Track1EpcrError("source-bound EPCR assignment failed") from None


def assign_sharded_epcr(
    selection,
    *,
    policy_bytes: bytes,
    policy_sha256: str,
    derivation_bytes: bytes,
    expected_derivation_sha256: str,
    expected_result: EpcrAssignment | None = None,
) -> EpcrAssignment:
    """Assign EPCR to an authenticated bounded selection manifest.

    Complete selection-audit authentication remains the preflight consumer's
    responsibility; this numerical step consumes only the at-most-ten selected
    groups already bound by the manifest identity.
    """

    try:
        from genosapiens.wgs_scientific_products import ScientificSelectionProduct

        if type(selection) is not ScientificSelectionProduct:
            raise ValueError
        return _assign_epcr_groups(
            selection.sha256,
            selection.selected,
            policy_bytes=policy_bytes,
            policy_sha256=policy_sha256,
            derivation_bytes=derivation_bytes,
            expected_derivation_sha256=expected_derivation_sha256,
            expected_result=expected_result,
        )
    except _ERRORS:
        raise Track1EpcrError("source-bound sharded EPCR assignment failed") from None
