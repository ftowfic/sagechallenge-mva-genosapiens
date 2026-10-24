"""Reference filesystem layout and measured capacity guardrails."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from decimal import Decimal, InvalidOperation
import json
import os
from pathlib import Path
import shutil
from typing import Any


DECIMAL_GB = 1_000_000_000
DEFAULT_RESERVE_BYTES = 100 * DECIMAL_GB

PUBLIC_LAYOUT = (
    "manifests",
    "releases",
    "releases/grch38",
    "releases/opentargets",
    "releases/gnomad",
    "releases/ensembl",
    "releases/clinvar",
    "releases/ontologies",
    "releases/giab",
    "derived",
    "derived/gnomad",
    "derived/opentargets",
    "derived/phenotype",
    "scratch",
    "scratch/public-build",
)
CONTROLLED_LAYOUT = (
    "controlled-work",
    "scratch/controlled",
)


class StorageError(RuntimeError):
    """Raised when storage initialization or capacity checks fail."""


class InsufficientCapacity(StorageError):
    """Raised when an operation would consume the emergency reserve."""


@dataclass(frozen=True, slots=True)
class CapacityPlan:
    operation: str
    source_bytes: int = 0
    output_bytes: int = 0
    scratch_bytes: int = 0
    reserve_bytes: int = DEFAULT_RESERVE_BYTES

    def __post_init__(self) -> None:
        if not self.operation.strip():
            raise StorageError("capacity plan requires an operation name")
        for field_name in (
            "source_bytes",
            "output_bytes",
            "scratch_bytes",
            "reserve_bytes",
        ):
            value = getattr(self, field_name)
            if not isinstance(value, int) or isinstance(value, bool) or value < 0:
                raise StorageError(f"{field_name} must be a non-negative integer")
        if self.reserve_bytes < DEFAULT_RESERVE_BYTES:
            raise StorageError("capacity plans must preserve at least 100 GB")

    @property
    def peak_bytes(self) -> int:
        return self.source_bytes + self.output_bytes + self.scratch_bytes


@dataclass(frozen=True, slots=True)
class CapacityReport:
    filesystem_root: str
    operation: str
    total_bytes: int
    used_bytes: int
    available_bytes: int
    source_bytes: int
    output_bytes: int
    scratch_bytes: int
    projected_free_bytes: int
    reserve_bytes: int
    allowed: bool
    additional_bytes_required: int


def gb_to_bytes(value: str | int | float | Decimal) -> int:
    try:
        amount = Decimal(str(value))
    except InvalidOperation as error:
        raise StorageError(f"invalid GB value: {value!r}") from error
    if not amount.is_finite() or amount < 0:
        raise StorageError("GB values must be finite and non-negative")
    result = amount * DECIMAL_GB
    if result != result.to_integral_value():
        raise StorageError("GB value must resolve to a whole number of bytes")
    return int(result)


def initialize_layout(root: str | Path) -> tuple[Path, ...]:
    reference_root = Path(root)
    if not reference_root.is_absolute():
        raise StorageError("reference root must be absolute")
    try:
        reference_root.mkdir(parents=True, exist_ok=True)
        if not reference_root.is_dir():
            raise StorageError(f"reference root is not a directory: {reference_root}")
        initialized: list[Path] = []
        for relative in PUBLIC_LAYOUT:
            path = reference_root / relative
            path.mkdir(parents=True, exist_ok=True, mode=0o755)
            if not path.is_dir():
                raise StorageError(f"layout path is not a directory: {path}")
            initialized.append(path)
        for relative in CONTROLLED_LAYOUT:
            path = reference_root / relative
            path.mkdir(parents=True, exist_ok=True, mode=0o700)
            if not path.is_dir():
                raise StorageError(f"layout path is not a directory: {path}")
            os.chmod(path, 0o700)
            initialized.append(path)
        return tuple(initialized)
    except OSError as error:
        raise StorageError(f"cannot initialize {reference_root}: {error}") from error


def capacity_report(root: str | Path, plan: CapacityPlan) -> CapacityReport:
    reference_root = Path(root)
    try:
        usage = shutil.disk_usage(reference_root)
    except OSError as error:
        raise StorageError(f"cannot measure {reference_root}: {error}") from error
    projected_free = usage.free - plan.peak_bytes
    additional_required = max(0, plan.reserve_bytes - projected_free)
    return CapacityReport(
        filesystem_root=str(reference_root.resolve()),
        operation=plan.operation,
        total_bytes=usage.total,
        used_bytes=usage.used,
        available_bytes=usage.free,
        source_bytes=plan.source_bytes,
        output_bytes=plan.output_bytes,
        scratch_bytes=plan.scratch_bytes,
        projected_free_bytes=projected_free,
        reserve_bytes=plan.reserve_bytes,
        allowed=additional_required == 0,
        additional_bytes_required=additional_required,
    )


def require_capacity(root: str | Path, plan: CapacityPlan) -> CapacityReport:
    report = capacity_report(root, plan)
    if not report.allowed:
        raise InsufficientCapacity(
            f"{plan.operation} requires {report.additional_bytes_required} additional bytes "
            f"to preserve the {report.reserve_bytes}-byte reserve"
        )
    return report


def report_as_json(report: CapacityReport) -> str:
    return json.dumps(asdict(report), indent=2, sort_keys=True)
