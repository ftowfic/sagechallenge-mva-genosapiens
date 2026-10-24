"""Explicit scientific execution choices and persistent operator budgets."""
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
import time

from genosapiens.track1_science_scope import ScienceScopeMode, ScienceScopePolicy
from genosapiens.wgs_ledger import canonical
from genosapiens.wgs_science_checkpoints import _path, _require, _sha
from genosapiens.wgs_science_runtime import ScienceBudget, ScienceClock


@dataclass(frozen=True, slots=True)
class ScienceRunOptions:
    mode: str = 'exhaustive'
    target_contexts: int = 1500
    budget_seconds: int = 18000
    report_reserve_seconds: int = 1800
    pilot: bool = False
    pilot_minutes: int = 20
    fresh_trial: bool = False
    no_budget: bool = False
    workers: int = 1
    publication_readers: int = 1

    def __post_init__(self):
        _require(type(self.workers) is int and self.workers in (1, 2, 4, 8)
            and type(self.publication_readers) is int and 1 <= self.publication_readers <= 8)
        ScienceScopePolicy(ScienceScopeMode(self.mode), self.target_contexts)
        _require(type(self.mode) is str and type(self.budget_seconds) is int
            and type(self.report_reserve_seconds) is int
            and 0 < self.report_reserve_seconds < self.budget_seconds <= 30 * 86400
            and type(self.pilot) is bool and type(self.pilot_minutes) is int
            and self.pilot_minutes in (20, 40) and (self.pilot or self.pilot_minutes == 20)
            and type(self.fresh_trial) is bool and not (self.fresh_trial and self.pilot)
            and type(self.no_budget) is bool and not (self.no_budget and (self.pilot or self.fresh_trial)))

    @property
    def policy(self):
        return ScienceScopePolicy(ScienceScopeMode(self.mode), self.target_contexts)

    def as_dict(self):
        return asdict(self)

    def bind_contracts(self, contracts):
        result = {key: dict(value) for key, value in contracts.items()}
        result['science']['scientific_scope'] = self.policy.sha256
        return result

    def open_budget(self, metadata_root, *, run_sha256, started_ns=None):
        root = _path(metadata_root)
        _require(root.is_dir() and _sha(run_sha256))
        now = time.time_ns()
        start = now if started_ns is None else started_ns
        _require(type(start) is int and 0 < start <= now)
        if self.no_budget:
            return ScienceClock(run_sha256, start)
        if self.fresh_trial:
            # Explicit operator authority: one new 30-minute timer for this
            # invocation. Historical cumulative records remain untouched.
            return ScienceBudget(run_sha256, start, start + 1800 * 10**9, 120 * 10**9)
        # Existing deadlines remain authoritative. Fresh creation uses the
        # invocation start, so intake and assembly consume the same budget.
        path = root / 'science-budget.json'
        overall = ScienceBudget.open(path, run_sha256=run_sha256, seconds=self.budget_seconds,
            reserve_seconds=self.report_reserve_seconds, now_ns=now if path.exists() else start)
        if not self.pilot:
            return overall
        path = root / 'science-pilot-budget.json'
        pilot = ScienceBudget.open(path, run_sha256=run_sha256, seconds=1200, reserve_seconds=120,
            now_ns=now if path.exists() else start)
        if self.pilot_minutes == 40:
            # Explicit extension retains the original 20-minute record. Its
            # timestamp anchors the longer pilot, including time between runs.
            from genosapiens.wgs_run_manifest import atomic_current
            original = pilot
            path = root / 'science-pilot-40-budget.json'
            pilot = ScienceBudget.open(path, run_sha256=run_sha256, seconds=2400, reserve_seconds=120,
                now_ns=now if path.exists() else original.started_ns)
            _require(pilot.started_ns == original.started_ns)
            receipt_path = root / 'science-pilot-40-authorization.json'
            receipt = dict(schema='genosapiens.explicit-pilot-extension/v1',
                reason='explicit_operator_40_minute_pilot', original_budget_sha256=original.sha256,
                extended_budget_sha256=pilot.sha256, overall_budget_sha256=overall.sha256,
                original_minutes=20, extended_minutes=40, overall_deadline_unchanged=True)
            if receipt_path.exists():
                from genosapiens.bounded_json import read_bounded_bytes
                _require(read_bounded_bytes(_path(receipt_path), maximum_bytes=4096) == canonical(receipt))
            else:
                atomic_current(receipt_path, receipt)
        final = min(overall.deadline_ns, pilot.deadline_ns)
        stop = min(overall.scheduling_deadline_ns, pilot.scheduling_deadline_ns)
        return ScienceBudget(run_sha256, min(overall.started_ns, pilot.started_ns), final, final - stop)


    def supervisor_seconds(self, metadata_root, *, elapsed_seconds):
        """Bound the process owner even if a child cannot deliver its alarm."""
        from genosapiens.bounded_json import read_bounded_bytes
        from genosapiens.wgs_science_runtime import SCHEMA
        root = _path(metadata_root)
        if self.no_budget:
            return None
        if self.fresh_trial:
            remaining = 1800 - elapsed_seconds
            _require(remaining > 0)
            return remaining
        limits = [self.budget_seconds - elapsed_seconds]
        paths = [('science-budget.json', self.budget_seconds, self.report_reserve_seconds)]
        if self.pilot:
            limits.append(self.pilot_minutes * 60 - elapsed_seconds)
            paths.append(('science-pilot-budget.json', 1200, 120))
        for name, seconds, reserve in paths:
            path = _path(root / name)
            if path.exists():
                raw = read_bounded_bytes(path, maximum_bytes=4096)
                value = json.loads(raw)
                _require(value.pop('schema') == SCHEMA)
                budget = ScienceBudget(**value)
                _require(budget.as_bytes() == raw
                    and budget.deadline_ns - budget.started_ns == seconds * 10**9
                    and budget.reserve_ns == reserve * 10**9)
                if name == 'science-pilot-budget.json' and self.pilot_minutes == 40:
                    # The child records the extension after obtaining case identity.
                    remaining = max(0, budget.deadline_ns + 1200 * 10**9 - time.time_ns()) / 10**9
                else:
                    remaining = budget.remaining(reporting=True)
                limits.append(remaining)
        remaining = min(limits)
        _require(remaining > 0)
        return remaining


def computation_contract(contracts):
    """Scope changes invalidate case aggregation, while direct science stays bound."""
    return {key: value for key, value in contracts['science'].items() if key != 'scientific_scope'}


def run_identity(*, input_sha256, hpo_sha256, sample):
    _require(_sha(input_sha256) and _sha(hpo_sha256) and type(sample) is str and 0 < len(sample) <= 1024)
    return hashlib.sha256(canonical(dict(input_sha256=input_sha256, hpo_sha256=hpo_sha256, sample=sample))).hexdigest()
