"""Persistent scientific deadlines and bounded synchronous cancellation.

The caller holds the case output lock. The deadline belongs to the run, not an
attempt or scope choice. A timeout never supplies a scientific exclusion.
"""
from contextlib import contextmanager
from dataclasses import asdict, dataclass
import hashlib
import json
import math
import os
from pathlib import Path
import signal
import threading
import time
import uuid

from genosapiens.bounded_json import read_bounded_bytes
from genosapiens.wgs_ledger import canonical
from genosapiens.wgs_science_checkpoints import _path, _require, _sha, _sync_directory, _write

SCHEMA = 'genosapiens.science-runtime-budget/v1'
MAX_BUDGET_SECONDS = 30 * 24 * 3600
_ACTIVE_ALARMS = []


class ScienceDeadlineExpired(TimeoutError):
    """An operation reached its fixed deadline. Committed evidence remains valid."""


@dataclass(frozen=True, slots=True)
class ScienceBudget:
    run_sha256: str
    started_ns: int
    deadline_ns: int
    reserve_ns: int

    def __post_init__(self):
        _require(_sha(self.run_sha256)
            and all(type(v) is int for v in (self.started_ns, self.deadline_ns, self.reserve_ns))
            and 0 < self.started_ns < self.deadline_ns
            and 0 < self.reserve_ns < self.deadline_ns - self.started_ns
            and self.deadline_ns - self.started_ns <= MAX_BUDGET_SECONDS * 10**9)

    @property
    def scheduling_deadline_ns(self):
        return self.deadline_ns - self.reserve_ns

    def remaining(self, *, reporting=False, now_ns=None):
        now_ns = time.time_ns() if now_ns is None else now_ns
        _require(type(now_ns) is int and now_ns >= self.started_ns)
        deadline = self.deadline_ns if reporting else self.scheduling_deadline_ns
        return max(0, deadline - now_ns) / 10**9

    @property
    def sha256(self):
        return hashlib.sha256(self.as_bytes()).hexdigest()

    def as_bytes(self):
        return canonical(dict(schema=SCHEMA, **asdict(self)))

    @classmethod
    def open(cls, path, *, run_sha256, seconds, reserve_seconds, now_ns=None):
        """Create once or validate the same budget. Never replace a prior deadline."""
        _require(type(seconds) is int and type(reserve_seconds) is int
            and 0 < reserve_seconds < seconds <= MAX_BUDGET_SECONDS and _sha(run_sha256))
        path = _path(path)
        _require(path.parent.is_dir())
        now_ns = time.time_ns() if now_ns is None else now_ns
        if not path.exists():
            budget = cls(run_sha256, now_ns, now_ns + seconds * 10**9, reserve_seconds * 10**9)
            pending = path.with_name('.' + path.name + '-' + uuid.uuid4().hex)
            _write(pending, budget.as_bytes())
            try:
                os.link(pending, path, follow_symlinks=False)
                _sync_directory(path.parent)
            except FileExistsError:
                # A concurrent creator wins. Its deadline must satisfy the same
                # contract below, and its timestamp remains authoritative.
                pass
            finally:
                pending.unlink()
        raw = read_bounded_bytes(_path(path), maximum_bytes=4096)
        value = json.loads(raw)
        _require(type(value) is dict and set(value) == {'schema', 'run_sha256', 'started_ns', 'deadline_ns', 'reserve_ns'}
            and value.pop('schema') == SCHEMA)
        budget = cls(**value)
        _require(budget.as_bytes() == raw and budget.run_sha256 == run_sha256
            and budget.deadline_ns - budget.started_ns == seconds * 10**9
            and budget.reserve_ns == reserve_seconds * 10**9 and now_ns >= budget.started_ns)
        return budget


CLOCK_SCHEMA = 'genosapiens.science-elapsed-clock/v1'


@dataclass(frozen=True, slots=True)
class ScienceClock:
    """Invocation timing without a cumulative or absolute runtime deadline."""
    run_sha256: str
    started_ns: int

    def __post_init__(self):
        _require(_sha(self.run_sha256) and type(self.started_ns) is int and self.started_ns > 0)

    def remaining(self, *, reporting=False, now_ns=None):
        return math.inf

    @property
    def sha256(self):
        return hashlib.sha256(self.as_bytes()).hexdigest()

    def as_bytes(self):
        return canonical(dict(schema=CLOCK_SCHEMA, **asdict(self)))


@contextmanager
def bounded_operation(budget, *, maximum_seconds, reporting=False):
    """Cancel main-thread Python work with a real-time alarm and restore it.

    Worker processes and external commands still need their own bounded process
    owner. This scope neither creates nor kills unrelated processes.
    """
    _require(type(budget) in (ScienceBudget, ScienceClock) and type(maximum_seconds) in (int, float)
        and (math.isfinite(maximum_seconds) or type(budget) is ScienceClock and maximum_seconds == math.inf)
        and maximum_seconds > 0
        and threading.current_thread() is threading.main_thread())
    remaining = min(maximum_seconds, budget.remaining(reporting=reporting))
    if remaining == math.inf:
        yield
        return
    if remaining <= 0:
        raise ScienceDeadlineExpired('scientific runtime deadline expired')
    previous = signal.getsignal(signal.SIGALRM)
    timer = signal.getitimer(signal.ITIMER_REAL)
    # Do not steal an outer owner's cancellation mechanism.
    _require(timer == (0.0, 0.0) or bool(_ACTIVE_ALARMS) and previous is _ACTIVE_ALARMS[-1])
    outer_stop = None if timer[0] == 0 else time.monotonic() + timer[0]
    if outer_stop is not None:
        remaining = min(remaining, timer[0])
    def expired(_signum, _frame):
        raise ScienceDeadlineExpired('scientific operation deadline expired')
    signal.signal(signal.SIGALRM, expired)
    _ACTIVE_ALARMS.append(expired)
    try:
        signal.setitimer(signal.ITIMER_REAL, remaining)
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        _ACTIVE_ALARMS.pop()
        signal.signal(signal.SIGALRM, previous)
        if outer_stop is not None:
            signal.setitimer(signal.ITIMER_REAL, max(0.000001, outer_stop - time.monotonic()), timer[1])
