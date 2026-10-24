"""Small cumulative byte/deadline guards for authenticated public lookups."""

from contextlib import contextmanager
from contextvars import ContextVar
import sqlite3
import threading
import time

from genosapiens.bounded_wgs_config import PublicQueryLimits
from genosapiens.wgs_failures import WgsFailure


_CANCELLATION = ContextVar('wgs_query_cancellation',default=None)


@contextmanager
def _cancellable_queries(event):
    token=_CANCELLATION.set(event)
    try:
        yield
    finally:
        _CANCELLATION.reset(token)


def _query_cancelled():
    event=_CANCELLATION.get()
    return event is not None and event.is_set()


def _check_query_cancellation(origin='population_lookup',operation='population_query'):
    if _query_cancelled():
        raise WgsFailure(origin,'interrupted',operation=operation,exception_kind='interrupted')


class QueryBudget:
    """Count canonical record bytes (including newline) before decoding/appending.

    Each materialized representation has its own cumulative byte counter; a
    single deadline applies to the complete batch, including cache replay.
    Raw transport bytes have an independent counter with the same ceiling.
    """

    def __init__(self, limits: PublicQueryLimits, *, max_row_bytes: int, origin: str):
        if type(limits) is not PublicQueryLimits or type(max_row_bytes) is not int or max_row_bytes <= 0:
            raise ValueError("invalid public query limits")
        self.limits = limits
        self.max_row_bytes = min(max_row_bytes, limits.max_response_bytes)
        self.origin = origin
        self.deadline = time.monotonic() + limits.timeout_seconds
        self.bytes = {}

    def check(self):
        _check_query_cancellation(self.origin,self._operation())
        if time.monotonic() >= self.deadline:
            raise WgsFailure(
                self.origin, "timeout", operation=self._operation(),
                exception_kind="project_failure",
            )

    def charge(self, size: int, representation="records"):
        self.check()
        used = self.bytes.get(representation, 0)
        if type(size) is not int or size < 0 or size > self.max_row_bytes or used + size > self.limits.max_response_bytes:
            raise WgsFailure(
                self.origin, "response_capacity", operation=self._operation(),
                exception_kind="project_failure",
            )
        self.bytes[representation] = used + size

    def batch(self, size):
        self.check()
        if size > self.limits.batch_rows:
            raise WgsFailure(
                self.origin, "response_capacity", operation=self._operation(),
                exception_kind="project_failure",
            )

    def _operation(self):
        return {
            "population_lookup": "population_query",
            "clinical_index": "clinical_index_query",
            "clinical_replay": "clinical_response_bind",
            "clinical_adaptation": "clinical_response_adapt",
            "annotation_child": "annotation_query",
        }[self.origin]

    @contextmanager
    def at_origin(self, origin):
        """Share the batch accounting while identifying its active boundary."""
        previous = self.origin
        self.origin = origin
        try:
            yield
        finally:
            self.origin = previous

    @contextmanager
    def sqlite(self, connection):
        self.check()
        connection.set_progress_handler(lambda: int(time.monotonic() >= self.deadline or _query_cancelled()), 1000)
        try:
            yield
            self.check()
        except sqlite3.OperationalError:
            self.check()
            raise
        finally:
            connection.set_progress_handler(None, 0)

    @contextmanager
    def duckdb(self, connection):
        self.check()
        timer = threading.Timer(max(0, self.deadline - time.monotonic()), connection.interrupt)
        timer.daemon = True
        timer.start()
        try:
            yield
            self.check()
        except Exception:
            self.check()
            raise
        finally:
            timer.cancel()
            timer.join()
