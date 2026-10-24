"""Bounded I/O helpers for untrusted remote responses."""

from __future__ import annotations

from queue import Empty, Queue
import threading
import time
from typing import Any, Callable, TypeVar, cast


Result = TypeVar("Result")


def call_before_deadline(operation: Callable[[], Result], *, deadline: float) -> Result:
    """Run blocking I/O without allowing it to exceed a caller's wall-clock budget."""

    result: Queue[tuple[Result | None, BaseException | None]] = Queue(maxsize=1)

    def call() -> None:
        try:
            result.put((operation(), None))
        except BaseException as error:
            result.put((None, error))

    remaining = deadline - time.monotonic()
    if remaining <= 0:
        raise TimeoutError("operation deadline exceeded")
    threading.Thread(target=call, daemon=True).start()
    try:
        value, error = result.get(timeout=remaining)
    except Empty as error:
        raise TimeoutError("operation deadline exceeded") from error
    if error is not None:
        raise error
    return cast(Result, value)


def read_response_before_deadline(
    response: Any,
    *,
    maximum_bytes: int,
    deadline: float,
) -> bytes:
    """Read at most ``maximum_bytes`` while enforcing a wall-clock deadline."""

    if maximum_bytes < 0:
        raise ValueError("maximum_bytes must be non-negative")

    def read_to_eof() -> bytes:
        content = bytearray()
        while len(content) < maximum_bytes:
            chunk = response.read(min(65_536, maximum_bytes - len(content)))
            if not isinstance(chunk, bytes):
                raise TypeError("response body is not bytes")
            if not chunk:
                break
            content.extend(chunk)
        return bytes(content)

    try:
        content = call_before_deadline(read_to_eof, deadline=deadline)
    except TimeoutError:
        try:
            response.close()
        except Exception:
            pass
        raise
    return content
