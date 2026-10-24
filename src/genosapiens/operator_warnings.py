"""Readable, narrowly scoped warnings for the local operator terminal."""

from contextlib import contextmanager
import sys
import warnings


class OperatorInputWarning(UserWarning):
    """An actionable input notice whose text contains no case values."""


class HpoBlankLineWarning(OperatorInputWarning):
    pass


@contextmanager
def friendly_operator_warnings(*, suppress_hpo=False):
    with warnings.catch_warnings():
        previous = warnings.showwarning

        def show(message, category, filename, lineno, file=None, line=None):
            if issubclass(category, OperatorInputWarning):
                print(f"Warning: {message}", file=file or sys.stderr, flush=True)
            else:
                previous(message, category, filename, lineno, file=file, line=line)

        warnings.showwarning = show
        if suppress_hpo:
            warnings.simplefilter("ignore", HpoBlankLineWarning)
        yield
