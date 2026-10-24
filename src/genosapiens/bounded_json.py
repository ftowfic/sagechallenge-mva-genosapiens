"""Bounded, duplicate-key-rejecting JSON input helpers."""

from __future__ import annotations

import json
import os
from pathlib import Path
import stat
from typing import Any


MAX_GOVERNING_JSON_BYTES = 2 * 1024 * 1024


class BoundedJSONError(ValueError):
    """Raised when a governing JSON document cannot be read safely."""


def read_bounded_bytes(
    path: str | Path,
    *,
    maximum_bytes: int = MAX_GOVERNING_JSON_BYTES,
    reference_root: str | Path | None = None,
    expected_size_bytes: int | None = None,
) -> bytes:
    source = Path(path)
    directory_descriptor = -1
    descriptor = -1
    try:
        if (
            not isinstance(maximum_bytes, int)
            or isinstance(maximum_bytes, bool)
            or maximum_bytes < 0
            or expected_size_bytes is not None
            and (
                not isinstance(expected_size_bytes, int)
                or isinstance(expected_size_bytes, bool)
                or expected_size_bytes < 0
            )
        ):
            raise BoundedJSONError("JSON input read limit is invalid")
        file_flags = os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK
        if reference_root is None:
            descriptor = os.open(source, file_flags)
        else:
            root = Path(reference_root)
            if (
                not source.is_absolute()
                or not root.is_absolute()
                or ".." in source.parts
                or ".." in root.parts
            ):
                raise BoundedJSONError("JSON input cannot be read")
            relative = source.relative_to(root)
            if not relative.parts:
                raise BoundedJSONError("JSON input cannot be read")
            directory_flags = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW
            directory_descriptor = os.open(root.anchor, directory_flags)
            for component in (*root.parts[1:], *relative.parts[:-1]):
                child = os.open(
                    component,
                    directory_flags,
                    dir_fd=directory_descriptor,
                )
                os.close(directory_descriptor)
                directory_descriptor = child
            descriptor = os.open(
                relative.name,
                file_flags,
                dir_fd=directory_descriptor,
            )
        metadata = os.fstat(descriptor)
        if (
            not stat.S_ISREG(metadata.st_mode)
            or metadata.st_size > maximum_bytes
            or expected_size_bytes is not None
            and metadata.st_size != expected_size_bytes
        ):
            raise BoundedJSONError("JSON input exceeds its safe read limit")
        with os.fdopen(descriptor, "rb") as handle:
            descriptor = -1
            content = handle.read(maximum_bytes + 1)
    except BoundedJSONError:
        raise
    except (AttributeError, OSError, ValueError) as error:
        raise BoundedJSONError("JSON input cannot be read") from error
    finally:
        if descriptor >= 0:
            os.close(descriptor)
        if directory_descriptor >= 0:
            os.close(directory_descriptor)
    if len(content) > maximum_bytes or len(content) != metadata.st_size:
        raise BoundedJSONError("JSON input exceeds its safe read limit")
    return content


def load_bounded_json(
    path: str | Path,
    *,
    maximum_bytes: int = MAX_GOVERNING_JSON_BYTES,
    reference_root: str | Path | None = None,
    expected_size_bytes: int | None = None,
) -> tuple[Any, bytes]:
    content = read_bounded_bytes(
        path,
        maximum_bytes=maximum_bytes,
        reference_root=reference_root,
        expected_size_bytes=expected_size_bytes,
    )

    def reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        value: dict[str, Any] = {}
        for key, child in pairs:
            if key in value:
                raise BoundedJSONError("JSON input contains duplicate fields")
            value[key] = child
        return value

    try:
        return json.loads(content, object_pairs_hook=reject_duplicate_keys), content
    except (UnicodeDecodeError, json.JSONDecodeError, RecursionError) as error:
        raise BoundedJSONError("JSON input is invalid") from error
