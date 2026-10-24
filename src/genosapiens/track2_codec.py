"""Strict decoding of statically declared scientific records."""

from collections.abc import Mapping
from dataclasses import fields, is_dataclass
import math
from pathlib import Path
from types import UnionType
from typing import Any, Literal, Union, get_args, get_origin, get_type_hints


def json_value(value: Any) -> Any:
    """Copy JSON containers without serializing and duplicating large strings.

    Scientific records use string keys and finite numbers. Convert tuples to
    arrays, as the canonical JSON round trip does, and preserve immutable text.
    """
    if isinstance(value, dict):
        if any(not isinstance(key, str) for key in value):
            raise ValueError("scientific JSON keys must be strings")
        return {key: json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_value(item) for item in value]
    if value is None or type(value) in (str, bool, int):
        return value
    if type(value) is float and math.isfinite(value):
        return value
    raise ValueError("invalid scientific JSON value")


def decode(model: Any, value: Any) -> Any:
    origin, args = get_origin(model), get_args(model)
    if model in (Any, object):
        return value
    if origin in (Union, UnionType):
        for option in args:
            try:
                return decode(option, value)
            except (TypeError, ValueError):
                pass
        raise ValueError("invalid scientific record union")
    if origin is Literal:
        if value not in args or not any(type(value) is type(item) for item in args):
            raise ValueError("invalid scientific record literal")
        return value
    if origin is tuple:
        if not isinstance(value, list):
            raise ValueError("scientific sequence must be a JSON array")
        if len(args) == 2 and args[1] is Ellipsis:
            return tuple(decode(args[0], item) for item in value)
        if len(value) != len(args):
            raise ValueError("scientific sequence has invalid cardinality")
        return tuple(decode(kind, item) for kind, item in zip(args, value, strict=True))
    if origin in (dict, Mapping):
        if not isinstance(value, dict):
            raise ValueError("scientific mapping must be an object")
        return {decode(args[0], key): decode(args[1], item) for key, item in value.items()}
    if model is Path:
        if not isinstance(value, str):
            raise ValueError("path must be a string")
        return Path(value)
    if isinstance(model, type) and is_dataclass(model):
        if not isinstance(value, dict) or set(value) != {item.name for item in fields(model)}:
            raise ValueError("scientific record fields differ from the current contract")
        hints = get_type_hints(model)
        return model(**{name: decode(hints[name], item) for name, item in value.items()})
    if model is float and type(value) in (int, float) and math.isfinite(value):
        return float(value)
    if type(value) is not model:
        raise ValueError("scientific record value has invalid type")
    return value
