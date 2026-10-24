"""Small Draft 2020-12 JSON Schema evaluator for committed contracts."""

from __future__ import annotations

import re
from typing import Any, Mapping


class JsonSchemaError(ValueError):
    """Raised when a value violates a supported JSON Schema keyword."""


def _matches_type(value: Any, expected: str) -> bool:
    return {
        "array": lambda: isinstance(value, list),
        "boolean": lambda: isinstance(value, bool),
        "integer": lambda: isinstance(value, int) and not isinstance(value, bool),
        "null": lambda: value is None,
        "number": lambda: isinstance(value, (int, float)) and not isinstance(value, bool),
        "object": lambda: isinstance(value, Mapping),
        "string": lambda: isinstance(value, str),
    }.get(expected, lambda: False)()


def _resolve_ref(root: Mapping[str, Any], reference: str) -> Mapping[str, Any]:
    if not reference.startswith("#/"):
        raise JsonSchemaError(f"unsupported external schema reference: {reference}")
    value: Any = root
    for part in reference[2:].split("/"):
        key = part.replace("~1", "/").replace("~0", "~")
        if not isinstance(value, Mapping) or key not in value:
            raise JsonSchemaError(f"unresolved schema reference: {reference}")
        value = value[key]
    if not isinstance(value, Mapping):
        raise JsonSchemaError(f"schema reference is not an object: {reference}")
    return value


def validate_json_schema(
    value: Any,
    schema: Mapping[str, Any],
    *,
    root_schema: Mapping[str, Any] | None = None,
    path: str = "$",
) -> None:
    """Validate the JSON Schema subset used by project-owned contracts."""
    root = root_schema or schema
    if "anyOf" in schema:
        alternatives = schema["anyOf"]
        if not isinstance(alternatives, list) or not alternatives:
            raise JsonSchemaError(f"{path}: anyOf must be a non-empty list")
        matches = 0
        for alternative in alternatives:
            if not isinstance(alternative, Mapping):
                raise JsonSchemaError(f"{path}: anyOf entries must be objects")
            try:
                validate_json_schema(
                    value,
                    alternative,
                    root_schema=root,
                    path=path,
                )
            except JsonSchemaError:
                continue
            matches += 1
        if matches == 0:
            raise JsonSchemaError(f"{path}: value does not match any allowed schema")
    if "$ref" in schema:
        reference = schema["$ref"]
        if not isinstance(reference, str):
            raise JsonSchemaError(f"{path}: $ref must be a string")
        validate_json_schema(
            value,
            _resolve_ref(root, reference),
            root_schema=root,
            path=path,
        )
        return

    if "const" in schema and value != schema["const"]:
        raise JsonSchemaError(f"{path}: value must equal {schema['const']!r}")
    if "enum" in schema:
        choices = schema["enum"]
        if not isinstance(choices, list) or value not in choices:
            raise JsonSchemaError(f"{path}: value is outside the allowed enum")

    expected = schema.get("type")
    if expected is not None:
        expected_types = [expected] if isinstance(expected, str) else expected
        if (
            not isinstance(expected_types, list)
            or any(not isinstance(item, str) for item in expected_types)
            or not any(_matches_type(value, item) for item in expected_types)
        ):
            raise JsonSchemaError(f"{path}: value has the wrong JSON type")

    if isinstance(value, Mapping):
        required = schema.get("required", [])
        if not isinstance(required, list):
            raise JsonSchemaError(f"{path}: schema required must be a list")
        missing = [key for key in required if key not in value]
        if missing:
            raise JsonSchemaError(f"{path}: missing required fields: {', '.join(missing)}")
        properties = schema.get("properties", {})
        if not isinstance(properties, Mapping):
            raise JsonSchemaError(f"{path}: schema properties must be an object")
        if schema.get("additionalProperties") is False:
            unknown = sorted(set(value) - set(properties))
            if unknown:
                raise JsonSchemaError(f"{path}: unknown fields: {', '.join(unknown)}")
        for key, item in value.items():
            child_schema = properties.get(key)
            if child_schema is not None:
                if not isinstance(child_schema, Mapping):
                    raise JsonSchemaError(f"{path}.{key}: schema must be an object")
                validate_json_schema(
                    item,
                    child_schema,
                    root_schema=root,
                    path=f"{path}.{key}",
                )

    if isinstance(value, list):
        minimum = schema.get("minItems")
        if isinstance(minimum, int) and len(value) < minimum:
            raise JsonSchemaError(f"{path}: array has fewer than {minimum} items")
        if schema.get("uniqueItems") is True:
            representations = [repr(item) for item in value]
            if len(representations) != len(set(representations)):
                raise JsonSchemaError(f"{path}: array items must be unique")
        item_schema = schema.get("items")
        if item_schema is not None:
            if not isinstance(item_schema, Mapping):
                raise JsonSchemaError(f"{path}: items schema must be an object")
            for index, item in enumerate(value):
                validate_json_schema(
                    item,
                    item_schema,
                    root_schema=root,
                    path=f"{path}[{index}]",
                )

    if isinstance(value, str):
        minimum = schema.get("minLength")
        if isinstance(minimum, int) and len(value) < minimum:
            raise JsonSchemaError(f"{path}: string is shorter than {minimum}")
        pattern = schema.get("pattern")
        if isinstance(pattern, str) and re.search(pattern, value) is None:
            raise JsonSchemaError(f"{path}: string does not match {pattern!r}")

    if isinstance(value, (int, float)) and not isinstance(value, bool):
        minimum = schema.get("minimum")
        maximum = schema.get("maximum")
        if isinstance(minimum, (int, float)) and value < minimum:
            raise JsonSchemaError(f"{path}: number is below {minimum}")
        if isinstance(maximum, (int, float)) and value > maximum:
            raise JsonSchemaError(f"{path}: number is above {maximum}")
