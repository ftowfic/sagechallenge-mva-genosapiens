"""Typed current configuration for reference and case workflows."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Any, Mapping, TypeVar


REFERENCE_SCHEMA = "genosapiens.reference-builder"
TRACK1_SCHEMA = "genosapiens.track1-pipeline"


class ConfigError(ValueError):
    """Raised when configuration does not satisfy its declared schema."""


def _required_string(data: Mapping[str, Any], key: str) -> str:
    value = data.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ConfigError(f"{key} must be a non-empty string")
    return value


def _absolute_path(data: Mapping[str, Any], key: str) -> Path:
    value = Path(_required_string(data, key))
    if not value.is_absolute():
        raise ConfigError(f"{key} must be an absolute container path")
    return value


def _reject_unknown(data: Mapping[str, Any], allowed: set[str]) -> None:
    unknown = sorted(set(data) - allowed)
    if unknown:
        raise ConfigError(f"unknown configuration fields: {', '.join(unknown)}")


@dataclass(frozen=True, slots=True)
class ReferenceBuilderConfig:
    schema: str
    reference_root: Path
    scratch_root: Path
    minimum_free_gb: int = 100

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "ReferenceBuilderConfig":
        _reject_unknown(
            data,
            {"schema", "reference_root", "scratch_root", "minimum_free_gb"},
        )
        schema = _required_string(data, "schema")
        if schema != REFERENCE_SCHEMA:
            raise ConfigError(f"schema must be {REFERENCE_SCHEMA!r}")
        minimum_free_gb = data.get("minimum_free_gb", 100)
        if not isinstance(minimum_free_gb, int) or isinstance(minimum_free_gb, bool):
            raise ConfigError("minimum_free_gb must be an integer")
        if minimum_free_gb < 100:
            raise ConfigError("minimum_free_gb must preserve at least 100 GB")
        reference_root = _absolute_path(data, "reference_root")
        scratch_root = _absolute_path(data, "scratch_root")
        if reference_root not in scratch_root.parents:
            raise ConfigError("scratch_root must be inside reference_root")
        return cls(schema, reference_root, scratch_root, minimum_free_gb)


@dataclass(frozen=True, slots=True)
class Track1PipelineConfig:
    schema: str
    reference_root: Path
    work_root: Path
    output_root: Path
    submission_contig_prefix: str = "chr"
    maximum_submission_rows: int = 10

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "Track1PipelineConfig":
        _reject_unknown(
            data,
            {
                "schema",
                "reference_root",
                "work_root",
                "output_root",
                "submission_contig_prefix",
                "maximum_submission_rows",
            },
        )
        schema = _required_string(data, "schema")
        if schema != TRACK1_SCHEMA:
            raise ConfigError(f"schema must be {TRACK1_SCHEMA!r}")
        prefix = data.get("submission_contig_prefix", "chr")
        if prefix not in {"", "chr"}:
            raise ConfigError("submission_contig_prefix must be empty or 'chr'")
        maximum_rows = data.get("maximum_submission_rows", 10)
        if not isinstance(maximum_rows, int) or isinstance(maximum_rows, bool):
            raise ConfigError("maximum_submission_rows must be an integer")
        if not 1 <= maximum_rows <= 10:
            raise ConfigError("maximum_submission_rows must be between 1 and 10")
        return cls(
            schema=schema,
            reference_root=_absolute_path(data, "reference_root"),
            work_root=_absolute_path(data, "work_root"),
            output_root=_absolute_path(data, "output_root"),
            submission_contig_prefix=prefix,
            maximum_submission_rows=maximum_rows,
        )




Config = ReferenceBuilderConfig | Track1PipelineConfig
T = TypeVar("T", bound=Config)


def load_config(path: str | Path, config_type: type[T]) -> T:
    config_path = Path(path)
    try:
        data = json.loads(config_path.read_text(encoding="utf-8"))
    except OSError as error:
        raise ConfigError(f"cannot read configuration {config_path}: {error}") from error
    except json.JSONDecodeError as error:
        raise ConfigError(f"invalid JSON in {config_path}: {error}") from error
    if not isinstance(data, Mapping):
        raise ConfigError("configuration root must be an object")
    return config_type.from_mapping(data)


def config_as_json(config: Config) -> str:
    def encode(value: Any) -> Any:
        if isinstance(value, Path):
            return str(value)
        raise TypeError(f"cannot encode {type(value).__name__}")

    return json.dumps(asdict(config), default=encode, indent=2, sort_keys=True)
