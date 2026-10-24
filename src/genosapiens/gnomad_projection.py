"""Schema-driven projection of the gnomAD v4.1 joint frequency table."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
import json
from pathlib import Path
import re
from typing import Any


SCHEMA = "genosapiens.gnomad-joint-projection"
DATASET_NAMES = ("exomes", "genomes", "joint")
DATASET_VCF_PREFIXES = {"exomes": "EX", "genomes": "GN", "joint": "JT"}
EXPECTED_REGION_FLAGS = (
    "fail_interval_qc",
    "outside_broad_capture_region",
    "outside_ukb_capture_region",
    "outside_broad_calling_region",
    "outside_ukb_calling_region",
    "not_called_in_exomes",
    "not_called_in_genomes",
)
IDENTIFIER = re.compile(r"^[a-z0-9][a-z0-9._-]*$")
ANCESTRY = re.compile(r"^[a-z][a-z0-9_]*$")


class GnomadProjectionError(RuntimeError):
    """Raised when the projection contract or source schema is invalid."""


def _mapping(value: Any, context: str) -> Mapping[str, Any]:
    if isinstance(value, Mapping):
        return value
    try:
        converted = dict(value)
    except (TypeError, ValueError) as error:
        raise GnomadProjectionError(f"{context} must be an object") from error
    if not all(isinstance(key, str) for key in converted):
        raise GnomadProjectionError(f"{context} keys must be strings")
    return converted


def _exact_fields(data: Mapping[str, Any], expected: set[str], context: str) -> None:
    if set(data) != expected:
        raise GnomadProjectionError(f"{context} fields are invalid")


def _string(data: Mapping[str, Any], name: str, context: str) -> str:
    value = data.get(name)
    if not isinstance(value, str) or not value:
        raise GnomadProjectionError(f"{context}.{name} must be a non-empty string")
    return value


def _positive_int(data: Mapping[str, Any], name: str, context: str) -> int:
    value = data.get(name)
    if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
        raise GnomadProjectionError(f"{context}.{name} must be a positive integer")
    return value


def _strings(value: Any, context: str) -> tuple[str, ...]:
    if (
        not isinstance(value, list)
        or not value
        or any(not isinstance(item, str) or not item for item in value)
        or len(value) != len(set(value))
    ):
        raise GnomadProjectionError(f"{context} must be a unique string array")
    return tuple(value)


@dataclass(frozen=True, slots=True)
class GnomadSourceContract:
    resource_id: str
    release: str
    assembly: str
    hail_path: str
    root_metadata_url: str
    root_metadata_object_generation: str
    root_metadata_md5_base64: str
    written_hail_version: str
    expected_rows: int
    expected_partitions: int
    license: str
    license_url: str

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any]) -> "GnomadSourceContract":
        expected = {
            "resource_id",
            "release",
            "assembly",
            "hail_path",
            "root_metadata_url",
            "root_metadata_object_generation",
            "root_metadata_md5_base64",
            "written_hail_version",
            "expected_rows",
            "expected_partitions",
            "license",
            "license_url",
        }
        _exact_fields(value, expected, "source")
        source = cls(
            resource_id=_string(value, "resource_id", "source"),
            release=_string(value, "release", "source"),
            assembly=_string(value, "assembly", "source"),
            hail_path=_string(value, "hail_path", "source"),
            root_metadata_url=_string(value, "root_metadata_url", "source"),
            root_metadata_object_generation=_string(
                value, "root_metadata_object_generation", "source"
            ),
            root_metadata_md5_base64=_string(
                value, "root_metadata_md5_base64", "source"
            ),
            written_hail_version=_string(value, "written_hail_version", "source"),
            expected_rows=_positive_int(value, "expected_rows", "source"),
            expected_partitions=_positive_int(
                value, "expected_partitions", "source"
            ),
            license=_string(value, "license", "source"),
            license_url=_string(value, "license_url", "source"),
        )
        if (
            not IDENTIFIER.fullmatch(source.resource_id)
            or not source.hail_path.startswith("gs://")
            or not source.root_metadata_url.startswith("https://")
            or not source.license_url.startswith("https://")
        ):
            raise GnomadProjectionError("source identifiers or URLs are invalid")
        return source


@dataclass(frozen=True, slots=True)
class DatasetContract:
    name: str
    row_field: str
    globals_field: str
    expected_frequency_entries: int
    expected_faf_entries: int
    frequency_ancestries: tuple[str, ...]
    faf_ancestries: tuple[str, ...]
    has_filters: bool

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any]) -> "DatasetContract":
        expected = {
            "name",
            "row_field",
            "globals_field",
            "expected_frequency_entries",
            "expected_faf_entries",
            "frequency_ancestries",
            "faf_ancestries",
            "has_filters",
        }
        _exact_fields(value, expected, "dataset")
        has_filters = value.get("has_filters")
        if not isinstance(has_filters, bool):
            raise GnomadProjectionError("dataset.has_filters must be a Boolean")
        dataset = cls(
            name=_string(value, "name", "dataset"),
            row_field=_string(value, "row_field", "dataset"),
            globals_field=_string(value, "globals_field", "dataset"),
            expected_frequency_entries=_positive_int(
                value, "expected_frequency_entries", "dataset"
            ),
            expected_faf_entries=_positive_int(
                value, "expected_faf_entries", "dataset"
            ),
            frequency_ancestries=_strings(
                value.get("frequency_ancestries"), "dataset.frequency_ancestries"
            ),
            faf_ancestries=_strings(
                value.get("faf_ancestries"), "dataset.faf_ancestries"
            ),
            has_filters=has_filters,
        )
        if dataset.name not in DATASET_NAMES:
            raise GnomadProjectionError(f"unsupported gnomAD dataset: {dataset.name}")
        if dataset.row_field != dataset.name:
            raise GnomadProjectionError("dataset row field must match its name")
        if dataset.globals_field != f"{dataset.name}_globals":
            raise GnomadProjectionError("dataset globals field is invalid")
        if any(not ANCESTRY.fullmatch(item) for item in dataset.frequency_ancestries):
            raise GnomadProjectionError("frequency ancestry identifier is invalid")
        if not set(dataset.faf_ancestries).issubset(dataset.frequency_ancestries):
            raise GnomadProjectionError("FAF ancestries must have frequency entries")
        return dataset


@dataclass(frozen=True, slots=True)
class GnomadProjectionContract:
    projection_id: str
    source: GnomadSourceContract
    major_genetic_ancestries: tuple[str, ...]
    region_flags: tuple[str, ...]
    datasets: tuple[DatasetContract, ...]
    contingency_index_dataset: str
    discrepant_p_value_threshold: float

    def dataset(self, name: str) -> DatasetContract:
        for dataset in self.datasets:
            if dataset.name == name:
                return dataset
        raise GnomadProjectionError(f"projection dataset not found: {name}")


def load_gnomad_projection_contract(
    path: str | Path,
) -> GnomadProjectionContract:
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise GnomadProjectionError(f"cannot load gnomAD projection contract: {error}") from error
    data = _mapping(data, "projection contract")
    expected = {
        "schema",
        "projection_id",
        "source",
        "major_genetic_ancestries",
        "region_flags",
        "datasets",
        "discordance",
    }
    _exact_fields(data, expected, "projection contract")
    if data.get("schema") != SCHEMA:
        raise GnomadProjectionError(f"projection schema must be {SCHEMA!r}")
    source = GnomadSourceContract.from_mapping(_mapping(data.get("source"), "source"))
    major = _strings(
        data.get("major_genetic_ancestries"), "major_genetic_ancestries"
    )
    if any(not ANCESTRY.fullmatch(item) for item in major):
        raise GnomadProjectionError("major ancestry identifier is invalid")
    region_flags = _strings(data.get("region_flags"), "region_flags")
    if region_flags != EXPECTED_REGION_FLAGS:
        raise GnomadProjectionError("region_flags do not match the v4.1 source schema")
    dataset_values = data.get("datasets")
    if not isinstance(dataset_values, list):
        raise GnomadProjectionError("datasets must be an array")
    datasets = tuple(
        DatasetContract.from_mapping(_mapping(item, "dataset"))
        for item in dataset_values
    )
    if tuple(dataset.name for dataset in datasets) != DATASET_NAMES:
        raise GnomadProjectionError("datasets must be exomes, genomes, and joint")
    if any(not set(major).issubset(dataset.frequency_ancestries) for dataset in datasets):
        raise GnomadProjectionError("all major ancestries require frequency entries")
    discordance = _mapping(data.get("discordance"), "discordance")
    _exact_fields(
        discordance,
        {"contingency_index_dataset", "discrepant_p_value_threshold"},
        "discordance",
    )
    index_dataset = _string(
        discordance, "contingency_index_dataset", "discordance"
    )
    threshold = discordance.get("discrepant_p_value_threshold")
    if index_dataset != "joint":
        raise GnomadProjectionError(
            "contingency statistics must use joint frequency indices"
        )
    if (
        not isinstance(threshold, (int, float))
        or isinstance(threshold, bool)
        or not 0 < float(threshold) < 1
    ):
        raise GnomadProjectionError("discordance threshold must be between zero and one")
    return GnomadProjectionContract(
        projection_id=_string(data, "projection_id", "projection contract"),
        source=source,
        major_genetic_ancestries=major,
        region_flags=region_flags,
        datasets=datasets,
        contingency_index_dataset=index_dataset,
        discrepant_p_value_threshold=float(threshold),
    )


@dataclass(frozen=True, slots=True)
class AncestryProjectionIndices:
    genetic_ancestry: str
    adjusted_frequency: int
    xy_adjusted_frequency: int
    adjusted_faf: int | None


@dataclass(frozen=True, slots=True)
class DatasetProjectionIndices:
    name: str
    row_field: str
    adjusted_frequency: int
    xy_adjusted_frequency: int
    adjusted_faf: int
    ancestries: tuple[AncestryProjectionIndices, ...]
    has_filters: bool


@dataclass(frozen=True, slots=True)
class GnomadProjectionLayout:
    datasets: tuple[DatasetProjectionIndices, ...]
    region_flags: tuple[str, ...]
    discrepant_p_value_threshold: float

    def dataset(self, name: str) -> DatasetProjectionIndices:
        for dataset in self.datasets:
            if dataset.name == name:
                return dataset
        raise GnomadProjectionError(f"resolved projection dataset not found: {name}")


def _metadata_array(value: Any, expected_size: int, context: str) -> tuple[Mapping[str, str], ...]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise GnomadProjectionError(f"{context} must be an array")
    if len(value) != expected_size:
        raise GnomadProjectionError(
            f"{context} has {len(value)} entries; expected {expected_size}"
        )
    result = []
    for offset, item in enumerate(value):
        mapping = _mapping(item, f"{context}[{offset}]")
        if not mapping or any(
            not isinstance(key, str)
            or not key
            or not isinstance(entry, str)
            or not entry
            for key, entry in mapping.items()
        ):
            raise GnomadProjectionError(f"{context}[{offset}] is invalid")
        result.append(dict(mapping))
    return tuple(result)


def _metadata_key(value: Mapping[str, str], context: str) -> str:
    if not set(value).issubset({"gen_anc", "sex", "group"}):
        raise GnomadProjectionError(f"{context} contains unsupported dimensions")
    return "_".join(
        item
        for item in (value.get("gen_anc"), value.get("sex"), value.get("group"))
        if item is not None
    )


def _index_dictionary(
    value: Any,
    metadata: tuple[Mapping[str, str], ...],
    context: str,
) -> Mapping[str, int]:
    mapping = _mapping(value, context)
    size = len(metadata)
    if len(mapping) != size:
        raise GnomadProjectionError(
            f"{context} has {len(mapping)} entries; expected {size}"
        )
    if any(
        not isinstance(key, str)
        or not key
        or not isinstance(index, int)
        or isinstance(index, bool)
        or index < 0
        or index >= size
        for key, index in mapping.items()
    ):
        raise GnomadProjectionError(f"{context} contains an invalid index")
    if set(mapping.values()) != set(range(size)):
        raise GnomadProjectionError(f"{context} must map each array index exactly once")
    if any(
        key != _metadata_key(metadata[index], f"{context}[{key!r}]")
        for key, index in mapping.items()
    ):
        raise GnomadProjectionError(f"{context} key-to-metadata mapping is invalid")
    return dict(mapping)


def _find_index(
    metadata: tuple[Mapping[str, str], ...],
    index_dictionary: Mapping[str, int],
    target: Mapping[str, str],
    context: str,
) -> int:
    key = _metadata_key(target, context)
    index = index_dictionary.get(key)
    if index is None:
        raise GnomadProjectionError(
            f"{context} source-global index is missing for {key!r}"
        )
    if metadata[index] != target:
        raise GnomadProjectionError(
            f"{context} source-global index {key!r} does not match its metadata"
        )
    return index


def _ancestries(metadata: tuple[Mapping[str, str], ...]) -> tuple[str, ...]:
    return tuple(
        sorted(
            item["gen_anc"]
            for item in metadata
            if set(item) == {"group", "gen_anc"} and item["group"] == "adj"
        )
    )


def resolve_gnomad_projection_layout(
    contract: GnomadProjectionContract,
    source_globals: Mapping[str, Any],
) -> GnomadProjectionLayout:
    source_globals = _mapping(source_globals, "source globals")
    resolved = []
    for dataset in contract.datasets:
        values = _mapping(
            source_globals.get(dataset.globals_field), dataset.globals_field
        )
        freq_meta = _metadata_array(
            values.get("freq_meta"),
            dataset.expected_frequency_entries,
            f"{dataset.globals_field}.freq_meta",
        )
        freq_index = _index_dictionary(
            values.get("freq_index_dict"),
            freq_meta,
            f"{dataset.globals_field}.freq_index_dict",
        )
        sample_counts = values.get("freq_meta_sample_count")
        if (
            not isinstance(sample_counts, Sequence)
            or isinstance(sample_counts, (str, bytes))
            or len(sample_counts) != dataset.expected_frequency_entries
        ):
            raise GnomadProjectionError(
                f"{dataset.globals_field}.freq_meta_sample_count has invalid length"
            )
        faf_meta = _metadata_array(
            values.get("faf_meta"),
            dataset.expected_faf_entries,
            f"{dataset.globals_field}.faf_meta",
        )
        faf_index = _index_dictionary(
            values.get("faf_index_dict"),
            faf_meta,
            f"{dataset.globals_field}.faf_index_dict",
        )
        if _ancestries(freq_meta) != tuple(sorted(dataset.frequency_ancestries)):
            raise GnomadProjectionError(
                f"{dataset.globals_field}.freq_meta ancestry set changed"
            )
        if _ancestries(faf_meta) != tuple(sorted(dataset.faf_ancestries)):
            raise GnomadProjectionError(
                f"{dataset.globals_field}.faf_meta ancestry set changed"
            )
        _find_index(freq_meta, freq_index, {"group": "raw"}, dataset.name)
        _find_index(
            freq_meta,
            freq_index,
            {"sex": "XX", "group": "adj"},
            dataset.name,
        )
        adjusted = _find_index(
            freq_meta, freq_index, {"group": "adj"}, dataset.name
        )
        xy_adjusted = _find_index(
            freq_meta,
            freq_index,
            {"sex": "XY", "group": "adj"},
            dataset.name,
        )
        adjusted_faf = _find_index(
            faf_meta, faf_index, {"group": "adj"}, dataset.name
        )
        ancestry_indices = []
        for ancestry in contract.major_genetic_ancestries:
            ancestry_frequency = _find_index(
                freq_meta,
                freq_index,
                {"gen_anc": ancestry, "group": "adj"},
                dataset.name,
            )
            ancestry_xy = _find_index(
                freq_meta,
                freq_index,
                {"gen_anc": ancestry, "sex": "XY", "group": "adj"},
                dataset.name,
            )
            ancestry_faf = None
            if ancestry in dataset.faf_ancestries:
                ancestry_faf = _find_index(
                    faf_meta,
                    faf_index,
                    {"gen_anc": ancestry, "group": "adj"},
                    dataset.name,
                )
            ancestry_indices.append(
                AncestryProjectionIndices(
                    genetic_ancestry=ancestry,
                    adjusted_frequency=ancestry_frequency,
                    xy_adjusted_frequency=ancestry_xy,
                    adjusted_faf=ancestry_faf,
                )
            )
        resolved.append(
            DatasetProjectionIndices(
                name=dataset.name,
                row_field=dataset.row_field,
                adjusted_frequency=adjusted,
                xy_adjusted_frequency=xy_adjusted,
                adjusted_faf=adjusted_faf,
                ancestries=tuple(ancestry_indices),
                has_filters=dataset.has_filters,
            )
        )
    return GnomadProjectionLayout(
        datasets=tuple(resolved),
        region_flags=contract.region_flags,
        discrepant_p_value_threshold=contract.discrepant_p_value_threshold,
    )


def _array_item(value: Any, index: int, context: str) -> Any:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise GnomadProjectionError(f"{context} must be an array")
    if index >= len(value):
        raise GnomadProjectionError(f"{context} is shorter than its global metadata")
    return value[index]


def _nullable_fields(value: Any, fields: tuple[str, ...], context: str) -> dict[str, Any]:
    if value is None:
        return {field: None for field in fields}
    value = _mapping(value, context)
    return {field: value.get(field) for field in fields}


def _count(
    source: Mapping[str, Any] | None,
    adjusted_index: int,
    xy_index: int,
    *,
    is_autosome_or_par: bool,
    context: str,
) -> dict[str, Any]:
    if source is None:
        adjusted = None
        xy_adjusted = None
    else:
        frequencies = source.get("freq")
        if frequencies is None:
            adjusted = None
            xy_adjusted = None
        else:
            adjusted = _array_item(frequencies, adjusted_index, f"{context}.freq")
            xy_adjusted = _array_item(frequencies, xy_index, f"{context}.freq")
    count = _nullable_fields(
        adjusted, ("AC", "AN", "homozygote_count"), f"{context}.adjusted"
    )
    if is_autosome_or_par:
        hemizygote_count = (
            0 if any(value is not None for value in count.values()) else None
        )
    else:
        hemizygote_count = _nullable_fields(
            xy_adjusted, ("AC",), f"{context}.xy_adjusted"
        )["AC"]
    return {
        "ac": count["AC"],
        "an": count["AN"],
        "homozygote_count": count["homozygote_count"],
        "hemizygote_count": hemizygote_count,
    }


def _faf(source: Mapping[str, Any] | None, index: int | None, context: str) -> dict[str, Any]:
    value = None
    if source is not None and index is not None and source.get("faf") is not None:
        value = _array_item(source.get("faf"), index, f"{context}.faf")
    return _nullable_fields(value, ("faf95", "faf99"), context)


def _project_dataset(
    row: Mapping[str, Any],
    indices: DatasetProjectionIndices,
    *,
    is_autosome_or_par: bool,
) -> dict[str, Any]:
    raw_source = row.get(indices.row_field)
    source = None if raw_source is None else _mapping(raw_source, indices.row_field)
    result = {
        "adjusted": _count(
            source,
            indices.adjusted_frequency,
            indices.xy_adjusted_frequency,
            is_autosome_or_par=is_autosome_or_par,
            context=indices.name,
        ),
        "ancestries": {
            ancestry.genetic_ancestry: {
                **_count(
                    source,
                    ancestry.adjusted_frequency,
                    ancestry.xy_adjusted_frequency,
                    is_autosome_or_par=is_autosome_or_par,
                    context=f"{indices.name}.{ancestry.genetic_ancestry}",
                ),
                **_faf(
                    source,
                    ancestry.adjusted_faf,
                    f"{indices.name}.{ancestry.genetic_ancestry}",
                ),
            }
            for ancestry in indices.ancestries
        },
        "faf": _faf(source, indices.adjusted_faf, f"{indices.name}.adjusted"),
        "grpmax": _nullable_fields(
            None if source is None else source.get("grpmax"),
            ("AC", "AF", "AN", "homozygote_count", "gen_anc"),
            f"{indices.name}.grpmax",
        ),
        "fafmax": _nullable_fields(
            None if source is None else source.get("fafmax"),
            (
                "faf95_max",
                "faf95_max_gen_anc",
                "faf99_max",
                "faf99_max_gen_anc",
            ),
            f"{indices.name}.fafmax",
        ),
    }
    if indices.has_filters:
        filters = None if source is None else source.get("filters")
        if filters is not None:
            if not isinstance(filters, (set, frozenset, list, tuple)) or any(
                not isinstance(item, str) for item in filters
            ):
                raise GnomadProjectionError(f"{indices.name}.filters is invalid")
            filters = sorted(filters)
        result["filters"] = filters
    return result


def _test_stat(value: Any, context: str) -> dict[str, Any]:
    return _nullable_fields(value, ("p_value", "odds_ratio"), context)


def project_gnomad_joint_row(
    row: Mapping[str, Any],
    layout: GnomadProjectionLayout,
    *,
    is_autosome_or_par: bool,
) -> dict[str, Any]:
    if not isinstance(is_autosome_or_par, bool):
        raise GnomadProjectionError("is_autosome_or_par must be a Boolean")
    row = _mapping(row, "source row")
    joint = layout.dataset("joint")
    stats_value = row.get("freq_comparison_stats")
    stats = None if stats_value is None else _mapping(stats_value, "freq_comparison_stats")
    contingency = None if stats is None else stats.get("contingency_table_test")
    overall_test = None
    ancestry_tests: dict[str, dict[str, Any]] = {}
    if contingency is not None:
        overall_test = _array_item(
            contingency, joint.adjusted_frequency, "contingency_table_test"
        )
    for ancestry in joint.ancestries:
        test = None
        if contingency is not None:
            test = _array_item(
                contingency,
                ancestry.adjusted_frequency,
                "contingency_table_test",
            )
        ancestry_tests[ancestry.genetic_ancestry] = _test_stat(
            test, f"discordance.{ancestry.genetic_ancestry}"
        )
    stat_union = _nullable_fields(
        None if stats is None else stats.get("stat_union"),
        ("p_value", "stat_test_name", "gen_ancs"),
        "discordance.stat_union",
    )
    p_value = stat_union["p_value"]
    discordant = None if p_value is None else p_value < layout.discrepant_p_value_threshold
    region_source_value = row.get("region_flags")
    region_source = (
        {} if region_source_value is None else _mapping(region_source_value, "region_flags")
    )
    return {
        "locus": row.get("locus"),
        "alleles": row.get("alleles"),
        **{
            dataset.name: _project_dataset(
                row, dataset, is_autosome_or_par=is_autosome_or_par
            )
            for dataset in layout.datasets
        },
        "region_flags": {
            flag: region_source.get(flag) for flag in layout.region_flags
        },
        "discordance": {
            "overall": _test_stat(overall_test, "discordance.overall"),
            "ancestries": ancestry_tests,
            "cochran_mantel_haenszel_test": _nullable_fields(
                None if stats is None else stats.get("cochran_mantel_haenszel_test"),
                ("p_value", "chisq"),
                "discordance.cochran_mantel_haenszel_test",
            ),
            "stat_union": stat_union,
            "is_discrepant": discordant,
        },
    }


def project_hail_joint_table(
    table: Any,
    hail_module: Any,
    contract: GnomadProjectionContract,
) -> Any:
    """Build the scalable Hail projection after evaluating source globals once."""
    globals_value = hail_module.eval(table.index_globals())
    layout = resolve_gnomad_projection_layout(contract, globals_value)
    return project_hail_joint_table_with_layout(table, hail_module, layout)


def project_hail_joint_table_with_layout(
    table: Any,
    hail_module: Any,
    layout: GnomadProjectionLayout,
) -> Any:
    """Build the Hail projection with a previously validated global layout."""
    in_autosome_or_par = table.locus.in_autosome_or_par()

    def selected_struct(source: Any, fields: tuple[str, ...]) -> Any:
        selected = source.select(*fields)
        missing = hail_module.struct(
            **{field: hail_module.missing(source[field].dtype) for field in fields}
        )
        return hail_module.or_else(selected, missing)

    def count_expr(source: Any, adjusted: int, xy_adjusted: int) -> Any:
        adjusted_frequency = source.freq[adjusted]
        frequency = selected_struct(
            adjusted_frequency, ("AC", "AN", "homozygote_count")
        )
        xy_frequency = selected_struct(source.freq[xy_adjusted], ("AC",))
        adjusted_count_is_defined = (
            hail_module.is_defined(frequency.AC)
            | hail_module.is_defined(frequency.AN)
            | hail_module.is_defined(frequency.homozygote_count)
        )
        return hail_module.struct(
            ac=frequency.AC,
            an=frequency.AN,
            homozygote_count=frequency.homozygote_count,
            hemizygote_count=hail_module.if_else(
                in_autosome_or_par,
                hail_module.or_missing(adjusted_count_is_defined, 0),
                xy_frequency.AC,
            ),
        )

    def faf_expr(source: Any, index: int | None) -> Any:
        if index is None:
            return hail_module.struct(
                faf95=hail_module.missing(hail_module.tfloat64),
                faf99=hail_module.missing(hail_module.tfloat64),
            )
        return selected_struct(source.faf[index], ("faf95", "faf99"))

    def dataset_expr(indices: DatasetProjectionIndices) -> Any:
        source = table[indices.row_field]
        ancestry_expressions = {}
        for ancestry in indices.ancestries:
            count = count_expr(
                source,
                ancestry.adjusted_frequency,
                ancestry.xy_adjusted_frequency,
            )
            faf = faf_expr(source, ancestry.adjusted_faf)
            ancestry_expressions[ancestry.genetic_ancestry] = hail_module.struct(
                ac=count.ac,
                an=count.an,
                homozygote_count=count.homozygote_count,
                hemizygote_count=count.hemizygote_count,
                faf95=faf.faf95,
                faf99=faf.faf99,
            )
        fields = {
            "adjusted": count_expr(
                source, indices.adjusted_frequency, indices.xy_adjusted_frequency
            ),
            "ancestries": hail_module.struct(**ancestry_expressions),
            "faf": faf_expr(source, indices.adjusted_faf),
            "grpmax": selected_struct(
                source.grpmax,
                ("AC", "AF", "AN", "homozygote_count", "gen_anc"),
            ),
            "fafmax": selected_struct(
                source.fafmax,
                (
                    "faf95_max",
                    "faf95_max_gen_anc",
                    "faf99_max",
                    "faf99_max_gen_anc",
                ),
            ),
        }
        if indices.has_filters:
            fields["filters"] = hail_module.sorted(hail_module.array(source.filters))
        return hail_module.struct(**fields)

    joint = layout.dataset("joint")
    stats = table.freq_comparison_stats
    stat_union_p = stats.stat_union.p_value
    ancestry_discordance = {}
    for ancestry in joint.ancestries:
        ancestry_discordance[ancestry.genetic_ancestry] = selected_struct(
            stats.contingency_table_test[ancestry.adjusted_frequency],
            ("p_value", "odds_ratio"),
        )
    overall_test = selected_struct(
        stats.contingency_table_test[joint.adjusted_frequency],
        ("p_value", "odds_ratio"),
    )
    cmh_test = selected_struct(
        stats.cochran_mantel_haenszel_test, ("p_value", "chisq")
    )
    stat_union = selected_struct(
        stats.stat_union, ("p_value", "stat_test_name", "gen_ancs")
    )
    return table.select(
        **{dataset.name: dataset_expr(dataset) for dataset in layout.datasets},
        region_flags=hail_module.struct(
            **{flag: table.region_flags[flag] for flag in layout.region_flags}
        ),
        discordance=hail_module.struct(
            overall=hail_module.struct(
                p_value=overall_test.p_value,
                odds_ratio=overall_test.odds_ratio,
            ),
            ancestries=hail_module.struct(**ancestry_discordance),
            cochran_mantel_haenszel_test=hail_module.struct(
                p_value=cmh_test.p_value,
                chisq=cmh_test.chisq,
            ),
            stat_union=hail_module.struct(
                p_value=stat_union.p_value,
                stat_test_name=stat_union.stat_test_name,
                gen_ancs=stat_union.gen_ancs,
            ),
            is_discrepant=hail_module.or_missing(
                hail_module.is_defined(stat_union_p),
                stat_union_p < layout.discrepant_p_value_threshold,
            ),
        ),
    )


def gnomad_vcf_metadata(
    contract: GnomadProjectionContract,
) -> dict[str, dict[str, dict[str, str]]]:
    """Return stable metadata for the BCF-oriented VCF INFO projection."""
    ancestry_order = ",".join(contract.major_genetic_ancestries)
    info: dict[str, dict[str, str]] = {}

    def add(
        identifier: str,
        number: str,
        field_type: str,
        description: str,
    ) -> None:
        info[identifier] = {
            "Description": description,
            "Number": number,
            "Type": field_type,
        }

    for dataset in contract.datasets:
        prefix = DATASET_VCF_PREFIXES[dataset.name]
        label = dataset.name.removesuffix("s")
        for suffix, description in (
            ("AC", "adjusted alternate allele count"),
            ("AN", "adjusted allele number"),
            ("NHOMALT", "adjusted homozygous alternate count"),
            ("NHEMI", "adjusted hemizygous alternate count"),
        ):
            add(f"{prefix}_{suffix}", "1", "Integer", f"gnomAD {label} {description}")
            add(
                f"{prefix}_{suffix}_ANC",
                str(len(contract.major_genetic_ancestries)),
                "Integer",
                f"gnomAD {label} {description} by ancestry ({ancestry_order})",
            )
        for suffix, description in (
            ("FAF95", "95% filtering allele frequency"),
            ("FAF99", "99% filtering allele frequency"),
        ):
            add(f"{prefix}_{suffix}", "1", "Float", f"gnomAD {label} {description}")
            add(
                f"{prefix}_{suffix}_ANC",
                str(len(contract.major_genetic_ancestries)),
                "Float",
                f"gnomAD {label} {description} by ancestry ({ancestry_order})",
            )
        for suffix, field_type, description in (
            ("GRPMAX_AC", "Integer", "grpmax alternate allele count"),
            ("GRPMAX_AF", "Float", "grpmax allele frequency"),
            ("GRPMAX_AN", "Integer", "grpmax allele number"),
            ("GRPMAX_NHOMALT", "Integer", "grpmax homozygous alternate count"),
            ("GRPMAX_ANC", "String", "grpmax ancestry"),
            ("FAFMAX95", "Float", "maximum ancestry 95% filtering allele frequency"),
            ("FAFMAX95_ANC", "String", "ancestry for maximum 95% filtering allele frequency"),
            ("FAFMAX99", "Float", "maximum ancestry 99% filtering allele frequency"),
            ("FAFMAX99_ANC", "String", "ancestry for maximum 99% filtering allele frequency"),
        ):
            add(f"{prefix}_{suffix}", "1", field_type, f"gnomAD {label} {description}")
        if dataset.has_filters:
            add(
                f"{prefix}_FILTERS",
                ".",
                "String",
                f"gnomAD {label} source filters",
            )
            add(
                f"{prefix}_FILTERS_DEFINED",
                "1",
                "Integer",
                f"gnomAD {label} source-filter availability (0 missing, 1 defined)",
            )

    for flag in contract.region_flags:
        add(
            f"REGION_{flag.upper()}",
            "1",
            "Integer",
            f"gnomAD region flag {flag} (0 false, 1 true, missing unknown)",
        )
    for suffix, number, field_type, description in (
        ("DISC_OVERALL_P", "1", "Float", "overall exome-genome test p-value"),
        ("DISC_OVERALL_OR", "1", "Float", "overall exome-genome odds ratio"),
        (
            "DISC_ANC_P",
            str(len(contract.major_genetic_ancestries)),
            "Float",
            f"exome-genome test p-value by ancestry ({ancestry_order})",
        ),
        (
            "DISC_ANC_OR",
            str(len(contract.major_genetic_ancestries)),
            "Float",
            f"exome-genome odds ratio by ancestry ({ancestry_order})",
        ),
        ("DISC_CMH_P", "1", "Float", "Cochran-Mantel-Haenszel p-value"),
        ("DISC_CMH_CHISQ", "1", "Float", "Cochran-Mantel-Haenszel chi-square"),
        ("DISC_UNION_P", "1", "Float", "source stat-union p-value"),
        ("DISC_UNION_TEST", "1", "String", "source stat-union test name"),
        ("DISC_UNION_ANC", ".", "String", "source stat-union ancestries"),
        (
            "DISC",
            "1",
            "Integer",
            "release-intent exome-genome discordance (0 false, 1 true, missing unknown)",
        ),
    ):
        add(suffix, number, field_type, description)
    return {"info": info}


def _nullable_boolean_as_integer(value: Any, hail_module: Any) -> Any:
    return hail_module.or_missing(
        hail_module.is_defined(value),
        hail_module.if_else(value, 1, 0),
    )


def project_hail_projection_to_vcf(
    table: Any,
    hail_module: Any,
    contract: GnomadProjectionContract,
) -> Any:
    """Encode every retained logical field in VCF-compatible INFO expressions."""
    info: dict[str, Any] = {}
    ancestries = contract.major_genetic_ancestries

    def ancestry_array(dataset: Any, field: str) -> Any:
        return hail_module.array(
            [dataset.ancestries[ancestry][field] for ancestry in ancestries]
        )

    for dataset_contract in contract.datasets:
        dataset = table[dataset_contract.name]
        prefix = DATASET_VCF_PREFIXES[dataset_contract.name]
        info.update(
            {
                f"{prefix}_AC": dataset.adjusted.ac,
                f"{prefix}_AN": dataset.adjusted.an,
                f"{prefix}_NHOMALT": dataset.adjusted.homozygote_count,
                f"{prefix}_NHEMI": dataset.adjusted.hemizygote_count,
                f"{prefix}_AC_ANC": ancestry_array(dataset, "ac"),
                f"{prefix}_AN_ANC": ancestry_array(dataset, "an"),
                f"{prefix}_NHOMALT_ANC": ancestry_array(
                    dataset, "homozygote_count"
                ),
                f"{prefix}_NHEMI_ANC": ancestry_array(dataset, "hemizygote_count"),
                f"{prefix}_FAF95": dataset.faf.faf95,
                f"{prefix}_FAF99": dataset.faf.faf99,
                f"{prefix}_FAF95_ANC": ancestry_array(dataset, "faf95"),
                f"{prefix}_FAF99_ANC": ancestry_array(dataset, "faf99"),
                f"{prefix}_GRPMAX_AC": dataset.grpmax.AC,
                f"{prefix}_GRPMAX_AF": dataset.grpmax.AF,
                f"{prefix}_GRPMAX_AN": dataset.grpmax.AN,
                f"{prefix}_GRPMAX_NHOMALT": dataset.grpmax.homozygote_count,
                f"{prefix}_GRPMAX_ANC": dataset.grpmax.gen_anc,
                f"{prefix}_FAFMAX95": dataset.fafmax.faf95_max,
                f"{prefix}_FAFMAX95_ANC": dataset.fafmax.faf95_max_gen_anc,
                f"{prefix}_FAFMAX99": dataset.fafmax.faf99_max,
                f"{prefix}_FAFMAX99_ANC": dataset.fafmax.faf99_max_gen_anc,
            }
        )
        if dataset_contract.has_filters:
            info[f"{prefix}_FILTERS"] = dataset.filters
            info[f"{prefix}_FILTERS_DEFINED"] = hail_module.if_else(
                hail_module.is_defined(dataset.filters), 1, 0
            )

    for flag in contract.region_flags:
        info[f"REGION_{flag.upper()}"] = _nullable_boolean_as_integer(
            table.region_flags[flag], hail_module
        )
    discordance = table.discordance
    info.update(
        {
            "DISC_OVERALL_P": discordance.overall.p_value,
            "DISC_OVERALL_OR": discordance.overall.odds_ratio,
            "DISC_ANC_P": hail_module.array(
                [discordance.ancestries[item].p_value for item in ancestries]
            ),
            "DISC_ANC_OR": hail_module.array(
                [discordance.ancestries[item].odds_ratio for item in ancestries]
            ),
            "DISC_CMH_P": discordance.cochran_mantel_haenszel_test.p_value,
            "DISC_CMH_CHISQ": discordance.cochran_mantel_haenszel_test.chisq,
            "DISC_UNION_P": discordance.stat_union.p_value,
            "DISC_UNION_TEST": discordance.stat_union.stat_test_name,
            "DISC_UNION_ANC": hail_module.sorted(
                hail_module.array(discordance.stat_union.gen_ancs)
            ),
            "DISC": _nullable_boolean_as_integer(
                discordance.is_discrepant, hail_module
            ),
        }
    )
    return table.select(
        filters=hail_module.empty_set(hail_module.tstr),
        info=hail_module.struct(**info),
    )
