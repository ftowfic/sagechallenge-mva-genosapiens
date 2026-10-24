"""Qualified gnomAD joint FAF95 field mapping and float32 uncertainty."""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import re
import struct
from typing import Any, Literal, Mapping

from genosapiens.evidence_model import PopulationEvidence
from genosapiens.gnomad_projection import (
    gnomad_vcf_metadata, load_gnomad_projection_contract,
)
from genosapiens.gnomad_publication import load_gnomad_publication_contract


CONTRACT_SCHEMA = "genosapiens.gnomad-faf95-precision-contract"
REPORT_SCHEMA = "genosapiens.gnomad-faf95-precision-qualification"
QUALIFICATION_ID = "gnomad-joint-v4.1-faf95-float32-precision-v1"
THRESHOLD = Fraction(1, 20)
_SHA = re.compile(r"[0-9a-f]{64}")
_ANCESTRIES = ("afr", "amr", "eas", "mid", "nfe", "sas")
_ERRORS = (ValueError, TypeError, KeyError, AttributeError, ArithmeticError,
           UnicodeError, json.JSONDecodeError)


class GnomadFaf95AdapterError(ValueError):
    """Static diagnostics only; no allele or resource paths."""


def _fail() -> None:
    raise GnomadFaf95AdapterError("invalid gnomAD FAF95 precision evidence")


def _canonical(value: object) -> bytes:
    try:
        return (json.dumps(value, sort_keys=True, separators=(",", ":"),
                           ensure_ascii=True, allow_nan=False) + "\n").encode("ascii")
    except _ERRORS:
        _fail()


def _sha(value: object) -> str:
    if type(value) is not str or _SHA.fullmatch(value) is None:
        _fail()
    return value


def _mapping(value: object, fields: set[str]) -> Mapping[str, Any]:
    if type(value) is not dict or set(value) != fields:
        _fail()
    return value


def _text(value: object) -> str:
    if type(value) is not str or not value or value.strip() != value or any(
            character in value for character in "\r\n\t"):
        _fail()
    return value


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as source:
            while chunk := source.read(1024 * 1024):
                digest.update(chunk)
    except OSError:
        _fail()
    return digest.hexdigest()


def _bound_file(root: Path, relative: str, sha256: str) -> Path:
    try:
        path = root / relative
        if (Path(relative).is_absolute() or not path.resolve(strict=True).is_relative_to(root)
                or path.is_symlink() or not path.is_file() or _file_sha256(path) != sha256):
            _fail()
        return path
    except OSError:
        _fail()


@dataclass(frozen=True, slots=True)
class BoundarySourceValue:
    label: str
    value: str
    expected_bcf_token: str


@dataclass(frozen=True, slots=True)
class GnomadFaf95PrecisionContract:
    descriptor_sha256: str
    projection_contract_path: str
    projection_contract_sha256: str
    publication_contract_path: str
    publication_contract_sha256: str
    manifest_path: str
    manifest_sha256: str
    converter_path: str
    converter_sha256: str
    ancestry_order: tuple[str, ...]
    boundary_values: tuple[BoundarySourceValue, ...]
    boundary_ancestry: str
    boundary_ac: int
    boundary_an: int


@dataclass(frozen=True, slots=True)
class Float32Interval:
    raw_token: str
    stored_value: Fraction
    lower_bound: Fraction
    upper_bound: Fraction
    stored_hex: str

    @property
    def threshold_relation(self) -> Literal[
        "upper_bound_below_threshold", "overlaps_or_touches_threshold",
        "lower_bound_above_threshold",
    ]:
        if self.lower_bound > THRESHOLD:
            return "lower_bound_above_threshold"
        if self.upper_bound < THRESHOLD:
            return "upper_bound_below_threshold"
        return "overlaps_or_touches_threshold"


@dataclass(frozen=True, slots=True)
class GnomadFaf95PrecisionQualification:
    report_sha256: str
    descriptor_sha256: str
    projection_contract_sha256: str
    publication_manifest_sha256: str
    converter_sha256: str
    ancestry_order: tuple[str, ...]
    runtime_image: str


@dataclass(frozen=True, slots=True)
class JointFaf95Facts:
    maximum_raw_token: str | None
    maximum_ancestry: str | None
    matching_ancestry_index: int | None
    matching_ancestry_ac: int | None
    matching_ancestry_an: int | None
    precision_interval: Float32Interval | None
    precision_status: Literal[
        "qualified", "missing_qualification", "unassessable_numeric",
    ]
    projection_filter: str
    projection_filter_is_source_pass: Literal[False]
    maximum_record: PopulationEvidence
    matching_ancestry_record: PopulationEvidence | None


def load_gnomad_faf95_precision_contract(
    content: bytes, expected_sha256: str,
) -> GnomadFaf95PrecisionContract:
    try:
        if (type(content) is not bytes or _sha(expected_sha256) != hashlib.sha256(content).hexdigest()):
            _fail()
        value = json.loads(content.decode("ascii"))
        if _canonical(value) != content:
            _fail()
        top = _mapping(value, {
            "schema", "version", "qualification_id", "classification",
            "controlled_use_authorized", "projection", "converter", "mapping",
            "storage", "boundary_fixture",
        })
        if (top["schema"] != CONTRACT_SCHEMA or top["version"] != 1
                or top["qualification_id"] != QUALIFICATION_ID
                or top["classification"] != "public"
                or top["controlled_use_authorized"] is not False):
            _fail()
        projection = _mapping(top["projection"], {
            "contract_path", "contract_sha256", "publication_contract_path",
            "publication_contract_sha256", "manifest_path", "manifest_sha256",
        })
        converter = _mapping(top["converter"], {"path", "sha256", "symbol"})
        mapping = _mapping(top["mapping"], {
            "maximum_frequency_field", "maximum_ancestry_field", "ancestry_an_field",
            "ancestry_ac_field", "ancestry_order", "upstream_frequency_path",
            "upstream_ancestry_path", "projection_filter_is_source_pass",
        })
        storage = _mapping(top["storage"], {
            "float_encoding", "rounding", "text_export", "comparison",
            "text_significant_decimal_digits",
        })
        boundary = _mapping(top["boundary_fixture"], {
            "source_values", "ancestry", "ancestry_ac", "ancestry_an",
        })
        if (converter["symbol"] != "project_hail_projection_to_vcf"
                or mapping != {
                    "maximum_frequency_field": "JT_FAFMAX95",
                    "maximum_ancestry_field": "JT_FAFMAX95_ANC",
                    "ancestry_an_field": "JT_AN_ANC", "ancestry_ac_field": "JT_AC_ANC",
                    "ancestry_order": list(_ANCESTRIES),
                    "upstream_frequency_path": "row.joint.fafmax.faf95_max",
                    "upstream_ancestry_path": "row.joint.fafmax.faf95_max_gen_anc",
                    "projection_filter_is_source_pass": False,
                } or storage != {
                    "float_encoding": "ieee754-binary32",
                    "rounding": "round_to_nearest_ties_to_even",
                    "text_export": "six_significant_decimal_digits_then_binary32",
                    "text_significant_decimal_digits": 6,
                    "comparison": "lower_bound_strictly_greater_than_1_over_20",
                }):
            _fail()
        source_values = tuple(BoundarySourceValue(
            _text(item["label"]), _text(item["value"]), _text(item["expected_bcf_token"])
        ) for item in boundary["source_values"] if type(item) is dict
            and set(item) == {"label", "value", "expected_bcf_token"})
        if (len(source_values) != 6
                or tuple(item.label for item in source_values) != (
                    "immediately_below", "equal", "immediately_above",
                    "inside_text_rounding_bin", "outside_text_rounding_bin",
                    "separably_above")
                or boundary["ancestry"] != "nfe"
                or type(boundary["ancestry_ac"]) is not int
                or type(boundary["ancestry_an"]) is not int
                or not 0 < boundary["ancestry_ac"] <= boundary["ancestry_an"]):
            _fail()
        return GnomadFaf95PrecisionContract(
            expected_sha256, _text(projection["contract_path"]),
            _sha(projection["contract_sha256"]),
            _text(projection["publication_contract_path"]),
            _sha(projection["publication_contract_sha256"]),
            _text(projection["manifest_path"]), _sha(projection["manifest_sha256"]),
            _text(converter["path"]), _sha(converter["sha256"]), _ANCESTRIES,
            source_values, boundary["ancestry"], boundary["ancestry_ac"],
            boundary["ancestry_an"],
        )
    except GnomadFaf95AdapterError:
        raise
    except _ERRORS:
        _fail()


def authenticate_gnomad_faf95_precision_contract(
    contract: GnomadFaf95PrecisionContract, *, repository_root: str | Path,
    reference_root: str | Path,
) -> None:
    try:
        if type(contract) is not GnomadFaf95PrecisionContract:
            _fail()
        repository = Path(repository_root).resolve(strict=True)
        references = Path(reference_root).resolve(strict=True)
        projection_path = _bound_file(
            repository, contract.projection_contract_path, contract.projection_contract_sha256
        )
        publication_path = _bound_file(
            repository, contract.publication_contract_path,
            contract.publication_contract_sha256,
        )
        _bound_file(repository, contract.converter_path, contract.converter_sha256)
        _bound_file(references, contract.manifest_path, contract.manifest_sha256)
        projection = load_gnomad_projection_contract(projection_path)
        publication = load_gnomad_publication_contract(publication_path)
        metadata = gnomad_vcf_metadata(projection)["info"]
        if (projection.major_genetic_ancestries != contract.ancestry_order
                or publication.projection_id != projection.projection_id
                or publication.manifest_path != contract.manifest_path
                or metadata["JT_FAFMAX95"] != {
                    "Description": "gnomAD joint maximum ancestry 95% filtering allele frequency",
                    "Number": "1", "Type": "Float",
                }
                or metadata["JT_FAFMAX95_ANC"]["Type"] != "String"
                or metadata["JT_AN_ANC"]["Number"] != "6"
                or metadata["JT_AN_ANC"]["Type"] != "Integer"
                or metadata["JT_AC_ANC"]["Number"] != "6"
                or metadata["JT_AC_ANC"]["Type"] != "Integer"):
            _fail()
    except GnomadFaf95AdapterError:
        raise
    except _ERRORS + (OSError,):
        _fail()


def float32_interval(raw_token: str) -> Float32Interval:
    """Return the full source interval that can round to one BCF Float value."""
    try:
        token = _text(raw_token)
        decimal = Decimal(token)
        if not decimal.is_finite() or not Decimal(0) <= decimal <= Decimal(1):
            _fail()
        packed = struct.pack("!f", float(decimal))
        bits = struct.unpack("!I", packed)[0]
        value = struct.unpack("!f", packed)[0]
        if not math.isfinite(value) or value < 0 or value > 1:
            _fail()
        current = Fraction(*value.as_integer_ratio())
        previous_value = 0.0 if bits == 0 else struct.unpack("!f", struct.pack("!I", bits - 1))[0]
        next_value = struct.unpack("!f", struct.pack("!I", bits + 1))[0]
        previous = Fraction(*previous_value.as_integer_ratio())
        following = Fraction(*next_value.as_integer_ratio())
        stored_lower = (previous + current) / 2
        stored_upper = (current + following) / 2
        # Hail 0.2.139's VCF exporter renders six significant decimal digits
        # before bcftools stores an IEEE-754 binary32. Account for both lossy
        # steps. Suppressed trailing zeros do not narrow this fixed precision.
        decimal_exponent = decimal.adjusted() - 6 + 1
        text_half_unit = Fraction(5, 10) * (
            Fraction(10**decimal_exponent, 1)
            if decimal_exponent >= 0 else Fraction(1, 10 ** -decimal_exponent)
        )
        decimal_fraction = Fraction(decimal)
        return Float32Interval(
            token, current, min(stored_lower, decimal_fraction - text_half_unit),
            max(stored_upper, decimal_fraction + text_half_unit), packed.hex(),
        )
    except GnomadFaf95AdapterError:
        raise
    except (InvalidOperation, OverflowError, struct.error, ValueError):
        _fail()


def _ratio(value: Fraction) -> Mapping[str, int]:
    return {"numerator": value.numerator, "denominator": value.denominator}


def _raw_integer_at(value: object, index: int) -> int | None:
    if value is None or value == ".":
        return None
    text = _text(value)
    items = text.split(",")
    if len(items) != len(_ANCESTRIES):
        _fail()
    item = items[index]
    if item == ".":
        return None
    if re.fullmatch(r"0|[1-9][0-9]*", item) is None:
        _fail()
    return int(item)


def qualification_report(
    contract: GnomadFaf95PrecisionContract, *, boundary_rows: tuple[Mapping[str, str], ...],
    runtime_image: str, hail_version: str, bcftools_version: str,
) -> bytes:
    """Encode results produced by the tiny real Hail-to-BCF boundary fixture."""
    try:
        if (type(contract) is not GnomadFaf95PrecisionContract
                or type(boundary_rows) is not tuple
                or len(boundary_rows) != len(contract.boundary_values)
                or re.fullmatch(r"[^\s]+@sha256:[0-9a-f]{64}", runtime_image) is None):
            _fail()
        cases = []
        for expected, observed in zip(contract.boundary_values, boundary_rows, strict=True):
            row = _mapping(observed, {
                "label", "source_value", "bcf_token", "filter", "ancestry",
                "ancestry_ac", "ancestry_an",
            })
            if (row["label"] != expected.label or row["source_value"] != expected.value
                    or row["bcf_token"] != expected.expected_bcf_token
                    or row["filter"] != "PASS" or row["ancestry"] != contract.boundary_ancestry
                    or row["ancestry_ac"] != str(contract.boundary_ac)
                    or row["ancestry_an"] != str(contract.boundary_an)):
                _fail()
            interval = float32_interval(row["bcf_token"])
            cases.append({
                "label": expected.label, "source_value": expected.value,
                "bcf_token": row["bcf_token"], "stored_float32_hex": interval.stored_hex,
                "stored_value": _ratio(interval.stored_value),
                "lower_bound": _ratio(interval.lower_bound),
                "upper_bound": _ratio(interval.upper_bound),
                "threshold_relation": interval.threshold_relation,
            })
        if ([case["threshold_relation"] for case in cases[:4]]
                != ["overlaps_or_touches_threshold"] * 4
                or [case["threshold_relation"] for case in cases[4:]]
                != ["lower_bound_above_threshold"] * 2):
            _fail()
        return _canonical({
            "schema": REPORT_SCHEMA, "version": 1, "qualification_id": QUALIFICATION_ID,
            "classification": "public", "controlled_use_authorized": False,
            "status": "qualified", "descriptor_sha256": contract.descriptor_sha256,
            "bindings": {
                "projection_contract_sha256": contract.projection_contract_sha256,
                "publication_manifest_sha256": contract.manifest_sha256,
                "converter_sha256": contract.converter_sha256,
            },
            "mapping": {
                "maximum_frequency_field": "JT_FAFMAX95",
                "maximum_ancestry_field": "JT_FAFMAX95_ANC",
                "matching_ancestry_an_field": "JT_AN_ANC",
                "matching_ancestry_ac_field": "JT_AC_ANC",
                "ancestry_order": list(contract.ancestry_order),
                "projection_filter_is_source_pass": False,
            },
            "runtime": {"image": runtime_image, "hail_version": _text(hail_version),
                        "bcftools_version": _text(bcftools_version)},
            "boundary_cases": cases,
        })
    except GnomadFaf95AdapterError:
        raise
    except _ERRORS:
        _fail()


def verify_gnomad_faf95_precision_qualification(
    content: bytes, expected_sha256: str, contract: GnomadFaf95PrecisionContract,
) -> GnomadFaf95PrecisionQualification:
    try:
        if (type(content) is not bytes or type(contract) is not GnomadFaf95PrecisionContract
                or _sha(expected_sha256) != hashlib.sha256(content).hexdigest()):
            _fail()
        report = json.loads(content.decode("ascii"))
        if _canonical(report) != content:
            _fail()
        top = _mapping(report, {
            "schema", "version", "qualification_id", "classification",
            "controlled_use_authorized", "status", "descriptor_sha256", "bindings",
            "mapping", "runtime", "boundary_cases",
        })
        if (top["schema"] != REPORT_SCHEMA or top["version"] != 1
                or top["qualification_id"] != QUALIFICATION_ID
                or top["classification"] != "public"
                or top["controlled_use_authorized"] is not False
                or top["status"] != "qualified"
                or top["descriptor_sha256"] != contract.descriptor_sha256):
            _fail()
        bindings = _mapping(top["bindings"], {
            "projection_contract_sha256", "publication_manifest_sha256", "converter_sha256",
        })
        if bindings != {
            "projection_contract_sha256": contract.projection_contract_sha256,
            "publication_manifest_sha256": contract.manifest_sha256,
            "converter_sha256": contract.converter_sha256,
        }:
            _fail()
        mapping = _mapping(top["mapping"], {
            "maximum_frequency_field", "maximum_ancestry_field",
            "matching_ancestry_an_field", "matching_ancestry_ac_field",
            "ancestry_order", "projection_filter_is_source_pass",
        })
        if mapping != {
            "maximum_frequency_field": "JT_FAFMAX95",
            "maximum_ancestry_field": "JT_FAFMAX95_ANC",
            "matching_ancestry_an_field": "JT_AN_ANC",
            "matching_ancestry_ac_field": "JT_AC_ANC",
            "ancestry_order": list(contract.ancestry_order),
            "projection_filter_is_source_pass": False,
        }:
            _fail()
        runtime = _mapping(top["runtime"], {"image", "hail_version", "bcftools_version"})
        rows = []
        for case in top["boundary_cases"]:
            current = _mapping(case, {
                "label", "source_value", "bcf_token", "stored_float32_hex",
                "stored_value", "lower_bound", "upper_bound", "threshold_relation",
            })
            interval = float32_interval(current["bcf_token"])
            if (current["stored_float32_hex"] != interval.stored_hex
                    or current["stored_value"] != _ratio(interval.stored_value)
                    or current["lower_bound"] != _ratio(interval.lower_bound)
                    or current["upper_bound"] != _ratio(interval.upper_bound)
                    or current["threshold_relation"] != interval.threshold_relation):
                _fail()
            rows.append(current)
        expected = contract.boundary_values
        if (tuple((row["label"], row["source_value"], row["bcf_token"]) for row in rows)
                != tuple((item.label, item.value, item.expected_bcf_token) for item in expected)
                or [row["threshold_relation"] for row in rows[:4]]
                != ["overlaps_or_touches_threshold"] * 4
                or [row["threshold_relation"] for row in rows[4:]]
                != ["lower_bound_above_threshold"] * 2):
            _fail()
        return GnomadFaf95PrecisionQualification(
            expected_sha256, contract.descriptor_sha256,
            contract.projection_contract_sha256, contract.manifest_sha256,
            contract.converter_sha256, contract.ancestry_order, _text(runtime["image"]),
        )
    except GnomadFaf95AdapterError:
        raise
    except _ERRORS:
        _fail()


def adapt_joint_faf95(
    records: tuple[PopulationEvidence, ...],
    qualification: GnomadFaf95PrecisionQualification | None,
) -> JointFaf95Facts:
    """Select the publisher maximum and its exact ancestry counts without fallback."""
    try:
        if (type(records) is not tuple
                or any(type(record) is not PopulationEvidence for record in records)
                or qualification is not None
                and type(qualification) is not GnomadFaf95PrecisionQualification):
            _fail()
        maxima = tuple(record for record in records
                       if record.subset == "joint" and record.stratum == "fafmax95")
        if len(maxima) != 1:
            _fail()
        maximum = maxima[0]
        fields = maximum.source.source_fields
        if (not isinstance(fields, Mapping)
                or fields.get("schema") != "genosapiens.population-annotation"):
            _fail()
        order = tuple(fields.get("ancestry_order", ()))
        if order != _ANCESTRIES or qualification is not None and order != qualification.ancestry_order:
            _fail()
        raw = fields.get("raw_fields")
        if (not isinstance(raw, Mapping)
                or set(raw) != {"JT_FAFMAX95", "JT_FAFMAX95_ANC"}):
            _fail()
        token = raw["JT_FAFMAX95"]
        ancestry = maximum.ancestry
        raw_ancestry = raw["JT_FAFMAX95_ANC"]
        if ((ancestry is None and raw_ancestry not in (None, "."))
                or ancestry is not None and raw_ancestry != ancestry):
            _fail()
        if token not in (None, "."):
            decimal_token = Decimal(_text(token))
            parsed = float(decimal_token)
            if (decimal_token.is_finite()
                    and maximum.filtering_allele_frequency != parsed
                    or not decimal_token.is_finite()
                    and maximum.filtering_allele_frequency is not None):
                _fail()
        ancestry_record = None
        index = ac = an = None
        if ancestry is not None:
            if type(ancestry) is not str or ancestry not in order:
                _fail()
            index = order.index(ancestry)
            matches = tuple(record for record in records if
                            record.subset == "joint" and record.stratum == "ancestry"
                            and record.ancestry == ancestry)
            if len(matches) != 1:
                _fail()
            ancestry_record = matches[0]
            ancestry_fields = ancestry_record.source.source_fields
            ancestry_raw = (ancestry_fields.get("raw_fields")
                            if isinstance(ancestry_fields, Mapping) else None)
            if (not isinstance(ancestry_fields, Mapping)
                    or ancestry_fields.get("ancestry_index") != index
                    or not isinstance(ancestry_raw, Mapping)
                    or "JT_AC_ANC" not in ancestry_raw or "JT_AN_ANC" not in ancestry_raw
                    or ancestry_record.source.source_sha256 != maximum.source.source_sha256
                    or ancestry_record.source.public_manifest_sha256
                    != maximum.source.public_manifest_sha256):
                _fail()
            ac, an = ancestry_record.allele_count, ancestry_record.allele_number
            if (_raw_integer_at(ancestry_raw["JT_AC_ANC"], index) != ac
                    or _raw_integer_at(ancestry_raw["JT_AN_ANC"], index) != an):
                _fail()
        interval = None
        status: Literal["qualified", "missing_qualification", "unassessable_numeric"]
        if token is None or maximum.filtering_allele_frequency is None:
            status = "unassessable_numeric"
        elif qualification is None:
            _text(token)
            status = "missing_qualification"
        else:
            interval = float32_interval(token)
            status = "qualified"
        projection_filter = fields.get("projection_filter")
        if type(projection_filter) is not str or not projection_filter:
            _fail()
        return JointFaf95Facts(
            token, ancestry, index, ac, an, interval, status, projection_filter, False,
            maximum, ancestry_record,
        )
    except GnomadFaf95AdapterError:
        raise
    except _ERRORS:
        _fail()
