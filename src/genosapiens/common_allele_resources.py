"""Authentication for public common-allele exception source snapshots."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import re
from typing import Any


DESCRIPTOR_SCHEMA = "genosapiens.track1-common-allele-exception-resource"
SNAPSHOT_SCHEMA = "genosapiens.track1-common-allele-exception-source-snapshot"
SHA256 = re.compile(r"^[0-9a-f]{64}$")
SAFE_PATH = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._/-]*$")
SVID = re.compile(r"^GN[0-9]{3}$")


class CommonAlleleResourceError(ValueError):
    """Static public-resource contract error."""


def _fail() -> None:
    raise CommonAlleleResourceError("invalid common-allele public resource contract")


def _pairs(items: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in items:
        if key in result:
            _fail()
        result[key] = value
    return result


def _constant(_value: str) -> None:
    _fail()


def _canonical(value: object) -> bytes:
    try:
        return (json.dumps(value, sort_keys=True, separators=(",", ":"),
                           ensure_ascii=True, allow_nan=False) + "\n").encode("ascii")
    except (ValueError, TypeError, UnicodeError, RecursionError, OverflowError):
        _fail()


def _object(value: object, fields: set[str]) -> dict[str, Any]:
    if type(value) is not dict or set(value) != fields:
        _fail()
    return value


def _sha(value: object) -> str:
    if type(value) is not str or SHA256.fullmatch(value) is None:
        _fail()
    return value


def _count(value: object, *, positive: bool = False) -> int:
    if type(value) is not int or value < (1 if positive else 0) or value > 2**63 - 1:
        _fail()
    return value


def _relative(value: object) -> str:
    if (type(value) is not str or SAFE_PATH.fullmatch(value) is None
            or Path(value).is_absolute() or ".." in Path(value).parts):
        _fail()
    return value


@dataclass(frozen=True, slots=True)
class CommonAlleleExceptionDescriptor:
    content: bytes
    sha256: str
    snapshot_path: str
    snapshot_sha256: str
    snapshot_size_bytes: int
    retrieved_date: str
    listing_count: int
    released_document_count: int
    source_identities: tuple[tuple[str, str, int, str, str], ...]


@dataclass(frozen=True, slots=True)
class CommonAlleleSourceSnapshot:
    content: bytes
    sha256: str
    document_ids: tuple[str, ...]
    document_sha256s: tuple[str, ...]
    missing_gene_document_ids: tuple[str, ...]


def load_common_allele_exception_descriptor(
    content: bytes, expected_sha256: str
) -> CommonAlleleExceptionDescriptor:
    try:
        _sha(expected_sha256)
        if type(content) is not bytes or hashlib.sha256(content).hexdigest() != expected_sha256:
            _fail()
        data = json.loads(content.decode("ascii"), object_pairs_hook=_pairs,
                          parse_constant=_constant)
        if _canonical(data) != content:
            _fail()
        top = _object(data, {
            "schema", "version", "resource_id", "assembly", "retrieved_date", "license",
            "license_url", "source_snapshot", "upstream_bindings", "curation_scope",
            "output", "scientific_input_prohibition", "qualification",
        })
        if (top["schema"] != DESCRIPTOR_SCHEMA or top["version"] != 1
                or top["resource_id"] != "track1-common-allele-exceptions-v1"
                or top["assembly"] != "GRCh38" or top["license"] != "CC0-1.0"
                or top["license_url"] != "https://creativecommons.org/publicdomain/zero/1.0/"
                or re.fullmatch(r"20[0-9]{2}-[0-9]{2}-[0-9]{2}", top["retrieved_date"]) is None):
            _fail()
        snapshot = _object(top["source_snapshot"], {
            "path", "sha256", "size_bytes", "retrieved_date", "listing_count",
            "released_document_count", "document_count", "sources",
        })
        if (snapshot["retrieved_date"] != top["retrieved_date"]
                or _count(snapshot["document_count"], positive=True)
                != _count(snapshot["released_document_count"], positive=True)):
            _fail()
        sources = snapshot["sources"]
        if type(sources) is not list or len(sources) != 3:
            _fail()
        identities = []
        for source in sources:
            row = _object(source, {"role", "requested_url", "resolved_url", "size_bytes", "sha256"})
            if (type(row["role"]) is not str or type(row["requested_url"]) is not str
                    or type(row["resolved_url"]) is not str
                    or not row["requested_url"].startswith("https://")
                    or not row["resolved_url"].startswith("https://")):
                _fail()
            identities.append((row["role"], _sha(row["sha256"]),
                               _count(row["size_bytes"], positive=True),
                               row["requested_url"], row["resolved_url"]))
        if tuple(item[0] for item in identities) != (
            "clingen_ba1_exception_list", "clingen_terms", "cspec_registry_listing",
        ):
            _fail()
        prohibited = top["scientific_input_prohibition"]
        if (type(prohibited) is not list
                or not {"challenge_case", "diagnosis", "hpo", "previous_campaign_outcome",
                        "selected_rows", "validation_truth"}.issubset(prohibited)):
            _fail()
        qualification = _object(top["qualification"], {
            "source_snapshot_status", "projection_status", "controlled_use_authorized",
        })
        if qualification != {"source_snapshot_status": "complete",
                              "projection_status": "unqualified",
                              "controlled_use_authorized": False}:
            _fail()
        downstream = _object(top["curation_scope"], {
            "ba1_exception_list", "cspec_registry", "downstream_public_curation",
        })["downstream_public_curation"]
        _object(downstream, {"empty_source_declaration", "founder_records",
                             "reduced_penetrance_records", "hypomorphic_records", "policy_sha256"})
        if (downstream["empty_source_declaration"] is not True
                or any(downstream[key] != 0 for key in (
                    "founder_records", "reduced_penetrance_records", "hypomorphic_records"))
                or _sha(downstream["policy_sha256"]) is None):
            _fail()
        bindings = top["upstream_bindings"]
        if type(bindings) is not list or len(bindings) != 4:
            _fail()
        roles = []
        for binding in bindings:
            row = _object(binding, {"role", "path", "sha256"})
            roles.append(row["role"])
            _relative(row["path"])
            _sha(row["sha256"])
        if roles != ["clinvar", "open_targets", "reference_identity", "selected_downstream_policy"]:
            _fail()
        return CommonAlleleExceptionDescriptor(
            content, expected_sha256, _relative(snapshot["path"]), _sha(snapshot["sha256"]),
            _count(snapshot["size_bytes"], positive=True), top["retrieved_date"],
            _count(snapshot["listing_count"], positive=True),
            _count(snapshot["released_document_count"], positive=True), tuple(identities),
        )
    except CommonAlleleResourceError:
        raise
    except (ValueError, TypeError, KeyError, UnicodeError, AttributeError, RecursionError):
        _fail()


def verify_common_allele_source_snapshot(
    descriptor: CommonAlleleExceptionDescriptor, reference_root: str | Path
) -> CommonAlleleSourceSnapshot:
    try:
        if type(descriptor) is not CommonAlleleExceptionDescriptor:
            _fail()
        root = Path(reference_root).resolve(strict=True)
        path = root / descriptor.snapshot_path
        if not path.is_relative_to(root) or path.is_symlink() or not path.is_file():
            _fail()
        content = path.read_bytes()
        if (len(content) != descriptor.snapshot_size_bytes
                or hashlib.sha256(content).hexdigest() != descriptor.snapshot_sha256):
            _fail()
        data = json.loads(content.decode("ascii"), object_pairs_hook=_pairs,
                          parse_constant=_constant)
        if _canonical(data) != content:
            _fail()
        top = _object(data, {
            "schema", "version", "classification", "snapshot_id", "retrieved_date",
            "license", "license_url", "terms_url", "challenge_or_truth_inputs_used",
            "registry_listing_count", "registry_released_count", "sources",
            "registry_documents", "downstream_public_curation",
        })
        if (top["schema"] != SNAPSHOT_SCHEMA or top["version"] != 1
                or top["classification"] != "public"
                or top["retrieved_date"] != descriptor.retrieved_date
                or top["license"] != "CC0-1.0"
                or top["challenge_or_truth_inputs_used"] is not False
                or top["registry_listing_count"] != descriptor.listing_count
                or top["registry_released_count"] != descriptor.released_document_count):
            _fail()
        manifest_sources = top["sources"]
        if type(manifest_sources) is not list or len(manifest_sources) != len(descriptor.source_identities):
            _fail()
        for source, expected in zip(manifest_sources, descriptor.source_identities, strict=True):
            role, sha256, size, requested, resolved = expected
            if (source.get("role"), source.get("sha256"), source.get("size_bytes"),
                    source.get("requested_url"), source.get("resolved_url")) != expected:
                _fail()
            child = root / _relative(source.get("path"))
            if (not child.is_relative_to(root) or child.is_symlink() or not child.is_file()
                    or child.stat().st_size != size
                    or hashlib.sha256(child.read_bytes()).hexdigest() != sha256):
                _fail()
        documents = top["registry_documents"]
        if type(documents) is not list or len(documents) != descriptor.released_document_count:
            _fail()
        ids, hashes, missing_genes = [], [], []
        for row in documents:
            svi_id = row.get("svi_id")
            if type(svi_id) is not str or SVID.fullmatch(svi_id) is None:
                _fail()
            child = root / _relative(row.get("path"))
            sha256, size = _sha(row.get("sha256")), _count(row.get("size_bytes"), positive=True)
            if (not child.is_relative_to(root) or child.is_symlink() or not child.is_file()
                    or child.stat().st_size != size
                    or hashlib.sha256(child.read_bytes()).hexdigest() != sha256
                    or row.get("ba1_coverage") != "present"
                    or type(row.get("ba1_records")) is not list or not row["ba1_records"]
                    or hashlib.sha256(_canonical(row["ba1_records"])).hexdigest()
                    != _sha(row.get("ba1_records_sha256"))):
                _fail()
            for record in row["ba1_records"]:
                if (type(record) is not dict or type(record.get("rule_set_id")) is not str
                        or type(record.get("criteria_id")) is not str
                        or "applicability" not in record or "description" not in record):
                    _fail()
            if row.get("gene_coverage") == "missing_declared":
                if row.get("genes") != []:
                    _fail()
                missing_genes.append(svi_id)
            elif row.get("gene_coverage") != "present" or not row.get("genes"):
                _fail()
            ids.append(svi_id)
            hashes.append(sha256)
        if ids != sorted(set(ids)):
            _fail()
        downstream = top["downstream_public_curation"]
        if (type(downstream) is not dict or downstream.get("empty_source_declaration") is not True
                or any(downstream.get(key) != [] for key in (
                    "founder_records", "reduced_penetrance_records", "hypomorphic_records"))):
            _fail()
        return CommonAlleleSourceSnapshot(
            content, descriptor.snapshot_sha256, tuple(ids), tuple(hashes), tuple(missing_genes)
        )
    except CommonAlleleResourceError:
        raise
    except (ValueError, TypeError, KeyError, UnicodeError, AttributeError, OSError,
            RecursionError):
        _fail()
