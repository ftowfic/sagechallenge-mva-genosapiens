"""Bounded header-only sanitization, not a privacy or reference certifier."""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import os
from pathlib import Path
import re
from typing import Mapping
import uuid

from genosapiens.bounded_json import BoundedJSONError, read_bounded_bytes


MAX_HEADER_BYTES = 2 * 1024 * 1024
PUBLIC_ASSEMBLY_LABEL = "GCA_000001405.15"
NEUTRAL_DESCRIPTION = "Technical declaration"
PUBLIC_ACCESSION = re.compile(r"(?:GCA|GCF)_[0-9]{9}\.[1-9][0-9]*")
FIELD_ID = re.compile(r"[A-Za-z_][A-Za-z0-9_.]*")
# hts-specs 510c107ba89aa13113818c9e4fa3071283cb7f7b, VCF 4.3-4.5.
CONTIG_ID = re.compile(r"[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*")
TOKEN = re.compile(r'[^\s,<>"\\]+', re.ASCII)
ATTRIBUTE = re.compile(
    r'([A-Za-z_][A-Za-z0-9_.]*)=("(?:[^"\\\x00-\x1f\x7f]|\\["\\])*"|[^\s,<>"\\]+)(?:,|$)'
)
COLUMNS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
ATTRIBUTES = {
    "contig": ("ID", "length"),
    "ALT": ("ID", "Description"),
    "FILTER": ("ID", "Description"),
    "INFO": ("ID", "Number", "Type", "Description"),
    "FORMAT": ("ID", "Number", "Type", "Description"),
}
# Only interpreted/generated fields have semantic contracts. This is not an ID allowlist.
FIELD_CONTRACTS = {
    ("FORMAT", "GT"): {("1", "String")},
    ("FORMAT", "GQ"): {("1", "Integer")},
    ("FORMAT", "DP"): {("1", "Integer")},
    ("FORMAT", "AD"): {("R", "Integer")},
    ("FORMAT", "PL"): {("G", "Integer")},
    ("FORMAT", "PGT"): {("1", "String")},
    ("FORMAT", "PID"): {("1", "String")},
    ("FORMAT", "PS"): {("1", "Integer")},
    ("FORMAT", "CN"): {("1", "Integer"), ("1", "Float")},
    ("FORMAT", "REPCN"): {("1", "String")},
    ("INFO", "AC"): {("A", "Integer")},
    ("INFO", "AF"): {("A", "Float")},
    ("INFO", "AN"): {("1", "Integer")},
    ("INFO", "DP"): {("1", "Integer")},
    ("INFO", "END"): {("1", "Integer")},
    ("INFO", "SVTYPE"): {("1", "String")},
    ("INFO", "RU"): {("1", "String")},
}


class VcfHeaderError(RuntimeError):
    """An error code only; never include input lines, identifiers or file paths."""


@dataclass(frozen=True, slots=True, repr=False)
class HeaderDeclaration:
    category: str
    attributes: tuple[tuple[str, str], ...]

    @property
    def identifier(self) -> str:
        return dict(self.attributes)["ID"]


@dataclass(frozen=True, slots=True, repr=False)
class SchemaHeader:
    content: bytes
    declarations: tuple[HeaderDeclaration, ...]
    contigs: tuple[tuple[str, int], ...]
    assembly_label: str
    sha256: str
    classification: str = field(default="controlled", init=False)
    llm_policy: str = field(default="no-llm", init=False)
    retention: str = field(default="delete", init=False)


def _parse_header(content: bytes, *, sanitized: bool) -> SchemaHeader:
    if not isinstance(content, bytes) or len(content) > MAX_HEADER_BYTES:
        raise VcfHeaderError("header-size-invalid")
    try:
        text = content.decode("utf-8", errors="strict")
    except UnicodeDecodeError:
        raise VcfHeaderError("header-encoding-invalid") from None
    if not text.endswith("\n") or any(
        (ord(char) < 32 and char not in "\t\r\n") or ord(char) == 127
        for char in text
    ):
        raise VcfHeaderError("header-text-invalid")
    lines = text.removesuffix("\n").split("\n")
    lines = [line.removesuffix("\r") for line in lines]
    if any("\r" in line for line in lines):
        raise VcfHeaderError("header-text-invalid")
    if not re.fullmatch(r"##fileformat=VCFv4\.[1-5]", lines[0]):
        raise VcfHeaderError("header-fileformat-invalid")
    version = int(lines[0][-1])
    columns = lines[-1].split("\t")
    if (
        len(columns) != 10 or "\t".join(columns[:9]) != COLUMNS
        or not columns[9] or any(char.isspace() for char in columns[9])
        or (sanitized and columns[9] != "SYNTHETIC")
    ):
        raise VcfHeaderError("header-single-sample-required")

    declarations: list[HeaderDeclaration] = []
    contigs: list[tuple[str, int]] = []
    seen: set[tuple[str, str]] = set()
    assembly: str | None = None
    output = [lines[0]]
    for line in lines[1:-1]:
        if not line.startswith("##") or "\t" in line:
            raise VcfHeaderError("header-record-or-order-invalid")
        category, separator, value = line[2:].partition("=")
        if not separator or not value or not FIELD_ID.fullmatch(category):
            raise VcfHeaderError("header-metadata-invalid")
        if category == "fileformat":
            raise VcfHeaderError("header-fileformat-invalid")
        if category == "reference":
            if assembly is not None or (sanitized and not PUBLIC_ACCESSION.fullmatch(value)):
                raise VcfHeaderError("header-reference-invalid")
            assembly = value if sanitized else PUBLIC_ASSEMBLY_LABEL
            output.append(f"##reference={assembly}")
            continue
        if category not in ATTRIBUTES:
            if sanitized:
                raise VcfHeaderError("header-category-forbidden")
            continue
        if not value.startswith("<") or not value.endswith(">"):
            raise VcfHeaderError("header-declaration-invalid")
        inner = value[1:-1]
        attributes: dict[str, str] = {}
        offset = 0
        while offset < len(inner):
            match = ATTRIBUTE.match(inner, offset)
            if match is None or match[1] in attributes:
                raise VcfHeaderError("header-attributes-invalid")
            attributes[match[1]] = match[2]
            offset = match.end()
        required = ATTRIBUTES[category]
        if (
            not inner or inner.endswith(",") or not set(required) <= attributes.keys()
            or (sanitized and set(attributes) != set(required))
        ):
            raise VcfHeaderError("header-attributes-invalid")
        identifier = attributes["ID"]
        syntax = TOKEN
        if category == "contig" and version >= 3:
            syntax = CONTIG_ID
        elif category in {"INFO", "FORMAT"} and version >= 3:
            syntax = FIELD_ID
        special_id = (category, identifier) == ("INFO", "1000G")
        if (
            (not syntax.fullmatch(identifier) and not special_id)
            or (category != "ALT" and identifier == ".")
            or (category in {"FILTER", "INFO"} and ";" in identifier)
            or (category == "INFO" and "=" in identifier)
            or (category == "FORMAT" and (":" in identifier or "=" in identifier))
            or not identifier.isascii()
            or (category == "FILTER" and identifier == "0")
            or (category, identifier) in seen
        ):
            raise VcfHeaderError("header-identifier-invalid")
        seen.add((category, identifier))
        if category == "contig":
            if not re.fullmatch(r"[0-9]{1,19}", attributes["length"]) or int(attributes["length"]) <= 0:
                raise VcfHeaderError("header-contig-length-invalid")
            contigs.append((identifier, int(attributes["length"])))
        else:
            description = attributes["Description"]
            if not description.startswith('"') or not description.endswith('"'):
                raise VcfHeaderError("header-description-invalid")
            # Existing technical descriptions remain acceptable; arbitrary prose does not.
            neutral = {f'"{NEUTRAL_DESCRIPTION}"'}
            if (category, identifier) == ("FORMAT", "GT"):
                neutral.add('"Genotype"')
            if (category, identifier) == ("FILTER", "PASS"):
                neutral.add('"All filters passed"')
            if sanitized and description not in neutral:
                raise VcfHeaderError("header-description-invalid")
            attributes["Description"] = f'"{NEUTRAL_DESCRIPTION}"'
        if category in {"INFO", "FORMAT"}:
            number, kind = attributes["Number"], attributes["Type"]
            numbers = {"A", "G", "."}
            if version >= 2:
                numbers.add("R")
            if category == "FORMAT" and version >= 4:
                numbers.add("P")
            if category == "FORMAT" and version >= 5:
                numbers.update({"LA", "LR", "LG", "M"})
            semantic_number = str(int(number)) if re.fullmatch(r"[0-9]{1,19}", number) else number
            if (
                (number not in numbers and not re.fullmatch(r"[0-9]{1,19}", number))
                or kind not in {"Integer", "Float", "Flag", "Character", "String"}
                or (kind == "Flag" and (category != "INFO" or semantic_number != "0"))
                or (number == "M" and identifier[-1] not in "ACGTUN")
                or ((category, identifier) in FIELD_CONTRACTS
                    and (semantic_number, kind) not in FIELD_CONTRACTS[(category, identifier)])
            ):
                raise VcfHeaderError("header-field-contract-invalid")
        retained = tuple((key, attributes[key]) for key in required)
        declarations.append(HeaderDeclaration(category, retained))
        output.append(f"##{category}=<" + ",".join(f"{key}={item}" for key, item in retained) + ">")
    if assembly is None or not contigs or not any(item.category == "FORMAT" for item in declarations):
        raise VcfHeaderError("header-required-declaration-missing")
    output.append(COLUMNS + "\tSYNTHETIC")
    result = content if sanitized else ("\n".join(output) + "\n").encode("utf-8")
    return SchemaHeader(result, tuple(declarations), tuple(contigs), assembly, hashlib.sha256(result).hexdigest())


def validate_schema_header(content: bytes) -> SchemaHeader:
    """Check the entire artifact without rewriting bytes or granting LLM access."""
    return _parse_header(content, sanitized=True)


def sanitize_header(content: bytes) -> SchemaHeader:
    """Reject bodies/multiple samples before returning any controlled output."""
    return validate_schema_header(_parse_header(content, sanitized=False).content)


def read_schema_header(path: str | Path, *, sanitize: bool = False) -> SchemaHeader:
    try:
        content = read_bounded_bytes(path, maximum_bytes=MAX_HEADER_BYTES)
    except BoundedJSONError:
        raise VcfHeaderError("header-read-invalid") from None
    return sanitize_header(content) if sanitize else validate_schema_header(content)


def sanitize_header_file(
    input_path: str | Path, output_path: str | Path, *, controlled_root: str | Path
) -> SchemaHeader:
    """Publish a new private header only after full validation; never overwrite approval."""
    source, destination = Path(input_path), Path(output_path)
    try:
        root = Path(controlled_root).resolve(strict=True)
        if (
            not destination.parent.resolve().is_relative_to(root)
            or destination.resolve() == source.resolve()
            or destination.exists() or destination.is_symlink()
        ):
            raise VcfHeaderError("header-output-invalid")
    except (OSError, RuntimeError):
        raise VcfHeaderError("header-output-invalid") from None
    header = read_schema_header(source, sanitize=True)
    parent_fd = -1
    published = False
    temporary_created = False
    temporary = f".header-{uuid.uuid4().hex}.partial"
    try:
        destination.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
        parent_fd = os.open(destination.parent, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
        if not Path(f"/proc/self/fd/{parent_fd}").resolve(strict=True).is_relative_to(root):
            raise VcfHeaderError("header-output-invalid")
        descriptor = os.open(temporary, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600, dir_fd=parent_fd)
        temporary_created = True
        with os.fdopen(descriptor, "wb") as output:
            output.write(header.content)
            output.flush()
            os.fsync(output.fileno())
        os.link(temporary, destination.name, src_dir_fd=parent_fd, dst_dir_fd=parent_fd)
        published = True
        os.fsync(parent_fd)
    except (OSError, RuntimeError):
        if published:
            try:
                current = os.stat(destination.name, dir_fd=parent_fd, follow_symlinks=False)
                original = os.stat(temporary, dir_fd=parent_fd, follow_symlinks=False)
                if (current.st_dev, current.st_ino) == (original.st_dev, original.st_ino):
                    os.unlink(destination.name, dir_fd=parent_fd)
            except OSError:
                raise VcfHeaderError("header-cleanup-incomplete") from None
        raise VcfHeaderError("header-output-invalid") from None
    finally:
        if parent_fd >= 0:
            try:
                if temporary_created:
                    os.unlink(temporary, dir_fd=parent_fd)
            except OSError:
                raise VcfHeaderError("header-cleanup-incomplete") from None
            finally:
                os.close(parent_fd)
    return header


@dataclass(frozen=True, slots=True, repr=False)
class SkeletonCompatibility:
    header_sha256: str
    reference_manifest_sha256: str
    status: str
    missing_contigs: tuple[str, ...]
    extra_contigs: tuple[str, ...]
    length_mismatches: tuple[str, ...]
    schema: str = field(default="genosapiens.skeleton-compatibility", init=False)
    classification: str = field(default="controlled", init=False)
    llm_policy: str = field(default="no-llm", init=False)
    retention: str = field(default="delete", init=False)
    reference_accepted: bool = field(default=False, init=False)
    primary_rescue_authorized: bool = field(default=False, init=False)


def compare_skeleton_dictionary(
    header: SchemaHeader,
    expected_contigs: Mapping[str, int],
    reference_manifest_sha256: str,
) -> SkeletonCompatibility:
    """Compare complete dictionaries; order-only differences still block compatibility."""
    header = validate_schema_header(header.content)
    if (
        not expected_contigs or not re.fullmatch(r"[0-9a-f]{64}", reference_manifest_sha256)
        or any(not CONTIG_ID.fullmatch(name) or type(length) is not int or length <= 0
               for name, length in expected_contigs.items())
    ):
        raise VcfHeaderError("skeleton-reference-binding-invalid")
    actual = dict(header.contigs)
    missing = tuple(name for name in expected_contigs if name not in actual)
    extra = tuple(name for name in actual if name not in expected_contigs)
    lengths = tuple(name for name in actual if name in expected_contigs and actual[name] != expected_contigs[name])
    status = "incompatible" if missing or extra or lengths else (
        "compatible" if header.contigs == tuple(expected_contigs.items()) else "order-adjudication-required"
    )
    return SkeletonCompatibility(header.sha256, reference_manifest_sha256, status, missing, extra, lengths)
