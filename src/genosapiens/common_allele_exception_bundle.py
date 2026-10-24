"""Build the immutable public common-allele exception projection."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import tempfile
from typing import Any, Iterable, Mapping, Sequence

from genosapiens.admission_links import parse_clinvar_vcf_row
from genosapiens.common_allele_resources import (
    load_common_allele_exception_descriptor,
    verify_common_allele_source_snapshot,
)
from genosapiens.tabix import TabixIndex


CURATION_SCHEMA = "genosapiens.track1-common-allele-exception-curation"
RECORD_SCHEMA = "genosapiens.track1-common-allele-exception-record"
IDENTITY_SCHEMA = "genosapiens.track1-common-allele-exception-bundle-identity"
MANIFEST_SCHEMA = "genosapiens.track1-common-allele-exception-bundle-manifest"
SHA256 = re.compile(r"^[0-9a-f]{64}$")
CAID = re.compile(r"^CA[0-9]{6}$")
CHROMOSOMES = {str(value) for value in range(1, 23)} | {"X", "Y", "MT"}


class CommonAlleleBundleError(RuntimeError):
    """Raised when public exception qualification cannot close safely."""


def _fail() -> None:
    raise CommonAlleleBundleError("common-allele exception bundle qualification failed")


def canonical(value: object) -> bytes:
    try:
        return (json.dumps(value, sort_keys=True, separators=(",", ":"),
                           ensure_ascii=True, allow_nan=False) + "\n").encode("ascii")
    except (TypeError, ValueError, UnicodeError, OverflowError, RecursionError):
        _fail()


def _pairs(items: list[tuple[str, Any]]) -> dict[str, Any]:
    value: dict[str, Any] = {}
    for key, item in items:
        if key in value:
            _fail()
        value[key] = item
    return value


def _json(content: bytes, *, require_canonical: bool = False) -> dict[str, Any]:
    try:
        value = json.loads(content.decode("utf-8"), object_pairs_hook=_pairs,
                           parse_constant=lambda _value: _fail())
    except CommonAlleleBundleError:
        raise
    except (ValueError, TypeError, UnicodeError, RecursionError):
        _fail()
    if type(value) is not dict or (require_canonical and canonical(value) != content):
        _fail()
    return value


def _sha(value: object) -> str:
    if type(value) is not str or SHA256.fullmatch(value) is None:
        _fail()
    return value


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as stream:
            while chunk := stream.read(1024 * 1024):
                digest.update(chunk)
    except OSError:
        _fail()
    return digest.hexdigest()


def _regular_file(root: Path, relative: str, sha256: str, size: int | None = None) -> Path:
    path = root / relative
    try:
        if (not path.is_relative_to(root) or path.is_symlink() or not path.is_file()
                or (size is not None and path.stat().st_size != size)
                or _file_sha256(path) != sha256):
            _fail()
    except (OSError, ValueError):
        _fail()
    return path


@dataclass(frozen=True, slots=True)
class ExactAlleleMapping:
    gene: str
    source_table_row: int
    source_hgvs: str
    source_classification: str
    source_chromosome: str
    source_position: int
    source_reference: str
    source_alternate: str
    clingen_allele_registry_id: str
    clinvar_variation_id: str
    clinvar_allele_id: str
    chromosome: str
    position: int
    reference: str
    alternate: str


@dataclass(frozen=True, slots=True)
class ExceptionCuration:
    content: bytes
    sha256: str
    source_snapshot_sha256: str
    clinvar_manifest_path: str
    clinvar_manifest_sha256: str
    open_targets_manifest_path: str
    open_targets_manifest_sha256: str
    mappings: tuple[ExactAlleleMapping, ...]
    mitochondrial_document_id: str


@dataclass(frozen=True, slots=True)
class BuiltExceptionBundle:
    bundle_sha256: str
    bundle_path: str
    bundle_file_sha256: str
    manifest_path: str
    manifest_sha256: str
    record_count: int
    exact_allele_count: int
    gene_interval_count: int
    native_scope_count: int


@dataclass(frozen=True, slots=True)
class VerifiedExceptionBundle:
    bundle_sha256: str
    manifest_sha256: str
    record_sha256: str
    records: tuple[Mapping[str, Any], ...]


def _allele(value: object) -> tuple[str, int, str, str]:
    if type(value) is not dict or set(value) != {
        "chromosome", "position", "reference", "alternate"
    }:
        _fail()
    chrom, pos, ref, alt = (value[key] for key in (
        "chromosome", "position", "reference", "alternate"
    ))
    if (type(chrom) is not str or chrom not in CHROMOSOMES or type(pos) is not int
            or isinstance(pos, bool) or pos <= 0 or type(ref) is not str
            or type(alt) is not str or not ref or not alt
            or re.fullmatch(r"[ACGTN]+", ref) is None
            or re.fullmatch(r"[ACGTN]+", alt) is None or ref == alt):
        _fail()
    return chrom, pos, ref, alt


def load_exception_curation(content: bytes, expected_sha256: str) -> ExceptionCuration:
    try:
        if (type(content) is not bytes or _sha(expected_sha256) != hashlib.sha256(content).hexdigest()):
            _fail()
        value = _json(content, require_canonical=True)
        if set(value) != {
            "schema", "version", "assembly", "source_build", "source_snapshot_sha256",
            "clinvar_manifest", "open_targets_manifest", "exact_allele_mappings",
            "mitochondrial_resolution", "qualification",
        } or value["schema"] != CURATION_SCHEMA or value["version"] != 1 \
                or value["assembly"] != "GRCh38" or value["source_build"] != "GRCh37":
            _fail()
        qualification = value["qualification"]
        if qualification != {"ambiguous_mapping": "block", "challenge_or_truth_inputs_used": False,
                              "missing_mapping": "block", "source_build_mismatch": "block"}:
            _fail()
        def manifest(row: object) -> tuple[str, str]:
            if type(row) is not dict or set(row) != {"path", "sha256"}:
                _fail()
            path = row["path"]
            if (type(path) is not str or Path(path).is_absolute() or ".." in Path(path).parts
                    or not path.startswith("manifests/")):
                _fail()
            return path, _sha(row["sha256"])
        clinvar_manifest = manifest(value["clinvar_manifest"])
        open_targets_manifest = manifest(value["open_targets_manifest"])
        rows = value["exact_allele_mappings"]
        if type(rows) is not list or len(rows) != 9:
            _fail()
        mappings = []
        for index, row in enumerate(rows, 1):
            if type(row) is not dict or set(row) != {
                "gene", "source_table_row", "source_hgvs", "source_classification",
                "source_grch37", "clingen_allele_registry_id", "clinvar_variation_id",
                "expected_clinvar_allele_id", "expected_grch38",
            } or row["source_table_row"] != index:
                _fail()
            source = _allele(row["source_grch37"])
            target = _allele(row["expected_grch38"])
            if (type(row["gene"]) is not str or not row["gene"]
                    or type(row["source_hgvs"]) is not str or not row["source_hgvs"]
                    or row["source_classification"] not in {"Pathogenic", "VUS"}
                    or type(row["clinvar_variation_id"]) is not str
                    or not row["clinvar_variation_id"].isdecimal()
                    or type(row["expected_clinvar_allele_id"]) is not str
                    or not row["expected_clinvar_allele_id"].isdecimal()
                    or type(row["clingen_allele_registry_id"]) is not str
                    or CAID.fullmatch(row["clingen_allele_registry_id"]) is None):
                _fail()
            mappings.append(ExactAlleleMapping(
                row["gene"], index, row["source_hgvs"], row["source_classification"],
                *source, row["clingen_allele_registry_id"], row["clinvar_variation_id"],
                row["expected_clinvar_allele_id"], *target,
            ))
        if len({item.clingen_allele_registry_id for item in mappings}) != 9 \
                or len({(item.chromosome, item.position, item.reference, item.alternate)
                        for item in mappings}) != 9:
            _fail()
        mitochondrial = value["mitochondrial_resolution"]
        if mitochondrial != {
            "document_id": "GN015",
            "protection": "all_native_mitochondrial_occurrences_bypass_common_allele_exclusion",
            "reason": "released Mitochondrial Diseases VCEP BA1 specification has native mitochondrial haplogroup scope and no gene list",
        }:
            _fail()
        return ExceptionCuration(
            content, expected_sha256, _sha(value["source_snapshot_sha256"]),
            *clinvar_manifest, *open_targets_manifest, tuple(mappings), "GN015",
        )
    except CommonAlleleBundleError:
        raise
    except (KeyError, TypeError, ValueError, AttributeError, RecursionError):
        _fail()


def select_clinvar_mapping(
    rows: Sequence[Mapping[str, object]], mapping: ExactAlleleMapping
) -> Mapping[str, object]:
    """Require one exact GRCh38 row with both ClinVar and ClinGen identities."""

    matches = []
    for parsed in rows:
        columns, info = parsed.get("columns"), parsed.get("info")
        if not isinstance(columns, tuple) or len(columns) != 8 or not isinstance(info, Mapping):
            _fail()
        if (columns[0], columns[1], columns[2], columns[3], columns[4]) != (
            mapping.chromosome, str(mapping.position), mapping.clinvar_variation_id,
            mapping.reference, mapping.alternate,
        ):
            continue
        sources = info.get("CLNVI")
        if (info.get("ALLELEID") == mapping.clinvar_allele_id and isinstance(sources, str)
                and f"ClinGen:{mapping.clingen_allele_registry_id}" in sources.split("|") + sources.split(",")):
            matches.append(parsed)
    if len(matches) != 1:
        _fail()
    return matches[0]


def select_target_mapping(rows: Sequence[Sequence[object]], symbol: str) -> tuple[str, str, int, int, int]:
    """Require a unique exact-symbol Open Targets target span."""

    matches = []
    for row in rows:
        if len(row) != 6 or row[1] != symbol:
            continue
        target, _symbol, chrom, start, end, strand = row
        if (type(target) is not str or re.fullmatch(r"ENSG[0-9]{11}", target) is None
                or type(chrom) is not str or chrom not in CHROMOSOMES
                or type(start) is not int or type(end) is not int or type(strand) is not int
                or start <= 0 or end < start or strand not in {-1, 1}):
            _fail()
        matches.append((target, chrom, start, end, strand))
    if len(matches) != 1:
        _fail()
    return matches[0]


def _atomic(path: Path, content: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.is_symlink() or path.parent.is_symlink():
        _fail()
    descriptor, temporary = tempfile.mkstemp(prefix=".partial-", dir=path.parent)
    try:
        os.fchmod(descriptor, 0o444)
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
        try:
            os.link(temporary, path)
        except FileExistsError:
            if path.is_symlink() or not path.is_file() or path.read_bytes() != content:
                _fail()
        directory = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
    finally:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass


def _gzip(content: bytes) -> bytes:
    temporary = tempfile.SpooledTemporaryFile(max_size=16 * 1024 * 1024)
    with gzip.GzipFile(filename="", mode="wb", fileobj=temporary, mtime=0, compresslevel=9) as stream:
        stream.write(content)
    temporary.seek(0)
    return temporary.read()


def _manifest_outputs(reference_root: Path, path: str, sha256: str) -> dict[str, dict[str, Any]]:
    content = _regular_file(reference_root, path, sha256).read_bytes()
    value = _json(content)
    outputs = value.get("outputs")
    if (value.get("classification") != "public" or value.get("assembly") != "GRCh38"
            or type(outputs) is not list):
        _fail()
    result = {}
    for row in outputs:
        if type(row) is not dict or set(row) != {"path", "retention", "sha256", "size_bytes"}:
            _fail()
        relative = row["path"]
        if type(relative) is not str or relative in result:
            _fail()
        result[relative] = row
    return result


def _repository_bindings(repository: Path, descriptor_value: Mapping[str, object]) -> dict[str, str]:
    result = {}
    rows = descriptor_value.get("upstream_bindings")
    if type(rows) is not list:
        _fail()
    for row in rows:
        if type(row) is not dict:
            _fail()
        role, relative, sha256 = row.get("role"), row.get("path"), row.get("sha256")
        if type(role) is not str or type(relative) is not str or role in result:
            _fail()
        _regular_file(repository, relative, _sha(sha256))
        result[role] = sha256
    if set(result) != {"clinvar", "open_targets", "reference_identity", "selected_downstream_policy"}:
        _fail()
    return result


def _source_documents(snapshot_value: Mapping[str, object]) -> tuple[list[dict[str, Any]], dict[str, list[dict[str, Any]]]]:
    documents = snapshot_value.get("registry_documents")
    if type(documents) is not list:
        _fail()
    genes: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for document in documents:
        if type(document) is not dict:
            _fail()
        for gene in document.get("genes", []):
            if type(gene) is not dict or type(gene.get("symbol")) is not str:
                _fail()
            genes[gene["symbol"]].append({
                "document_id": document["svi_id"],
                "document_version": document["version"],
                "document_sha256": document["sha256"],
                "document_path": document["path"],
                "criteria": document["ba1_records"],
            })
    return documents, genes


def build_exception_bundle(
    *, repository_root: str | Path, reference_root: str | Path,
    descriptor_path: str | Path, curation_path: str | Path,
) -> BuiltExceptionBundle:
    repository = Path(repository_root).resolve(strict=True)
    references = Path(reference_root).resolve(strict=True)
    descriptor_file = Path(descriptor_path).resolve(strict=True)
    curation_file = Path(curation_path).resolve(strict=True)
    if (not descriptor_file.is_relative_to(repository) or not curation_file.is_relative_to(repository)
            or descriptor_file.is_symlink() or curation_file.is_symlink()):
        _fail()
    descriptor_content = descriptor_file.read_bytes()
    descriptor_sha = hashlib.sha256(descriptor_content).hexdigest()
    descriptor = load_common_allele_exception_descriptor(descriptor_content, descriptor_sha)
    snapshot = verify_common_allele_source_snapshot(descriptor, references)
    descriptor_value = _json(descriptor_content, require_canonical=True)
    bindings = _repository_bindings(repository, descriptor_value)
    curation_content = curation_file.read_bytes()
    curation = load_exception_curation(curation_content, hashlib.sha256(curation_content).hexdigest())
    if curation.source_snapshot_sha256 != snapshot.sha256:
        _fail()
    snapshot_value = _json(snapshot.content, require_canonical=True)
    documents, gene_sources = _source_documents(snapshot_value)
    if set(snapshot.missing_gene_document_ids) != {curation.mitochondrial_document_id}:
        _fail()

    clinvar_outputs = _manifest_outputs(
        references, curation.clinvar_manifest_path, curation.clinvar_manifest_sha256
    )
    clinvar_contract = _json((repository / descriptor_value["upstream_bindings"][0]["path"]).read_bytes())
    artifacts = {row["role"]: row for row in clinvar_contract.get("artifacts", [])}
    if set(artifacts) != {"vcf", "tbi"}:
        _fail()
    clinvar_paths = {}
    for role, row in artifacts.items():
        output = clinvar_outputs.get(row["path"])
        if (output is None or output["sha256"] != row["sha256"]
                or output["size_bytes"] != row["size_bytes"]):
            _fail()
        clinvar_paths[role] = _regular_file(references, row["path"], row["sha256"], row["size_bytes"])
    tabix = TabixIndex.load(clinvar_paths["tbi"])
    tabix.require_vcf_layout()

    records: list[dict[str, Any]] = []
    for mapping in curation.mappings:
        candidates = tuple(
            parse_clinvar_vcf_row(row) for row in tabix.query(
                clinvar_paths["vcf"], mapping.chromosome, mapping.position,
                mapping.position, max_records=100,
            )
        )
        selected = select_clinvar_mapping(candidates, mapping)
        info = selected["info"]
        records.append({
            "schema": RECORD_SCHEMA, "version": 1,
            "record_id": f"exact:{mapping.chromosome}:{mapping.position}:{mapping.reference}:{mapping.alternate}",
            "record_type": "exact_allele", "assembly": "GRCh38",
            "chromosome": mapping.chromosome, "position": mapping.position,
            "reference": mapping.reference, "alternate": mapping.alternate,
            "protection_reason": "clingen_ba1_exception_list",
            "curated_vus_exception": mapping.source_classification == "VUS",
            "source": {
                "source_table": "ClinGen SVI BA1 Exception List 2018-07-30",
                "source_table_row": mapping.source_table_row,
                "source_build": "GRCh37", "source_hgvs": mapping.source_hgvs,
                "source_classification": mapping.source_classification,
                "source_allele": {"chromosome": mapping.source_chromosome,
                                  "position": mapping.source_position,
                                  "reference": mapping.source_reference,
                                  "alternate": mapping.source_alternate},
                "clingen_allele_registry_id": mapping.clingen_allele_registry_id,
                "clinvar_variation_id": mapping.clinvar_variation_id,
                "clinvar_allele_id": mapping.clinvar_allele_id,
                "clinvar_clinical_significance": info.get("CLNSIG"),
                "clinvar_conflicting_significance": info.get("CLNSIGCONF"),
                "clinvar_review_status": info.get("CLNREVSTAT"),
            },
        })

    open_targets_outputs = _manifest_outputs(
        references, curation.open_targets_manifest_path, curation.open_targets_manifest_sha256
    )
    target_outputs = sorted(
        (row for path, row in open_targets_outputs.items()
         if path.startswith("releases/open-targets/26.06/output/target/")
         and path.endswith(".parquet")), key=lambda row: row["path"]
    )
    if len(target_outputs) != 3:
        _fail()
    target_paths = [
        _regular_file(references, row["path"], row["sha256"], row["size_bytes"])
        for row in target_outputs
    ]
    try:
        import duckdb
        connection = duckdb.connect(":memory:")
        connection.execute("SET threads=1")
        connection.execute("SET memory_limit='512MB'")
        connection.execute("SET preserve_insertion_order=false")
        symbols = sorted(gene_sources)
        placeholders = ",".join("?" for _ in symbols)
        rows = connection.execute(
            f"SELECT id, approvedSymbol, genomicLocation.chromosome, "
            f"genomicLocation.start, genomicLocation.end, genomicLocation.strand "
            f"FROM read_parquet(?, union_by_name=true) WHERE approvedSymbol IN ({placeholders}) "
            f"ORDER BY approvedSymbol, id",
            [[str(path) for path in target_paths], *symbols],
        ).fetchall()
        connection.close()
    except Exception as error:
        raise CommonAlleleBundleError("common-allele Open Targets projection failed") from error
    by_symbol: dict[str, list[Sequence[object]]] = defaultdict(list)
    for row in rows:
        if len(row) == 6 and type(row[1]) is str:
            by_symbol[row[1]].append(row)
    for symbol in symbols:
        target, chrom, start, end, strand = select_target_mapping(by_symbol.get(symbol, []), symbol)
        records.append({
            "schema": RECORD_SCHEMA, "version": 1,
            "record_id": f"gene:{target}", "record_type": "gene_interval",
            "assembly": "GRCh38", "chromosome": chrom, "start": start, "end": end,
            "coordinate_system": "one_based_inclusive", "strand": strand,
            "gene_symbol": symbol, "target_id": target,
            "protection_reason": "released_cspec_ba1_modification_or_not_applicable_rule",
            "source_documents": sorted(gene_sources[symbol], key=lambda row: row["document_id"]),
        })

    gn015 = [row for row in documents if row.get("svi_id") == curation.mitochondrial_document_id]
    if len(gn015) != 1 or gn015[0].get("gene_coverage") != "missing_declared" \
            or not gn015[0].get("ba1_records"):
        _fail()
    records.append({
        "schema": RECORD_SCHEMA, "version": 1, "record_id": "scope:native_mitochondrial",
        "record_type": "native_scope", "assembly": "GRCh38", "scope": "mitochondrial",
        "protection_reason": "released_mitochondrial_vcep_ba1_specification",
        "source_documents": [{
            "document_id": gn015[0]["svi_id"], "document_version": gn015[0]["version"],
            "document_sha256": gn015[0]["sha256"], "document_path": gn015[0]["path"],
            "criteria": gn015[0]["ba1_records"],
        }],
    })

    records.sort(key=lambda row: row["record_id"])
    if len({row["record_id"] for row in records}) != len(records):
        _fail()
    covered_documents = {
        source["document_id"] for row in records if row["record_type"] != "exact_allele"
        for source in row["source_documents"]
    }
    if covered_documents != set(snapshot.document_ids):
        _fail()
    record_content = b"".join(canonical(row) for row in records)
    record_sha = hashlib.sha256(record_content).hexdigest()
    counts = {
        "exact_alleles": sum(row["record_type"] == "exact_allele" for row in records),
        "gene_intervals": sum(row["record_type"] == "gene_interval" for row in records),
        "native_scopes": sum(row["record_type"] == "native_scope" for row in records),
        "records": len(records),
        "curated_vus_exact_alleles": sum(
            row.get("curated_vus_exception") is True for row in records
        ),
        "source_registry_documents": len(documents),
        "source_ba1_criteria": sum(len(row["ba1_records"]) for row in documents),
    }
    if counts != {"exact_alleles": 9, "gene_intervals": 134, "native_scopes": 1,
                  "records": 144, "curated_vus_exact_alleles": 5,
                  "source_registry_documents": 123, "source_ba1_criteria": 126}:
        _fail()
    identity = {
        "schema": IDENTITY_SCHEMA, "version": 1, "classification": "public",
        "assembly": "GRCh38", "source_descriptor_sha256": descriptor_sha,
        "curation_sha256": curation.sha256, "source_snapshot_sha256": snapshot.sha256,
        "upstream_contract_sha256s": bindings,
        "upstream_manifest_sha256s": {
            "clinvar": curation.clinvar_manifest_sha256,
            "open_targets": curation.open_targets_manifest_sha256,
        },
        "output_identity": {
            "format": "canonical-jsonl-gzip", "gzip_mtime": 0,
            "uncompressed_sha256": record_sha, "counts": counts,
        },
    }
    bundle_sha = hashlib.sha256(canonical(identity)).hexdigest()
    bundle_relative = f"derived/track1-common-allele-exceptions/{bundle_sha}/exceptions.jsonl.gz"
    manifest_relative = f"manifests/track1-common-allele-exceptions/{bundle_sha}.json"
    bundle_content = _gzip(record_content)
    bundle_file_sha = hashlib.sha256(bundle_content).hexdigest()
    manifest = {
        "schema": MANIFEST_SCHEMA, "version": 1, "classification": "public",
        "qualification_status": "qualified", "controlled_use_authorized": False,
        "challenge_or_truth_inputs_used": False, "bundle_sha256": bundle_sha,
        "identity": identity,
        "output": {"path": bundle_relative, "size_bytes": len(bundle_content),
                   "sha256": bundle_file_sha, "uncompressed_size_bytes": len(record_content),
                   "uncompressed_sha256": record_sha},
        "coverage": {"document_ids": list(snapshot.document_ids), **counts},
        "residual_scope": "public snapshots through 2026-09-19; no claim of exhaustive future clinical knowledge",
    }
    manifest_content = canonical(manifest)
    _atomic(references / bundle_relative, bundle_content)
    _atomic(references / manifest_relative, manifest_content)
    if gzip.decompress((references / bundle_relative).read_bytes()) != record_content \
            or (references / manifest_relative).read_bytes() != manifest_content:
        _fail()
    return BuiltExceptionBundle(
        bundle_sha, bundle_relative, bundle_file_sha, manifest_relative,
        hashlib.sha256(manifest_content).hexdigest(), counts["records"],
        counts["exact_alleles"], counts["gene_intervals"], counts["native_scopes"],
    )


def verify_exception_bundle(
    *, reference_root: str | Path, manifest_path: str,
    expected_manifest_sha256: str,
) -> VerifiedExceptionBundle:
    """Authenticate and reconcile a qualified canonical bundle from scratch."""

    references = Path(reference_root).resolve(strict=True)
    manifest_file = _regular_file(references, manifest_path, _sha(expected_manifest_sha256))
    content = manifest_file.read_bytes()
    manifest = _json(content, require_canonical=True)
    if set(manifest) != {
        "schema", "version", "classification", "qualification_status",
        "controlled_use_authorized", "challenge_or_truth_inputs_used", "bundle_sha256",
        "identity", "output", "coverage", "residual_scope",
    } or manifest["schema"] != MANIFEST_SCHEMA or manifest["version"] != 1 \
            or manifest["classification"] != "public" \
            or manifest["qualification_status"] != "qualified" \
            or manifest["controlled_use_authorized"] is not False \
            or manifest["challenge_or_truth_inputs_used"] is not False:
        _fail()
    bundle_sha = _sha(manifest["bundle_sha256"])
    identity = manifest["identity"]
    if type(identity) is not dict or identity.get("schema") != IDENTITY_SCHEMA \
            or hashlib.sha256(canonical(identity)).hexdigest() != bundle_sha:
        _fail()
    output = manifest["output"]
    expected_path = f"derived/track1-common-allele-exceptions/{bundle_sha}/exceptions.jsonl.gz"
    if (type(output) is not dict or set(output) != {
            "path", "size_bytes", "sha256", "uncompressed_size_bytes", "uncompressed_sha256"
        } or output["path"] != expected_path):
        _fail()
    bundle_file = _regular_file(
        references, expected_path, _sha(output["sha256"]), output["size_bytes"]
    )
    try:
        uncompressed = gzip.decompress(bundle_file.read_bytes())
    except (OSError, gzip.BadGzipFile):
        _fail()
    if (len(uncompressed) != output["uncompressed_size_bytes"]
            or hashlib.sha256(uncompressed).hexdigest() != _sha(output["uncompressed_sha256"])
            or identity.get("output_identity", {}).get("uncompressed_sha256")
            != output["uncompressed_sha256"] or not uncompressed.endswith(b"\n")):
        _fail()
    records = []
    for line in uncompressed.splitlines(keepends=True):
        row = _json(line, require_canonical=True)
        if row.get("schema") != RECORD_SCHEMA or row.get("version") != 1:
            _fail()
        records.append(row)
    ids = [row.get("record_id") for row in records]
    if (ids != sorted(ids) or len(ids) != len(set(ids))
            or not all(type(item) is str for item in ids)):
        _fail()
    counts = {
        "exact_alleles": sum(row.get("record_type") == "exact_allele" for row in records),
        "gene_intervals": sum(row.get("record_type") == "gene_interval" for row in records),
        "native_scopes": sum(row.get("record_type") == "native_scope" for row in records),
        "records": len(records),
        "curated_vus_exact_alleles": sum(row.get("curated_vus_exception") is True for row in records),
    }
    coverage = manifest["coverage"]
    identity_counts = identity.get("output_identity", {}).get("counts")
    if (type(coverage) is not dict or type(identity_counts) is not dict
            or any(coverage.get(key) != value or identity_counts.get(key) != value
                   for key, value in counts.items())
            or coverage.get("source_registry_documents") != 123
            or coverage.get("source_ba1_criteria") != 126
            or identity_counts.get("source_registry_documents") != 123
            or identity_counts.get("source_ba1_criteria") != 126):
        _fail()
    covered = {
        source.get("document_id")
        for row in records if row.get("record_type") in {"gene_interval", "native_scope"}
        for source in row.get("source_documents", []) if type(source) is dict
    }
    document_ids = coverage.get("document_ids")
    if type(document_ids) is not list or covered != set(document_ids) or len(covered) != 123:
        _fail()
    return VerifiedExceptionBundle(
        bundle_sha, expected_manifest_sha256, output["uncompressed_sha256"], tuple(records)
    )
