"""Bootstrap direct rare-disease sources and build provenance-aware projections."""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping, Sequence
import csv
from dataclasses import asdict, dataclass
from datetime import date, datetime
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import tempfile
from typing import Any
import uuid
import xml.etree.ElementTree as ET

from genosapiens.download import DownloadRequest, DownloadResult, download_atomic
from genosapiens.manifest import (
    ManifestError,
    PublicResourceManifest,
    ResourceOutput,
    RetentionPolicy,
    SourceIdentity,
    Transformation,
    publish_manifest,
)
from genosapiens.storage import DEFAULT_RESERVE_BYTES, CapacityPlan, require_capacity


SCHEMA = "genosapiens.rare-disease-evidence-contract"
SOURCE_TRANSFORMATION_SCHEMA = "genosapiens.rare-disease-source-bootstrap"
BUNDLE_MANIFEST_SCHEMA = "genosapiens.rare-disease-evidence-bundle"
DIGEST = re.compile(r"^sha256:[0-9a-f]{64}$")
SHA256 = re.compile(r"^[0-9a-f]{64}$")
IDENTIFIER = re.compile(r"^[a-z0-9][a-z0-9._-]*$")
ORPHA = re.compile(r"^ORPHA:[0-9]+$")
HPO = re.compile(r"^HP:[0-9]{7}$")
SGC = re.compile(r"^SGC-[0-9]+$")
PMID = re.compile(r"([0-9]+)\s*\[PMID\]", re.IGNORECASE)
EXPECTED_ROLES = {
    "orphadata_diseases",
    "orphadata_phenotypes",
    "orphadata_genes",
    "orphadata_natural_history",
    "orphadata_prevalence",
    "orpha_nomenclature",
    "gencc_submissions",
}
EXPECTED_TABLES = {
    "disease",
    "disease_alias",
    "disease_relation",
    "phenotype_assertion",
    "gene_assertion",
    "disease_onset",
    "disease_inheritance",
    "disease_prevalence",
}
EXPECTED_VALIDATIONS = {
    "hpo_noncanonical_source_ids",
    "open_targets_hpo_evidence",
    "open_targets_non_hpo_excluded",
    "open_targets_gene_evidence",
    "open_targets_gene_aggregates",
    "unresolved_gencc_genes",
    "unresolved_gencc_diseases",
    "unresolved_orpha_relation_endpoints",
    "preferred_phenotype_assertions",
    "preferred_gene_assertions",
    "shared_phenotype_lineages",
    "shared_gene_lineages",
    "phenotype_direct_dependency_matches",
    "phenotype_dependencies_with_multiple_mapped_diseases",
    "orphanet_direct_gene_shared_lineages",
    "integrated_preferred_over_direct_dependencies",
    "phenotype_lineages_with_multiple_independent",
    "gene_lineages_with_multiple_independent",
    "assertions_spanning_mapped_diseases",
}
HPOA_COLUMNS = (
    "database_id",
    "disease_name",
    "qualifier",
    "hpo_id",
    "reference",
    "evidence",
    "onset",
    "frequency",
    "sex",
    "modifier",
    "aspect",
    "biocuration",
)
OT_EVIDENCE_DATASETS = {
    "evidence_orphanet": "orphanet",
    "evidence_clingen": "clingen",
    "evidence_gene2phenotype": "gene2phenotype",
    "evidence_genomics_england": "genomics_england",
}


class RareDiseaseError(RuntimeError):
    """Raised when rare-disease acquisition or projection validation fails."""


def _required_string(data: Mapping[str, Any], name: str, context: str) -> str:
    value = data.get(name)
    if not isinstance(value, str) or not value:
        raise RareDiseaseError(f"{context}.{name} must be a non-empty string")
    return value


def _positive_int(data: Mapping[str, Any], name: str, context: str) -> int:
    value = data.get(name)
    if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
        raise RareDiseaseError(f"{context}.{name} must be a positive integer")
    return value


def _safe_path(value: str, context: str) -> PurePosixPath:
    path = PurePosixPath(value)
    if path.is_absolute() or path == PurePosixPath(".") or ".." in path.parts:
        raise RareDiseaseError(f"{context} must be a safe relative path")
    return path


def _string_list(data: Mapping[str, Any], name: str, context: str) -> tuple[str, ...]:
    value = data.get(name)
    if (
        not isinstance(value, list)
        or not value
        or any(not isinstance(item, str) or not item for item in value)
        or len(value) != len(set(value))
    ):
        raise RareDiseaseError(f"{context}.{name} must be a unique string array")
    return tuple(value)


@dataclass(frozen=True, slots=True)
class RareDiseaseArtifact:
    role: str
    resource_id: str
    source_family: str
    release: str
    extraction_date: str
    source_revision: str
    url: str
    path: str
    size_bytes: int
    sha256: str
    license: str
    license_url: str
    expected_records: int

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "RareDiseaseArtifact":
        expected = {
            "role",
            "resource_id",
            "source_family",
            "release",
            "extraction_date",
            "source_revision",
            "url",
            "path",
            "size_bytes",
            "sha256",
            "license",
            "license_url",
            "expected_records",
        }
        if set(data) != expected:
            raise RareDiseaseError("rare-disease artifact fields are invalid")
        artifact = cls(
            role=_required_string(data, "role", "artifact"),
            resource_id=_required_string(data, "resource_id", "artifact"),
            source_family=_required_string(data, "source_family", "artifact"),
            release=_required_string(data, "release", "artifact"),
            extraction_date=_required_string(data, "extraction_date", "artifact"),
            source_revision=_required_string(data, "source_revision", "artifact"),
            url=_required_string(data, "url", "artifact"),
            path=_required_string(data, "path", "artifact"),
            size_bytes=_positive_int(data, "size_bytes", "artifact"),
            sha256=_required_string(data, "sha256", "artifact"),
            license=_required_string(data, "license", "artifact"),
            license_url=_required_string(data, "license_url", "artifact"),
            expected_records=_positive_int(data, "expected_records", "artifact"),
        )
        try:
            datetime.fromisoformat(artifact.extraction_date)
        except ValueError as error:
            raise RareDiseaseError("artifact extraction_date must be ISO-like") from error
        path = _safe_path(artifact.path, "artifact.path")
        if (
            not IDENTIFIER.fullmatch(artifact.resource_id)
            or not artifact.url.startswith("https://")
            or not artifact.license_url.startswith("https://")
            or not SHA256.fullmatch(artifact.sha256)
            or path.parts[0] != "releases"
        ):
            raise RareDiseaseError(f"invalid rare-disease artifact: {artifact.role}")
        return artifact


@dataclass(frozen=True, slots=True)
class HpoInput:
    release: str
    ontology_path: str
    ontology_sha256: str
    annotations_path: str
    annotations_sha256: str
    expected_annotation_records: int

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "HpoInput":
        expected = {
            "release",
            "ontology_path",
            "ontology_sha256",
            "annotations_path",
            "annotations_sha256",
            "expected_annotation_records",
        }
        if set(data) != expected:
            raise RareDiseaseError("HPO input fields are invalid")
        value = cls(
            release=_required_string(data, "release", "hpo"),
            ontology_path=_required_string(data, "ontology_path", "hpo"),
            ontology_sha256=_required_string(data, "ontology_sha256", "hpo"),
            annotations_path=_required_string(data, "annotations_path", "hpo"),
            annotations_sha256=_required_string(data, "annotations_sha256", "hpo"),
            expected_annotation_records=_positive_int(
                data, "expected_annotation_records", "hpo"
            ),
        )
        for path_value in (value.ontology_path, value.annotations_path):
            if _safe_path(path_value, "HPO input path").parts[0] != "releases":
                raise RareDiseaseError("HPO inputs must be release artifacts")
        if not SHA256.fullmatch(value.ontology_sha256) or not SHA256.fullmatch(
            value.annotations_sha256
        ):
            raise RareDiseaseError("HPO input checksum is invalid")
        return value


@dataclass(frozen=True, slots=True)
class OpenTargetsInput:
    release: str
    manifest_path: str
    croissant_sha256: str
    dataset_rows: Mapping[str, int]

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "OpenTargetsInput":
        expected = {"release", "manifest_path", "croissant_sha256", "dataset_rows"}
        rows = data.get("dataset_rows")
        required = set(OT_EVIDENCE_DATASETS) | {
            "association_by_datasource_direct",
            "disease",
            "disease_phenotype",
            "target",
        }
        if (
            set(data) != expected
            or not isinstance(rows, Mapping)
            or set(rows) != required
            or any(
                not isinstance(value, int) or isinstance(value, bool) or value <= 0
                for value in rows.values()
            )
        ):
            raise RareDiseaseError("Open Targets input fields are invalid")
        value = cls(
            release=_required_string(data, "release", "open_targets"),
            manifest_path=_required_string(data, "manifest_path", "open_targets"),
            croissant_sha256=_required_string(
                data, "croissant_sha256", "open_targets"
            ),
            dataset_rows=dict(rows),
        )
        if (
            _safe_path(value.manifest_path, "Open Targets manifest path").parts[0]
            != "manifests"
            or not SHA256.fullmatch(value.croissant_sha256)
        ):
            raise RareDiseaseError("Open Targets input identity is invalid")
        return value


@dataclass(frozen=True, slots=True)
class ProjectionContract:
    schema_version: str
    output_prefix: str
    expected_row_counts: Mapping[str, int]

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "ProjectionContract":
        expected = {"schema_version", "output_prefix", "expected_row_counts"}
        counts = data.get("expected_row_counts")
        if (
            set(data) != expected
            or not isinstance(counts, Mapping)
            or set(counts) != EXPECTED_TABLES
            or any(
                not isinstance(value, int) or isinstance(value, bool) or value <= 0
                for value in counts.values()
            )
        ):
            raise RareDiseaseError("rare-disease projection fields are invalid")
        value = cls(
            schema_version=_required_string(data, "schema_version", "projection"),
            output_prefix=_required_string(data, "output_prefix", "projection"),
            expected_row_counts=dict(counts),
        )
        if _safe_path(value.output_prefix, "projection.output_prefix").parts[:2] != (
            "derived",
            "rare-disease-evidence",
        ):
            raise RareDiseaseError("projection output must be derived rare-disease data")
        return value


@dataclass(frozen=True, slots=True)
class RareDiseaseContract:
    bundle_id: str
    projection_id: str
    assembly: str
    projected_output_bytes: int
    artifacts: tuple[RareDiseaseArtifact, ...]
    gencc_columns: tuple[str, ...]
    hpo: HpoInput
    open_targets: OpenTargetsInput
    validation_expectations: Mapping[str, int]
    projection: ProjectionContract
    schema: str = SCHEMA

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "RareDiseaseContract":
        expected = {
            "schema",
            "bundle_id",
            "projection_id",
            "assembly",
            "projected_output_bytes",
            "artifacts",
            "gencc_columns",
            "hpo",
            "open_targets",
            "validation_expectations",
            "projection",
        }
        artifacts = data.get("artifacts")
        validation_expectations = data.get("validation_expectations")
        if (
            set(data) != expected
            or data.get("schema") != SCHEMA
            or not isinstance(artifacts, list)
            or not all(isinstance(item, Mapping) for item in artifacts)
            or not isinstance(data.get("hpo"), Mapping)
            or not isinstance(data.get("open_targets"), Mapping)
            or not isinstance(validation_expectations, Mapping)
            or set(validation_expectations) != EXPECTED_VALIDATIONS
            or any(
                not isinstance(value, int) or isinstance(value, bool) or value < 0
                for value in validation_expectations.values()
            )
            or not isinstance(data.get("projection"), Mapping)
        ):
            raise RareDiseaseError("rare-disease contract fields or schema are invalid")
        contract = cls(
            schema=SCHEMA,
            bundle_id=_required_string(data, "bundle_id", "contract"),
            projection_id=_required_string(data, "projection_id", "contract"),
            assembly=_required_string(data, "assembly", "contract"),
            projected_output_bytes=_positive_int(
                data, "projected_output_bytes", "contract"
            ),
            artifacts=tuple(RareDiseaseArtifact.from_mapping(item) for item in artifacts),
            gencc_columns=_string_list(data, "gencc_columns", "contract"),
            hpo=HpoInput.from_mapping(data["hpo"]),
            open_targets=OpenTargetsInput.from_mapping(data["open_targets"]),
            validation_expectations=dict(validation_expectations),
            projection=ProjectionContract.from_mapping(data["projection"]),
        )
        roles = [item.role for item in contract.artifacts]
        resource_ids = [item.resource_id for item in contract.artifacts]
        paths = [item.path for item in contract.artifacts]
        if (
            not IDENTIFIER.fullmatch(contract.bundle_id)
            or contract.assembly != "not-applicable"
            or set(roles) != EXPECTED_ROLES
            or len(roles) != len(set(roles))
            or len(resource_ids) != len(set(resource_ids))
            or len(paths) != len(set(paths))
            or PurePosixPath(contract.projection.output_prefix).name
            != contract.projection_id
        ):
            raise RareDiseaseError("rare-disease contract identities are inconsistent")
        return contract

    def artifact(self, role: str) -> RareDiseaseArtifact:
        for artifact in self.artifacts:
            if artifact.role == role:
                return artifact
        raise RareDiseaseError(f"rare-disease contract lacks artifact role {role}")


def load_rare_disease_contract(path: str | Path) -> RareDiseaseContract:
    source = Path(path)
    try:
        value = json.loads(source.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise RareDiseaseError(f"cannot read rare-disease contract {source}: {error}") from error
    if not isinstance(value, Mapping):
        raise RareDiseaseError("rare-disease contract root must be an object")
    return RareDiseaseContract.from_mapping(value)


@dataclass(frozen=True, slots=True)
class RareDiseaseValidation:
    row_counts: Mapping[str, int]
    preferred_phenotype_assertions: int
    preferred_gene_assertions: int
    shared_phenotype_lineages: int
    shared_gene_lineages: int
    hpo_noncanonical_source_ids: int
    open_targets_hpo_evidence: int
    open_targets_non_hpo_excluded: int
    unresolved_gencc_genes: int
    unresolved_gencc_diseases: int
    unresolved_orpha_relation_endpoints: int


@dataclass(frozen=True, slots=True)
class RareDiseaseBootstrapResult:
    bundle_id: str
    projection_id: str
    bundle_manifest_path: str
    bundle_manifest_sha256: str
    source_manifest_paths: tuple[str, ...]
    source_manifest_sha256: tuple[str, ...]
    validation: RareDiseaseValidation
    reused: bool


Downloader = Callable[[DownloadRequest, str | Path], DownloadResult]
CapacityChecker = Callable[[str | Path, CapacityPlan], Any]


def _file_identity(path: Path) -> tuple[str, int]:
    digest = hashlib.sha256()
    size = 0
    try:
        with path.open("rb") as source:
            while chunk := source.read(1024 * 1024):
                digest.update(chunk)
                size += len(chunk)
    except OSError as error:
        raise RareDiseaseError(f"cannot read rare-disease artifact {path}: {error}") from error
    return digest.hexdigest(), size


def _verify_file(path: Path, sha256: str, size_bytes: int | None = None) -> None:
    if path.is_symlink() or not path.is_file():
        raise RareDiseaseError(f"rare-disease input is not a regular file: {path}")
    observed_sha256, observed_size = _file_identity(path)
    if observed_sha256 != sha256 or (
        size_bytes is not None and observed_size != size_bytes
    ):
        raise RareDiseaseError(f"rare-disease input checksum mismatch: {path}")


def _require_within_root(root: Path, path: Path, context: str) -> None:
    if not path.is_absolute() or not path.is_relative_to(root):
        raise RareDiseaseError(f"{context} escapes reference_root: {path}")
    current = root
    for part in path.relative_to(root).parts:
        current /= part
        if current.is_symlink():
            raise RareDiseaseError(f"{context} contains a symlink: {current}")


def _canonical_json(value: Mapping[str, Any]) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def _publish_json(value: Mapping[str, Any], path: Path, root: Path) -> str:
    _require_within_root(root, path, "rare-disease bundle manifest")
    content = _canonical_json(value)
    digest = hashlib.sha256(content).hexdigest()
    if path.exists() or path.is_symlink():
        if path.is_symlink() or not path.is_file():
            raise RareDiseaseError("rare-disease bundle manifest is not a regular file")
        if path.read_bytes() != content:
            raise RareDiseaseError("completed rare-disease bundle manifest is immutable")
        return digest
    path.parent.mkdir(parents=True, exist_ok=True)
    partial_root = path.parent / ".partial"
    _require_within_root(root, partial_root, "rare-disease manifest partial path")
    partial_root.mkdir(mode=0o700, exist_ok=True)
    partial = partial_root / f"{path.name}.{uuid.uuid4().hex}.partial"
    try:
        with partial.open("xb") as output:
            os.chmod(partial, 0o600)
            output.write(content)
            output.flush()
            os.fsync(output.fileno())
        os.chmod(partial, 0o444)
        try:
            os.link(partial, path)
        except FileExistsError:
            if path.is_symlink() or not path.is_file() or path.read_bytes() != content:
                raise RareDiseaseError("completed rare-disease bundle manifest is immutable")
        directory_fd = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    finally:
        partial.unlink(missing_ok=True)
        try:
            partial_root.rmdir()
        except OSError:
            pass
    return digest


def _text(element: ET.Element, path: str) -> str | None:
    value = element.findtext(path)
    if value is None:
        return None
    stripped = value.strip()
    return stripped or None


def _controlled_value(
    element: ET.Element, path: str
) -> tuple[str | None, str | None]:
    target = element.find(path)
    if target is None:
        return None, None
    return target.get("id"), _text(target, "Name")


def _normalise_curie(value: str | None) -> str | None:
    if value is None:
        return None
    result = value.replace("\u00a0", "").strip(" \t\n\r\v\f")
    if not result:
        return None
    match = re.fullmatch(r"(MONDO|EFO|HP|SO|Orphanet)_([A-Za-z0-9.]+)", result)
    if match:
        prefix = "ORPHA" if match.group(1) == "Orphanet" else match.group(1)
        return f"{prefix}:{match.group(2)}"
    if result.startswith("Orphanet:"):
        return "ORPHA:" + result.partition(":")[2]
    return result


def _register_normalise_curie(connection: Any) -> None:
    connection.execute(
        """
        CREATE MACRO clean_curie(value) AS CASE
          WHEN value IS NULL THEN NULL
          ELSE nullif(regexp_replace(replace(value, chr(160), ''),
            '^[[:space:]]+|[[:space:]]+$', '', 'g'), '')
        END;
        CREATE MACRO norm_curie(value) AS CASE
          WHEN clean_curie(value) IS NULL THEN NULL
          WHEN regexp_full_match(clean_curie(value), 'MONDO_[A-Za-z0-9.]+')
            THEN 'MONDO:' || substr(clean_curie(value), 7)
          WHEN regexp_full_match(clean_curie(value), 'EFO_[A-Za-z0-9.]+')
            THEN 'EFO:' || substr(clean_curie(value), 5)
          WHEN regexp_full_match(clean_curie(value), 'HP_[A-Za-z0-9.]+')
            THEN 'HP:' || substr(clean_curie(value), 4)
          WHEN regexp_full_match(clean_curie(value), 'SO_[A-Za-z0-9.]+')
            THEN 'SO:' || substr(clean_curie(value), 4)
          WHEN regexp_full_match(clean_curie(value), 'Orphanet_[A-Za-z0-9.]+')
            THEN 'ORPHA:' || substr(clean_curie(value), 10)
          WHEN starts_with(clean_curie(value), 'Orphanet:')
            THEN 'ORPHA:' || substr(clean_curie(value), 10)
          ELSE clean_curie(value)
        END
        """
    )


def _lineage_key(*values: str | None) -> str:
    content = "|".join(value or "" for value in values).encode()
    return hashlib.sha256(content).hexdigest()


def _list_key(values: Sequence[str]) -> str:
    return json.dumps(sorted(set(values)), separators=(",", ":"))


def _dependency_key(lineage: str, *values: str | None) -> str:
    return _lineage_key("dependency-v2", lineage, *values)


def _assertion_key(dependency: str, mapped_disease_id: str | None) -> str:
    return _lineage_key("assertion-v2", dependency, mapped_disease_id)


def _pmids(value: str | None) -> list[str]:
    if not value:
        return []
    seen: set[str] = set()
    result = []
    for match in PMID.finditer(value):
        pmid = match.group(1)
        if pmid not in seen:
            seen.add(pmid)
            result.append(pmid)
    return result


def _split_field(value: str | None) -> list[str]:
    if not value:
        return []
    return [item.strip() for item in value.split(";") if item.strip()]


def _specificity(classification: str | None) -> str:
    if classification == "Subtype of disorder":
        return "subtype"
    if classification == "Group of disorders":
        return "broad"
    if classification == "Disorder":
        return "disorder"
    return "unresolved"


def _load_xml(path: Path, artifact: RareDiseaseArtifact, list_tag: str) -> ET.Element:
    try:
        root = ET.parse(path).getroot()
    except (OSError, ET.ParseError) as error:
        raise RareDiseaseError(f"cannot parse {artifact.role}: {error}") from error
    extraction = root.get("ExtractionDate") or root.get("date")
    license_id = _text(root, "Availability/Licence/ShortIdentifier")
    if extraction != artifact.extraction_date or license_id != artifact.license:
        raise RareDiseaseError(f"{artifact.role} release or license differs from contract")
    for element in root.iter():
        declared = element.get("count")
        if declared is None:
            continue
        try:
            expected = int(declared)
        except ValueError as error:
            raise RareDiseaseError(f"{artifact.role} has a non-integer count") from error
        if expected != len(element):
            raise RareDiseaseError(
                f"{artifact.role} count mismatch for {element.tag}: "
                f"expected {expected}, observed {len(element)}"
            )
    record_list = root.find(list_tag)
    if record_list is None or int(record_list.get("count", "-1")) != artifact.expected_records:
        raise RareDiseaseError(f"{artifact.role} primary record count differs from contract")
    return root


@dataclass(frozen=True, slots=True)
class HpoTerms:
    names: Mapping[str, str]
    aliases: Mapping[str, str]
    obsolete: Mapping[str, tuple[str, ...]]

    def resolve(self, term_id: str) -> tuple[str, str, bool]:
        if term_id in self.names:
            return term_id, self.names[term_id], False
        alias = self.aliases.get(term_id)
        if alias is not None:
            return alias, self.names[alias], False
        replacements = self.obsolete.get(term_id, ())
        resolved = [
            self.aliases.get(candidate, candidate)
            for candidate in replacements
            if self.aliases.get(candidate, candidate) in self.names
        ]
        if not resolved:
            raise RareDiseaseError(f"Orphadata HPO term does not resolve: {term_id}")
        canonical = sorted(set(resolved))[0]
        return canonical, self.names[canonical], True


def _load_hpo_terms(path: Path) -> HpoTerms:
    names: dict[str, str] = {}
    aliases: dict[str, str] = {}
    obsolete: dict[str, tuple[str, ...]] = {}
    current: dict[str, list[str]] | None = None

    def finish() -> None:
        nonlocal current
        if current is None:
            return
        ids = current.get("id", [])
        term_names = current.get("name", [])
        if len(ids) != 1 or len(term_names) != 1:
            raise RareDiseaseError("HPO term lacks one ID and name")
        term_id = ids[0]
        if current.get("is_obsolete") == ["true"]:
            replacements = tuple(
                item.split(" !", 1)[0].strip()
                for item in current.get("replaced_by", []) + current.get("consider", [])
            )
            obsolete[term_id] = replacements
        else:
            names[term_id] = term_names[0]
            for alias in current.get("alt_id", []):
                if alias in aliases and aliases[alias] != term_id:
                    raise RareDiseaseError(f"HPO alias is ambiguous: {alias}")
                aliases[alias] = term_id
        current = None

    try:
        with path.open("rt", encoding="utf-8") as source:
            for raw_line in source:
                line = raw_line.rstrip("\r\n")
                if line == "[Term]":
                    finish()
                    current = {}
                    continue
                if line.startswith("[") and line.endswith("]"):
                    finish()
                    current = None
                    continue
                if current is None or ": " not in line:
                    continue
                key, value = line.split(": ", 1)
                current.setdefault(key, []).append(value)
            finish()
    except (OSError, UnicodeDecodeError) as error:
        raise RareDiseaseError(f"cannot parse HPO ontology: {error}") from error
    if not names or "HP:0000001" not in names:
        raise RareDiseaseError("HPO ontology is empty or lacks its root")
    return HpoTerms(names, aliases, obsolete)


def _insert_rows(
    connection: Any, table: str, rows: Iterable[Sequence[Any]], width: int
) -> int:
    count = 0
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="",
            prefix=f"{table}-",
            suffix=".csv",
            delete=False,
        ) as temporary:
            temporary_path = Path(temporary.name)
            writer = csv.writer(temporary, lineterminator="\n")
            for row in rows:
                if len(row) != width:
                    raise RareDiseaseError(f"internal {table} row width is invalid")
                writer.writerow([r"\N" if value is None else value for value in row])
                count += 1
        if count:
            escaped = str(temporary_path).replace("'", "''")
            connection.execute(
                f"COPY {table} FROM '{escaped}' "
                "(FORMAT CSV, HEADER false, NULL '\\N', STRICT_MODE true)"
            )
        return count
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


TABLE_SCHEMAS = {
    "disease": """
        disease_id VARCHAR, internal_id VARCHAR, name VARCHAR,
        disorder_type_id VARCHAR, disorder_type VARCHAR,
        classification_level_id VARCHAR, classification_level VARCHAR,
        status VARCHAR, flag_value BIGINT, aggregation_status VARCHAR,
        aggregation_target_id VARCHAR, aggregation_target_name VARCHAR,
        expert_link VARCHAR, synonyms VARCHAR[], definition VARCHAR,
        source_dataset VARCHAR, source_release VARCHAR
    """,
    "disease_alias": """
        source_record_id VARCHAR, disease_id VARCHAR, namespace VARCHAR,
        reference VARCHAR, normalised_reference VARCHAR,
        mapping_relation_id VARCHAR, mapping_relation VARCHAR,
        mapping_validation_status VARCHAR, icd_mapping_relation VARCHAR,
        reference_url VARCHAR, source_dataset VARCHAR, source_release VARCHAR
    """,
    "disease_relation": """
        root_disease_id VARCHAR, target_disease_id VARCHAR,
        relation_type_id VARCHAR, relation_type VARCHAR,
        root_internal_id VARCHAR, target_internal_id VARCHAR,
        unresolved_endpoint BOOLEAN, source_release VARCHAR
    """,
    "phenotype_assertion": """
        source_family VARCHAR, source_dataset VARCHAR, source_release VARCHAR,
        source_record_id VARCHAR, lineage_key VARCHAR, dependency_key VARCHAR,
        assertion_key VARCHAR,
        representation_kind VARCHAR, association_scope VARCHAR,
        is_ontology_propagated BOOLEAN, representation_priority INTEGER,
        preferred_representation BOOLEAN, counts_as_independent BOOLEAN,
        source_disease_id VARCHAR, mapped_disease_id VARCHAR,
        disease_specificity VARCHAR, source_hpo_id VARCHAR,
        source_hpo_label VARCHAR, canonical_hpo_id VARCHAR,
        canonical_hpo_label VARCHAR, frequency_id VARCHAR,
        frequency_name VARCHAR, qualifier VARCHAR, qualifier_not BOOLEAN,
        aspect VARCHAR, onset VARCHAR[], sex VARCHAR, modifiers VARCHAR[],
        diagnostic_status_id VARCHAR, diagnostic_status VARCHAR,
        evidence_type VARCHAR, "references" VARCHAR[], biocuration VARCHAR,
        provenance_raw VARCHAR, validation_status VARCHAR, validation_date VARCHAR
    """,
    "gene_assertion": """
        source_family VARCHAR, source_dataset VARCHAR, source_release VARCHAR,
        source_record_id VARCHAR, lineage_key VARCHAR, dependency_key VARCHAR,
        assertion_key VARCHAR,
        representation_kind VARCHAR, association_scope VARCHAR,
        is_ontology_propagated BOOLEAN, representation_priority INTEGER,
        preferred_representation BOOLEAN, counts_as_independent BOOLEAN,
        source_disease_id VARCHAR, mapped_disease_id VARCHAR,
        source_disease_name VARCHAR, mapped_disease_name VARCHAR,
        disease_specificity VARCHAR, gene_id VARCHAR, ensembl_gene_id VARCHAR,
        hgnc_id VARCHAR, gene_symbol VARCHAR, internal_gene_id VARCHAR,
        gene_type_id VARCHAR, gene_type VARCHAR, association_type_id VARCHAR,
        association_type VARCHAR, validity_id VARCHAR, validity VARCHAR,
        inheritance_id VARCHAR, inheritance VARCHAR,
        allelic_requirements VARCHAR[], consequence_requirement VARCHAR,
        confidence VARCHAR, evidence_date VARCHAR, publication_date VARCHAR,
        source_run_date VARCHAR, publications VARCHAR[], study_id VARCHAR, score DOUBLE,
        evidence_count BIGINT, public_report_url VARCHAR,
        assertion_criteria_url VARCHAR, notes VARCHAR, provenance_raw VARCHAR
    """,
    "disease_onset": """
        disease_id VARCHAR, onset_id VARCHAR, onset VARCHAR,
        disease_specificity VARCHAR, source_record_id VARCHAR,
        source_dataset VARCHAR, source_release VARCHAR, provenance_raw VARCHAR
    """,
    "disease_inheritance": """
        disease_id VARCHAR, inheritance_id VARCHAR, inheritance VARCHAR,
        disease_specificity VARCHAR, source_record_id VARCHAR,
        source_dataset VARCHAR, source_release VARCHAR, provenance_raw VARCHAR
    """,
    "disease_prevalence": """
        source_record_id VARCHAR, disease_id VARCHAR, disease_specificity VARCHAR,
        prevalence_type_id VARCHAR, prevalence_type VARCHAR,
        qualification_id VARCHAR, qualification VARCHAR,
        prevalence_class_id VARCHAR, prevalence_class VARCHAR,
        mean_value DOUBLE, geography_id VARCHAR, geography VARCHAR,
        validation_status_id VARCHAR, validation_status VARCHAR,
        provenance_raw VARCHAR, source_dataset VARCHAR, source_release VARCHAR
    """,
}
TABLE_COLUMNS = {
    "disease": (
        "disease_id",
        "internal_id",
        "name",
        "disorder_type_id",
        "disorder_type",
        "classification_level_id",
        "classification_level",
        "status",
        "flag_value",
        "aggregation_status",
        "aggregation_target_id",
        "aggregation_target_name",
        "expert_link",
        "synonyms",
        "definition",
        "source_dataset",
        "source_release",
    ),
    "disease_alias": (
        "source_record_id",
        "disease_id",
        "namespace",
        "reference",
        "normalised_reference",
        "mapping_relation_id",
        "mapping_relation",
        "mapping_validation_status",
        "icd_mapping_relation",
        "reference_url",
        "source_dataset",
        "source_release",
    ),
    "disease_relation": (
        "root_disease_id",
        "target_disease_id",
        "relation_type_id",
        "relation_type",
        "root_internal_id",
        "target_internal_id",
        "unresolved_endpoint",
        "source_release",
    ),
    "phenotype_assertion": (
        "source_family",
        "source_dataset",
        "source_release",
        "source_record_id",
        "lineage_key",
        "dependency_key",
        "assertion_key",
        "representation_kind",
        "association_scope",
        "is_ontology_propagated",
        "representation_priority",
        "preferred_representation",
        "counts_as_independent",
        "source_disease_id",
        "mapped_disease_id",
        "disease_specificity",
        "source_hpo_id",
        "source_hpo_label",
        "canonical_hpo_id",
        "canonical_hpo_label",
        "frequency_id",
        "frequency_name",
        "qualifier",
        "qualifier_not",
        "aspect",
        "onset",
        "sex",
        "modifiers",
        "diagnostic_status_id",
        "diagnostic_status",
        "evidence_type",
        "references",
        "biocuration",
        "provenance_raw",
        "validation_status",
        "validation_date",
    ),
    "gene_assertion": (
        "source_family",
        "source_dataset",
        "source_release",
        "source_record_id",
        "lineage_key",
        "dependency_key",
        "assertion_key",
        "representation_kind",
        "association_scope",
        "is_ontology_propagated",
        "representation_priority",
        "preferred_representation",
        "counts_as_independent",
        "source_disease_id",
        "mapped_disease_id",
        "source_disease_name",
        "mapped_disease_name",
        "disease_specificity",
        "gene_id",
        "ensembl_gene_id",
        "hgnc_id",
        "gene_symbol",
        "internal_gene_id",
        "gene_type_id",
        "gene_type",
        "association_type_id",
        "association_type",
        "validity_id",
        "validity",
        "inheritance_id",
        "inheritance",
        "allelic_requirements",
        "consequence_requirement",
        "confidence",
        "evidence_date",
        "publication_date",
        "source_run_date",
        "publications",
        "study_id",
        "score",
        "evidence_count",
        "public_report_url",
        "assertion_criteria_url",
        "notes",
        "provenance_raw",
    ),
    "disease_onset": (
        "disease_id",
        "onset_id",
        "onset",
        "disease_specificity",
        "source_record_id",
        "source_dataset",
        "source_release",
        "provenance_raw",
    ),
    "disease_inheritance": (
        "disease_id",
        "inheritance_id",
        "inheritance",
        "disease_specificity",
        "source_record_id",
        "source_dataset",
        "source_release",
        "provenance_raw",
    ),
    "disease_prevalence": (
        "source_record_id",
        "disease_id",
        "disease_specificity",
        "prevalence_type_id",
        "prevalence_type",
        "qualification_id",
        "qualification",
        "prevalence_class_id",
        "prevalence_class",
        "mean_value",
        "geography_id",
        "geography",
        "validation_status_id",
        "validation_status",
        "provenance_raw",
        "source_dataset",
        "source_release",
    ),
}
TABLE_SORTS = {
    "disease": "disease_id",
    "disease_alias": "disease_id, namespace, reference, source_record_id",
    "disease_relation": "root_disease_id, target_disease_id, relation_type_id",
    "phenotype_assertion": (
        "lineage_key, representation_priority, source_dataset, source_record_id, "
        "dependency_key, assertion_key"
    ),
    "gene_assertion": (
        "lineage_key, representation_priority, source_dataset, source_record_id, "
        "dependency_key, assertion_key"
    ),
    "disease_onset": "disease_id, onset_id",
    "disease_inheritance": "disease_id, inheritance_id",
    "disease_prevalence": "disease_id, source_record_id",
}


def _create_tables(connection: Any) -> None:
    for table, schema in TABLE_SCHEMAS.items():
        connection.execute(f"CREATE TABLE {table} ({schema})")


def _parse_nomenclature(
    connection: Any,
    path: Path,
    artifact: RareDiseaseArtifact,
) -> tuple[dict[str, str], int]:
    root = _load_xml(path, artifact, "DisorderList")
    disorders = root.find("DisorderList")
    assert disorders is not None
    specificity: dict[str, str] = {}
    internal_to_orpha: dict[str, str] = {}
    disease_rows = []
    raw_relations: set[tuple[str, str, str, str, str, str]] = set()
    for disorder in disorders:
        code = _text(disorder, "OrphaCode")
        if code is None:
            raise RareDiseaseError("nomenclature disorder lacks an ORPHA code")
        disease_id = f"ORPHA:{code}"
        internal_id = disorder.get("id")
        if not internal_id or disease_id in specificity:
            raise RareDiseaseError("nomenclature disease identity is duplicated")
        internal_to_orpha[internal_id] = disease_id
        type_id, type_name = _controlled_value(disorder, "DisorderType")
        level_id, level_name = _controlled_value(disorder, "ClassificationLevel")
        specificity[disease_id] = _specificity(level_name)
        section = disorder.find("AggregationLevelSection")
        aggregation = (
            section.find("AggregationLevelList/AggregationLevel")
            if section is not None
            else None
        )
        aggregation_status = (
            _text(aggregation, "AggregationLevelStatus")
            if aggregation is not None
            else _text(section, "AggregationLevelStatus") if section is not None else None
        )
        aggregation_code = _text(aggregation, "OrphaCode") if aggregation is not None else None
        aggregation_name = (
            _text(aggregation, "PreferredTerm") if aggregation is not None else None
        )
        definition = None
        for text_section in disorder.findall(
            "SummaryInformationList/SummaryInformation/TextSectionList/TextSection"
        ):
            if _text(text_section, "TextSectionType/Name") == "Definition":
                definition = _text(text_section, "Contents")
                break
        disease_rows.append(
            (
                disease_id,
                internal_id,
                _text(disorder, "Name"),
                type_id,
                type_name,
                level_id,
                level_name,
                _text(disorder, "Totalstatus"),
                int(_text(disorder, "FlagValue") or "0"),
                aggregation_status or "missing",
                f"ORPHA:{aggregation_code}" if aggregation_code else None,
                aggregation_name,
                _text(disorder, "ExpertLink"),
                [
                    synonym.text.strip()
                    for synonym in disorder.findall("SynonymList/Synonym")
                    if synonym.text and synonym.text.strip()
                ],
                definition,
                "orpha_nomenclature",
                artifact.release,
            )
        )
        for association in disorder.findall(
            "DisorderDisorderAssociationList/DisorderDisorderAssociation"
        ):
            root_disorder = association.find("RootDisorder")
            target_disorder = association.find("TargetDisorder")
            if root_disorder is None or target_disorder is None:
                raise RareDiseaseError("nomenclature relation lacks an endpoint")
            root_code = _text(root_disorder, "OrphaCode")
            target_code = _text(target_disorder, "OrphaCode")
            if root_disorder.get("cycle") == "true":
                root_code = code
            if target_disorder.get("cycle") == "true":
                target_code = code
            if root_code is None or target_code is None:
                raise RareDiseaseError("nomenclature relation endpoint lacks an ORPHA code")
            relation_id, relation_name = _controlled_value(
                association, "DisorderDisorderAssociationType"
            )
            if relation_id is None or relation_name is None:
                raise RareDiseaseError("nomenclature relation lacks a type")
            raw_relations.add(
                (
                    f"ORPHA:{root_code}",
                    f"ORPHA:{target_code}",
                    relation_id,
                    relation_name,
                    root_disorder.get("id") or internal_id,
                    target_disorder.get("id") or internal_id,
                )
            )
    _insert_rows(connection, "disease", disease_rows, len(TABLE_COLUMNS["disease"]))
    relation_rows = [
        (*row, row[0] not in specificity or row[1] not in specificity, artifact.release)
        for row in sorted(raw_relations)
    ]
    _insert_rows(
        connection,
        "disease_relation",
        relation_rows,
        len(TABLE_COLUMNS["disease_relation"]),
    )
    return specificity, sum(row[6] for row in relation_rows)


def _normalise_reference(namespace: str, reference: str) -> str:
    prefixes = {
        "MONDO": "MONDO",
        "OMIM": "OMIM",
        "UMLS": "UMLS",
        "MeSH": "MESH",
        "MedDRA": "MEDDRA",
        "GARD": "GARD",
        "ICD-10": "ICD10",
        "ICD-11": "ICD11",
    }
    prefix = prefixes.get(namespace, namespace.replace(" ", "_"))
    return f"{prefix}:{reference}"


def _parse_disease_aliases(
    connection: Any,
    path: Path,
    artifact: RareDiseaseArtifact,
    disease_ids: set[str],
) -> int:
    root = _load_xml(path, artifact, "DisorderList")
    disorders = root.find("DisorderList")
    assert disorders is not None
    observed: set[str] = set()
    rows = []
    for disorder in disorders:
        disease_id = f"ORPHA:{_text(disorder, 'OrphaCode')}"
        observed.add(disease_id)
        for reference in disorder.findall("ExternalReferenceList/ExternalReference"):
            namespace = _text(reference, "Source")
            value = _text(reference, "Reference")
            if namespace is None or value is None or not reference.get("id"):
                raise RareDiseaseError("Orphadata external reference lacks identity")
            relation_id, relation = _controlled_value(
                reference, "DisorderMappingRelation"
            )
            rows.append(
                (
                    reference.get("id"),
                    disease_id,
                    namespace,
                    value,
                    _normalise_reference(namespace, value),
                    relation_id,
                    relation,
                    _text(reference, "DisorderMappingValidationStatus/Name"),
                    _text(reference, "DisorderMappingICDRelation/Name"),
                    _text(reference, "DisorderMappingICDRefUrl"),
                    "orphadata_product1",
                    artifact.release,
                )
            )
    if observed != disease_ids:
        raise RareDiseaseError("Orphadata disease and nomenclature identifiers differ")
    return _insert_rows(
        connection, "disease_alias", rows, len(TABLE_COLUMNS["disease_alias"])
    )


FREQUENCY_HPO = {
    "28405": "HP:0040280",
    "28412": "HP:0040281",
    "28419": "HP:0040282",
    "28426": "HP:0040283",
    "28433": "HP:0040284",
    "28440": "HP:0040285",
}


def _parse_orphadata_phenotypes(
    connection: Any,
    path: Path,
    artifact: RareDiseaseArtifact,
    disease_specificity: Mapping[str, str],
    hpo: HpoTerms,
) -> tuple[int, set[str]]:
    root = _load_xml(path, artifact, "HPODisorderSetStatusList")
    sets = root.find("HPODisorderSetStatusList")
    assert sets is not None
    rows = []
    noncanonical_ids: set[str] = set()
    for status_set in sets:
        disorder = status_set.find("Disorder")
        if disorder is None:
            raise RareDiseaseError("Orphadata phenotype set lacks a disease")
        disease_id = f"ORPHA:{_text(disorder, 'OrphaCode')}"
        if disease_id not in disease_specificity:
            raise RareDiseaseError(f"Orphadata phenotype disease is unknown: {disease_id}")
        raw_source = _text(status_set, "Source")
        references = _pmids(raw_source)
        validation_status = _text(status_set, "ValidationStatus")
        validation_date = _text(status_set, "ValidationDate")
        for association in disorder.findall(
            "HPODisorderAssociationList/HPODisorderAssociation"
        ):
            source_hpo_id = _text(association, "HPO/HPOId")
            source_hpo_label = _text(association, "HPO/HPOTerm")
            if source_hpo_id is None or not HPO.fullmatch(source_hpo_id):
                raise RareDiseaseError("Orphadata phenotype assertion lacks an HPO ID")
            canonical_hpo_id, canonical_hpo_label, _was_obsolete = hpo.resolve(
                source_hpo_id
            )
            if source_hpo_id != canonical_hpo_id:
                noncanonical_ids.add(source_hpo_id)
            frequency_id, frequency_name = _controlled_value(
                association, "HPOFrequency"
            )
            normalised_frequency = FREQUENCY_HPO.get(frequency_id or "", frequency_id)
            qualifier_not = frequency_id == "28440"
            qualifier = "NOT" if qualifier_not else ""
            diagnostic_id, diagnostic = _controlled_value(
                association, "DiagnosticCriteria"
            )
            lineage = _lineage_key(
                "phenotype", "orphanet", disease_id, canonical_hpo_id
            )
            dependency = _dependency_key(
                lineage,
                normalised_frequency,
                qualifier,
                "P",
                _list_key([]),
                None,
                _list_key([]),
                diagnostic_id,
                "phenotype",
                _list_key(references),
            )
            assertion = _assertion_key(dependency, disease_id)
            provenance = json.dumps(
                {
                    "assertion_xml": ET.tostring(association, encoding="unicode"),
                    "source": raw_source,
                    "validation_date": validation_date,
                    "validation_status": validation_status,
                },
                sort_keys=True,
                separators=(",", ":"),
            )
            rows.append(
                (
                    "orphanet",
                    "orphadata_product4",
                    artifact.release,
                    association.get("id"),
                    lineage,
                    dependency,
                    assertion,
                    "source_direct",
                    "source_direct",
                    False,
                    0,
                    False,
                    False,
                    disease_id,
                    disease_id,
                    disease_specificity[disease_id],
                    source_hpo_id,
                    source_hpo_label,
                    canonical_hpo_id,
                    canonical_hpo_label,
                    normalised_frequency,
                    frequency_name,
                    qualifier,
                    qualifier_not,
                    "P",
                    [],
                    None,
                    [],
                    diagnostic_id,
                    diagnostic,
                    "phenotype",
                    references,
                    None,
                    provenance,
                    validation_status,
                    validation_date,
                )
            )
    count = _insert_rows(
        connection,
        "phenotype_assertion",
        rows,
        len(TABLE_COLUMNS["phenotype_assertion"]),
    )
    return count, noncanonical_ids


def _gene_external_references(gene: ET.Element) -> dict[str, str]:
    references: dict[str, str] = {}
    for item in gene.findall("ExternalReferenceList/ExternalReference"):
        source = _text(item, "Source")
        value = _text(item, "Reference")
        if source and value and source not in references:
            references[source] = value
    return references


def _parse_orphadata_genes(
    connection: Any,
    path: Path,
    artifact: RareDiseaseArtifact,
    disease_specificity: Mapping[str, str],
) -> tuple[int, dict[str, str]]:
    root = _load_xml(path, artifact, "DisorderList")
    disorders = root.find("DisorderList")
    assert disorders is not None
    rows = []
    hgnc_to_ensembl: dict[str, str] = {}
    source_keys: set[tuple[str, str, str]] = set()
    for disorder in disorders:
        disease_id = f"ORPHA:{_text(disorder, 'OrphaCode')}"
        if disease_id not in disease_specificity:
            raise RareDiseaseError(f"Orphadata gene disease is unknown: {disease_id}")
        disease_name = _text(disorder, "Name")
        for association in disorder.findall(
            "DisorderGeneAssociationList/DisorderGeneAssociation"
        ):
            gene = association.find("Gene")
            if gene is None or not gene.get("id"):
                raise RareDiseaseError("Orphadata gene assertion lacks a gene")
            references = _gene_external_references(gene)
            ensembl = references.get("Ensembl")
            hgnc = f"HGNC:{references['HGNC']}" if references.get("HGNC") else None
            if hgnc and ensembl:
                previous = hgnc_to_ensembl.setdefault(hgnc, ensembl)
                if previous != ensembl:
                    raise RareDiseaseError(f"HGNC maps to multiple Ensembl genes: {hgnc}")
            gene_id = ensembl or hgnc or f"ORPHA-GENE:{gene.get('id')}"
            association_id, association_name = _controlled_value(
                association, "DisorderGeneAssociationType"
            )
            status_id, status = _controlled_value(
                association, "DisorderGeneAssociationStatus"
            )
            if association_id is None or association_name is None:
                raise RareDiseaseError("Orphadata gene assertion lacks an association type")
            source_key = (disease_id, gene.get("id") or "", association_id)
            if source_key in source_keys:
                raise RareDiseaseError("Orphadata gene assertion key is duplicated")
            source_keys.add(source_key)
            raw_source = _text(association, "SourceOfValidation")
            publications = _pmids(raw_source)
            provenance = json.dumps(
                {
                    "assertion_xml": ET.tostring(association, encoding="unicode"),
                    "source_of_validation": raw_source,
                },
                sort_keys=True,
                separators=(",", ":"),
            )
            gene_type_id, gene_type = _controlled_value(gene, "GeneType")
            record_id = "|".join(source_key)
            lineage = _lineage_key("gene", "orphanet", disease_id, gene_id)
            dependency = _dependency_key(
                lineage,
                association_id,
                status_id,
                None,
                _list_key([]),
                None,
                status,
                _list_key(publications),
                None,
            )
            assertion = _assertion_key(dependency, disease_id)
            rows.append(
                (
                    "orphanet",
                    "orphadata_product6",
                    artifact.release,
                    record_id,
                    lineage,
                    dependency,
                    assertion,
                    "source_direct",
                    "source_direct",
                    False,
                    0,
                    False,
                    False,
                    disease_id,
                    disease_id,
                    disease_name,
                    disease_name,
                    disease_specificity[disease_id],
                    gene_id,
                    ensembl,
                    hgnc,
                    _text(gene, "Symbol"),
                    gene.get("id"),
                    gene_type_id,
                    gene_type,
                    association_id,
                    association_name,
                    status_id,
                    status,
                    None,
                    None,
                    [],
                    None,
                    status,
                    None,
                    None,
                    None,
                    publications,
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                    provenance,
                )
            )
    count = _insert_rows(
        connection, "gene_assertion", rows, len(TABLE_COLUMNS["gene_assertion"])
    )
    return count, hgnc_to_ensembl


def _parse_natural_history(
    connection: Any,
    path: Path,
    artifact: RareDiseaseArtifact,
    disease_specificity: Mapping[str, str],
) -> tuple[int, int]:
    root = _load_xml(path, artifact, "DisorderList")
    disorders = root.find("DisorderList")
    assert disorders is not None
    onset_rows = []
    inheritance_rows = []
    onset_keys: set[tuple[str, str]] = set()
    inheritance_keys: set[tuple[str, str]] = set()
    for disorder in disorders:
        disease_id = f"ORPHA:{_text(disorder, 'OrphaCode')}"
        if disease_id not in disease_specificity:
            raise RareDiseaseError(f"natural-history disease is unknown: {disease_id}")
        for onset in disorder.findall("AverageAgeOfOnsetList/AverageAgeOfOnset"):
            onset_id = onset.get("id")
            onset_name = _text(onset, "Name")
            if not onset_id or not onset_name or (disease_id, onset_id) in onset_keys:
                raise RareDiseaseError("natural-history onset identity is invalid")
            onset_keys.add((disease_id, onset_id))
            onset_rows.append(
                (
                    disease_id,
                    onset_id,
                    onset_name,
                    disease_specificity[disease_id],
                    f"{disease_id}|{onset_id}",
                    "orphadata_product9_ages",
                    artifact.release,
                    onset_name,
                )
            )
        for inheritance in disorder.findall(
            "TypeOfInheritanceList/TypeOfInheritance"
        ):
            inheritance_id = inheritance.get("id")
            inheritance_name = _text(inheritance, "Name")
            if (
                not inheritance_id
                or not inheritance_name
                or (disease_id, inheritance_id) in inheritance_keys
            ):
                raise RareDiseaseError("natural-history inheritance identity is invalid")
            inheritance_keys.add((disease_id, inheritance_id))
            inheritance_rows.append(
                (
                    disease_id,
                    inheritance_id,
                    inheritance_name,
                    disease_specificity[disease_id],
                    f"{disease_id}|{inheritance_id}",
                    "orphadata_product9_ages",
                    artifact.release,
                    inheritance_name,
                )
            )
    onset_count = _insert_rows(
        connection, "disease_onset", onset_rows, len(TABLE_COLUMNS["disease_onset"])
    )
    inheritance_count = _insert_rows(
        connection,
        "disease_inheritance",
        inheritance_rows,
        len(TABLE_COLUMNS["disease_inheritance"]),
    )
    return onset_count, inheritance_count


def _parse_prevalence(
    connection: Any,
    path: Path,
    artifact: RareDiseaseArtifact,
    disease_specificity: Mapping[str, str],
) -> int:
    root = _load_xml(path, artifact, "DisorderList")
    disorders = root.find("DisorderList")
    assert disorders is not None
    rows = []
    record_ids: set[str] = set()
    for disorder in disorders:
        disease_id = f"ORPHA:{_text(disorder, 'OrphaCode')}"
        if disease_id not in disease_specificity:
            raise RareDiseaseError(f"prevalence disease is unknown: {disease_id}")
        for prevalence in disorder.findall("PrevalenceList/Prevalence"):
            record_id = prevalence.get("id")
            if not record_id or record_id in record_ids:
                raise RareDiseaseError("Orphadata prevalence ID is invalid or duplicated")
            record_ids.add(record_id)
            type_id, type_name = _controlled_value(prevalence, "PrevalenceType")
            qualification_id, qualification = _controlled_value(
                prevalence, "PrevalenceQualification"
            )
            class_id, class_name = _controlled_value(prevalence, "PrevalenceClass")
            geography_id, geography = _controlled_value(
                prevalence, "PrevalenceGeographic"
            )
            validation_id, validation = _controlled_value(
                prevalence, "PrevalenceValidationStatus"
            )
            raw_value = _text(prevalence, "ValMoy")
            raw_source = _text(prevalence, "Source")
            try:
                mean_value = float(raw_value or "")
            except ValueError as error:
                raise RareDiseaseError("Orphadata prevalence value is invalid") from error
            rows.append(
                (
                    record_id,
                    disease_id,
                    disease_specificity[disease_id],
                    type_id,
                    type_name,
                    qualification_id,
                    qualification,
                    class_id,
                    class_name,
                    mean_value,
                    geography_id,
                    geography,
                    validation_id,
                    validation,
                    raw_source,
                    "orphadata_product9_prevalence",
                    artifact.release,
                )
            )
    return _insert_rows(
        connection,
        "disease_prevalence",
        rows,
        len(TABLE_COLUMNS["disease_prevalence"]),
    )


def _parse_hpo_annotations(
    connection: Any,
    path: Path,
    release: str,
    expected_records: int,
    disease_specificity: Mapping[str, str],
    hpo: HpoTerms,
) -> int:
    try:
        source = path.open("rt", encoding="utf-8", newline="")
    except OSError as error:
        raise RareDiseaseError(f"cannot open HPO annotations: {error}") from error
    rows = []
    try:
        reader = csv.DictReader(
            (line for line in source if not line.startswith("#")), delimiter="\t"
        )
        if tuple(reader.fieldnames or ()) != HPOA_COLUMNS:
            raise RareDiseaseError("HPO annotation columns differ from the pinned schema")
        for index, row in enumerate(reader, start=1):
            if None in row or any(value is None for value in row.values()):
                raise RareDiseaseError("HPO annotation row width is invalid")
            source_disease = _normalise_curie(row["database_id"])
            source_hpo = row["hpo_id"]
            if source_disease is None or not HPO.fullmatch(source_hpo):
                raise RareDiseaseError("HPO annotation identity is invalid")
            canonical_hpo, canonical_name, _obsolete = hpo.resolve(source_hpo)
            source_family = "orphanet" if source_disease.startswith("ORPHA:") else "hpo"
            frequency = _normalise_curie(row["frequency"]) or row["frequency"] or None
            onset = _split_field(row["onset"])
            modifiers = _split_field(row["modifier"])
            references = _split_field(row["reference"])
            lineage = _lineage_key(
                "phenotype", source_family, source_disease, canonical_hpo
            )
            dependency = _dependency_key(
                lineage,
                frequency,
                row["qualifier"],
                row["aspect"],
                _list_key(onset),
                row["sex"] or None,
                _list_key(modifiers),
                None,
                row["evidence"],
                _list_key(references),
            )
            assertion = _assertion_key(dependency, source_disease)
            rows.append(
                (
                    source_family,
                    "hpo_disease_annotations",
                    release,
                    f"hpoa:{index}",
                    lineage,
                    dependency,
                    assertion,
                    "source_direct",
                    "source_direct",
                    False,
                    1,
                    False,
                    False,
                    source_disease,
                    source_disease,
                    disease_specificity.get(source_disease, "unresolved"),
                    source_hpo,
                    hpo.names.get(source_hpo),
                    canonical_hpo,
                    canonical_name,
                    frequency,
                    None,
                    row["qualifier"],
                    row["qualifier"] == "NOT",
                    row["aspect"],
                    onset,
                    row["sex"] or None,
                    modifiers,
                    None,
                    None,
                    row["evidence"],
                    references,
                    row["biocuration"],
                    json.dumps(row, sort_keys=True, separators=(",", ":")),
                    None,
                    None,
                )
            )
    except (csv.Error, OSError, UnicodeDecodeError) as error:
        raise RareDiseaseError(f"cannot parse HPO annotations: {error}") from error
    finally:
        source.close()
    if len(rows) != expected_records:
        raise RareDiseaseError(
            f"HPO annotation count differs: expected {expected_records}, observed {len(rows)}"
        )
    return _insert_rows(
        connection,
        "phenotype_assertion",
        rows,
        len(TABLE_COLUMNS["phenotype_assertion"]),
    )


def _gencc_source_family(submitter_title: str, submitter_curie: str) -> str:
    known = {
        "Orphanet": "orphanet",
        "ClinGen": "clingen",
        "G2P": "gene2phenotype",
        "Genomics England PanelApp": "genomics_england",
    }
    return known.get(submitter_title, f"gencc:{submitter_curie.lower()}")


def _gencc_publications(value: str | None) -> list[str]:
    if not value:
        return []
    return list(dict.fromkeys(re.findall(r"[0-9]+", value)))


def _parse_gencc(
    connection: Any,
    path: Path,
    artifact: RareDiseaseArtifact,
    columns: Sequence[str],
    disease_specificity: Mapping[str, str],
    hgnc_to_ensembl: Mapping[str, str],
    ot_disease_specificity: Mapping[str, str],
) -> tuple[int, int, int]:
    try:
        source = path.open("rt", encoding="utf-8-sig", newline="")
    except OSError as error:
        raise RareDiseaseError(f"cannot open GenCC submissions: {error}") from error
    rows = []
    sgc_ids: set[str] = set()
    unresolved_genes: set[str] = set()
    unresolved_diseases: set[str] = set()
    try:
        reader = csv.DictReader(source)
        if tuple(reader.fieldnames or ()) != tuple(columns):
            raise RareDiseaseError("GenCC columns differ from the pinned schema")
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise RareDiseaseError("GenCC row width is invalid")
            sgc_id = row["sgc_id"].strip()
            if not SGC.fullmatch(sgc_id) or sgc_id in sgc_ids:
                raise RareDiseaseError("GenCC SGC identity is invalid or duplicated")
            sgc_ids.add(sgc_id)
            hgnc = _normalise_curie(row["gene_curie"])
            mapped_disease = _normalise_curie(row["disease_curie"])
            source_disease = _normalise_curie(row["disease_original_curie"])
            if (
                hgnc is None
                or mapped_disease is None
                or source_disease is None
                or not row["classification_curie"].strip()
                or not row["classification_title"].strip()
                or not row["moi_curie"].strip()
                or not row["moi_title"].strip()
            ):
                raise RareDiseaseError("GenCC assertion lacks required validity or inheritance")
            ensembl = hgnc_to_ensembl.get(hgnc)
            if ensembl is None:
                unresolved_genes.add(hgnc)
            ot_disease_id = mapped_disease.replace(":", "_", 1)
            if ot_disease_id not in ot_disease_specificity:
                unresolved_diseases.add(mapped_disease)
            gene_id = ensembl or hgnc
            source_family = _gencc_source_family(
                row["submitter_title"].strip(), row["submitter_curie"].strip()
            )
            lineage = _lineage_key("gene", source_family, source_disease, gene_id)
            publications = _gencc_publications(row["submitted_as_pmids"])
            dependency = _dependency_key(
                lineage,
                "gene_disease_validity",
                row["classification_curie"].strip(),
                row["moi_curie"].strip(),
                _list_key([row["moi_title"].strip()]),
                None,
                row["classification_title"].strip(),
                _list_key(publications),
                row["submitted_as_submission_id"].strip() or None,
            )
            assertion = _assertion_key(dependency, mapped_disease)
            submitted = dict(row)
            rows.append(
                (
                    source_family,
                    "gencc_submissions",
                    artifact.release,
                    sgc_id,
                    lineage,
                    dependency,
                    assertion,
                    "source_submission",
                    "source_direct",
                    False,
                    1 if source_family == "orphanet" else 0,
                    False,
                    False,
                    source_disease,
                    mapped_disease,
                    row["disease_original_title"].strip() or None,
                    row["disease_title"].strip() or None,
                    disease_specificity.get(
                        source_disease,
                        ot_disease_specificity.get(ot_disease_id, "unresolved"),
                    ),
                    gene_id,
                    ensembl,
                    hgnc,
                    row["gene_symbol"].strip(),
                    None,
                    None,
                    None,
                    "gene_disease_validity",
                    "Gene-disease validity assertion",
                    row["classification_curie"].strip(),
                    row["classification_title"].strip(),
                    row["moi_curie"].strip(),
                    row["moi_title"].strip(),
                    [row["moi_title"].strip()],
                    None,
                    row["classification_title"].strip(),
                    row["submitted_as_date"].strip() or None,
                    None,
                    row["submitted_run_date"].strip() or None,
                    publications,
                    row["submitted_as_submission_id"].strip() or None,
                    None,
                    None,
                    row["submitted_as_public_report_url"].strip() or None,
                    row["submitted_as_assertion_criteria_url"].strip() or None,
                    row["submitted_as_notes"],
                    json.dumps(submitted, sort_keys=True, separators=(",", ":")),
                )
            )
    except (csv.Error, OSError, UnicodeDecodeError) as error:
        raise RareDiseaseError(f"cannot parse GenCC submissions: {error}") from error
    finally:
        source.close()
    if len(rows) != artifact.expected_records:
        raise RareDiseaseError(
            f"GenCC row count differs: expected {artifact.expected_records}, observed {len(rows)}"
        )
    count = _insert_rows(
        connection, "gene_assertion", rows, len(TABLE_COLUMNS["gene_assertion"])
    )
    return count, len(unresolved_genes), len(unresolved_diseases)


def _load_public_manifest(path: Path, context: str) -> PublicResourceManifest:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise RareDiseaseError(f"cannot read {context} manifest: {error}") from error
    if not isinstance(value, Mapping):
        raise RareDiseaseError(f"{context} manifest root must be an object")
    try:
        return PublicResourceManifest.from_mapping(value)
    except ManifestError as error:
        raise RareDiseaseError(f"invalid {context} manifest: {error}") from error


def _open_targets_files(
    root: Path, contract: OpenTargetsInput
) -> tuple[dict[str, tuple[Path, ...]], str]:
    manifest_path = root / contract.manifest_path
    _require_within_root(root, manifest_path, "Open Targets input manifest")
    manifest = _load_public_manifest(manifest_path, "Open Targets")
    if (
        manifest.release != contract.release
        or manifest.source.sha256 != contract.croissant_sha256
    ):
        raise RareDiseaseError("Open Targets input release differs from the contract")
    by_path = {output.path: output for output in manifest.outputs}
    croissant_path = f"releases/open-targets/{contract.release}/croissant.json"
    croissant = by_path.get(croissant_path)
    if croissant is None or croissant.sha256 != contract.croissant_sha256:
        raise RareDiseaseError("Open Targets Croissant output is missing")
    required_files: dict[str, tuple[Path, ...]] = {}
    for dataset in contract.dataset_rows:
        prefix = f"releases/open-targets/{contract.release}/output/{dataset}/"
        outputs = sorted(
            (output for path, output in by_path.items() if path.startswith(prefix)),
            key=lambda item: item.path,
        )
        if not outputs:
            raise RareDiseaseError(f"Open Targets input dataset is missing: {dataset}")
        paths = []
        for output in outputs:
            path = root / output.path
            _require_within_root(root, path, f"Open Targets {dataset} input")
            _verify_file(path, output.sha256, output.size_bytes)
            paths.append(path)
        required_files[dataset] = tuple(paths)
    _verify_file(root / croissant_path, croissant.sha256, croissant.size_bytes)
    return required_files, _file_identity(manifest_path)[0]


def _sql_paths(paths: Sequence[Path]) -> str:
    values = ["'" + str(path).replace("'", "''") + "'" for path in paths]
    return "[" + ",".join(values) + "]"


def _prepare_open_targets(
    connection: Any,
    files: Mapping[str, Sequence[Path]],
    contract: OpenTargetsInput,
) -> tuple[dict[str, str], dict[str, str]]:
    _register_normalise_curie(connection)
    for dataset, paths in files.items():
        path_list = _sql_paths(paths)
        observed = int(
            connection.sql(
                f"SELECT count(*) FROM read_parquet({path_list}, union_by_name=true)"
            ).fetchone()[0]
        )
        expected = contract.dataset_rows[dataset]
        if observed != expected:
            raise RareDiseaseError(
                f"Open Targets row count differs for {dataset}: "
                f"expected {expected}, observed {observed}"
            )
        connection.execute(
            f"CREATE VIEW ot_{dataset} AS "
            f"SELECT * FROM read_parquet({path_list}, union_by_name=true)"
        )
    connection.execute(
        """
        CREATE TEMP TABLE ot_disease_map AS
        SELECT id, name,
               CASE WHEN coalesce(len(children), 0) > 0 THEN 'broad'
                    ELSE 'disorder' END AS specificity
        FROM ot_disease
        """
    )
    connection.execute(
        """
        CREATE TEMP TABLE ot_gene_map AS
        SELECT target.id AS ensembl_gene_id,
               target.approvedSymbol AS gene_symbol,
               'HGNC:' || x.id AS hgnc_id
        FROM ot_target AS target, UNNEST(target.dbXrefs) AS item(x)
        WHERE x.source = 'HGNC'
        """
    )
    duplicate_hgnc = int(
        connection.sql(
            """
            SELECT count(*) FROM (
              SELECT hgnc_id FROM ot_gene_map
              GROUP BY hgnc_id HAVING count(DISTINCT ensembl_gene_id) > 1
            )
            """
        ).fetchone()[0]
    )
    if duplicate_hgnc:
        raise RareDiseaseError("Open Targets maps one HGNC ID to multiple targets")
    hgnc_to_ensembl = dict(
        connection.sql(
            "SELECT hgnc_id, min(ensembl_gene_id) FROM ot_gene_map GROUP BY hgnc_id"
        ).fetchall()
    )
    disease_specificity = dict(
        connection.sql("SELECT id, specificity FROM ot_disease_map").fetchall()
    )
    return hgnc_to_ensembl, disease_specificity


def _insert_open_targets_phenotypes(
    connection: Any,
    release: str,
    hpo: HpoTerms,
) -> tuple[int, int]:
    source_ids = [
        _normalise_curie(row[0])
        for row in connection.sql(
            "SELECT DISTINCT phenotype FROM ot_disease_phenotype "
            "WHERE starts_with(phenotype, 'HP_')"
        ).fetchall()
    ]
    resolutions = []
    for source_id in source_ids:
        if source_id is None:
            raise RareDiseaseError("Open Targets HPO identifier is null")
        canonical_id, canonical_name, _obsolete = hpo.resolve(source_id)
        resolutions.append((source_id, canonical_id, canonical_name))
    connection.execute(
        "CREATE TEMP TABLE hpo_resolution "
        "(source_id VARCHAR, canonical_id VARCHAR, canonical_name VARCHAR)"
    )
    connection.executemany(
        "INSERT INTO hpo_resolution VALUES (?, ?, ?)", resolutions
    )
    non_hpo = int(
        connection.sql(
            "SELECT count(*) FROM ot_disease_phenotype "
            "WHERE NOT starts_with(phenotype, 'HP_')"
        ).fetchone()[0]
    )
    expected_hpo_rows = int(
        connection.sql(
            "SELECT coalesce(sum(len(evidence)), 0) FROM ot_disease_phenotype "
            "WHERE starts_with(phenotype, 'HP_')"
        ).fetchone()[0]
    )
    connection.execute(
        f"""
        INSERT INTO phenotype_assertion
        WITH flattened AS (
          SELECT source.disease, source.phenotype, evidence_item.e,
                 evidence_item.evidence_index
          FROM ot_disease_phenotype AS source
          CROSS JOIN UNNEST(source.evidence) WITH ORDINALITY
            AS evidence_item(e, evidence_index)
          WHERE starts_with(source.phenotype, 'HP_')
        ), keyed AS (
          SELECT flattened.*,
                 CASE WHEN starts_with(norm_curie(e.diseaseFromSourceId), 'ORPHA:')
                      THEN 'orphanet'
                      ELSE lower(coalesce(e.resource, 'hpo')) END AS source_family,
                 norm_curie(e.diseaseFromSourceId) AS source_disease_id,
                 norm_curie(disease) AS mapped_disease_id,
                 coalesce(
                   norm_curie(e.diseaseFromSourceId),
                   'record:' || sha256(concat_ws('|',
                     lower(coalesce(e.resource, 'hpo')), norm_curie(phenotype), to_json(e)))
                 ) AS source_identity
          FROM flattened
        ), lineaged AS (
          SELECT keyed.*,
                 sha256(concat_ws('|', 'phenotype', source_family,
                   source_identity, resolution.canonical_id))
                   AS lineage_key
          FROM keyed
          JOIN hpo_resolution AS resolution
            ON resolution.source_id = norm_curie(keyed.phenotype)
        ), dependencies AS (
          SELECT lineaged.*,
                 sha256(concat_ws('|', 'dependency-v2', lineage_key,
                   coalesce(norm_curie(e.frequency), ''),
                   coalesce(e.qualifier, ''), coalesce(e.aspect, ''),
                   to_json(list_sort(list_distinct(coalesce(e.onset, []::VARCHAR[])))),
                   coalesce(e.sex, ''),
                   to_json(list_sort(list_distinct(coalesce(e.modifiers, []::VARCHAR[])))),
                   '', coalesce(e.evidenceType, ''),
                   to_json(list_sort(list_distinct(
                     coalesce(e."references", []::VARCHAR[])))))) AS dependency_key
          FROM lineaged
        )
        SELECT
          p.source_family,
          'open_targets_disease_phenotype',
          '{release}',
          sha256(concat_ws('|', p.disease, p.phenotype, cast(evidence_index AS VARCHAR), to_json(e))),
          p.lineage_key,
          p.dependency_key,
          sha256(concat_ws('|', 'assertion-v2', p.dependency_key,
            coalesce(p.mapped_disease_id, ''))),
          'integrated_copy', 'ot_mapped_evidence', false, 2, false, false,
          p.source_disease_id, p.mapped_disease_id,
          coalesce(orpha.specificity, mapped.specificity, 'unresolved'),
          p.phenotype, NULL, resolution.canonical_id, resolution.canonical_name,
          norm_curie(e.frequency), NULL, coalesce(e.qualifier, ''),
          coalesce(e.qualifierNot, false), e.aspect,
          coalesce(e.onset, []::VARCHAR[]), e.sex,
          coalesce(e.modifiers, []::VARCHAR[]), NULL, NULL,
          e.evidenceType, coalesce(e."references", []::VARCHAR[]), e.bioCuration,
          to_json(e), NULL, NULL
        FROM dependencies AS p
        JOIN hpo_resolution AS resolution ON resolution.source_id = norm_curie(p.phenotype)
        LEFT JOIN ot_disease_map AS mapped ON mapped.id = p.disease
        LEFT JOIN (
          SELECT disease_id, min(specificity) AS specificity
          FROM orpha_specificity GROUP BY disease_id
        ) AS orpha ON orpha.disease_id = norm_curie(e.diseaseFromSourceId)
        """
    )
    inserted = int(
        connection.sql(
            "SELECT count(*) FROM phenotype_assertion "
            "WHERE source_dataset = 'open_targets_disease_phenotype'"
        ).fetchone()[0]
    )
    if inserted != expected_hpo_rows:
        raise RareDiseaseError("Open Targets HPO evidence flattening lost records")
    return inserted, non_hpo


def _ot_evidence_insert_sql(dataset: str, source_family: str, release: str) -> str:
    fields = {
        "evidence_orphanet": {
            "inheritance": "[]::VARCHAR[]",
            "consequence": "evidence.variantFunctionalConsequenceId",
            "publication_date": "evidence.publicationDate",
            "publications": "coalesce(evidence.literature, []::VARCHAR[])",
            "study": "NULL",
            "report": "NULL",
            "notes": "CASE WHEN len(evidence.qualityControls) > 0 "
            "THEN to_json(evidence.qualityControls) END",
        },
        "evidence_clingen": {
            "inheritance": "coalesce(evidence.allelicRequirements, []::VARCHAR[])",
            "consequence": "NULL",
            "publication_date": "evidence.releaseDate",
            "publications": "[]::VARCHAR[]",
            "study": "evidence.studyId",
            "report": "CASE WHEN len(evidence.urls) > 0 THEN evidence.urls[1].url END",
            "notes": "CASE WHEN len(evidence.qualityControls) > 0 "
            "THEN to_json(evidence.qualityControls) END",
        },
        "evidence_gene2phenotype": {
            "inheritance": "coalesce(evidence.allelicRequirements, []::VARCHAR[])",
            "consequence": "evidence.variantFunctionalConsequenceId",
            "publication_date": "evidence.publicationDate",
            "publications": "coalesce(evidence.literature, []::VARCHAR[])",
            "study": "evidence.studyId",
            "report": "NULL",
            "notes": "CASE WHEN len(evidence.qualityControls) > 0 "
            "THEN to_json(evidence.qualityControls) END",
        },
        "evidence_genomics_england": {
            "inheritance": "coalesce(evidence.allelicRequirements, []::VARCHAR[])",
            "consequence": "NULL",
            "publication_date": "evidence.publicationDate",
            "publications": "coalesce(evidence.literature, []::VARCHAR[])",
            "study": "evidence.studyId",
            "report": "NULL",
            "notes": "evidence.studyOverview",
        },
    }[dataset]
    return f"""
        INSERT INTO gene_assertion
        WITH keyed AS (
          SELECT source.*,
                 coalesce(norm_curie(source.diseaseFromSourceId),
                          'record:' || source.id) AS source_identity,
                 to_json(source) AS source_payload
          FROM ot_{dataset} AS source
        ), lineaged AS (
          SELECT keyed.*,
                 sha256(concat_ws('|', 'gene', '{source_family}',
                   source_identity, coalesce(targetId, '')))
                   AS lineage_key
          FROM keyed
        ), dependencies AS (
          SELECT evidence.*,
                 sha256(concat_ws('|', 'dependency-v2', lineage_key,
                   coalesce(datasourceId, ''), coalesce(confidence, ''), '',
                   to_json(list_sort(list_distinct({fields['inheritance']}))),
                   coalesce({fields['consequence']}, ''), coalesce(confidence, ''),
                   to_json(list_sort(list_distinct({fields['publications']}))),
                   coalesce({fields['study']}, ''))) AS dependency_key
          FROM lineaged AS evidence
        )
        SELECT
          '{source_family}', '{dataset}', '{release}', evidence.id,
          evidence.lineage_key,
          evidence.dependency_key,
          sha256(concat_ws('|', 'assertion-v2', evidence.dependency_key,
            coalesce(norm_curie(evidence.diseaseId), ''))),
          'integrated_copy', 'ot_mapped_evidence', false, 2, false, false,
          norm_curie(evidence.diseaseFromSourceId), norm_curie(evidence.diseaseId),
          evidence.diseaseFromSource, mapped.name,
          coalesce(orpha.specificity, mapped.specificity, 'unresolved'),
          evidence.targetId, evidence.targetId, genes.hgnc_id,
          coalesce(evidence.targetFromSourceId, genes.gene_symbol), NULL, NULL, NULL,
          evidence.datasourceId, evidence.datasourceId, NULL, evidence.confidence,
          NULL, NULL, {fields['inheritance']}, {fields['consequence']},
          evidence.confidence, evidence.evidenceDate, {fields['publication_date']},
          NULL, {fields['publications']}, {fields['study']}, evidence.score,
          NULL, {fields['report']}, NULL, {fields['notes']}, evidence.source_payload
        FROM dependencies AS evidence
        LEFT JOIN ot_disease_map AS mapped ON mapped.id = evidence.diseaseId
        LEFT JOIN ot_gene_map AS genes ON genes.ensembl_gene_id = evidence.targetId
        LEFT JOIN (
          SELECT disease_id, min(specificity) AS specificity
          FROM orpha_specificity GROUP BY disease_id
        ) AS orpha ON orpha.disease_id = norm_curie(evidence.diseaseFromSourceId)
    """


def _insert_open_targets_genes(
    connection: Any,
    release: str,
) -> tuple[int, int]:
    evidence_count = 0
    for dataset, source_family in OT_EVIDENCE_DATASETS.items():
        connection.execute(_ot_evidence_insert_sql(dataset, source_family, release))
        evidence_count += int(
            connection.sql(
                "SELECT count(*) FROM gene_assertion WHERE source_dataset = ?",
                params=[dataset],
            ).fetchone()[0]
        )
    values = ", ".join(f"'{value}'" for value in OT_EVIDENCE_DATASETS.values())
    connection.execute(
        f"""
        INSERT INTO gene_assertion
        WITH lineaged AS (
          SELECT association.*,
                 sha256(concat_ws('|', 'gene', aggregationValue,
                   norm_curie(diseaseId), targetId)) AS lineage_key
          FROM ot_association_by_datasource_direct AS association
          WHERE aggregationType = 'datasourceId'
            AND aggregationValue IN ({values})
        ), dependencies AS (
          SELECT lineaged.*,
                 sha256(concat_ws('|', 'dependency-v2', lineage_key,
                   aggregationType, aggregationValue,
                   cast(associationScore AS VARCHAR),
                   cast(evidenceCount AS VARCHAR))) AS dependency_key
          FROM lineaged
        )
        SELECT
          aggregationValue, 'association_by_datasource_direct', '{release}',
          sha256(concat_ws('|', diseaseId, targetId, aggregationValue)),
          association.lineage_key, association.dependency_key,
          sha256(concat_ws('|', 'assertion-v2', association.dependency_key,
            coalesce(norm_curie(diseaseId), ''))),
          'integrated_aggregate', 'ot_direct_aggregate', false, 3, false, false,
          norm_curie(diseaseId), norm_curie(diseaseId), mapped.name, mapped.name,
          coalesce(mapped.specificity, 'unresolved'), targetId, targetId,
          genes.hgnc_id, genes.gene_symbol, NULL, NULL, NULL,
          aggregationType, aggregationValue, NULL, NULL, NULL, NULL,
          []::VARCHAR[], NULL, NULL, NULL, NULL, NULL, []::VARCHAR[], NULL,
          associationScore, evidenceCount, NULL, NULL, NULL,
          to_json(association)
        FROM dependencies AS association
        LEFT JOIN ot_disease_map AS mapped ON mapped.id = association.diseaseId
        LEFT JOIN ot_gene_map AS genes ON genes.ensembl_gene_id = association.targetId
        """
    )
    aggregate_count = int(
        connection.sql(
            "SELECT count(*) FROM gene_assertion "
            "WHERE source_dataset = 'association_by_datasource_direct'"
        ).fetchone()[0]
    )
    return evidence_count, aggregate_count


def _set_representation_preferences(connection: Any) -> None:
    for table in ("phenotype_assertion", "gene_assertion"):
        connection.execute(
            f"""
            CREATE TEMP TABLE {table}_choice AS
            SELECT source_dataset, source_record_id, assertion_key
            FROM (
              SELECT source_dataset, source_record_id, assertion_key,
                     row_number() OVER (
                       PARTITION BY dependency_key
                       ORDER BY representation_priority, source_dataset,
                                source_record_id, assertion_key
                     ) AS preference_rank
              FROM {table}
            )
            WHERE preference_rank = 1
            """
        )
        connection.execute(
            f"""
            UPDATE {table} AS assertion
            SET preferred_representation = true
            FROM {table}_choice AS choice
            WHERE assertion.source_dataset = choice.source_dataset
              AND assertion.source_record_id = choice.source_record_id
              AND assertion.assertion_key = choice.assertion_key
            """
        )
        connection.execute(
            f"""
            CREATE TEMP TABLE {table}_independent_choice AS
            SELECT source_dataset, source_record_id, assertion_key
            FROM (
              SELECT source_dataset, source_record_id, assertion_key,
                     row_number() OVER (
                       PARTITION BY lineage_key
                       ORDER BY representation_priority, source_dataset,
                                source_record_id, assertion_key
                     ) AS independence_rank
              FROM {table}
              WHERE preferred_representation
                AND association_scope = 'source_direct'
            )
            WHERE independence_rank = 1
            """
        )
        connection.execute(
            f"""
            UPDATE {table} AS assertion
            SET counts_as_independent = true
            FROM {table}_independent_choice AS choice
            WHERE assertion.source_dataset = choice.source_dataset
              AND assertion.source_record_id = choice.source_record_id
              AND assertion.assertion_key = choice.assertion_key
            """
        )


def _semantic_validation_observations(connection: Any) -> dict[str, int]:
    def scalar(query: str) -> int:
        return int(connection.sql(query).fetchone()[0])

    return {
        "preferred_phenotype_assertions": scalar(
            "SELECT count(*) FROM phenotype_assertion WHERE preferred_representation"
        ),
        "preferred_gene_assertions": scalar(
            "SELECT count(*) FROM gene_assertion WHERE preferred_representation"
        ),
        "shared_phenotype_lineages": scalar(
            "SELECT count(*) FROM (SELECT lineage_key FROM phenotype_assertion "
            "GROUP BY lineage_key HAVING count(DISTINCT source_dataset) > 1)"
        ),
        "shared_gene_lineages": scalar(
            "SELECT count(*) FROM (SELECT lineage_key FROM gene_assertion "
            "GROUP BY lineage_key HAVING count(DISTINCT source_dataset) > 1)"
        ),
        "phenotype_direct_dependency_matches": scalar(
            "SELECT count(*) FROM phenotype_assertion AS integrated WHERE "
            "integrated.source_dataset = 'open_targets_disease_phenotype' AND EXISTS "
            "(SELECT 1 FROM phenotype_assertion AS direct WHERE "
            "direct.dependency_key = integrated.dependency_key AND "
            "direct.association_scope = 'source_direct')"
        ),
        "phenotype_dependencies_with_multiple_mapped_diseases": scalar(
            "SELECT count(*) FROM (SELECT dependency_key FROM phenotype_assertion "
            "GROUP BY dependency_key HAVING count(DISTINCT mapped_disease_id) > 1)"
        ),
        "orphanet_direct_gene_shared_lineages": scalar(
            "SELECT count(*) FROM (SELECT lineage_key FROM gene_assertion WHERE "
            "source_dataset IN ('orphadata_product6', 'gencc_submissions') "
            "GROUP BY lineage_key HAVING count(*) FILTER "
            "(WHERE source_dataset = 'orphadata_product6') > 0 AND count(*) FILTER "
            "(WHERE source_dataset = 'gencc_submissions') > 0)"
        ),
        "integrated_preferred_over_direct_dependencies": scalar(
            "SELECT (SELECT count(*) FROM phenotype_assertion AS integrated WHERE "
            "integrated.preferred_representation AND integrated.association_scope "
            "<> 'source_direct' AND EXISTS (SELECT 1 FROM phenotype_assertion AS direct "
            "WHERE direct.dependency_key = integrated.dependency_key AND "
            "direct.association_scope = 'source_direct')) + "
            "(SELECT count(*) FROM gene_assertion AS integrated WHERE "
            "integrated.preferred_representation AND integrated.association_scope "
            "<> 'source_direct' AND EXISTS (SELECT 1 FROM gene_assertion AS direct "
            "WHERE direct.dependency_key = integrated.dependency_key AND "
            "direct.association_scope = 'source_direct'))"
        ),
        "phenotype_lineages_with_multiple_independent": scalar(
            "SELECT count(*) FROM (SELECT lineage_key FROM phenotype_assertion "
            "GROUP BY lineage_key HAVING count(*) FILTER "
            "(WHERE counts_as_independent) > 1)"
        ),
        "gene_lineages_with_multiple_independent": scalar(
            "SELECT count(*) FROM (SELECT lineage_key FROM gene_assertion "
            "GROUP BY lineage_key HAVING count(*) FILTER "
            "(WHERE counts_as_independent) > 1)"
        ),
        "assertions_spanning_mapped_diseases": scalar(
            "SELECT (SELECT count(*) FROM (SELECT assertion_key FROM "
            "phenotype_assertion GROUP BY assertion_key HAVING "
            "count(DISTINCT mapped_disease_id) > 1)) + (SELECT count(*) FROM "
            "(SELECT assertion_key FROM gene_assertion GROUP BY assertion_key HAVING "
            "count(DISTINCT mapped_disease_id) > 1))"
        ),
    }


def _export_tables(connection: Any, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True, mode=0o700)
    for table in sorted(TABLE_SCHEMAS):
        path = destination / f"{table}.parquet"
        escaped = str(path).replace("'", "''")
        connection.execute(
            f"COPY (SELECT * FROM {table} ORDER BY {TABLE_SORTS[table]}) "
            f"TO '{escaped}' (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 122880)"
        )
        os.chmod(path, 0o444)


def _build_projection(
    *,
    contract: RareDiseaseContract,
    source_paths: Mapping[str, Path],
    hpo_ontology_path: Path,
    hpo_annotations_path: Path,
    open_targets_files: Mapping[str, Sequence[Path]],
    destination: Path,
) -> None:
    try:
        import duckdb
    except ImportError as error:
        raise RareDiseaseError("DuckDB is required to build rare-disease projections") from error
    database_path = destination.parent / f".{destination.name}.duckdb"
    database_path.unlink(missing_ok=True)
    connection = duckdb.connect(str(database_path))
    try:
        connection.execute("PRAGMA disable_progress_bar")
        connection.execute("SET threads = 4")
        connection.execute(
            "SET temp_directory = ?", [str(destination.parent / ".duckdb-tmp")]
        )
        _create_tables(connection)
        connection.execute("BEGIN TRANSACTION")
        ot_hgnc, ot_disease_specificity = _prepare_open_targets(
            connection, open_targets_files, contract.open_targets
        )
        nomenclature = contract.artifact("orpha_nomenclature")
        specificity, unresolved_relations = _parse_nomenclature(
            connection, source_paths[nomenclature.role], nomenclature
        )
        disease_ids = set(specificity)
        connection.execute(
            "CREATE TEMP TABLE orpha_specificity "
            "(disease_id VARCHAR, specificity VARCHAR)"
        )
        connection.executemany(
            "INSERT INTO orpha_specificity VALUES (?, ?)", sorted(specificity.items())
        )
        alias_count = _parse_disease_aliases(
            connection,
            source_paths["orphadata_diseases"],
            contract.artifact("orphadata_diseases"),
            disease_ids,
        )
        hpo = _load_hpo_terms(hpo_ontology_path)
        orphadata_phenotypes, noncanonical_hpo = _parse_orphadata_phenotypes(
            connection,
            source_paths["orphadata_phenotypes"],
            contract.artifact("orphadata_phenotypes"),
            specificity,
            hpo,
        )
        orphadata_genes, direct_hgnc = _parse_orphadata_genes(
            connection,
            source_paths["orphadata_genes"],
            contract.artifact("orphadata_genes"),
            specificity,
        )
        unresolved_direct_targets = int(
            connection.sql(
                """
                SELECT count(*) FROM gene_assertion AS assertion
                LEFT JOIN ot_gene_map AS target
                  ON target.ensembl_gene_id = assertion.ensembl_gene_id
                WHERE assertion.source_dataset = 'orphadata_product6'
                  AND assertion.ensembl_gene_id IS NOT NULL
                  AND target.ensembl_gene_id IS NULL
                """
            ).fetchone()[0]
        )
        if unresolved_direct_targets:
            raise RareDiseaseError("Orphadata Ensembl genes do not resolve in Open Targets")
        onset_count, inheritance_count = _parse_natural_history(
            connection,
            source_paths["orphadata_natural_history"],
            contract.artifact("orphadata_natural_history"),
            specificity,
        )
        prevalence_count = _parse_prevalence(
            connection,
            source_paths["orphadata_prevalence"],
            contract.artifact("orphadata_prevalence"),
            specificity,
        )
        hpo_annotation_count = _parse_hpo_annotations(
            connection,
            hpo_annotations_path,
            contract.hpo.release,
            contract.hpo.expected_annotation_records,
            specificity,
            hpo,
        )
        hgnc_to_ensembl = dict(ot_hgnc)
        for hgnc, ensembl in direct_hgnc.items():
            hgnc_to_ensembl.setdefault(hgnc, ensembl)
        gencc_count, unresolved_genes, unresolved_diseases = _parse_gencc(
            connection,
            source_paths["gencc_submissions"],
            contract.artifact("gencc_submissions"),
            contract.gencc_columns,
            specificity,
            hgnc_to_ensembl,
            ot_disease_specificity,
        )
        ot_hpo_count, ot_non_hpo = _insert_open_targets_phenotypes(
            connection, contract.open_targets.release, hpo
        )
        ot_gene_count, ot_aggregate_count = _insert_open_targets_genes(
            connection, contract.open_targets.release
        )
        _set_representation_preferences(connection)

        source_counts = {
            "disease": len(disease_ids),
            "disease_alias": alias_count,
            "disease_relation": int(
                connection.sql("SELECT count(*) FROM disease_relation").fetchone()[0]
            ),
            "phenotype_assertion": (
                orphadata_phenotypes + hpo_annotation_count + ot_hpo_count
            ),
            "gene_assertion": (
                orphadata_genes + gencc_count + ot_gene_count + ot_aggregate_count
            ),
            "disease_onset": onset_count,
            "disease_inheritance": inheritance_count,
            "disease_prevalence": prevalence_count,
        }
        if source_counts != dict(contract.projection.expected_row_counts):
            raise RareDiseaseError(
                "rare-disease projection row counts differ from the contract: "
                f"{source_counts}"
            )
        observations = {
            "hpo_noncanonical_source_ids": len(noncanonical_hpo),
            "open_targets_hpo_evidence": ot_hpo_count,
            "open_targets_non_hpo_excluded": ot_non_hpo,
            "open_targets_gene_evidence": ot_gene_count,
            "open_targets_gene_aggregates": ot_aggregate_count,
            "unresolved_gencc_genes": unresolved_genes,
            "unresolved_gencc_diseases": unresolved_diseases,
            "unresolved_orpha_relation_endpoints": unresolved_relations,
        }
        observations.update(_semantic_validation_observations(connection))
        if observations != dict(contract.validation_expectations):
            raise RareDiseaseError(
                "rare-disease validation observations differ from the contract: "
                f"{observations}"
            )
        _export_tables(connection, destination)
        connection.execute("COMMIT")
    finally:
        connection.close()
        database_path.unlink(missing_ok=True)
        shutil.rmtree(destination.parent / ".duckdb-tmp", ignore_errors=True)


def _verify_projection_tree(path: Path) -> None:
    if path.is_symlink() or not path.is_dir():
        raise RareDiseaseError(f"rare-disease projection is not a directory: {path}")
    expected = {f"{table}.parquet" for table in TABLE_SCHEMAS}
    observed: set[str] = set()
    for item in path.iterdir():
        if item.is_symlink() or not item.is_file():
            raise RareDiseaseError(f"rare-disease projection has an invalid entry: {item}")
        observed.add(item.name)
    if observed != expected:
        raise RareDiseaseError("rare-disease projection tree is incomplete or contaminated")


def _validate_projection(
    path: Path,
    contract: RareDiseaseContract,
    open_targets_files: Mapping[str, Sequence[Path]],
) -> RareDiseaseValidation:
    try:
        import duckdb
    except ImportError as error:
        raise RareDiseaseError("DuckDB is required to validate rare-disease projections") from error
    _verify_projection_tree(path)
    connection = duckdb.connect(database=":memory:")
    try:
        _register_normalise_curie(connection)
        row_counts: dict[str, int] = {}
        for table, columns in TABLE_COLUMNS.items():
            parquet = path / f"{table}.parquet"
            relation = connection.from_parquet(str(parquet))
            if tuple(relation.columns) != columns:
                raise RareDiseaseError(f"rare-disease projection schema differs: {table}")
            row_counts[table] = int(relation.aggregate("count(*)").fetchone()[0])
            escaped = str(parquet).replace("'", "''")
            connection.execute(
                f"CREATE VIEW {table} AS SELECT * FROM read_parquet('{escaped}')"
            )
        if row_counts != dict(contract.projection.expected_row_counts):
            raise RareDiseaseError("rare-disease projection row counts differ from contract")
        checks = {
            "duplicate diseases": "SELECT count(*) - count(DISTINCT disease_id) FROM disease",
            "orphan disease aliases": (
                "SELECT count(*) FROM disease_alias AS alias "
                "LEFT JOIN disease USING (disease_id) WHERE disease.disease_id IS NULL"
            ),
            "orphan phenotype diseases": (
                "SELECT count(*) FROM phenotype_assertion WHERE "
                "starts_with(source_disease_id, 'ORPHA:') AND source_disease_id "
                "NOT IN (SELECT disease_id FROM disease)"
            ),
            "orphan gene diseases": (
                "SELECT count(*) FROM gene_assertion WHERE "
                "starts_with(source_disease_id, 'ORPHA:') AND source_disease_id "
                "NOT IN (SELECT disease_id FROM disease)"
            ),
            "invalid HPO projection": (
                "SELECT count(*) FROM phenotype_assertion "
                "WHERE canonical_hpo_id IS NULL OR NOT starts_with(canonical_hpo_id, 'HP:')"
            ),
            "missing disease subtype identity": (
                "SELECT count(*) FROM disease WHERE classification_level_id IS NULL "
                "OR classification_level IS NULL"
            ),
            "missing natural-history identity": (
                "SELECT (SELECT count(*) FROM disease_onset WHERE onset_id IS NULL "
                "OR onset IS NULL) + (SELECT count(*) FROM disease_inheritance "
                "WHERE inheritance_id IS NULL OR inheritance IS NULL)"
            ),
            "missing direct validity": (
                "SELECT count(*) FROM gene_assertion WHERE source_dataset IN "
                "('orphadata_product6', 'gencc_submissions') AND "
                "(validity_id IS NULL OR validity IS NULL)"
            ),
            "missing GenCC inheritance": (
                "SELECT count(*) FROM gene_assertion WHERE "
                "source_dataset = 'gencc_submissions' AND "
                "(inheritance_id IS NULL OR inheritance IS NULL)"
            ),
            "misclassified GenCC run date": (
                "SELECT count(*) FROM gene_assertion WHERE "
                "source_dataset = 'gencc_submissions' AND "
                "(publication_date IS NOT NULL OR source_run_date IS NULL)"
            ),
            "lossy Open Targets gene provenance": (
                "SELECT count(*) FROM gene_assertion WHERE source_dataset IN "
                f"({','.join(repr(item) for item in OT_EVIDENCE_DATASETS)}) AND "
                "(provenance_raw IS NULL OR NOT json_valid(provenance_raw) OR "
                "try(json_extract_string(provenance_raw, '$.id')) <> source_record_id)"
            ),
            "lossy Orphadata assertion provenance": (
                "SELECT (SELECT count(*) FROM phenotype_assertion WHERE "
                "source_dataset = 'orphadata_product4' AND "
                "try(json_extract_string(provenance_raw, '$.assertion_xml')) IS NULL) + "
                "(SELECT count(*) FROM gene_assertion WHERE "
                "source_dataset = 'orphadata_product6' AND "
                "try(json_extract_string(provenance_raw, '$.assertion_xml')) IS NULL)"
            ),
            "lossy HPO annotation provenance": (
                "SELECT count(*) FROM phenotype_assertion WHERE "
                "source_dataset = 'hpo_disease_annotations' AND "
                "(norm_curie(try(json_extract_string(provenance_raw, '$.database_id'))) "
                "IS DISTINCT FROM source_disease_id OR "
                "try(json_extract_string(provenance_raw, '$.hpo_id')) "
                "IS DISTINCT FROM source_hpo_id)"
            ),
            "lossy Open Targets phenotype provenance": (
                "SELECT count(*) FROM phenotype_assertion WHERE "
                "source_dataset = 'open_targets_disease_phenotype' AND "
                "norm_curie(try(json_extract_string(provenance_raw, "
                "'$.diseaseFromSourceId'))) IS DISTINCT FROM source_disease_id"
            ),
            "lossy Open Targets aggregate provenance": (
                "SELECT count(*) FROM gene_assertion WHERE "
                "source_dataset = 'association_by_datasource_direct' AND "
                "try(json_extract_string(provenance_raw, '$.aggregationValue')) "
                "IS DISTINCT FROM source_family"
            ),
            "duplicate assertion identity": (
                "SELECT (SELECT count(*) FROM (SELECT source_dataset, source_record_id, "
                "assertion_key FROM phenotype_assertion GROUP BY ALL HAVING count(*) > 1)) + "
                "(SELECT count(*) FROM (SELECT source_dataset, source_record_id, "
                "assertion_key FROM gene_assertion GROUP BY ALL HAVING count(*) > 1))"
            ),
            "phenotype preference cardinality": (
                "SELECT count(*) FROM (SELECT dependency_key FROM phenotype_assertion "
                "GROUP BY dependency_key HAVING count(*) FILTER "
                "(WHERE preferred_representation) <> 1)"
            ),
            "gene preference cardinality": (
                "SELECT count(*) FROM (SELECT dependency_key FROM gene_assertion "
                "GROUP BY dependency_key HAVING count(*) FILTER "
                "(WHERE preferred_representation) <> 1)"
            ),
            "phenotype independence cardinality": (
                "SELECT count(*) FROM (SELECT lineage_key FROM phenotype_assertion "
                "GROUP BY lineage_key HAVING count(*) FILTER (WHERE counts_as_independent) "
                "<> CASE WHEN count(*) FILTER (WHERE preferred_representation AND "
                "association_scope = 'source_direct') > 0 THEN 1 ELSE 0 END)"
            ),
            "gene independence cardinality": (
                "SELECT count(*) FROM (SELECT lineage_key FROM gene_assertion "
                "GROUP BY lineage_key HAVING count(*) FILTER (WHERE counts_as_independent) "
                "<> CASE WHEN count(*) FILTER (WHERE preferred_representation AND "
                "association_scope = 'source_direct') > 0 THEN 1 ELSE 0 END)"
            ),
            "indirect independent evidence": (
                "SELECT (SELECT count(*) FROM phenotype_assertion "
                "WHERE counts_as_independent AND association_scope <> 'source_direct') + "
                "(SELECT count(*) FROM gene_assertion WHERE counts_as_independent "
                "AND association_scope <> 'source_direct')"
            ),
            "integrated preference over direct evidence": (
                "SELECT (SELECT count(*) FROM phenotype_assertion AS preferred WHERE "
                "preferred.preferred_representation AND preferred.association_scope "
                "<> 'source_direct' AND EXISTS (SELECT 1 FROM phenotype_assertion AS direct "
                "WHERE direct.dependency_key = preferred.dependency_key AND "
                "direct.association_scope = 'source_direct')) + "
                "(SELECT count(*) FROM gene_assertion AS preferred WHERE "
                "preferred.preferred_representation AND preferred.association_scope "
                "<> 'source_direct' AND EXISTS (SELECT 1 FROM gene_assertion AS direct "
                "WHERE direct.dependency_key = preferred.dependency_key AND "
                "direct.association_scope = 'source_direct'))"
            ),
            "missing provenance": (
                "SELECT (SELECT count(*) FROM phenotype_assertion WHERE source_record_id IS NULL "
                "OR source_dataset IS NULL OR source_release IS NULL OR lineage_key IS NULL "
                "OR dependency_key IS NULL OR assertion_key IS NULL OR provenance_raw IS NULL "
                "OR NOT json_valid(provenance_raw)) + "
                "(SELECT count(*) FROM gene_assertion WHERE source_record_id IS NULL "
                "OR source_dataset IS NULL OR source_release IS NULL OR lineage_key IS NULL "
                "OR dependency_key IS NULL OR assertion_key IS NULL OR provenance_raw IS NULL "
                "OR NOT json_valid(provenance_raw))"
            ),
            "assertion spanning mapped diseases": (
                "SELECT (SELECT count(*) FROM (SELECT assertion_key FROM phenotype_assertion "
                "GROUP BY assertion_key HAVING count(DISTINCT mapped_disease_id) > 1)) + "
                "(SELECT count(*) FROM (SELECT assertion_key FROM gene_assertion "
                "GROUP BY assertion_key HAVING count(DISTINCT mapped_disease_id) > 1))"
            ),
        }
        gencc_missing_fields = " OR ".join(
            "try(json_extract(provenance_raw, '$.%s')) IS NULL" % column
            for column in contract.gencc_columns
        )
        checks["lossy GenCC provenance"] = (
            "SELECT count(*) FROM gene_assertion WHERE "
            "source_dataset = 'gencc_submissions' AND ("
            f"{gencc_missing_fields})"
        )
        for description, query in checks.items():
            if int(connection.sql(query).fetchone()[0]) != 0:
                raise RareDiseaseError(f"rare-disease projection has {description}")
        relation_unresolved = int(
            connection.sql(
                "SELECT count(*) FROM disease_relation WHERE unresolved_endpoint"
            ).fetchone()[0]
        )
        noncanonical_hpo = int(
            connection.sql(
                "SELECT count(DISTINCT source_hpo_id) FROM phenotype_assertion "
                "WHERE source_dataset = 'orphadata_product4' "
                "AND source_hpo_id <> canonical_hpo_id"
            ).fetchone()[0]
        )
        ot_hpo = int(
            connection.sql(
                "SELECT count(*) FROM phenotype_assertion "
                "WHERE source_dataset = 'open_targets_disease_phenotype'"
            ).fetchone()[0]
        )
        ot_non_hpo_paths = _sql_paths(open_targets_files["disease_phenotype"])
        ot_non_hpo = int(
            connection.sql(
                "SELECT count(*) FROM read_parquet("
                f"{ot_non_hpo_paths}, union_by_name=true) "
                "WHERE NOT starts_with(phenotype, 'HP_')"
            ).fetchone()[0]
        )
        ot_gene = int(
            connection.sql(
                "SELECT count(*) FROM gene_assertion "
                f"WHERE source_dataset IN ({','.join(repr(item) for item in OT_EVIDENCE_DATASETS)})"
            ).fetchone()[0]
        )
        ot_aggregate = int(
            connection.sql(
                "SELECT count(*) FROM gene_assertion "
                "WHERE source_dataset = 'association_by_datasource_direct'"
            ).fetchone()[0]
        )
        unresolved_genes = int(
            connection.sql(
                "SELECT count(DISTINCT hgnc_id) FROM gene_assertion "
                "WHERE source_dataset = 'gencc_submissions' AND ensembl_gene_id IS NULL"
            ).fetchone()[0]
        )
        disease_paths = _sql_paths(open_targets_files["disease"])
        unresolved_diseases = int(
            connection.sql(
                "SELECT count(DISTINCT mapped_disease_id) FROM gene_assertion "
                "WHERE source_dataset = 'gencc_submissions' AND replace(mapped_disease_id, ':', '_') "
                f"NOT IN (SELECT id FROM read_parquet({disease_paths}, union_by_name=true))"
            ).fetchone()[0]
        )
        mapped_without_specificity = int(
            connection.sql(
                "SELECT count(*) FROM gene_assertion WHERE "
                "source_dataset = 'gencc_submissions' AND disease_specificity = 'unresolved' "
                "AND replace(mapped_disease_id, ':', '_') IN "
                f"(SELECT id FROM read_parquet({disease_paths}, union_by_name=true))"
            ).fetchone()[0]
        )
        if mapped_without_specificity:
            raise RareDiseaseError(
                "mapped GenCC diseases lack Open Targets specificity"
            )
        observations = {
            "hpo_noncanonical_source_ids": noncanonical_hpo,
            "open_targets_hpo_evidence": ot_hpo,
            "open_targets_non_hpo_excluded": ot_non_hpo,
            "open_targets_gene_evidence": ot_gene,
            "open_targets_gene_aggregates": ot_aggregate,
            "unresolved_gencc_genes": unresolved_genes,
            "unresolved_gencc_diseases": unresolved_diseases,
            "unresolved_orpha_relation_endpoints": relation_unresolved,
        }
        semantic_observations = _semantic_validation_observations(connection)
        observations.update(semantic_observations)
        if observations != dict(contract.validation_expectations):
            raise RareDiseaseError(
                "rare-disease projection validation differs from contract: "
                f"{observations}"
            )
        preferred_phenotypes = semantic_observations["preferred_phenotype_assertions"]
        preferred_genes = semantic_observations["preferred_gene_assertions"]
        shared_phenotypes = semantic_observations["shared_phenotype_lineages"]
        shared_genes = semantic_observations["shared_gene_lineages"]
        scopes = {
            row[0]
            for row in connection.sql(
                "SELECT DISTINCT association_scope FROM gene_assertion"
            ).fetchall()
        }
        if not {"source_direct", "ot_mapped_evidence", "ot_direct_aggregate"}.issubset(scopes):
            raise RareDiseaseError("gene assertion scopes do not separate direct and integrated evidence")
        phenotype_scopes = {
            row[0]
            for row in connection.sql(
                "SELECT DISTINCT association_scope FROM phenotype_assertion"
            ).fetchall()
        }
        if not {"source_direct", "ot_mapped_evidence"}.issubset(phenotype_scopes):
            raise RareDiseaseError(
                "phenotype assertion scopes do not separate direct and integrated evidence"
            )
        return RareDiseaseValidation(
            row_counts=row_counts,
            preferred_phenotype_assertions=preferred_phenotypes,
            preferred_gene_assertions=preferred_genes,
            shared_phenotype_lineages=shared_phenotypes,
            shared_gene_lineages=shared_genes,
            hpo_noncanonical_source_ids=noncanonical_hpo,
            open_targets_hpo_evidence=ot_hpo,
            open_targets_non_hpo_excluded=ot_non_hpo,
            unresolved_gencc_genes=unresolved_genes,
            unresolved_gencc_diseases=unresolved_diseases,
            unresolved_orpha_relation_endpoints=relation_unresolved,
        )
    finally:
        connection.close()


def _source_manifest_path(
    root: Path,
    artifact: RareDiseaseArtifact,
    artifact_contract_sha256: str | None = None,
) -> Path:
    legacy = root / f"manifests/rare-disease-evidence/{artifact.resource_id}.v2.json"
    if artifact_contract_sha256 is None or not legacy.exists():
        return legacy
    try:
        document = json.loads(legacy.read_text(encoding="utf-8"))
        if (
            document.get("schema") == "genosapiens.public-resource"
            and document.get("transformation", {}).get("parameters_sha256")
            == artifact_contract_sha256
        ):
            return legacy
    except (OSError, json.JSONDecodeError, AttributeError):
        pass
    return legacy.with_name(
        f"{artifact.resource_id}.{artifact_contract_sha256[:12]}.json"
    )


def _artifact_contract_sha256(artifact: RareDiseaseArtifact) -> str:
    encoded = json.dumps(
        asdict(artifact), sort_keys=True, separators=(",", ":")
    ).encode()
    return hashlib.sha256(encoded).hexdigest()


def _source_manifest(
    artifact: RareDiseaseArtifact,
    artifact_contract_sha256: str,
    image_digest: str,
    retrieved_on: str,
) -> PublicResourceManifest:
    return PublicResourceManifest(
        resource_id=artifact.resource_id,
        release=artifact.release,
        source=SourceIdentity(artifact.url, sha256=artifact.sha256),
        license=artifact.license,
        assembly="not-applicable",
        retrieved_on=retrieved_on,
        transformation=Transformation(
            SOURCE_TRANSFORMATION_SCHEMA, artifact_contract_sha256
        ),
        image_digest=image_digest,
        outputs=(
            ResourceOutput(
                path=artifact.path,
                size_bytes=artifact.size_bytes,
                sha256=artifact.sha256,
                retention=RetentionPolicy("durable", None),
            ),
        ),
    )


def _verify_source_manifest(
    manifest: PublicResourceManifest,
    artifact: RareDiseaseArtifact,
    artifact_contract_sha256: str,
    root: Path,
) -> None:
    if (
        manifest.resource_id != artifact.resource_id
        or manifest.release != artifact.release
        or manifest.source != SourceIdentity(artifact.url, sha256=artifact.sha256)
        or manifest.license != artifact.license
        or manifest.assembly != "not-applicable"
        or manifest.transformation
        != Transformation(SOURCE_TRANSFORMATION_SCHEMA, artifact_contract_sha256)
        or len(manifest.outputs) != 1
        or manifest.outputs[0].path != artifact.path
        or manifest.outputs[0].size_bytes != artifact.size_bytes
        or manifest.outputs[0].sha256 != artifact.sha256
    ):
        raise RareDiseaseError(
            f"rare-disease source manifest differs from contract: {artifact.role}"
        )
    _verify_file(root / artifact.path, artifact.sha256, artifact.size_bytes)


def _publish_source_file(root: Path, artifact: RareDiseaseArtifact, source: Path) -> None:
    target = root / artifact.path
    _require_within_root(root, target, f"{artifact.role} release")
    if target.exists() or target.is_symlink():
        _verify_file(target, artifact.sha256, artifact.size_bytes)
        return
    target.parent.mkdir(parents=True, exist_ok=True)
    partial_root = target.parent / ".partial"
    _require_within_root(root, partial_root, f"{artifact.role} release partial path")
    partial_root.mkdir(mode=0o700, exist_ok=True)
    partial = partial_root / f"{target.name}.{uuid.uuid4().hex}.partial"
    try:
        with source.open("rb") as input_file, partial.open("xb") as output_file:
            os.chmod(partial, 0o600)
            shutil.copyfileobj(input_file, output_file, length=1024 * 1024)
            output_file.flush()
            os.fsync(output_file.fileno())
        _verify_file(partial, artifact.sha256, artifact.size_bytes)
        os.chmod(partial, 0o444)
        try:
            os.link(partial, target)
        except FileExistsError:
            _verify_file(target, artifact.sha256, artifact.size_bytes)
        directory_fd = os.open(target.parent, os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    finally:
        partial.unlink(missing_ok=True)
        try:
            partial_root.rmdir()
        except OSError:
            pass


def _source_publications(
    *,
    root: Path,
    contract: RareDiseaseContract,
    image_digest: str,
    retrieved_on: str,
    publish_missing: bool = True,
) -> tuple[tuple[Path, str], ...]:
    publications = []
    for artifact in contract.artifacts:
        artifact_contract_sha256 = _artifact_contract_sha256(artifact)
        path = _source_manifest_path(root, artifact, artifact_contract_sha256)
        _require_within_root(root, path, "rare-disease source manifest")
        if path.exists() or path.is_symlink():
            manifest = _load_public_manifest(path, artifact.role)
            _verify_source_manifest(
                manifest, artifact, artifact_contract_sha256, root
            )
            publications.append((path, _file_identity(path)[0]))
            continue
        if not publish_missing:
            raise RareDiseaseError(
                f"completed rare-disease bundle lacks source manifest: {artifact.role}"
            )
        publication = publish_manifest(
            _source_manifest(
                artifact, artifact_contract_sha256, image_digest, retrieved_on
            ),
            path,
            root,
        )
        publications.append((publication.path, publication.manifest_sha256))
    return tuple(publications)


def _projection_outputs(path: Path, contract: RareDiseaseContract) -> list[dict[str, Any]]:
    outputs = []
    for table in sorted(TABLE_SCHEMAS):
        file_path = path / f"{table}.parquet"
        sha256, size = _file_identity(file_path)
        outputs.append(
            {
                "table": table,
                "path": f"{contract.projection.output_prefix}/{table}.parquet",
                "size_bytes": size,
                "sha256": sha256,
                "row_count": contract.projection.expected_row_counts[table],
                "retention": "durable",
            }
        )
    return outputs


def _bundle_manifest_data(
    *,
    root: Path,
    contract: RareDiseaseContract,
    contract_sha256: str,
    image_digest: str,
    retrieved_on: str,
    source_publications: Sequence[tuple[Path, str]],
    open_targets_manifest_sha256: str,
    projection_path: Path,
    validation: RareDiseaseValidation,
) -> Mapping[str, Any]:
    sources = []
    for artifact, (manifest_path, manifest_sha256) in zip(
        contract.artifacts, source_publications, strict=True
    ):
        sources.append(
            {
                "role": artifact.role,
                "resource_id": artifact.resource_id,
                "release": artifact.release,
                "extraction_date": artifact.extraction_date,
                "source_revision": artifact.source_revision,
                "source_url": artifact.url,
                "source_sha256": artifact.sha256,
                "license": artifact.license,
                "license_url": artifact.license_url,
                "manifest_path": str(manifest_path.relative_to(root)),
                "manifest_sha256": manifest_sha256,
            }
        )
    return {
        "schema": BUNDLE_MANIFEST_SCHEMA,
        "classification": "public",
        "bundle_id": contract.bundle_id,
        "projection_id": contract.projection_id,
        "projection_schema": contract.projection.schema_version,
        "contract_sha256": contract_sha256,
        "assembly": contract.assembly,
        "retrieved_on": retrieved_on,
        "image_digest": image_digest,
        "sources": sources,
        "upstream": {
            "hpo_release": contract.hpo.release,
            "hpo_ontology_sha256": contract.hpo.ontology_sha256,
            "hpo_annotations_sha256": contract.hpo.annotations_sha256,
            "open_targets_release": contract.open_targets.release,
            "open_targets_croissant_sha256": contract.open_targets.croissant_sha256,
            "open_targets_manifest_path": contract.open_targets.manifest_path,
            "open_targets_manifest_sha256": open_targets_manifest_sha256,
        },
        "outputs": _projection_outputs(projection_path, contract),
        "validation": asdict(validation),
    }


def _load_bundle_manifest(path: Path) -> Mapping[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise RareDiseaseError(f"cannot read rare-disease bundle manifest: {error}") from error
    if not isinstance(value, Mapping):
        raise RareDiseaseError("rare-disease bundle manifest root must be an object")
    return value


def validate_bundle_manifest_metadata(value: Mapping[str, Any]) -> None:
    """Pure current-bundle metadata checks; output rows stay upstream-attested."""
    try:
        if (set(value) != {"schema", "classification", "bundle_id", "projection_id",
                "projection_schema", "contract_sha256", "assembly", "retrieved_on",
                "image_digest", "sources", "upstream", "outputs", "validation"}
                or value["schema"] != BUNDLE_MANIFEST_SCHEMA
                or value["classification"] != "public"
                or value["projection_schema"] != "v6"
                or any(not isinstance(value[key], str) or not value[key].strip()
                       for key in ("bundle_id", "projection_id", "assembly"))
                or SHA256.fullmatch(value["contract_sha256"]) is None
                or DIGEST.fullmatch(value["image_digest"]) is None):
            raise ValueError
        date.fromisoformat(value["retrieved_on"])
        upstream = value["upstream"]
        if set(upstream) != {"hpo_release", "hpo_ontology_sha256", "hpo_annotations_sha256",
                "open_targets_release", "open_targets_croissant_sha256",
                "open_targets_manifest_path", "open_targets_manifest_sha256"}:
            raise ValueError
        for key in ("hpo_ontology_sha256", "hpo_annotations_sha256",
                    "open_targets_croissant_sha256", "open_targets_manifest_sha256"):
            if SHA256.fullmatch(upstream[key]) is None:
                raise ValueError
        sources, outputs = value["sources"], value["outputs"]
        if (type(sources) is not list or {source["role"] for source in sources} != EXPECTED_ROLES
                or len(sources) != len(EXPECTED_ROLES) or type(outputs) is not list or not outputs):
            raise ValueError
        for source in sources:
            if (set(source) != {"role", "resource_id", "release", "extraction_date",
                    "source_revision", "source_url", "source_sha256", "license",
                    "license_url", "manifest_path", "manifest_sha256"}
                    or any(not isinstance(source[key], str) or not source[key].strip()
                           for key in ("resource_id", "release", "source_url", "license"))
                    or SHA256.fullmatch(source["source_sha256"]) is None
                    or SHA256.fullmatch(source["manifest_sha256"]) is None):
                raise ValueError
        row_counts = {}
        for output in outputs:
            if (set(output) != {"table", "path", "size_bytes", "sha256", "row_count", "retention"}
                    or output["table"] not in TABLE_SCHEMAS or output["table"] in row_counts
                    or type(output["size_bytes"]) is not int or output["size_bytes"] <= 0
                    or type(output["row_count"]) is not int or output["row_count"] < 0
                    or SHA256.fullmatch(output["sha256"]) is None or output["retention"] != "durable"):
                raise ValueError
            path = PurePosixPath(output["path"])
            if path.is_absolute() or ".." in path.parts or not output["path"].startswith("derived/"):
                raise ValueError
            row_counts[output["table"]] = output["row_count"]
        if set(row_counts) != set(TABLE_SCHEMAS) or value["validation"]["row_counts"] != row_counts:
            raise ValueError
    except (KeyError, TypeError, ValueError, AttributeError):
        raise RareDiseaseError("invalid current bundle provenance metadata") from None


def _verify_bundle_manifest(
    *,
    value: Mapping[str, Any],
    root: Path,
    contract: RareDiseaseContract,
    contract_sha256: str,
    source_publications: Sequence[tuple[Path, str]],
    open_targets_manifest_sha256: str,
    projection_path: Path,
    validation: RareDiseaseValidation,
) -> None:
    if (
        value.get("schema") != BUNDLE_MANIFEST_SCHEMA
        or value.get("classification") != "public"
        or value.get("bundle_id") != contract.bundle_id
        or value.get("projection_id") != contract.projection_id
        or value.get("projection_schema") != contract.projection.schema_version
        or value.get("contract_sha256") != contract_sha256
        or value.get("assembly") != contract.assembly
        or value.get("validation") != asdict(validation)
    ):
        raise RareDiseaseError("rare-disease bundle manifest identity is invalid")
    if not DIGEST.fullmatch(str(value.get("image_digest", ""))):
        raise RareDiseaseError("rare-disease bundle image digest is invalid")
    try:
        date.fromisoformat(str(value["retrieved_on"]))
    except (KeyError, ValueError) as error:
        raise RareDiseaseError("rare-disease retrieval date is invalid") from error
    expected_upstream = {
        "hpo_release": contract.hpo.release,
        "hpo_ontology_sha256": contract.hpo.ontology_sha256,
        "hpo_annotations_sha256": contract.hpo.annotations_sha256,
        "open_targets_release": contract.open_targets.release,
        "open_targets_croissant_sha256": contract.open_targets.croissant_sha256,
        "open_targets_manifest_path": contract.open_targets.manifest_path,
        "open_targets_manifest_sha256": open_targets_manifest_sha256,
    }
    if value.get("upstream") != expected_upstream:
        raise RareDiseaseError("rare-disease upstream manifest identities differ")
    sources = value.get("sources")
    if not isinstance(sources, list) or len(sources) != len(contract.artifacts):
        raise RareDiseaseError("rare-disease source manifest list is invalid")
    for artifact, source, (path, digest) in zip(
        contract.artifacts, sources, source_publications, strict=True
    ):
        if not isinstance(source, Mapping):
            raise RareDiseaseError("rare-disease source manifest entry is invalid")
        manifest_path, manifest_sha256 = str(path.relative_to(root)), digest
        expected_source = {
            "role": artifact.role,
            "resource_id": artifact.resource_id,
            "release": artifact.release,
            "extraction_date": artifact.extraction_date,
            "source_revision": artifact.source_revision,
            "source_url": artifact.url,
            "source_sha256": artifact.sha256,
            "license": artifact.license,
            "license_url": artifact.license_url,
            "manifest_path": manifest_path,
            "manifest_sha256": manifest_sha256,
        }
        if dict(source) != expected_source:
            raise RareDiseaseError(f"rare-disease source identity differs: {artifact.role}")
    expected_outputs = _projection_outputs(projection_path, contract)
    if value.get("outputs") != expected_outputs:
        raise RareDiseaseError("rare-disease projection output identity differs")


def _bundle_manifest_path(
    root: Path,
    contract: RareDiseaseContract,
    contract_sha256: str | None = None,
) -> Path:
    legacy = root / f"manifests/rare-disease-evidence/{contract.bundle_id}.bundle.json"
    if contract_sha256 is None or not legacy.exists():
        return legacy
    try:
        document = json.loads(legacy.read_text(encoding="utf-8"))
        if (
            document.get("schema") == BUNDLE_MANIFEST_SCHEMA
            and document.get("contract_sha256") == contract_sha256
        ):
            return legacy
    except (OSError, json.JSONDecodeError, AttributeError):
        pass
    return legacy.with_name(
        f"{contract.bundle_id}.{contract_sha256[:12]}.bundle.json"
    )


def _verify_historical_bundle_for_projection(
    *,
    value: Mapping[str, Any],
    root: Path,
    contract: RareDiseaseContract,
    projection_path: Path,
    validation: RareDiseaseValidation,
) -> None:
    """Authorize republishing an unchanged projection under a new schema contract."""
    expected_identity = {
        "schema": f"{BUNDLE_MANIFEST_SCHEMA}/v1",
        "classification": "public",
        "bundle_id": contract.bundle_id,
        "projection_id": contract.projection_id,
        "projection_schema": contract.projection.schema_version,
        "assembly": contract.assembly,
        "validation": asdict(validation),
    }
    if any(value.get(key) != expected for key, expected in expected_identity.items()):
        raise RareDiseaseError(
            "historical rare-disease bundle does not bind the existing projection"
        )
    expected_upstream = {
        "hpo_release": contract.hpo.release,
        "hpo_ontology_sha256": contract.hpo.ontology_sha256,
        "hpo_annotations_sha256": contract.hpo.annotations_sha256,
        "open_targets_release": contract.open_targets.release,
        "open_targets_croissant_sha256": contract.open_targets.croissant_sha256,
    }
    historical_upstream = value.get("upstream")
    if not isinstance(historical_upstream, Mapping) or any(
        historical_upstream.get(key) != expected
        for key, expected in expected_upstream.items()
    ):
        raise RareDiseaseError(
            "historical rare-disease bundle upstream identities differ"
        )
    historical_manifest_value = historical_upstream.get(
        "open_targets_manifest_path"
    )
    historical_manifest_sha256 = historical_upstream.get(
        "open_targets_manifest_sha256"
    )
    if (
        not isinstance(historical_manifest_value, str)
        or not isinstance(historical_manifest_sha256, str)
        or not SHA256.fullmatch(historical_manifest_sha256)
    ):
        raise RareDiseaseError(
            "historical Open Targets manifest identity is invalid"
        )
    historical_manifest_path = root / _safe_path(
        historical_manifest_value, "historical Open Targets manifest path"
    )
    _require_within_root(
        root, historical_manifest_path, "historical Open Targets manifest"
    )
    _verify_file(historical_manifest_path, historical_manifest_sha256)
    if value.get("outputs") != _projection_outputs(projection_path, contract):
        raise RareDiseaseError(
            "historical rare-disease bundle output identities differ"
        )
    sources = value.get("sources")
    if not isinstance(sources, list) or len(sources) != len(contract.artifacts):
        raise RareDiseaseError(
            "historical rare-disease bundle source identities are invalid"
        )
    for artifact, source in zip(contract.artifacts, sources, strict=True):
        if not isinstance(source, Mapping):
            raise RareDiseaseError(
                "historical rare-disease bundle source identity is invalid"
            )
        expected_source = {
            "role": artifact.role,
            "resource_id": artifact.resource_id,
            "release": artifact.release,
            "extraction_date": artifact.extraction_date,
            "source_revision": artifact.source_revision,
            "source_url": artifact.url,
            "source_sha256": artifact.sha256,
            "license": artifact.license,
            "license_url": artifact.license_url,
        }
        if any(source.get(key) != expected for key, expected in expected_source.items()):
            raise RareDiseaseError(
                f"historical rare-disease source identity differs: {artifact.role}"
            )
        manifest_path_value = source.get("manifest_path")
        if not isinstance(manifest_path_value, str):
            raise RareDiseaseError(
                f"historical rare-disease manifest path is invalid: {artifact.role}"
            )
        manifest_relative = _safe_path(
            manifest_path_value, f"historical {artifact.role} manifest path"
        )
        manifest_sha256 = source.get("manifest_sha256")
        if not isinstance(manifest_sha256, str) or not SHA256.fullmatch(
            manifest_sha256
        ):
            raise RareDiseaseError(
                f"historical rare-disease manifest digest is invalid: {artifact.role}"
            )
        manifest_path = root / manifest_relative
        _require_within_root(
            root, manifest_path, f"historical {artifact.role} manifest"
        )
        _verify_file(manifest_path, manifest_sha256)


def bootstrap_rare_disease_evidence(
    *,
    contract_path: str | Path,
    reference_root: str | Path,
    image_digest: str,
    retrieved_on: str | None = None,
    scratch_root: str | Path | None = None,
    reserve_bytes: int = DEFAULT_RESERVE_BYTES,
    downloader: Downloader = download_atomic,
    capacity_checker: CapacityChecker = require_capacity,
) -> RareDiseaseBootstrapResult:
    contract_file = Path(contract_path)
    root = Path(reference_root).resolve(strict=True)
    contract = load_rare_disease_contract(contract_file)
    if not DIGEST.fullmatch(image_digest):
        raise RareDiseaseError("rare-disease image_digest must be a SHA-256 digest")
    contract_sha256 = _file_identity(contract_file)[0]
    retrieval_date = retrieved_on or date.today().isoformat()
    try:
        date.fromisoformat(retrieval_date)
    except ValueError as error:
        raise RareDiseaseError("rare-disease retrieved_on must be an ISO date") from error
    hpo_ontology = root / contract.hpo.ontology_path
    hpo_annotations = root / contract.hpo.annotations_path
    _require_within_root(root, hpo_ontology, "HPO ontology input")
    _require_within_root(root, hpo_annotations, "HPO annotation input")
    _verify_file(hpo_ontology, contract.hpo.ontology_sha256)
    _verify_file(hpo_annotations, contract.hpo.annotations_sha256)
    open_targets_files, open_targets_manifest_sha256 = _open_targets_files(
        root, contract.open_targets
    )
    projection_path = root / contract.projection.output_prefix
    legacy_bundle_path = _bundle_manifest_path(root, contract)
    bundle_path = _bundle_manifest_path(root, contract, contract_sha256)
    _require_within_root(root, projection_path, "rare-disease projection")
    _require_within_root(root, bundle_path, "rare-disease bundle manifest")

    if bundle_path.exists() or bundle_path.is_symlink():
        if bundle_path.is_symlink() or not bundle_path.is_file():
            raise RareDiseaseError("rare-disease bundle manifest is not a regular file")
        for artifact in contract.artifacts:
            _verify_file(
                root / artifact.path, artifact.sha256, artifact.size_bytes
            )
        publications = _source_publications(
            root=root,
            contract=contract,
            image_digest=image_digest,
            retrieved_on=retrieval_date,
            publish_missing=False,
        )
        validation = _validate_projection(
            projection_path, contract, open_targets_files
        )
        _verify_bundle_manifest(
            value=_load_bundle_manifest(bundle_path),
            root=root,
            contract=contract,
            contract_sha256=contract_sha256,
            source_publications=publications,
            open_targets_manifest_sha256=open_targets_manifest_sha256,
            projection_path=projection_path,
            validation=validation,
        )
        return RareDiseaseBootstrapResult(
            bundle_id=contract.bundle_id,
            projection_id=contract.projection_id,
            bundle_manifest_path=str(bundle_path),
            bundle_manifest_sha256=_file_identity(bundle_path)[0],
            source_manifest_paths=tuple(str(path) for path, _digest in publications),
            source_manifest_sha256=tuple(digest for _path, digest in publications),
            validation=validation,
            reused=True,
        )

    if projection_path.exists() or projection_path.is_symlink():
        if (
            legacy_bundle_path == bundle_path
            or not legacy_bundle_path.exists()
            or legacy_bundle_path.is_symlink()
            or not legacy_bundle_path.is_file()
        ):
            raise RareDiseaseError(
                "rare-disease projection exists without its immutable completion manifest"
            )
        for artifact in contract.artifacts:
            _verify_file(root / artifact.path, artifact.sha256, artifact.size_bytes)
        validation = _validate_projection(
            projection_path, contract, open_targets_files
        )
        _verify_historical_bundle_for_projection(
            value=_load_bundle_manifest(legacy_bundle_path),
            root=root,
            contract=contract,
            projection_path=projection_path,
            validation=validation,
        )
        publications = _source_publications(
            root=root,
            contract=contract,
            image_digest=image_digest,
            retrieved_on=retrieval_date,
        )
        bundle_sha256 = _publish_json(
            _bundle_manifest_data(
                root=root,
                contract=contract,
                contract_sha256=contract_sha256,
                image_digest=image_digest,
                retrieved_on=retrieval_date,
                source_publications=publications,
                open_targets_manifest_sha256=open_targets_manifest_sha256,
                projection_path=projection_path,
                validation=validation,
            ),
            bundle_path,
            root,
        )
        return RareDiseaseBootstrapResult(
            bundle_id=contract.bundle_id,
            projection_id=contract.projection_id,
            bundle_manifest_path=str(bundle_path),
            bundle_manifest_sha256=bundle_sha256,
            source_manifest_paths=tuple(
                str(path) for path, _digest in publications
            ),
            source_manifest_sha256=tuple(
                digest for _path, digest in publications
            ),
            validation=validation,
            reused=True,
        )

    scratch_base = (
        Path(scratch_root) if scratch_root is not None else root / "scratch/public-build"
    )
    if not scratch_base.is_absolute():
        raise RareDiseaseError("rare-disease scratch_root must be absolute")
    _require_within_root(root, scratch_base, "rare-disease scratch_root")
    source_scratch = scratch_base / "rare-disease-evidence" / contract.bundle_id / "sources"
    _require_within_root(root, source_scratch, "rare-disease source scratch")
    source_scratch.mkdir(parents=True, exist_ok=True, mode=0o700)
    source_paths: dict[str, Path] = {}
    missing: list[tuple[RareDiseaseArtifact, Path]] = []
    for artifact in contract.artifacts:
        final = root / artifact.path
        _require_within_root(root, final, f"{artifact.role} release")
        if final.exists() or final.is_symlink():
            _verify_file(final, artifact.sha256, artifact.size_bytes)
            source_paths[artifact.role] = final
            continue
        staged = source_scratch / artifact.role
        _require_within_root(root, staged, f"{artifact.role} source scratch")
        source_paths[artifact.role] = staged
        if staged.exists() or staged.is_symlink():
            _verify_file(staged, artifact.sha256, artifact.size_bytes)
        else:
            missing.append((artifact, staged))
    missing_bytes = sum(artifact.size_bytes for artifact, _path in missing)
    capacity_checker(
        root,
        CapacityPlan(
            operation=f"bootstrap rare-disease evidence {contract.bundle_id}",
            source_bytes=missing_bytes,
            output_bytes=missing_bytes + contract.projected_output_bytes,
            scratch_bytes=missing_bytes + contract.projected_output_bytes,
            reserve_bytes=reserve_bytes,
        ),
    )
    for artifact, staged in missing:
        downloader(
            DownloadRequest(
                url=artifact.url,
                destination=staged,
                expected_size_bytes=artifact.size_bytes,
                expected_sha256=artifact.sha256,
                attempts=3,
                timeout_seconds=300,
            ),
            root,
        )
    for artifact in contract.artifacts:
        _verify_file(
            source_paths[artifact.role], artifact.sha256, artifact.size_bytes
        )

    build_root = scratch_base / "rare-disease-evidence" / contract.bundle_id
    staging = build_root / f".{contract.projection_id}.{uuid.uuid4().hex}.partial"
    _require_within_root(root, staging, "rare-disease projection staging")
    staging.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    _build_projection(
        contract=contract,
        source_paths=source_paths,
        hpo_ontology_path=hpo_ontology,
        hpo_annotations_path=hpo_annotations,
        open_targets_files=open_targets_files,
        destination=staging,
    )
    validation = _validate_projection(staging, contract, open_targets_files)
    os.chmod(staging, 0o755)
    projection_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        os.replace(staging, projection_path)
    except OSError:
        shutil.rmtree(staging, ignore_errors=True)
        raise RareDiseaseError(
            "rare-disease projection publication collided with another builder"
        )
    directory_fd = os.open(projection_path.parent, os.O_RDONLY)
    try:
        os.fsync(directory_fd)
    finally:
        os.close(directory_fd)

    for artifact in contract.artifacts:
        _publish_source_file(root, artifact, source_paths[artifact.role])
        source_paths[artifact.role] = root / artifact.path
    publications = _source_publications(
        root=root,
        contract=contract,
        image_digest=image_digest,
        retrieved_on=retrieval_date,
    )
    bundle_existed = bundle_path.exists() or bundle_path.is_symlink()
    if bundle_existed:
        if bundle_path.is_symlink() or not bundle_path.is_file():
            raise RareDiseaseError("rare-disease bundle manifest is not a regular file")
        _verify_bundle_manifest(
            value=_load_bundle_manifest(bundle_path),
            root=root,
            contract=contract,
            contract_sha256=contract_sha256,
            source_publications=publications,
            open_targets_manifest_sha256=open_targets_manifest_sha256,
            projection_path=projection_path,
            validation=validation,
        )
        bundle_sha256 = _file_identity(bundle_path)[0]
    else:
        bundle_data = _bundle_manifest_data(
            root=root,
            contract=contract,
            contract_sha256=contract_sha256,
            image_digest=image_digest,
            retrieved_on=retrieval_date,
            source_publications=publications,
            open_targets_manifest_sha256=open_targets_manifest_sha256,
            projection_path=projection_path,
            validation=validation,
        )
        bundle_sha256 = _publish_json(bundle_data, bundle_path, root)
    for artifact in contract.artifacts:
        staged = source_scratch / artifact.role
        if staged.exists() and staged != source_paths[artifact.role]:
            staged.unlink()
    return RareDiseaseBootstrapResult(
        bundle_id=contract.bundle_id,
        projection_id=contract.projection_id,
        bundle_manifest_path=str(bundle_path),
        bundle_manifest_sha256=bundle_sha256,
        source_manifest_paths=tuple(str(path) for path, _digest in publications),
        source_manifest_sha256=tuple(digest for _path, digest in publications),
        validation=validation,
        reused=bundle_existed,
    )


def rare_disease_result_as_json(result: RareDiseaseBootstrapResult) -> str:
    return json.dumps(asdict(result), indent=2, sort_keys=True) + "\n"
