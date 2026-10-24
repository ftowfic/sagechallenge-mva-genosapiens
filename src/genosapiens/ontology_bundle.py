"""Bootstrap and validate the pinned core ontology reference bundle."""

from __future__ import annotations

from collections.abc import Callable, Mapping
import csv
from dataclasses import asdict, dataclass, field
from datetime import date
import hashlib
import json
import os
import stat
from pathlib import Path, PurePosixPath
import re
import shutil
from types import MappingProxyType
from typing import Any
import uuid

from genosapiens.download import DownloadRequest, DownloadResult, download_atomic
from genosapiens.manifest import (
    ManifestError,
    PublicResourceManifest,
    ResourceOutput,
    RetentionPolicy,
    SourceIdentity,
    Transformation,
    publish_manifest,
    verify_outputs,
)
from genosapiens.storage import DEFAULT_RESERVE_BYTES, CapacityPlan, require_capacity


SCHEMA = "genosapiens.core-ontology-bundle"
TRANSFORMATION_SCHEMA = "genosapiens.core-ontology-bootstrap"
OBSOLETE_RESOLUTION_POLICY = "genosapiens.hpo-obsolete-resolution"
SHA256 = re.compile(r"^[0-9a-f]{64}$")
IDENTIFIER = re.compile(r"^[a-z0-9][a-z0-9._-]*$")
EXPECTED_ROLES = {
    "hpo_ontology",
    "hpo_disease_annotations",
    "mondo_efo_mappings",
    "sequence_ontology",
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
SSSOM_COLUMNS = (
    "subject_id",
    "subject_label",
    "predicate_id",
    "object_id",
    "mapping_justification",
)


class OntologyBundleError(RuntimeError):
    """Raised when ontology acquisition, provenance, or validation fails."""


@dataclass(frozen=True, slots=True)
class OntologyArtifact:
    role: str
    resource_id: str
    release: str
    url: str
    path: str
    size_bytes: int
    sha256: str
    license: str
    license_url: str

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "OntologyArtifact":
        expected = {
            "role",
            "resource_id",
            "release",
            "url",
            "path",
            "size_bytes",
            "sha256",
            "license",
            "license_url",
        }
        if set(data) != expected:
            raise OntologyBundleError("ontology artifact fields are invalid")
        try:
            artifact = cls(**data)
        except TypeError as error:
            raise OntologyBundleError(f"invalid ontology artifact: {error}") from error
        string_values = (
            artifact.role,
            artifact.resource_id,
            artifact.release,
            artifact.url,
            artifact.path,
            artifact.sha256,
            artifact.license,
            artifact.license_url,
        )
        if not all(isinstance(value, str) and value for value in string_values):
            raise OntologyBundleError("ontology artifact strings are invalid")
        path = PurePosixPath(artifact.path)
        if (
            not IDENTIFIER.fullmatch(artifact.resource_id)
            or not artifact.url.startswith("https://")
            or not artifact.license_url.startswith("https://")
            or path.is_absolute()
            or ".." in path.parts
            or path == PurePosixPath(".")
            or path.parts[:2] != ("releases", "ontologies")
            or not isinstance(artifact.size_bytes, int)
            or isinstance(artifact.size_bytes, bool)
            or artifact.size_bytes <= 0
            or not SHA256.fullmatch(artifact.sha256)
        ):
            raise OntologyBundleError(f"invalid ontology artifact: {artifact.role}")
        return artifact


@dataclass(frozen=True, slots=True)
class OntologyValidationContract:
    artifact_role: str
    ontology_id: str
    data_version: str
    identifier_prefix: str
    root_id: str
    sentinel_term_id: str
    sentinel_parent_id: str
    obsolete_term_id: str
    replaced_by: str

    @classmethod
    def from_mapping(
        cls, data: Mapping[str, Any]
    ) -> "OntologyValidationContract":
        expected = {
            "artifact_role",
            "ontology_id",
            "data_version",
            "identifier_prefix",
            "root_id",
            "sentinel_term_id",
            "sentinel_parent_id",
            "obsolete_term_id",
            "replaced_by",
        }
        if set(data) != expected:
            raise OntologyBundleError("ontology validation fields are invalid")
        try:
            value = cls(**data)
        except TypeError as error:
            raise OntologyBundleError(f"invalid ontology validation: {error}") from error
        if not all(
            isinstance(item, str) and item
            for item in asdict(value).values()
        ):
            raise OntologyBundleError("ontology validation strings must be non-empty")
        identifier = re.compile(rf"^{re.escape(value.identifier_prefix)}[0-9]{{7}}$")
        for term_id in (
            value.root_id,
            value.sentinel_term_id,
            value.sentinel_parent_id,
            value.obsolete_term_id,
            value.replaced_by,
        ):
            if not identifier.fullmatch(term_id):
                raise OntologyBundleError(
                    f"invalid {value.identifier_prefix} ontology sentinel: {term_id}"
                )
        return value


@dataclass(frozen=True, slots=True)
class HpoAnnotationValidationContract:
    artifact_role: str
    version: str
    hpo_version: str
    required_database_prefixes: tuple[str, ...]

    @classmethod
    def from_mapping(
        cls, data: Mapping[str, Any]
    ) -> "HpoAnnotationValidationContract":
        expected = {
            "artifact_role",
            "version",
            "hpo_version",
            "required_database_prefixes",
        }
        if set(data) != expected:
            raise OntologyBundleError("HPO annotation validation fields are invalid")
        prefixes = data.get("required_database_prefixes")
        if not isinstance(prefixes, list) or not all(
            isinstance(prefix, str) and prefix for prefix in prefixes
        ):
            raise OntologyBundleError("HPO annotation database prefixes are invalid")
        try:
            value = cls(
                artifact_role=data["artifact_role"],
                version=data["version"],
                hpo_version=data["hpo_version"],
                required_database_prefixes=tuple(prefixes),
            )
        except (KeyError, TypeError) as error:
            raise OntologyBundleError(
                f"invalid HPO annotation validation: {error}"
            ) from error
        if (
            not all(
                isinstance(item, str) and item
                for item in (value.artifact_role, value.version, value.hpo_version)
            )
            or not value.hpo_version.startswith("http://purl.obolibrary.org/obo/hp/")
            or len(value.required_database_prefixes)
            != len(set(value.required_database_prefixes))
        ):
            raise OntologyBundleError("HPO annotation validation is invalid")
        return value


@dataclass(frozen=True, slots=True)
class MappingValidationContract:
    artifact_role: str
    license_url: str
    mapping_set_id: str
    mapping_set_version: str
    subject_source: str
    object_source: str

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "MappingValidationContract":
        expected = {
            "artifact_role",
            "license_url",
            "mapping_set_id",
            "mapping_set_version",
            "subject_source",
            "object_source",
        }
        if set(data) != expected:
            raise OntologyBundleError("mapping validation fields are invalid")
        try:
            value = cls(**data)
        except TypeError as error:
            raise OntologyBundleError(f"invalid mapping validation: {error}") from error
        if not all(
            isinstance(item, str) and item
            for item in asdict(value).values()
        ) or not value.license_url.startswith("https://"):
            raise OntologyBundleError("mapping validation strings are invalid")
        return value


@dataclass(frozen=True, slots=True)
class InformationContentContract:
    source_artifact_role: str
    included_aspects: tuple[str, ...]
    excluded_qualifiers: tuple[str, ...]

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "InformationContentContract":
        expected = {
            "source_artifact_role",
            "included_aspects",
            "excluded_qualifiers",
        }
        if set(data) != expected:
            raise OntologyBundleError("information-content corpus fields are invalid")
        aspects = data.get("included_aspects")
        qualifiers = data.get("excluded_qualifiers")
        if (
            not isinstance(aspects, list)
            or not aspects
            or not all(isinstance(item, str) and item for item in aspects)
            or not isinstance(qualifiers, list)
            or not qualifiers
            or not all(isinstance(item, str) and item for item in qualifiers)
        ):
            raise OntologyBundleError("information-content corpus selectors are invalid")
        source_role = data["source_artifact_role"]
        if not isinstance(source_role, str) or not source_role:
            raise OntologyBundleError("information-content source role is invalid")
        return cls(source_role, tuple(aspects), tuple(qualifiers))


@dataclass(frozen=True, slots=True)
class OntologyBundleContract:
    bundle_id: str
    artifacts: tuple[OntologyArtifact, ...]
    hpo_ontology: OntologyValidationContract
    hpo_annotations: HpoAnnotationValidationContract
    mondo_efo_mappings: MappingValidationContract
    sequence_ontology: OntologyValidationContract
    information_content_corpus: InformationContentContract
    schema: str = SCHEMA

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "OntologyBundleContract":
        expected = {
            "schema",
            "bundle_id",
            "artifacts",
            "hpo_ontology",
            "hpo_annotations",
            "mondo_efo_mappings",
            "sequence_ontology",
            "information_content_corpus",
        }
        if set(data) != expected or data.get("schema") != SCHEMA:
            raise OntologyBundleError("ontology bundle fields or schema are invalid")
        artifacts = data.get("artifacts")
        nested = (
            data.get("hpo_ontology"),
            data.get("hpo_annotations"),
            data.get("mondo_efo_mappings"),
            data.get("sequence_ontology"),
            data.get("information_content_corpus"),
        )
        if (
            not isinstance(artifacts, list)
            or not all(isinstance(item, Mapping) for item in artifacts)
            or not all(isinstance(item, Mapping) for item in nested)
        ):
            raise OntologyBundleError("ontology bundle nested fields are invalid")
        try:
            contract = cls(
                schema=data["schema"],
                bundle_id=data["bundle_id"],
                artifacts=tuple(
                    OntologyArtifact.from_mapping(item) for item in artifacts
                ),
                hpo_ontology=OntologyValidationContract.from_mapping(
                    data["hpo_ontology"]
                ),
                hpo_annotations=HpoAnnotationValidationContract.from_mapping(
                    data["hpo_annotations"]
                ),
                mondo_efo_mappings=MappingValidationContract.from_mapping(
                    data["mondo_efo_mappings"]
                ),
                sequence_ontology=OntologyValidationContract.from_mapping(
                    data["sequence_ontology"]
                ),
                information_content_corpus=InformationContentContract.from_mapping(
                    data["information_content_corpus"]
                ),
            )
        except (KeyError, TypeError) as error:
            raise OntologyBundleError(f"invalid ontology bundle: {error}") from error
        if (
            not isinstance(contract.bundle_id, str)
            or not IDENTIFIER.fullmatch(contract.bundle_id)
        ):
            raise OntologyBundleError("ontology bundle_id is invalid")
        release_root = PurePosixPath("releases/ontologies") / contract.bundle_id
        if any(
            not PurePosixPath(artifact.path).is_relative_to(release_root)
            or PurePosixPath(artifact.path) == release_root
            for artifact in contract.artifacts
        ):
            raise OntologyBundleError(
                "ontology artifacts must be inside the atomic bundle directory"
            )
        roles = {artifact.role for artifact in contract.artifacts}
        if roles != EXPECTED_ROLES or len(roles) != len(contract.artifacts):
            raise OntologyBundleError("ontology bundle must define every artifact once")
        resource_ids = {artifact.resource_id for artifact in contract.artifacts}
        paths = {artifact.path for artifact in contract.artifacts}
        if len(resource_ids) != len(contract.artifacts) or len(paths) != len(
            contract.artifacts
        ):
            raise OntologyBundleError("ontology resource IDs and paths must be unique")
        validation_roles = {
            contract.hpo_ontology.artifact_role,
            contract.hpo_annotations.artifact_role,
            contract.mondo_efo_mappings.artifact_role,
            contract.sequence_ontology.artifact_role,
        }
        if validation_roles != EXPECTED_ROLES:
            raise OntologyBundleError("ontology validation roles do not match artifacts")
        if (
            contract.information_content_corpus.source_artifact_role
            != "hpo_disease_annotations"
            or contract.information_content_corpus.included_aspects != ("P",)
            or contract.information_content_corpus.excluded_qualifiers != ("NOT",)
        ):
            raise OntologyBundleError(
                "information-content corpus must use non-negated HPO phenotype rows"
            )
        if (
            contract.artifact("hpo_ontology").release
            != contract.artifact("hpo_disease_annotations").release
        ):
            raise OntologyBundleError("HPO ontology and annotations must share a release")
        return contract

    def artifact(self, role: str) -> OntologyArtifact:
        for artifact in self.artifacts:
            if artifact.role == role:
                return artifact
        raise OntologyBundleError(f"ontology contract lacks artifact role {role}")


def load_ontology_bundle_contract(path: str | Path) -> OntologyBundleContract:
    contract_path = Path(path)
    try:
        data = json.loads(contract_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise OntologyBundleError(
            f"cannot read ontology contract {contract_path}: {error}"
        ) from error
    if not isinstance(data, Mapping):
        raise OntologyBundleError("ontology contract root must be an object")
    return OntologyBundleContract.from_mapping(data)


@dataclass(frozen=True, slots=True)
class OntologyTerm:
    term_id: str
    name: str
    parents: tuple[str, ...]
    alt_ids: tuple[str, ...]
    obsolete: bool
    replaced_by: tuple[str, ...]
    consider: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class OntologyResolution:
    """Resolution under OBSOLETE_RESOLUTION_POLICY, never a guessed canonical ID.

    Chain includes the source, alias expansions, and attempted replacement IDs,
    including the repeated ID on a cycle. Mapping fields are sorted, distinct
    statements from encountered terms (all candidates at a conflicting ID).
    Original statements and their owning terms remain in ParsedOntology.terms.
    """

    source_id: str
    canonical_id: str | None
    status: str
    chain: tuple[str, ...]
    replaced_by: tuple[str, ...]
    consider: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class ParsedOntology:
    terms: Mapping[str, OntologyTerm]
    active_aliases: Mapping[str, str]
    parent_edges: int
    obsolete_mappings: int
    _alias_targets: Mapping[str, tuple[str, ...]] = field(
        init=False, repr=False, compare=False
    )

    def __post_init__(self) -> None:
        targets: dict[str, set[str]] = {}
        for term in self.terms.values():
            for alias in term.alt_ids:
                targets.setdefault(alias, set()).add(term.term_id)
        object.__setattr__(
            self,
            "_alias_targets",
            MappingProxyType({alias: tuple(sorted(owners)) for alias, owners in targets.items()}),
        )

    def resolve_active(self, term_id: str) -> str | None:
        term = self.terms.get(term_id)
        if term is not None and not term.obsolete:
            return term_id
        return self.active_aliases.get(term_id)

    def resolve_term(self, term_id: str) -> OntologyResolution:
        """Follow only unique aliases and unique explicit replaced_by edges.

        Success statuses: active, alias, replaced_by. Failure statuses: unknown,
        ambiguous_identifier, ambiguous_replacement, consider_only,
        no_replacement, missing_target, cycle. Only success has a canonical ID.
        Distinct replacement IDs are ambiguous even if they would converge.
        """
        current = term_id
        chain: list[str] = []
        seen: set[str] = set()
        replacements: set[str] = set()
        suggestions: set[str] = set()
        followed_replacement = False
        canonical_id = None
        while True:
            chain.append(current)
            if current in seen:
                status = "cycle"
                break
            seen.add(current)
            owners = set(self._alias_targets.get(current, ()))
            if current in self.terms:
                owners.add(current)
            if not owners:
                status = "missing_target" if followed_replacement else "unknown"
                break
            for owner in owners:
                term = self.terms[owner]
                replacements.update(term.replaced_by)
                suggestions.update(term.consider)
            if len(owners) != 1:
                status = "ambiguous_identifier"
                break
            owner = next(iter(owners))
            if owner != current:
                current = owner
                continue
            term = self.terms[current]
            if not term.obsolete:
                canonical_id = current
                status = "alias" if current != term_id else "active"
                if followed_replacement:
                    status = "replaced_by"
                break
            targets = set(term.replaced_by)
            if len(targets) > 1:
                status = "ambiguous_replacement"
                break
            if not targets:
                status = "consider_only" if term.consider else "no_replacement"
                break
            current = next(iter(targets))
            followed_replacement = True
        return OntologyResolution(
            source_id=term_id,
            canonical_id=canonical_id,
            status=status,
            chain=tuple(chain),
            replaced_by=tuple(sorted(replacements)),
            consider=tuple(sorted(suggestions)),
        )


@dataclass(frozen=True, slots=True)
class OntologyBundleValidation:
    hpo_terms: int
    hpo_parent_edges: int
    hpo_obsolete_mappings: int
    hpo_annotation_records: int
    hpo_diseases: int
    information_content_records: int
    mondo_efo_mappings: int
    sequence_ontology_terms: int
    sequence_ontology_parent_edges: int
    sequence_ontology_obsolete_mappings: int


@dataclass(frozen=True, slots=True)
class OntologyBootstrapResult:
    bundle_id: str
    bundle_manifest_path: str
    bundle_manifest_sha256: str
    manifest_paths: tuple[str, ...]
    manifest_sha256: tuple[str, ...]
    validation: OntologyBundleValidation
    reused: bool


Downloader = Callable[[DownloadRequest, str | Path], DownloadResult]
CapacityChecker = Callable[[str | Path, CapacityPlan], Any]


def _file_signature(path: Path) -> tuple[int, int, int, int, int]:
    value = path.lstat()
    if not stat.S_ISREG(value.st_mode):
        raise OntologyBundleError("resource is not a regular file")
    return value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns


def _file_identity(path: Path) -> tuple[str, int]:
    digest = hashlib.sha256()
    size = 0
    descriptor = -1
    try:
        descriptor = os.open(path, os.O_RDONLY | os.O_NONBLOCK | os.O_NOFOLLOW)
        if not stat.S_ISREG(os.fstat(descriptor).st_mode):
            raise OntologyBundleError("ontology artifact is not a regular file")
        with os.fdopen(descriptor, "rb") as source:
            descriptor = -1
            while chunk := source.read(1024 * 1024):
                digest.update(chunk)
                size += len(chunk)
    except OSError as error:
        raise OntologyBundleError(f"cannot read ontology artifact {path}: {error}") from error
    finally:
        if descriptor >= 0:
            os.close(descriptor)
    return digest.hexdigest(), size


def _require_within_root(root: Path, path: Path, description: str) -> None:
    if not path.is_absolute() or not path.is_relative_to(root):
        raise OntologyBundleError(f"{description} escapes reference_root: {path}")
    relative = path.relative_to(root)
    if ".." in relative.parts:
        raise OntologyBundleError(f"{description} escapes reference_root: {path}")
    current = root
    for part in relative.parts:
        current /= part
        if current.is_symlink():
            raise OntologyBundleError(f"{description} contains a symlink: {current}")
    try:
        resolved = path.resolve()
    except OSError as error:
        raise OntologyBundleError(f"cannot resolve {description} {path}: {error}") from error
    if not resolved.is_relative_to(root):
        raise OntologyBundleError(f"{description} escapes reference_root: {path}")


def _verify_artifact(path: Path, artifact: OntologyArtifact) -> None:
    digest, size = _file_identity(path)
    if size != artifact.size_bytes or digest != artifact.sha256:
        raise OntologyBundleError(
            f"ontology release checksum mismatch: {artifact.role}"
        )


def _obo_identifier(value: str) -> str:
    return value.split(" !", 1)[0].split(" {", 1)[0].strip()


def _parse_obo(
    path: Path,
    validation: OntologyValidationContract,
    *,
    retain_unresolved_obsolete: bool = False,
) -> ParsedOntology:
    """Validate the pinned ontology, optionally retaining unresolved mappings.

    Retention checks mapping/alias syntax without requiring obsolete targets to
    be immediately active. It does not relax source, sentinel, active-alias, or
    is_a graph validation, or rewrite any terms or mapping statements.
    """
    header: dict[str, list[str]] = {}
    terms: dict[str, OntologyTerm] = {}
    current: dict[str, list[str]] | None = None

    def finish_term() -> None:
        nonlocal current
        if current is None:
            return
        term_ids = current.get("id", [])
        names = current.get("name", [])
        if len(term_ids) != 1 or len(names) != 1:
            raise OntologyBundleError(f"ontology term lacks one ID and name: {path}")
        term_id = term_ids[0]
        if term_id in terms:
            raise OntologyBundleError(f"duplicate ontology term ID: {term_id}")
        terms[term_id] = OntologyTerm(
            term_id=term_id,
            name=names[0],
            parents=tuple(_obo_identifier(item) for item in current.get("is_a", [])),
            alt_ids=tuple(current.get("alt_id", [])),
            obsolete=current.get("is_obsolete", []) == ["true"],
            replaced_by=tuple(
                _obo_identifier(item) for item in current.get("replaced_by", [])
            ),
            consider=tuple(
                _obo_identifier(item) for item in current.get("consider", [])
            ),
        )
        current = None

    try:
        with path.open("rt", encoding="utf-8") as source:
            for raw_line in source:
                line = raw_line.rstrip("\r\n")
                if line == "[Term]":
                    finish_term()
                    current = {}
                    continue
                if line.startswith("[") and line.endswith("]"):
                    finish_term()
                    current = None
                    continue
                if ": " not in line:
                    continue
                key, value = line.split(": ", 1)
                target = current if current is not None else header
                target.setdefault(key, []).append(value)
            finish_term()
    except (OSError, UnicodeDecodeError) as error:
        raise OntologyBundleError(f"cannot parse ontology {path}: {error}") from error
    if header.get("format-version") != ["1.2"]:
        raise OntologyBundleError(f"ontology format-version is not 1.2: {path}")
    if header.get("ontology") != [validation.ontology_id]:
        raise OntologyBundleError(f"ontology identifier mismatch: {path}")
    if header.get("data-version") != [validation.data_version]:
        raise OntologyBundleError(f"ontology data version mismatch: {path}")
    if not terms:
        raise OntologyBundleError(f"ontology has no terms: {path}")

    identifier = re.compile(
        rf"^{re.escape(validation.identifier_prefix)}[0-9]{{7}}$"
    )
    if any(not identifier.fullmatch(term_id) for term_id in terms):
        raise OntologyBundleError(f"ontology contains an invalid term identifier: {path}")
    if retain_unresolved_obsolete:
        for term in terms.values():
            for target in term.alt_ids + term.replaced_by + term.consider + term.parents:
                if (
                    not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*:[A-Za-z0-9._-]+", target)
                    or (
                        target.startswith(validation.identifier_prefix)
                        and not identifier.fullmatch(target)
                    )
                ):
                    raise OntologyBundleError(
                        f"ontology contains an invalid mapping/alias identifier: {target}"
                    )
    active = {term_id for term_id, term in terms.items() if not term.obsolete}
    aliases: dict[str, str] = {}
    for term_id in active:
        for alias in terms[term_id].alt_ids:
            if (
                not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*:[A-Za-z0-9._-]+", alias)
                or alias in active
                or (alias in aliases and aliases[alias] != term_id)
            ):
                raise OntologyBundleError(f"ontology has an invalid active alias: {alias}")
            aliases[alias] = term_id

    graph: dict[str, tuple[str, ...]] = {}
    edge_count = 0
    for term_id in active:
        resolved_parents = []
        for parent in terms[term_id].parents:
            resolved = parent if parent in active else aliases.get(parent)
            if resolved is None:
                raise OntologyBundleError(
                    f"ontology parent does not resolve: {term_id} -> {parent}"
                )
            resolved_parents.append(resolved)
        graph[term_id] = tuple(resolved_parents)
        edge_count += len(resolved_parents)

    state: dict[str, int] = {}

    def visit(term_id: str) -> None:
        if state.get(term_id) == 1:
            raise OntologyBundleError(f"ontology parent graph contains a cycle: {term_id}")
        if state.get(term_id) == 2:
            return
        state[term_id] = 1
        for parent in graph[term_id]:
            visit(parent)
        state[term_id] = 2

    for term_id in graph:
        visit(term_id)

    mapping_count = 0
    for term in terms.values():
        if not term.obsolete:
            continue
        for target in term.replaced_by + term.consider:
            if (
                not retain_unresolved_obsolete
                and target not in active
                and target not in aliases
            ):
                raise OntologyBundleError(
                    f"obsolete ontology mapping does not resolve: {term.term_id} -> {target}"
                )
            mapping_count += 1
    root = terms.get(validation.root_id)
    sentinel = terms.get(validation.sentinel_term_id)
    obsolete = terms.get(validation.obsolete_term_id)
    if root is None or root.obsolete or root.parents:
        raise OntologyBundleError(f"ontology root sentinel is invalid: {validation.root_id}")
    if (
        sentinel is None
        or sentinel.obsolete
        or validation.sentinel_parent_id not in sentinel.parents
    ):
        raise OntologyBundleError(
            f"ontology parent sentinel is invalid: {validation.sentinel_term_id}"
        )
    if (
        obsolete is None
        or not obsolete.obsolete
        or validation.replaced_by not in obsolete.replaced_by
    ):
        raise OntologyBundleError(
            f"ontology obsolete sentinel is invalid: {validation.obsolete_term_id}"
        )
    if mapping_count == 0:
        raise OntologyBundleError(f"ontology has no obsolete-term mappings: {path}")
    return ParsedOntology(terms, aliases, edge_count, mapping_count)


def _parse_comment_table(
    path: Path,
) -> tuple[dict[str, str], csv.DictReader, Any]:
    try:
        source = path.open("rt", encoding="utf-8", newline="")
    except OSError as error:
        raise OntologyBundleError(f"cannot open ontology table {path}: {error}") from error
    metadata: dict[str, str] = {}
    try:
        while True:
            position = source.tell()
            line = source.readline()
            if not line.startswith("#"):
                source.seek(position)
                break
            content = line[1:].rstrip("\r\n")
            if content.startswith(" "):
                content = content[1:]
            if ": " in content and not content.startswith(" "):
                key, value = content.split(": ", 1)
                if key in metadata:
                    raise OntologyBundleError(f"duplicate table metadata field: {key}")
                metadata[key] = value
        return metadata, csv.DictReader(source, delimiter="\t"), source
    except Exception:
        source.close()
        raise


def _validate_hpo_annotations(
    path: Path,
    validation: HpoAnnotationValidationContract,
    corpus: InformationContentContract,
    hpo: ParsedOntology,
) -> tuple[int, int, int]:
    metadata, reader, source = _parse_comment_table(path)
    try:
        if tuple(reader.fieldnames or ()) != HPOA_COLUMNS:
            raise OntologyBundleError("HPO annotation columns are invalid")
        if metadata.get("version") != validation.version:
            raise OntologyBundleError("HPO annotation version mismatch")
        if metadata.get("hpo-version") != validation.hpo_version:
            raise OntologyBundleError("HPO annotation ontology provenance mismatch")
        if "HPO annotations for rare diseases" not in metadata.get("description", ""):
            raise OntologyBundleError("HPO annotation description provenance is missing")
        if not metadata.get("tracker", "").startswith("https://github.com/obophenotype/"):
            raise OntologyBundleError("HPO annotation tracker provenance is missing")
        records = 0
        ic_records = 0
        diseases: set[str] = set()
        prefixes: set[str] = set()
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise OntologyBundleError("HPO annotation row width is invalid")
            disease_id = row["database_id"]
            if ":" not in disease_id or not row["disease_name"]:
                raise OntologyBundleError("HPO annotation disease identity is invalid")
            prefix = disease_id.split(":", 1)[0]
            if prefix not in validation.required_database_prefixes:
                raise OntologyBundleError(
                    f"HPO annotation database prefix is not declared: {prefix}"
                )
            if hpo.resolve_active(row["hpo_id"]) is None:
                raise OntologyBundleError(
                    f"HPO annotation term does not resolve: {row['hpo_id']}"
                )
            if (
                not row["reference"]
                or not row["evidence"]
                or not row["aspect"]
                or not row["biocuration"]
                or row["qualifier"] not in {"", "NOT"}
            ):
                raise OntologyBundleError("HPO annotation row lacks direct provenance")
            records += 1
            diseases.add(disease_id)
            prefixes.add(prefix)
            if (
                row["aspect"] in corpus.included_aspects
                and row["qualifier"] not in corpus.excluded_qualifiers
            ):
                ic_records += 1
        if records == 0 or ic_records == 0:
            raise OntologyBundleError("HPO annotation corpus is empty")
        if prefixes != set(validation.required_database_prefixes):
            raise OntologyBundleError("HPO annotation database sources are incomplete")
        return records, len(diseases), ic_records
    except (csv.Error, OSError, UnicodeDecodeError) as error:
        raise OntologyBundleError(f"cannot parse HPO annotations {path}: {error}") from error
    finally:
        source.close()


def _validate_mondo_efo_mappings(
    path: Path, validation: MappingValidationContract
) -> int:
    metadata, reader, source = _parse_comment_table(path)
    try:
        if tuple(reader.fieldnames or ()) != SSSOM_COLUMNS:
            raise OntologyBundleError("MONDO/EFO mapping columns are invalid")
        expected_metadata = {
            "license": validation.license_url,
            "mapping_set_id": validation.mapping_set_id,
            "mapping_set_version": validation.mapping_set_version,
            "subject_source": validation.subject_source,
            "object_source": validation.object_source,
        }
        for key, value in expected_metadata.items():
            if metadata.get(key) != value:
                raise OntologyBundleError(f"MONDO/EFO mapping metadata mismatch: {key}")
        count = 0
        pairs: set[tuple[str, str]] = set()
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise OntologyBundleError("MONDO/EFO mapping row width is invalid")
            subject = row["subject_id"]
            object_id = row["object_id"]
            pair = (subject, object_id)
            if (
                not re.fullmatch(r"MONDO:[0-9]{7}", subject)
                or not re.fullmatch(r"EFO:[0-9]{7}", object_id)
                or row["predicate_id"] != "skos:exactMatch"
                or not row["subject_label"]
                or not row["mapping_justification"]
                or pair in pairs
            ):
                raise OntologyBundleError("MONDO/EFO exact mapping row is invalid")
            pairs.add(pair)
            count += 1
        if count == 0:
            raise OntologyBundleError("MONDO/EFO mapping set is empty")
        return count
    except (csv.Error, OSError, UnicodeDecodeError) as error:
        raise OntologyBundleError(f"cannot parse MONDO/EFO mappings {path}: {error}") from error
    finally:
        source.close()


def validate_ontology_bundle(
    contract: OntologyBundleContract,
    artifact_paths: Mapping[str, Path],
) -> OntologyBundleValidation:
    if set(artifact_paths) != EXPECTED_ROLES:
        raise OntologyBundleError("ontology validation paths are incomplete")
    for artifact in contract.artifacts:
        _verify_artifact(artifact_paths[artifact.role], artifact)
    hpo = _parse_obo(
        artifact_paths[contract.hpo_ontology.artifact_role], contract.hpo_ontology
    )
    annotations, diseases, ic_records = _validate_hpo_annotations(
        artifact_paths[contract.hpo_annotations.artifact_role],
        contract.hpo_annotations,
        contract.information_content_corpus,
        hpo,
    )
    mappings = _validate_mondo_efo_mappings(
        artifact_paths[contract.mondo_efo_mappings.artifact_role],
        contract.mondo_efo_mappings,
    )
    sequence = _parse_obo(
        artifact_paths[contract.sequence_ontology.artifact_role],
        contract.sequence_ontology,
    )
    return OntologyBundleValidation(
        hpo_terms=len(hpo.terms),
        hpo_parent_edges=hpo.parent_edges,
        hpo_obsolete_mappings=hpo.obsolete_mappings,
        hpo_annotation_records=annotations,
        hpo_diseases=diseases,
        information_content_records=ic_records,
        mondo_efo_mappings=mappings,
        sequence_ontology_terms=len(sequence.terms),
        sequence_ontology_parent_edges=sequence.parent_edges,
        sequence_ontology_obsolete_mappings=sequence.obsolete_mappings,
    )


def _manifest_relative_path(artifact: OntologyArtifact) -> str:
    return f"manifests/ontologies/{artifact.resource_id}.json"


def _load_manifest(path: Path) -> PublicResourceManifest:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise OntologyBundleError(f"cannot read ontology manifest {path}: {error}") from error
    if not isinstance(data, Mapping):
        raise OntologyBundleError("ontology manifest root must be an object")
    try:
        return PublicResourceManifest.from_mapping(data)
    except ManifestError as error:
        raise OntologyBundleError(f"invalid ontology manifest: {error}") from error


def _verify_manifest_contract(
    manifest: PublicResourceManifest,
    artifact: OntologyArtifact,
    contract_sha256: str,
) -> None:
    if (
        manifest.resource_id != artifact.resource_id
        or manifest.release != artifact.release
        or manifest.source != SourceIdentity(artifact.url, sha256=artifact.sha256)
        or manifest.license != artifact.license
        or manifest.assembly != "not-applicable"
        or manifest.transformation
        != Transformation(TRANSFORMATION_SCHEMA, contract_sha256)
        or len(manifest.outputs) != 1
        or manifest.outputs[0].path != artifact.path
        or manifest.outputs[0].size_bytes != artifact.size_bytes
        or manifest.outputs[0].sha256 != artifact.sha256
    ):
        raise OntologyBundleError(
            f"ontology manifest identity differs from contract: {artifact.role}"
        )


def _copy_file(source: Path, destination: Path) -> None:
    try:
        with source.open("rb") as input_file, destination.open("xb") as output_file:
            shutil.copyfileobj(input_file, output_file, length=1024 * 1024)
            output_file.flush()
            os.fsync(output_file.fileno())
    except OSError as error:
        raise OntologyBundleError(f"cannot stage ontology artifact {source}: {error}") from error


def _bundle_release_relative(contract: OntologyBundleContract) -> PurePosixPath:
    return PurePosixPath("releases/ontologies") / contract.bundle_id


def _validate_release_tree(root: Path, contract: OntologyBundleContract) -> None:
    relative_root = _bundle_release_relative(contract)
    target = root / str(relative_root)
    _require_within_root(root, target, "ontology release")
    if not target.is_dir() or target.is_symlink():
        raise OntologyBundleError(f"ontology release is not a directory: {target}")
    expected_files: set[str] = set()
    expected_directories: set[str] = set()
    for artifact in contract.artifacts:
        relative = PurePosixPath(artifact.path).relative_to(relative_root)
        expected_files.add(str(relative))
        parent = relative.parent
        while parent != PurePosixPath("."):
            expected_directories.add(str(parent))
            parent = parent.parent
    actual_files: set[str] = set()
    actual_directories: set[str] = set()
    for path in target.rglob("*"):
        if path.is_symlink():
            raise OntologyBundleError(f"ontology release contains a symlink: {path}")
        relative = str(path.relative_to(target))
        if path.is_file():
            actual_files.add(relative)
        elif path.is_dir():
            actual_directories.add(relative)
        else:
            raise OntologyBundleError(f"ontology release contains an invalid entry: {path}")
    if actual_files != expected_files or actual_directories != expected_directories:
        raise OntologyBundleError(
            f"ontology release directory is incomplete or contaminated: {target}"
        )
    for artifact in contract.artifacts:
        _verify_artifact(root / artifact.path, artifact)


def _publish_release_bundle(
    root: Path,
    contract: OntologyBundleContract,
    source_paths: Mapping[str, Path],
) -> None:
    relative_root = _bundle_release_relative(contract)
    for artifact in contract.artifacts:
        path = PurePosixPath(artifact.path)
        if not path.is_relative_to(relative_root) or path == relative_root:
            raise OntologyBundleError(
                f"ontology artifact is outside its atomic bundle: {artifact.role}"
            )
    target = root / str(relative_root)
    _require_within_root(root, target, "ontology release")
    if target.exists():
        _validate_release_tree(root, contract)
        return
    _require_within_root(root, target.parent, "ontology release parent")
    target.parent.mkdir(parents=True, exist_ok=True)
    partial_root = target.parent / ".partial"
    _require_within_root(root, partial_root, "ontology partial release")
    partial_root.mkdir(mode=0o700, exist_ok=True)
    staging = partial_root / f"{target.name}.{uuid.uuid4().hex}.partial"
    staging.mkdir(mode=0o700)
    try:
        for artifact in contract.artifacts:
            relative = PurePosixPath(artifact.path).relative_to(relative_root)
            destination = staging / str(relative)
            destination.parent.mkdir(parents=True, exist_ok=True)
            _copy_file(source_paths[artifact.role], destination)
            _verify_artifact(destination, artifact)
            os.chmod(destination, 0o444)
        for directory in (path for path in staging.rglob("*") if path.is_dir()):
            os.chmod(directory, 0o755)
        os.chmod(staging, 0o755)
        os.replace(staging, target)
        directory_fd = os.open(target.parent, os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    _validate_release_tree(root, contract)


def _bundle_manifest_path(root: Path, contract: OntologyBundleContract) -> Path:
    return root / f"manifests/ontologies/{contract.bundle_id}.bundle.json"


def _bundle_manifest_data(
    contract: OntologyBundleContract,
    contract_sha256: str,
    publications: tuple[tuple[Path, str], ...],
    validation: OntologyBundleValidation,
    root: Path,
) -> Mapping[str, Any]:
    by_resource = {
        Path(path).stem: (path, digest) for path, digest in publications
    }
    resources = []
    for artifact in contract.artifacts:
        path, digest = by_resource[artifact.resource_id]
        resources.append(
            {
                "role": artifact.role,
                "resource_id": artifact.resource_id,
                "release": artifact.release,
                "source_url": artifact.url,
                "source_sha256": artifact.sha256,
                "license": artifact.license,
                "license_url": artifact.license_url,
                "manifest_path": str(path.relative_to(root)),
                "manifest_sha256": digest,
            }
        )
    return {
        "schema": "genosapiens.ontology-bundle-manifest",
        "classification": "public",
        "bundle_id": contract.bundle_id,
        "contract_sha256": contract_sha256,
        "resources": resources,
        "validation": asdict(validation),
    }


def _publish_bundle_manifest(data: Mapping[str, Any], path: Path, root: Path) -> str:
    _require_within_root(root, path, "ontology bundle manifest")
    content = (
        json.dumps(data, sort_keys=True, separators=(",", ":")) + "\n"
    ).encode()
    digest = hashlib.sha256(content).hexdigest()
    if path.exists():
        if path.is_symlink() or not path.is_file():
            raise OntologyBundleError(
                f"ontology bundle manifest is not a regular file: {path}"
            )
        try:
            existing = path.read_bytes()
        except OSError as error:
            raise OntologyBundleError(
                f"cannot read ontology bundle manifest {path}: {error}"
            ) from error
        if existing != content:
            raise OntologyBundleError("completed ontology bundle manifest is immutable")
        return digest
    path.parent.mkdir(parents=True, exist_ok=True)
    partial_root = path.parent / ".partial"
    _require_within_root(root, partial_root, "ontology bundle manifest partial path")
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
            if path.is_symlink() or not path.is_file():
                raise OntologyBundleError(
                    f"ontology bundle manifest is not a regular file: {path}"
                )
            try:
                existing = path.read_bytes()
            except OSError as error:
                raise OntologyBundleError(
                    f"cannot read ontology bundle manifest {path}: {error}"
                ) from error
            if existing != content:
                raise OntologyBundleError(
                    "completed ontology bundle manifest is immutable"
                )
            return digest
        except OSError as error:
            raise OntologyBundleError(
                f"cannot publish ontology bundle manifest {path}: {error}"
            ) from error
        directory_fd = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    finally:
        partial.unlink(missing_ok=True)
    return digest


def bootstrap_ontology_bundle(
    *,
    contract_path: str | Path,
    reference_root: str | Path,
    image_digest: str,
    retrieved_on: str | None = None,
    scratch_root: str | Path | None = None,
    reserve_bytes: int = DEFAULT_RESERVE_BYTES,
    downloader: Downloader = download_atomic,
    capacity_checker: CapacityChecker = require_capacity,
) -> OntologyBootstrapResult:
    contract_file = Path(contract_path)
    root = Path(reference_root).resolve(strict=True)
    contract = load_ontology_bundle_contract(contract_file)
    contract_sha256 = _file_identity(contract_file)[0]
    scratch_base = (
        Path(scratch_root)
        if scratch_root is not None
        else root / "scratch/public-build"
    )
    if not scratch_base.is_absolute():
        raise OntologyBundleError("ontology scratch_root must be absolute")
    _require_within_root(root, scratch_base, "ontology scratch_root")
    bundle_path = _bundle_manifest_path(root, contract)
    _require_within_root(root, bundle_path, "ontology bundle manifest")
    manifest_partial_root = bundle_path.parent / ".partial"
    _require_within_root(
        root, manifest_partial_root, "ontology manifest partial path"
    )
    bundle_existed = bundle_path.exists()
    if bundle_existed and not bundle_path.is_file():
        raise OntologyBundleError(
            f"ontology bundle manifest is not a regular file: {bundle_path}"
        )
    for artifact in contract.artifacts:
        _require_within_root(root, root / artifact.path, "ontology artifact")
    release_target = root / str(_bundle_release_relative(contract))
    if release_target.exists():
        _validate_release_tree(root, contract)
    manifests: dict[str, PublicResourceManifest] = {}
    for artifact in contract.artifacts:
        path = root / _manifest_relative_path(artifact)
        _require_within_root(root, path, "ontology manifest")
        if not path.exists():
            if bundle_existed:
                raise OntologyBundleError(
                    "completed ontology bundle is missing a resource manifest: "
                    f"{artifact.role}"
                )
            continue
        manifest = _load_manifest(path)
        _verify_manifest_contract(manifest, artifact, contract_sha256)
        try:
            verify_outputs(manifest, root)
        except ManifestError as error:
            raise OntologyBundleError(
                f"completed ontology resource is invalid: {artifact.role}: {error}"
            ) from error
        manifests[artifact.role] = manifest

    if len(manifests) == len(contract.artifacts):
        paths = {artifact.role: root / artifact.path for artifact in contract.artifacts}
        validation = validate_ontology_bundle(contract, paths)
        publications = tuple(
            (
                root / _manifest_relative_path(artifact),
                _file_identity(root / _manifest_relative_path(artifact))[0],
            )
            for artifact in contract.artifacts
        )
        if not bundle_existed:
            scratch = scratch_base / f"ontologies/{contract.bundle_id}"
            _require_within_root(root, scratch, "ontology build scratch")
            for artifact in contract.artifacts:
                scratch_artifact = scratch / artifact.role
                _require_within_root(
                    root, scratch_artifact, "ontology scratch artifact"
                )
                scratch_artifact.unlink(missing_ok=True)
        bundle_digest = _publish_bundle_manifest(
            _bundle_manifest_data(
                contract, contract_sha256, publications, validation, root
            ),
            bundle_path,
            root,
        )
        return OntologyBootstrapResult(
            bundle_id=contract.bundle_id,
            bundle_manifest_path=str(bundle_path),
            bundle_manifest_sha256=bundle_digest,
            manifest_paths=tuple(str(path) for path, _digest in publications),
            manifest_sha256=tuple(digest for _path, digest in publications),
            validation=validation,
            reused=bundle_existed,
        )

    source_paths: dict[str, Path] = {}
    missing: list[OntologyArtifact] = []
    for artifact in contract.artifacts:
        installed = root / artifact.path
        if installed.exists():
            _verify_artifact(installed, artifact)
            source_paths[artifact.role] = installed
        else:
            missing.append(artifact)
    missing_size = sum(artifact.size_bytes for artifact in missing)
    if missing_size:
        capacity_checker(
            root,
            CapacityPlan(
                operation="bootstrap core ontology references",
                source_bytes=missing_size,
                output_bytes=missing_size,
                reserve_bytes=reserve_bytes,
            ),
        )
    scratch = scratch_base / f"ontologies/{contract.bundle_id}"
    _require_within_root(root, scratch, "ontology build scratch")
    _require_within_root(root, scratch / ".partial", "ontology download partial path")
    scratch.mkdir(parents=True, exist_ok=True, mode=0o700)
    for artifact in missing:
        destination = scratch / artifact.role
        _require_within_root(root, destination, "ontology scratch artifact")
        downloader(
            DownloadRequest(
                url=artifact.url,
                destination=destination,
                expected_size_bytes=artifact.size_bytes,
                expected_sha256=artifact.sha256,
                attempts=3,
                timeout_seconds=120,
            ),
            root,
        )
        source_paths[artifact.role] = destination

    validation = validate_ontology_bundle(contract, source_paths)
    _publish_release_bundle(root, contract, source_paths)
    final_paths = {artifact.role: root / artifact.path for artifact in contract.artifacts}
    validation = validate_ontology_bundle(contract, final_paths)
    publications = []
    for artifact in contract.artifacts:
        manifest_path = root / _manifest_relative_path(artifact)
        if artifact.role in manifests:
            publications.append((manifest_path, _file_identity(manifest_path)[0]))
            continue
        manifest = PublicResourceManifest(
            resource_id=artifact.resource_id,
            release=artifact.release,
            source=SourceIdentity(artifact.url, sha256=artifact.sha256),
            license=artifact.license,
            assembly="not-applicable",
            retrieved_on=retrieved_on or date.today().isoformat(),
            transformation=Transformation(TRANSFORMATION_SCHEMA, contract_sha256),
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
        publication = publish_manifest(manifest, manifest_path, root)
        publications.append((publication.path, publication.manifest_sha256))
    for artifact in contract.artifacts:
        scratch_artifact = scratch / artifact.role
        _require_within_root(root, scratch_artifact, "ontology scratch artifact")
        scratch_artifact.unlink(missing_ok=True)
    publication_tuple = tuple(publications)
    bundle_digest = _publish_bundle_manifest(
        _bundle_manifest_data(
            contract, contract_sha256, publication_tuple, validation, root
        ),
        bundle_path,
        root,
    )
    return OntologyBootstrapResult(
        bundle_id=contract.bundle_id,
        bundle_manifest_path=str(bundle_path),
        bundle_manifest_sha256=bundle_digest,
        manifest_paths=tuple(str(path) for path, _digest in publications),
        manifest_sha256=tuple(digest for _path, digest in publications),
        validation=validation,
        reused=False,
    )


def ontology_result_as_json(result: OntologyBootstrapResult) -> str:
    return json.dumps(asdict(result), indent=2, sort_keys=True) + "\n"
