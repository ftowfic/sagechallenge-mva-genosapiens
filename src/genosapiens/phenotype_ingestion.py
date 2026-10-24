"""Deterministic controlled ingestion for plain HPO identifier lists."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from enum import StrEnum
import hashlib
import json
from pathlib import Path
import re
import warnings

from genosapiens.operator_warnings import HpoBlankLineWarning

from genosapiens.bounded_json import BoundedJSONError, read_bounded_bytes
from genosapiens.evidence_model import Observation
from genosapiens.governance import ArtifactClassification
from genosapiens.llm_guard import LLMPolicy
from genosapiens.ontology_bundle import (
    OBSOLETE_RESOLUTION_POLICY,
    ParsedOntology,
)
from genosapiens.phenotype_scoring import (
    ObservationProvenance,
    ScoringObservation,
)


PHENOTYPE_INGESTION_SCHEMA = "genosapiens.phenotype-ingestion"
MAX_HPO_LIST_BYTES = 2 * 1024 * 1024
SHA256 = re.compile(r"^[0-9a-f]{64}$")
HPO_ID = re.compile(r"^HP:[0-9]{7}$")


class PhenotypeIngestionError(RuntimeError):
    """Raised without protected content when phenotype ingestion fails."""


class UnlistedTermStatus(StrEnum):
    UNKNOWN = "unknown"


class PhenotypeInputType(StrEnum):
    HPO_LIST = "hpo-list"


class InputConvention(StrEnum):
    LISTED_POSITIVE_DEFAULT_PROBAND = "listed-positive-default-proband"


@dataclass(frozen=True, slots=True)
class PhenotypeIngestionResult:
    input_sha256: str
    ontology_manifest_sha256: str
    resolution_policy: str
    observations: tuple[Observation, ...]
    proband_observation_ids: tuple[str, ...]
    perinatal_observation_ids: tuple[str, ...]
    family_history_observation_ids: tuple[str, ...]
    unresolved_source_hpo_ids: tuple[str, ...]
    unlisted_term_status: UnlistedTermStatus
    classification: ArtifactClassification
    llm_policy: LLMPolicy
    input_type: PhenotypeInputType = PhenotypeInputType.HPO_LIST
    input_convention: InputConvention = InputConvention.LISTED_POSITIVE_DEFAULT_PROBAND
    schema: str = PHENOTYPE_INGESTION_SCHEMA

    def __post_init__(self) -> None:
        identifiers = tuple(
            item.scoring.observation_id for item in self.observations
        )
        unresolved = tuple(
            sorted(
                {
                    item.scoring.source_hpo_id
                    for item in self.observations
                    if item.resolution.canonical_id is None
                }
            )
        )
        if (
            self.schema != PHENOTYPE_INGESTION_SCHEMA
            or SHA256.fullmatch(self.input_sha256) is None
            or SHA256.fullmatch(self.ontology_manifest_sha256) is None
            or self.resolution_policy != OBSOLETE_RESOLUTION_POLICY
            or not self.observations
            or len(identifiers) != len(set(identifiers))
            or self.proband_observation_ids != identifiers
            or self.perinatal_observation_ids
            or self.family_history_observation_ids
            or self.unresolved_source_hpo_ids != unresolved
            or any(
                item.ontology_manifest_sha256 != self.ontology_manifest_sha256
                or item.resolution_policy != self.resolution_policy
                or item.input_convention != self.input_convention.value
                or item.scoring.scope != "proband"
                or item.scoring.status != "present"
                or item.scoring.provenance.input_sha256 != self.input_sha256
                or item.scoring.provenance.input_type != self.input_type.value
                or re.fullmatch(
                    r"line:[1-9][0-9]*",
                    item.scoring.provenance.source_locator,
                )
                is None
                or any(
                    getattr(item.scoring, field) is not None
                    for field in (
                        "subject_id",
                        "onset",
                        "age_years",
                        "sex",
                        "domain",
                        "dependency_cluster",
                        "quality",
                    )
                )
                or item.measurements
                or item.notes
                or item.context_provenance
                or item.source_fields is not None
                for item in self.observations
            )
            or self.unlisted_term_status is not UnlistedTermStatus.UNKNOWN
            or self.classification is not ArtifactClassification.CONTROLLED
            or self.llm_policy is not LLMPolicy.NO_LLM
            or self.input_type is not PhenotypeInputType.HPO_LIST
            or self.input_convention
            is not InputConvention.LISTED_POSITIVE_DEFAULT_PROBAND
        ):
            raise PhenotypeIngestionError("phenotype ingestion result is invalid")

    @property
    def proband_vector(self) -> tuple[str, ...]:
        return tuple(
            item.resolution.canonical_id
            for item in self.observations
            if item.resolution.canonical_id is not None
        )


def _observation_id(
    source_hpo_id: str,
    line_number: int,
    *,
    input_sha256: str,
    resolution: object,
    ontology_manifest_sha256: str,
) -> str:
    semantic = {
        "source_hpo_id": source_hpo_id,
        "resolution": asdict(resolution),
        "ontology_manifest_sha256": ontology_manifest_sha256,
        "resolution_policy": OBSOLETE_RESOLUTION_POLICY,
        "scope": "proband",
        "status": "present",
        "input_sha256": input_sha256,
        "input_type": PhenotypeInputType.HPO_LIST.value,
        "source_locator": f"line:{line_number}",
        "input_convention": InputConvention.LISTED_POSITIVE_DEFAULT_PROBAND.value,
    }
    content = (
        json.dumps(semantic, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
        + "\n"
    ).encode()
    return f"observation-{hashlib.sha256(content).hexdigest()[:24]}"


def _hpo_lines(content: bytes) -> tuple[tuple[str, int], ...]:
    try:
        text = content.decode("utf-8", errors="strict")
    except UnicodeDecodeError as error:
        raise PhenotypeIngestionError("phenotype input is not valid UTF-8") from error
    if "\r" in text.replace("\r\n", ""):
        raise PhenotypeIngestionError("plain HPO lists require one HPO ID per line")
    lines = text.replace("\r\n", "\n").split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    rows = []
    blanks = False
    for number, line in enumerate(lines, start=1):
        if not line.strip(" \t"):
            blanks = True
            continue
        if HPO_ID.fullmatch(line) is None:
            raise PhenotypeIngestionError(
                f"plain HPO lists require one HP: identifier with seven digits per line; line {number}"
            )
        rows.append((line, number))
    if blanks:
        warnings.warn("Blank lines in the HPO file were ignored. Use one HPO identifier per nonblank line.",
                      HpoBlankLineWarning, stacklevel=2)
    if not rows:
        raise PhenotypeIngestionError("phenotype input contains no HPO identifiers")
    return tuple(rows)


def _make_observation(
    source_hpo_id: str,
    line_number: int,
    *,
    input_sha256: str,
    hpo_ontology: ParsedOntology,
    ontology_manifest_sha256: str,
) -> Observation:
    resolution = hpo_ontology.resolve_term(source_hpo_id)
    provenance = ObservationProvenance(
        input_sha256,
        f"line:{line_number}",
        PhenotypeInputType.HPO_LIST.value,
    )
    return Observation(
        scoring=ScoringObservation(
            observation_id=_observation_id(
                source_hpo_id,
                line_number,
                input_sha256=input_sha256,
                resolution=resolution,
                ontology_manifest_sha256=ontology_manifest_sha256,
            ),
            source_hpo_id=source_hpo_id,
            scope="proband",
            status="present",
            provenance=provenance,
        ),
        resolution=resolution,
        ontology_manifest_sha256=ontology_manifest_sha256,
        resolution_policy=OBSOLETE_RESOLUTION_POLICY,
        input_convention=InputConvention.LISTED_POSITIVE_DEFAULT_PROBAND.value,
    )


def ingest_phenotypes(
    path: str | Path,
    *,
    expected_input_sha256: str,
    hpo_ontology: ParsedOntology,
    ontology_manifest_sha256: str,
    reference_root: str | Path | None = None,
    expected_input_size_bytes: int | None = None,
) -> PhenotypeIngestionResult:
    if (
        SHA256.fullmatch(expected_input_sha256) is None
        or type(hpo_ontology) is not ParsedOntology
        or SHA256.fullmatch(ontology_manifest_sha256) is None
    ):
        raise PhenotypeIngestionError("phenotype input identity is invalid")
    try:
        content = read_bounded_bytes(
            path,
            maximum_bytes=MAX_HPO_LIST_BYTES,
            reference_root=reference_root,
            expected_size_bytes=expected_input_size_bytes,
        )
    except BoundedJSONError as error:
        raise PhenotypeIngestionError("phenotype input cannot be read safely") from error
    if hashlib.sha256(content).hexdigest() != expected_input_sha256:
        raise PhenotypeIngestionError("phenotype input identity changed")
    observations = tuple(
        _make_observation(
            source_hpo_id,
            line_number,
            input_sha256=expected_input_sha256,
            hpo_ontology=hpo_ontology,
            ontology_manifest_sha256=ontology_manifest_sha256,
        )
        for source_hpo_id, line_number in _hpo_lines(content)
    )
    return PhenotypeIngestionResult(
        input_sha256=expected_input_sha256,
        ontology_manifest_sha256=ontology_manifest_sha256,
        resolution_policy=OBSOLETE_RESOLUTION_POLICY,
        observations=observations,
        proband_observation_ids=tuple(
            item.scoring.observation_id for item in observations
        ),
        perinatal_observation_ids=(),
        family_history_observation_ids=(),
        unresolved_source_hpo_ids=tuple(
            sorted(
                {
                    item.scoring.source_hpo_id
                    for item in observations
                    if item.resolution.canonical_id is None
                }
            )
        ),
        unlisted_term_status=UnlistedTermStatus.UNKNOWN,
        classification=ArtifactClassification.CONTROLLED,
        llm_policy=LLMPolicy.NO_LLM,
    )


def phenotype_ingestion_as_json(result: PhenotypeIngestionResult) -> str:
    return json.dumps(asdict(result), indent=2, sort_keys=True) + "\n"
