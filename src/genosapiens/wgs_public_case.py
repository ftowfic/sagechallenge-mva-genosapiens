"""Persisted current public case preparation shared by initial science and replay.

Replaying this contract authenticates original public inputs and every prepared
source. It performs no preparation, annotation, or scientific invocation.
"""

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Literal

from genosapiens.case_manifest import CaseFileClass, case_manifest_as_bytes, load_case_manifest
from genosapiens.evidence_model import _sha, _typed
from genosapiens.wgs_admission import NativeSourceContract
from genosapiens.wgs_failures import WgsFailure
from genosapiens.wgs_ledger import canonical, publish_bytes
from genosapiens.wgs_provenance import execution_sha256
from genosapiens.wgs_read_rescue import (
    PreparedIntakeSource, MergedReadIntake, MergedIntakeIndex,
    load_prepared_intake_source, load_merged_intake,
)
from genosapiens.wgs_shards import read_bounded_metadata
from genosapiens.wgs_scientific import _bytes
from genosapiens.wgs_validation_authority import PublicValidationAuthority, load_public_validation_authority
from genosapiens.reference_concordance import ReferenceRoleBinding
from genosapiens.submission_rendering import load_submission_contig_renderer
from genosapiens.wgs_source_routing import SourceRoutingContext


def _relative(value):
    path = Path(value)
    if (not isinstance(value, str) or not value or not path.parts or path.is_absolute() or path.as_posix() != value
            or any(part in {"..", ".", "data", "controlled", "controlled-work"} for part in path.parts)):
        raise WgsFailure("prepared_authority", "malformed_evidence", operation="prepared_authority_validate", exception_kind="project_failure")


@dataclass(frozen=True, slots=True, repr=False)
class PublicPreparedSource:
    input_sha256: str
    preparation_directory: str
    prepared_sha256: str
    early_snapshot_sha256: str
    native_contract: NativeSourceContract

    def __post_init__(self):
        for value in (self.input_sha256, self.prepared_sha256, self.early_snapshot_sha256):
            _sha(value)
        _relative(self.preparation_directory)
        if type(self.native_contract) is not NativeSourceContract:
            raise WgsFailure("prepared_authority", "malformed_evidence", operation="prepared_authority_validate", exception_kind="project_failure")


@dataclass(frozen=True, slots=True, repr=False)
class PublicCaseManifest:
    authority_sha256: str
    case_manifest_sha256: str
    case_manifest_path: str
    execution_config_sha256: str
    sources: tuple[PublicPreparedSource, ...]
    intake_directory: str
    intake_sha256: str
    schema: Literal["genosapiens.wgs-public-prepared-case"] = "genosapiens.wgs-public-prepared-case"
    version: Literal[1] = 1

    def __post_init__(self):
        for value in (self.authority_sha256, self.case_manifest_sha256, self.execution_config_sha256, self.intake_sha256):
            _sha(value)
        _relative(self.case_manifest_path)
        _relative(self.intake_directory)
        if (self.schema != "genosapiens.wgs-public-prepared-case" or type(self.version) is not int or self.version != 1
                or not self.sources or tuple(row.input_sha256 for row in self.sources)
                    != tuple(sorted({row.input_sha256 for row in self.sources}))):
            raise WgsFailure("prepared_authority", "malformed_evidence", operation="prepared_authority_validate", exception_kind="project_failure")

    def as_bytes(self):
        content = _bytes({"classification": "controlled", "llm_policy": "no-llm", "manifest": self})
        if len(content) > 8 * 1024 * 1024:
            raise WgsFailure("prepared_authority", "response_capacity", operation="prepared_authority_validate", exception_kind="project_failure")
        return content

    @property
    def sha256(self):
        return hashlib.sha256(self.as_bytes()).hexdigest()


@dataclass(frozen=True, slots=True, repr=False)
class BoundedPublicCase:
    manifest: PublicCaseManifest
    authority: PublicValidationAuthority
    sources: tuple[PreparedIntakeSource, ...]
    intake: MergedReadIntake
    case_manifest_path: Path


def _source_routing_context(contract, resources):
    renderer = load_submission_contig_renderer(
        identity_path=resources.population_reference.identity_path,
        reference_manifest_path=resources.population_reference.manifest_path,
        reference_root=resources.population_reference.reference_root,
        reference_binding=ReferenceRoleBinding(
            contract.reference_role.value,
            contract.reference_identity_id,
        ),
    )
    return SourceRoutingContext(contract, renderer)


def load_public_case(
    content,
    *,
    expected_sha256,
    authority_bytes,
    execution_config,
    scratch_path,
    resources,
):
    """Authenticate saved preparation under its exact public-only authority."""
    try:
        if (type(content) is not bytes or len(content) > 8 * 1024 * 1024
                or hashlib.sha256(content).hexdigest() != expected_sha256):
            raise ValueError
        data = json.loads(content)
        if (set(data) != {"classification", "llm_policy", "manifest"}
                or data["classification"] != "controlled" or data["llm_policy"] != "no-llm"):
            raise ValueError
        manifest = _typed(data["manifest"], PublicCaseManifest, decode=True)
        if manifest.as_bytes() != content or manifest.execution_config_sha256 != execution_sha256(execution_config):
            raise ValueError
        authority = load_public_validation_authority(authority_bytes, expected_sha256=manifest.authority_sha256)
        root = Path(authority.output_root)
        inputs = {row.sha256: row for row in authority.inputs if row.kind == "vcf"}
        if set(inputs) != {row.input_sha256 for row in manifest.sources}:
            raise ValueError
        case_path = root / manifest.case_manifest_path
        case_content = read_bounded_metadata(case_path)
        if hashlib.sha256(case_content).hexdigest() != manifest.case_manifest_sha256:
            raise ValueError
        case = load_case_manifest(case_path)
        if (case_manifest_as_bytes(case) != case_content or {
                ("vcf" if row.file_class is CaseFileClass.VCF else "hpo_text", row.sha256, row.size_bytes)
                for row in case.inputs if row.file_class in {CaseFileClass.VCF, CaseFileClass.HPO_LIST}
            } != {(row.kind, row.sha256, row.size_bytes) for row in authority.inputs}
                or any(row.file_class not in {CaseFileClass.VCF, CaseFileClass.HPO_LIST} for row in case.inputs)):
            raise ValueError
        sources = tuple(load_prepared_intake_source(
            preparation_root=root / row.preparation_directory, original_input_path=Path(inputs[row.input_sha256].path),
            expected_prepared_sha256=row.prepared_sha256, expected_early_snapshot_sha256=row.early_snapshot_sha256,
            source_routing_context=_source_routing_context(row.native_contract, resources),
            execution_config=execution_config,
            base_release_bytes=None, public_validation_authority_bytes=authority_bytes) for row in manifest.sources)
        intake_root = root / manifest.intake_directory
        intake_manifest = load_merged_intake(read_bounded_metadata(intake_root / "merged-intake-manifest.json"),
                                             expected_sha256=manifest.intake_sha256)
        if (intake_manifest.deltas or intake_manifest.case_manifest_sha256 != manifest.case_manifest_sha256
                or intake_manifest.source_bindings != tuple(source.binding for source in sources)
                or intake_manifest.sample_id != case.sample_id
                or any(source.prepared.manifest.case_manifest_sha256 != manifest.case_manifest_sha256 for source in sources)):
            raise ValueError
        intake = MergedReadIntake(intake_manifest, intake_root)
        with MergedIntakeIndex(intake, execution_config=execution_config, scratch_path=scratch_path):
            pass
        return BoundedPublicCase(manifest, authority, sources, intake, case_path)
    except WgsFailure:
        raise
    except (ValueError, TypeError, RuntimeError, KeyError, IndexError, OSError):
        raise WgsFailure("prepared_authority", "malformed_evidence", operation="prepared_authority_validate", exception_kind="project_failure") from None


def publish_public_case(
    *, directory, authority, sources, intake, case_manifest_path, execution_config,
    resources,
):
    """Persist only exact prepared-source bindings, then independently read back."""
    root = Path(authority.output_root)
    case_content = read_bounded_metadata(case_manifest_path)
    manifest = PublicCaseManifest(authority.sha256, hashlib.sha256(case_content).hexdigest(),
        Path(case_manifest_path).relative_to(root).as_posix(), execution_sha256(execution_config),
        tuple(PublicPreparedSource(source.binding.input_sha256,
            source.prepared.root.parent.relative_to(root).as_posix(), source.prepared.manifest.sha256,
            source.prepared.early_snapshot.sha256, source.native_contract) for source in sources),
        intake.directory.relative_to(root).as_posix(), intake.manifest.sha256)
    directory = Path(directory)
    authority.require_execution(original_input_path=sources[0].original_input_path,
        hpo_input_sha256=sources[0].prepared.manifest.hpo_input_sha256, output_path=directory,
        executable_stack_sha256=authority.executable_stack_sha256, runtime_image_digest=authority.runtime_image_digest)
    directory.mkdir(mode=0o700, parents=True, exist_ok=True)
    result = load_public_case(manifest.as_bytes(), expected_sha256=manifest.sha256,
        authority_bytes=authority.as_bytes(), execution_config=execution_config,
        scratch_path=directory / "public-case-readback.sqlite",
        resources=resources)
    publish_bytes(directory, "public-case-manifest.json", manifest.as_bytes())
    return result
