"""Saved public preparation and fixed annotation for a separate replay worker."""

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Literal

from genosapiens.evidence_model import _sha
from genosapiens.wgs_annotation import load_fixed_annotation_source
from genosapiens.wgs_failures import WgsFailure
from genosapiens.wgs_ledger import publish_bytes
from genosapiens.wgs_public_case import _relative, load_public_case
from genosapiens.wgs_public_execution import load_public_execution
from genosapiens.wgs_scientific import _bytes, _load_manifest
from genosapiens.wgs_shards import read_bounded_metadata
from genosapiens.wgs_validation_authority import load_public_validation_authority


@dataclass(frozen=True, slots=True, repr=False)
class PublicReplayArtifact:
    path: str
    sha256: str

    def __post_init__(self):
        _relative(self.path)
        _sha(self.sha256)


@dataclass(frozen=True, slots=True, repr=False)
class PublicReplaySource:
    execution: PublicReplayArtifact
    authority: PublicReplayArtifact
    public_case: PublicReplayArtifact
    annotation: PublicReplayArtifact
    annotation_plan_directory: str
    annotation_output_directory: str
    schema: Literal["genosapiens.wgs-fixed-public-replay-source"] = "genosapiens.wgs-fixed-public-replay-source"
    version: Literal[2] = 2

    def __post_init__(self):
        _relative(self.annotation_plan_directory)
        _relative(self.annotation_output_directory)
        if (any(type(row) is not PublicReplayArtifact for row in (self.execution, self.authority, self.public_case, self.annotation))
                or self.schema != "genosapiens.wgs-fixed-public-replay-source" or type(self.version) is not int or self.version != 2):
            raise WgsFailure("prepared_authority", "malformed_evidence", operation="prepared_authority_validate", exception_kind="project_failure")

    def as_bytes(self):
        return _bytes({"classification": "controlled", "llm_policy": "no-llm", "manifest": self})

    @property
    def sha256(self):
        return hashlib.sha256(self.as_bytes()).hexdigest()


def capture_public_replay_source(execution, *, directory):
    root = Path(execution.public_case.authority.output_root)
    directory = Path(directory)
    source = execution.public_case.sources[0]
    authority = execution.public_case.authority
    authority.require_execution(original_input_path=source.original_input_path, hpo_input_sha256=source.prepared.manifest.hpo_input_sha256,
        output_path=directory, executable_stack_sha256=authority.executable_stack_sha256, runtime_image_digest=authority.runtime_image_digest)
    directory.mkdir(mode=0o700, parents=True, exist_ok=True)
    def artifact(name, content):
        value = publish_bytes(directory, name, content)
        return PublicReplayArtifact((directory / name).relative_to(root).as_posix(), value.sha256)
    annotation = execution.products.annotation_source
    result = PublicReplaySource(artifact("execution.json", execution.as_bytes()),
        artifact("authority.json", authority.as_bytes()), artifact("public-case.json", execution.public_case.manifest.as_bytes()),
        artifact("annotation.json", annotation.as_bytes()), annotation.plan_root.relative_to(root).as_posix(),
        annotation.output_root.relative_to(root).as_posix())
    publish_bytes(directory, "fixed-public-source.json", result.as_bytes())
    return result


def load_public_replay_source(content, *, expected_sha256, root, resources, policies):
    """Authenticate current source contracts; launch no preparation or annotation."""
    from genosapiens.bounded_wgs_config import load_bounded_wgs_execution_config
    from genosapiens.track1_public_validation import PreparedPublicCase, materialize_runner_view, COMMON_ALLELE_PUBLIC_INDEX_MANIFEST_PATH
    from genosapiens.track1_public_validation_coordinator import executable_stack_bytes
    from genosapiens.common_allele_public_indexes import verify_common_allele_public_indexes
    from genosapiens.wgs_public_analysis import (
        VEP_BUNDLE_SHA256, VEP_CONTRACT_PATH, _resource_session_sha256,
    )
    source = _load_manifest(content, expected_sha256=expected_sha256, result_type=PublicReplaySource, result_key="manifest")
    root = Path(root)
    def read(artifact):
        content = read_bounded_metadata(root / artifact.path)
        if hashlib.sha256(content).hexdigest() != artifact.sha256:
            raise WgsFailure("prepared_authority", "malformed_evidence", operation="prepared_authority_validate", exception_kind="project_failure")
        return content
    authority_content = read(source.authority)
    authority = load_public_validation_authority(authority_content, expected_sha256=source.authority.sha256)
    execution = load_public_execution(read(source.execution), expected_sha256=source.execution.sha256)
    stack = resources.executable_stack_content or executable_stack_bytes(resources.repository_root)
    if (Path(authority.output_root) != root or authority.executable_stack_sha256 != hashlib.sha256(stack).hexdigest()
            or authority.policy_family_sha256 != policies.policy_sha256 or execution.public_authority_sha256 != authority.sha256
            or execution.public_case_sha256 != source.public_case.sha256
            or execution.ranker_case_view_sha256 != authority.runner_view_sha256
            or dict(execution.output_bindings)["annotation_source"] != source.annotation.sha256
            or not any(row.candidate.bindings.sha256 == execution.candidate_stack_sha256 for row in policies.runtimes)):
        raise WgsFailure("prepared_authority", "malformed_evidence", operation="prepared_authority_validate", exception_kind="project_failure")
    case_content = read(source.public_case)
    config_sha = json.loads(case_content)["manifest"]["execution_config_sha256"]
    config = load_bounded_wgs_execution_config(read_bounded_metadata(resources.repository_root / "config/track1/bounded-wgs-execution-v5.json"), config_sha)
    scratch = resources.scratch_root / "fixed-source-readback"
    authority.require_execution(original_input_path=Path(next(row.path for row in authority.inputs if row.kind == "vcf")),
        hpo_input_sha256=next(row.sha256 for row in authority.inputs if row.kind == "hpo_text"), output_path=scratch,
        executable_stack_sha256=authority.executable_stack_sha256, runtime_image_digest=authority.runtime_image_digest)
    scratch.mkdir(mode=0o700, parents=True, exist_ok=True)
    public_case = load_public_case(case_content, expected_sha256=source.public_case.sha256, authority_bytes=authority_content,
        execution_config=config, scratch_path=scratch / "source-intake.sqlite",
        resources=resources)
    annotation = load_fixed_annotation_source(read(source.annotation), expected_sha256=source.annotation.sha256,
        plan_root=root / source.annotation_plan_directory, output_root=root / source.annotation_output_directory)
    original = public_case.sources[0]
    annotation.verify(expected_plan=annotation.plan, source_manifest=original.source, working_manifest=original.working,
        resource_session_sha256=_resource_session_sha256(resources, executable_stack_sha256=authority.executable_stack_sha256,
            runtime_image_digest=authority.runtime_image_digest,
            vep_contract_bytes=read_bounded_metadata(
                resources.repository_root / VEP_CONTRACT_PATH
            )), bundle_manifest_sha256=VEP_BUNDLE_SHA256,
        executable_stack_sha256=authority.executable_stack_sha256, runtime_image_digest=authority.runtime_image_digest,
        execution_config=config, intake_binding=public_case.intake.manifest.annotation_binding,
        prepared_index=original.prepared.prepared_index,
        prepared_state_root=original.prepared.root.parent / "state")
    indexes = verify_common_allele_public_indexes(reference_root=resources.reference_root,
        manifest_path=COMMON_ALLELE_PUBLIC_INDEX_MANIFEST_PATH,
        expected_manifest_sha256=authority.policy_bindings.common_allele_public_index_manifest_sha256,
        verify_sources=False, reuse_qualified_reconciliation=True)
    materialized = materialize_runner_view(authority.runner_view, resources)
    from genosapiens.wgs_public_analysis import load_analysis_screening_arguments
    screening=load_analysis_screening_arguments(resources,authority.policy_bindings,indexes,
        dependencies=(('merged_intake_manifest',str((public_case.intake.directory/'merged-intake-manifest.json').absolute()),
            public_case.intake.manifest.sha256,len(public_case.intake.manifest.as_bytes())),))
    return PreparedPublicCase(authority.case_id, authority.runner_view, materialized, public_case,
        indexes, "0", "0", screening, annotation)
