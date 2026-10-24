"""Construct a new bounded-WGS operator package from repository-owned sources."""

from __future__ import annotations

import ast
from dataclasses import asdict, dataclass
import hashlib
import json
import os
from pathlib import Path
import re
import stat

from genosapiens.bounded_wgs_config import load_bounded_wgs_execution_config
from genosapiens.common_allele_filter import load_common_allele_filter_policy
from genosapiens.track1_release import load_base_release
from genosapiens.track1_publication import READY_ROLES
from genosapiens.track1_public_validation_coordinator import executable_stack_bytes
from genosapiens.wgs_analysis_policy import current_interfaces
from genosapiens.wgs_preparation import PreparationRunOptions
from genosapiens.wgs_preparation import CompletedSourceArtifact, CompletedSourceProducer, completed_source_readonly_paths, load_completed_source_input, resolve_completed_source
from genosapiens.wgs_runtime import (
    WgsRuntimeError, load_stage_runtime_evidence,
    require_capacity_qualified_runtime,
)


CONTRACT_SCHEMA = "genosapiens.bounded-wgs-operator-contract"
BINDINGS_SCHEMA = "genosapiens.bounded-wgs-operator-package-bindings"
PRIMARY_REQUIRED_ROLES = tuple(sorted(("primary", *(str(role) for role in READY_ROLES))))
_SHA = re.compile(r"^[0-9a-f]{64}$")
MAX_PACKAGE_SOURCE_BYTES = 16 * 1024 * 1024
_DRIVERS = {
    "preparation": "prepare_bounded_case.py",
    "analysis": "analyze_bounded_case.py",
}
_RUNTIME_HELPERS = (
    "archive_bounded_resume.py",
    "archive_failed_preparation.py",
    "check_bounded_host_capacity.py",
    "check_stage_prerequisite.py",
    "finalize_bounded_runtime.py",
    "run_bounded_stage.py",
    "seal_bounded_stage_evidence.py",
)


class WgsOperatorPackageError(RuntimeError):
    """Static package-construction failure."""


def _digest(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def _sha(value: object) -> bool:
    return type(value) is str and _SHA.fullmatch(value) is not None


def _canonical(value: object) -> bytes:
    return (
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        )
        + "\n"
    ).encode("ascii")


def _private_directory(path: Path, *, create: bool = False) -> None:
    try:
        if create:
            path.mkdir(mode=0o700)
        metadata = path.lstat()
        if (
            not path.is_absolute()
            or path.is_symlink()
            or not stat.S_ISDIR(metadata.st_mode)
            or stat.S_IMODE(metadata.st_mode) != 0o700
            or metadata.st_uid != os.getuid()
        ):
            raise ValueError
    except (OSError, ValueError):
        raise WgsOperatorPackageError("operator package directory is unsafe") from None


def _source(path: Path) -> bytes:
    descriptor = -1
    try:
        descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
        before = os.fstat(descriptor)
        if (
            not stat.S_ISREG(before.st_mode)
            or before.st_size < 1
            or before.st_size > MAX_PACKAGE_SOURCE_BYTES
        ):
            raise ValueError
        with os.fdopen(descriptor, "rb", closefd=True) as stream:
            descriptor = -1
            content = stream.read(MAX_PACKAGE_SOURCE_BYTES + 1)
            after = os.fstat(stream.fileno())
        if (
            len(content) != before.st_size
            or before.st_dev != after.st_dev
            or before.st_ino != after.st_ino
            or before.st_size != after.st_size
            or before.st_mtime_ns != after.st_mtime_ns
            or before.st_ctime_ns != after.st_ctime_ns
        ):
            raise ValueError
        return content
    except (OSError, ValueError):
        raise WgsOperatorPackageError("operator package source is invalid") from None
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _write(path: Path, content: bytes, mode: int) -> None:
    descriptor = -1
    try:
        missing = []
        parent = path.parent
        while not parent.exists():
            missing.append(parent)
            parent = parent.parent
        for directory in reversed(missing):
            directory.mkdir(mode=0o700)
        descriptor = os.open(
            path,
            os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
            mode,
        )
        os.fchmod(descriptor, mode)
        view = memoryview(content)
        while view:
            written = os.write(descriptor, view)
            if written <= 0:
                raise OSError
            view = view[written:]
        os.fsync(descriptor)
    except OSError:
        raise WgsOperatorPackageError("operator package publication failed") from None
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _driver_contract(source: bytes, expected_import: str) -> frozenset[str]:
    try:
        text = source.decode("utf-8", errors="strict")
        tree = ast.parse(text)
        imports = {
            value.module
            for value in tree.body
            if isinstance(value, ast.ImportFrom) and value.module is not None
        }
        assignment = next(
            value for value in tree.body
            if isinstance(value, ast.Assign)
            and any(
                isinstance(target, ast.Name)
                and target.id == "SHARED_CONTRACT_FIELDS"
                for target in value.targets
            )
        )
        if not (
            isinstance(assignment.value, ast.Call)
            and isinstance(assignment.value.func, ast.Name)
            and assignment.value.func.id == "frozenset"
            and len(assignment.value.args) == 1
        ):
            raise ValueError
        fields = frozenset(ast.literal_eval(assignment.value.args[0]))
        if (
            expected_import not in imports
            or "controlled-analysis-components" in text
            or "controlled_analysis_components" in text
            or "operator-package-compatibility-repair" in text
            or "EvidenceDocument" in text
            or "SuppliedCallInputs" in text
            or not fields
        ):
            raise ValueError
        return fields
    except (UnicodeError, SyntaxError, StopIteration, ValueError, TypeError):
        raise WgsOperatorPackageError("bounded driver source contract is invalid") from None


@dataclass(frozen=True, slots=True)
class BoundedOperatorContract:
    base_release_sha256: str
    executable_stack_sha256: str
    runtime_image_digest: str
    execution_limits_sha256: str
    common_allele_policy_sha256: str
    population_manifest_sha256: str
    clinical_manifest_sha256: str
    exception_manifest_sha256: str
    precision_manifest_sha256: str
    common_allele_public_index_manifest_sha256: str
    ontology_manifest_sha256: str
    default_preparation_run_options_sha256: str
    audited_preparation_run_options_sha256: str
    preparation_runtime_evidence_sha256: str
    analysis_runtime_evidence_sha256: str
    shared_driver_fields: tuple[str, ...]
    evidence_interface_sha256: str
    canonical_handoff_schema_sha256: str
    analysis_eligibility_policy_sha256: str
    context_writer_profile_sha256: str
    analysis_profile: str
    analysis_profile_policy_sha256: str
    population_work_policy_sha256: str
    focused_work_contract_sha256: str
    qualification_scope: str = 'focused-fresh-subsample'
    full_wgs_runtime_validated: bool = False
    primary_required_roles: tuple[str, ...] = PRIMARY_REQUIRED_ROLES
    optional_analysis_required: bool = False
    completed_source_input: CompletedSourceArtifact | None = None
    source_producer: CompletedSourceProducer | None = None
    readonly_dependencies: tuple[str, ...] = ()
    operator_interface: str = 'track1-run'
    operator_runtime_identity: str = 'invoking-non-root-uid-gid'
    operator_default_deadline: str = 'none'
    schema: str = CONTRACT_SCHEMA
    version: int = 8
    classification: str = "public"
    controlled_use_authorized: bool = False

    def __post_init__(self) -> None:
        from genosapiens.track1_analysis_profiles import authenticate_analysis_profile
        authenticate_analysis_profile(self.analysis_profile,
            expected_policy_sha256=self.analysis_profile_policy_sha256)
        if self.completed_source_input is not None:
            value = load_completed_source_input(self.completed_source_input.read())
            if (self.source_producer != value.producer
                    or value.producer.source.executable_stack_sha256 != self.executable_stack_sha256
                    or value.producer.runtime_image_digest != self.runtime_image_digest
                    or value.producer.execution_config_sha256 != self.execution_limits_sha256
                    or self.readonly_dependencies != completed_source_readonly_paths(self.completed_source_input)):
                raise WgsOperatorPackageError("operator completed source dependency changed")
        elif self.source_producer is not None or self.readonly_dependencies:
            raise WgsOperatorPackageError("operator source dependency is incomplete")
        if (
            self.schema != CONTRACT_SCHEMA
            or type(self.version) is not int or self.version != 8
            or self.operator_interface != 'track1-run'
            or self.operator_runtime_identity != 'invoking-non-root-uid-gid'
            or self.operator_default_deadline != 'none'
            or self.qualification_scope != 'focused-fresh-subsample'
            or self.full_wgs_runtime_validated is not False
            or self.primary_required_roles != PRIMARY_REQUIRED_ROLES
            or self.optional_analysis_required is not False
            or self.classification != "public"
            or self.controlled_use_authorized is not False
            or any(not _sha(value) for value in (
                self.base_release_sha256,
                self.executable_stack_sha256,
                self.execution_limits_sha256,
                self.common_allele_policy_sha256,
                self.population_manifest_sha256,
                self.clinical_manifest_sha256,
                self.exception_manifest_sha256,
                self.precision_manifest_sha256,
                self.common_allele_public_index_manifest_sha256,
                self.ontology_manifest_sha256,
                self.default_preparation_run_options_sha256,
                self.audited_preparation_run_options_sha256,
                self.preparation_runtime_evidence_sha256,
                self.analysis_runtime_evidence_sha256,
                self.evidence_interface_sha256,
                self.canonical_handoff_schema_sha256,
                self.analysis_eligibility_policy_sha256,
                self.context_writer_profile_sha256,
                self.analysis_profile_policy_sha256,
                self.population_work_policy_sha256,
                self.focused_work_contract_sha256,
            ))
            or self.default_preparation_run_options_sha256
            != PreparationRunOptions(False).sha256
            or self.audited_preparation_run_options_sha256
            != PreparationRunOptions(True).sha256
            or re.fullmatch(r"sha256:[0-9a-f]{64}", self.runtime_image_digest)
            is None
            or self.shared_driver_fields
            != tuple(sorted(set(self.shared_driver_fields)))
            or not self.shared_driver_fields
        ):
            raise WgsOperatorPackageError("bounded operator contract is invalid")

    def as_bytes(self) -> bytes:
        return _canonical(asdict(self))


@dataclass(frozen=True, slots=True)
class PackageFileBinding:
    relative_path: str
    sha256: str
    size_bytes: int
    mode: int


def load_bounded_operator_contract(content: bytes, *, expected_sha256: str):
    """Read back the sole current package contract, without schema conversion."""
    from genosapiens.evidence_model import _typed
    try:
        if len(content) > 65536 or _digest(content) != expected_sha256:
            raise ValueError
        value = _typed(json.loads(content), BoundedOperatorContract, decode=True)
        if value.as_bytes() != content:
            raise ValueError
        return value
    except (ValueError, TypeError, KeyError, RuntimeError):
        raise WgsOperatorPackageError("operator contract is invalid") from None


@dataclass(frozen=True, slots=True)
class BoundedOperatorPackage:
    root: Path
    contract: BoundedOperatorContract
    bindings: tuple[PackageFileBinding, ...]
    bindings_sha256: str


def build_bounded_operator_package(
    *,
    repository_root: str | Path,
    destination_parent: str | Path,
    base_release_bytes: bytes,
    expected_base_release_sha256: str,
    execution_config_bytes: bytes,
    expected_execution_config_sha256: str,
    common_allele_policy_bytes: bytes,
    expected_common_allele_policy_sha256: str,
    population_manifest_sha256: str,
    clinical_manifest_sha256: str,
    ontology_manifest_sha256: str,
    completed_source_input: CompletedSourceArtifact | None = None,
    preparation_runtime_evidence_bytes: bytes | None = None,
    analysis_runtime_evidence_bytes: bytes | None = None,
    exhaustive: bool = False,
) -> BoundedOperatorPackage:
    """Publish a new package without reading or modifying predecessor packages."""

    from genosapiens.wgs_analysis_policy import select_analysis_profile
    from genosapiens.track1_analysis_profiles import authenticate_analysis_profile
    mode, scope_policy = select_analysis_profile(exhaustive=exhaustive,
        repository_root=repository_root)

    try:
        repository = Path(repository_root).resolve(strict=True)
        parent = Path(destination_parent).resolve(strict=True)
        _private_directory(parent)
        release = load_base_release(
            base_release_bytes,
            expected_sha256=expected_base_release_sha256,
        )
        release.require_analysis_profile(authenticate_analysis_profile(mode, policy=scope_policy))
        from genosapiens.wgs_work_contract import load_focused_work_contract
        work_contract = load_focused_work_contract(repository)
        release.require_focused_work(work_contract)
        bindings = release.policy_bindings
        if release.authorization_status != "accepted":
            raise ValueError
        stack_content = executable_stack_bytes(repository)
        stack = json.loads(stack_content)
        stack_files = {row['path']:row for row in stack['files']}
        def qualified_source(relative):
            row = stack_files[str(relative)]
            content = _source(repository/relative)
            if len(content) != row['size_bytes'] or _digest(content) != row['sha256']:
                raise ValueError("qualified operator source changed")
            return content
        interfaces = current_interfaces(repository_root=repository)
        if (release.executable_stack_sha256 != _digest(stack_content)
                or release.runtime_image_digest != stack['runtime_image'].split('@')[-1]
                or any(getattr(release,name) != digest for name,digest in interfaces.items())):
            raise ValueError("operator executable or interface changed")
        from genosapiens.wgs_analysis_eligibility import load_analysis_eligibility_policy
        from genosapiens.wgs_context_runtime import load_context_writer_profile
        selector_paths = ('config/track1/analysis-eligibility-v1.json',
            'config/track1/context-writer-runtime-v1.json')
        selector_bytes = tuple(_source(repository/path) for path in selector_paths)
        selector_sha = tuple(_digest(content) for content in selector_bytes)
        screening = load_analysis_eligibility_policy(selector_bytes[0],selector_sha[0])
        load_context_writer_profile(selector_bytes[1],selector_sha[1])
        if (screening.common_filter_policy_sha256 != expected_common_allele_policy_sha256
                or any(digest not in release.resource_manifest_sha256s for digest in selector_sha)):
            raise ValueError("operator current selectors are not qualified")
        load_bounded_wgs_execution_config(
            execution_config_bytes, expected_execution_config_sha256,
        )
        load_common_allele_filter_policy(
            common_allele_policy_bytes, expected_common_allele_policy_sha256,
        )
        if (preparation_runtime_evidence_bytes is None
                or analysis_runtime_evidence_bytes is None):
            raise ValueError("external runtime qualification is missing")
        preparation_runtime = load_stage_runtime_evidence(
            preparation_runtime_evidence_bytes, expected_stage="preparation",
            expected_execution_config_sha256=expected_execution_config_sha256,
        )
        analysis_runtime = load_stage_runtime_evidence(
            analysis_runtime_evidence_bytes, expected_stage="analysis",
            expected_execution_config_sha256=expected_execution_config_sha256,
        )
        completed = None
        if completed_source_input is not None:
            descriptor = load_completed_source_input(completed_source_input.read())
            if (descriptor.producer.source.executable_stack_sha256 != release.executable_stack_sha256
                    or descriptor.producer.source.reference_manifest_sha256 != release.reference_manifest_sha256
                    or descriptor.producer.runtime_image_digest != release.runtime_image_digest
                    or descriptor.producer.execution_config_sha256 != expected_execution_config_sha256):
                raise ValueError("operator source requires exact current producer bindings")
            completed = resolve_completed_source(descriptor, expected_source=descriptor.producer.source,
                                                 audit_requested=True)
        require_capacity_qualified_runtime(preparation_runtime, completed_source=completed, required_phases=(
            "input_validation", "source_accounting", "filter_index_commit",
            "preparation_terminal_publication",
        ))
        require_capacity_qualified_runtime(analysis_runtime, required_phases=(
            "prepared_authority_validation", "annotation_response_reconciliation",
            "admission_projection", "scientific_analysis", "final_publication",
            "analysis_terminal_publication",
        ))
        if (
            bindings.bounded_wgs_execution_sha256
            != expected_execution_config_sha256
            or bindings.common_allele_policy_sha256
            != expected_common_allele_policy_sha256
            or any(not _sha(value) for value in (
                population_manifest_sha256,
                clinical_manifest_sha256,
                ontology_manifest_sha256,
            ))
            or any(value not in release.resource_manifest_sha256s for value in (
                population_manifest_sha256,
                clinical_manifest_sha256,
                ontology_manifest_sha256,
                bindings.common_allele_exception_manifest_sha256,
                bindings.gnomad_projection_precision_sha256,
                bindings.common_allele_public_index_manifest_sha256,
            ))
        ):
            raise ValueError
        script_root = repository / "scripts" / "operator_compatibility_repair"
        driver_sources = {
            role: qualified_source((script_root / name).relative_to(repository))
            for role, name in _DRIVERS.items()
        }
        preparation_fields = _driver_contract(
            driver_sources["preparation"], "genosapiens.wgs_preparation",
        )
        analysis_fields = _driver_contract(
            driver_sources["analysis"], "genosapiens.wgs_analysis",
        )
        if preparation_fields != analysis_fields:
            raise ValueError
        contract = BoundedOperatorContract(
            expected_base_release_sha256,
            release.executable_stack_sha256,
            release.runtime_image_digest,
            expected_execution_config_sha256,
            expected_common_allele_policy_sha256,
            population_manifest_sha256,
            clinical_manifest_sha256,
            bindings.common_allele_exception_manifest_sha256,
            bindings.gnomad_projection_precision_sha256,
            bindings.common_allele_public_index_manifest_sha256,
            ontology_manifest_sha256,
            PreparationRunOptions(False).sha256,
            PreparationRunOptions(True).sha256,
            _digest(preparation_runtime_evidence_bytes),
            _digest(analysis_runtime_evidence_bytes),
            tuple(sorted(preparation_fields)),
            interfaces['evidence_interface_sha256'],
            interfaces['canonical_handoff_schema_sha256'],
            *selector_sha,
            analysis_profile=mode,
            analysis_profile_policy_sha256=scope_policy.sha256,
            population_work_policy_sha256=work_contract.policy.sha256,
            focused_work_contract_sha256=work_contract.sha256,
            completed_source_input=completed_source_input,
            source_producer=None if completed_source_input is None else load_completed_source_input(completed_source_input.read()).producer,
            readonly_dependencies=() if completed_source_input is None else completed_source_readonly_paths(completed_source_input),
        )
        package = parent / (
            "operator-package-bounded-wgs-" + expected_base_release_sha256[:16]
        )
        _private_directory(package, create=True)
        for role, name in _DRIVERS.items():
            _write(package / name, driver_sources[role], 0o500)
        for name in _RUNTIME_HELPERS:
            _write(package / name, qualified_source((script_root / name).relative_to(repository)), 0o500)
        source_root = repository / "src" / "genosapiens"
        for source in sorted(source_root.rglob("*.py")):
            relative = source.relative_to(source_root)
            _write(
                package / "python" / "genosapiens" / relative,
                qualified_source(source.relative_to(repository)),
                0o400,
            )
        # Retain qualified policy/schema bytes under their runtime-relative paths.
        # Large installed reference data remain separately bound external inputs.
        for row in stack['files']:
            relative = Path(row['path'])
            if relative.parts[0] not in {'config','schemas'}:
                continue
            content = qualified_source(relative)
            _write(package/relative,content,0o400)
        _write(package/'executable-stack.json',stack_content,0o400)
        _write(package / "operator-contract.json", contract.as_bytes(), 0o400)
        _write(package / "base-release.json", base_release_bytes, 0o400)
        _write(
            package / "bounded-wgs-execution.json", execution_config_bytes, 0o400,
        )
        _write(
            package / "common-allele-prefilter.json",
            common_allele_policy_bytes,
            0o400,
        )
        wrapper = (
            "#!/bin/sh\n"
            "set -eu\n"
            "package_dir=$(CDPATH= cd -- \"$(dirname -- \"$0\")\" && pwd)\n"
            "stage=${1:?preparation or analysis required}\n"
            "shift\n"
            "case \"$stage\" in\n"
            "  preparation) driver=prepare_bounded_case.py ;;\n"
            "  analysis) driver=analyze_bounded_case.py ;;\n"
            "  *) exit 2 ;;\n"
            "esac\n"
            "PYTHONPATH=\"$package_dir/python\" exec python \"$package_dir/$driver\" \"$@\"\n"
        ).encode("ascii")
        _write(package / "run-bounded-step.sh", wrapper, 0o500)
        _write(package / 'genosapiens', (
            '#!/bin/sh\nset -eu\n'
            'package_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)\n'
            'PYTHONPATH="$package_dir/python" exec python -m genosapiens "$@"\n'
        ).encode('ascii'), 0o500)

        file_bindings = []
        for path in sorted(value for value in package.rglob("*") if value.is_file()):
            relative = path.relative_to(package).as_posix()
            metadata = path.lstat()
            content = _source(path)
            file_bindings.append(PackageFileBinding(
                relative,
                _digest(content),
                len(content),
                stat.S_IMODE(metadata.st_mode),
            ))
        manifest = _canonical({
            "schema": BINDINGS_SCHEMA,
            "version": 1,
            "classification": "public",
            "controlled_use_authorized": False,
            "contract_sha256": _digest(contract.as_bytes()),
            "files": tuple(asdict(value) for value in file_bindings),
        })
        _write(package / "package-bindings.json", manifest, 0o400)
        return BoundedOperatorPackage(
            package,
            contract,
            tuple(file_bindings),
            _digest(manifest),
        )
    except WgsOperatorPackageError:
        raise
    except (ValueError, TypeError, KeyError, OSError, RuntimeError, WgsRuntimeError):
        raise WgsOperatorPackageError("bounded operator package build failed") from None
