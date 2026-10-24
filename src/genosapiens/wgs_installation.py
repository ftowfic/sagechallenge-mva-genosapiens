"""Prepare an ordinary Track 1 installation from code and public version metadata."""
from dataclasses import asdict, dataclass
import json
import os
from pathlib import Path
import re
import stat
from genosapiens.bounded_json import read_bounded_bytes
from genosapiens.track1_installation import InstalledVersions, canonical, digest, load_run_versions
from genosapiens.track1_publication import READY_ROLES

PRIMARY_REQUIRED_ROLES = tuple(sorted(("primary", *(str(role) for role in READY_ROLES))))


class WgsOperatorPackageError(RuntimeError):
    """Static installed-package integrity failure."""


def _write(path, content, mode):
    path.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, mode)
    with os.fdopen(descriptor, "wb") as stream:
        os.fchmod(stream.fileno(), mode)
        stream.write(content)
        stream.flush()
        os.fsync(stream.fileno())

SCHEMA = "genosapiens.track1-installed-package"
CONTRACT_SCHEMA = "genosapiens.track1-installed-contract"

@dataclass(frozen=True, slots=True)
class OperatorContract:
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
    scope_policy_sha256: str
    scope_manifest_sha256: str
    default_preparation_run_options_sha256: str
    audited_preparation_run_options_sha256: str
    shared_driver_fields: tuple[str, ...]
    evidence_interface_sha256: str
    canonical_handoff_schema_sha256: str
    analysis_eligibility_policy_sha256: str
    context_writer_profile_sha256: str
    analysis_profile: str
    analysis_profile_policy_sha256: str
    population_work_policy_sha256: str
    focused_work_contract_sha256: str
    primary_required_roles: tuple[str, ...] = PRIMARY_REQUIRED_ROLES
    optional_analysis_required: bool = False
    operator_interface: str = 'track1-run'
    operator_runtime_identity: str = 'invoking-non-root-uid-gid'
    operator_default_deadline: str = 'none'
    schema: str = CONTRACT_SCHEMA
    version: int = 2
    classification: str = "public"

    def __post_init__(self):
        from genosapiens.track1_analysis_profiles import authenticate_analysis_profile
        from genosapiens.wgs_preparation import PreparationRunOptions
        authenticate_analysis_profile(self.analysis_profile, expected_policy_sha256=self.analysis_profile_policy_sha256)
        if (self.schema != CONTRACT_SCHEMA or type(self.version) is not int or self.version != 2
                or self.classification != 'public' or self.operator_interface != 'track1-run'
                or self.operator_runtime_identity != 'invoking-non-root-uid-gid'
                or self.operator_default_deadline != 'none' or self.optional_analysis_required is not False
                or self.primary_required_roles != PRIMARY_REQUIRED_ROLES
                or self.default_preparation_run_options_sha256 != PreparationRunOptions(False).sha256
                or self.audited_preparation_run_options_sha256 != PreparationRunOptions(True).sha256
                or re.fullmatch(r'sha256:[0-9a-f]{64}', self.runtime_image_digest) is None):
            raise WgsOperatorPackageError('invalid installed contract')
        for name, value in asdict(self).items():
            if name.endswith('_sha256') and (type(value) is not str or re.fullmatch(r'[0-9a-f]{64}', value) is None):
                raise WgsOperatorPackageError('invalid installed dependency')

    def as_bytes(self):
        return canonical(asdict(self))

def _document(content):
    value = json.loads(content)
    if type(value) is not dict or canonical(value) != content:
        raise WgsOperatorPackageError("invalid canonical installed document")
    return value


def authenticate_package(root, expected):
    """Authenticate code and contracts without case or resource payload reads."""
    from genosapiens.evidence_model import _typed
    from genosapiens.track1_release import EVIDENCE_INTERFACE_FILES
    try:
        root = Path(root)
        content = read_bounded_bytes(root / "package-bindings.json", reference_root=root)
        manifest = _document(content)
        if (digest(content) != expected or set(manifest) != {"schema", "version", "classification", "contract_sha256", "files"}
                or manifest["schema"] != SCHEMA or type(manifest["version"]) is not int or manifest["version"] != 1
                or manifest["classification"] != "public"
                or not 0 < len(manifest["files"]) <= 4096):
            raise ValueError
        files = {}
        total = 0
        for row in manifest["files"]:
            name = row["relative_path"]
            path = Path(name)
            if (set(row) != {"relative_path", "sha256", "size_bytes", "mode"}
                    or path.is_absolute() or ".." in path.parts or str(path) != name or name in files
                    or name == "package-bindings.json" or row["mode"] not in (0o400, 0o500)
                    or type(row["size_bytes"]) is not int or not 0 < row["size_bytes"] <= 16 * 1024**2):
                raise ValueError
            total += row["size_bytes"]
            if total > 256 * 1024**2:
                raise ValueError
            data = read_bounded_bytes(root / path, reference_root=root, maximum_bytes=16 * 1024**2,
                                      expected_size_bytes=row["size_bytes"])
            if digest(data) != row["sha256"] or stat.S_IMODE((root / path).stat().st_mode) != row["mode"]:
                raise ValueError
            files[name] = row
        if tuple(files) != tuple(sorted(files)):
            raise ValueError
        def read(name):
            row = files[name]
            data = read_bounded_bytes(root / name, reference_root=root, maximum_bytes=16 * 1024**2,
                                      expected_size_bytes=row["size_bytes"])
            if digest(data) != row["sha256"]:
                raise ValueError
            return data
        raw = read("operator-contract.json")
        if digest(raw) != manifest["contract_sha256"]:
            raise ValueError
        contract = _typed(_document(raw), OperatorContract, decode=True)
        from genosapiens.installed_public_index import PIN_PATH, validate_index_pin
        if PIN_PATH in files:
            pin = validate_index_pin(_document(read(PIN_PATH)))
            if pin['manifest_sha256'] != contract.common_allele_public_index_manifest_sha256:
                raise ValueError
        from genosapiens.track1_early_scope import validate_scope_pin, load_early_scope_policy
        scope_pin = validate_scope_pin(_document(read('track1-early-scope-pin.json')))
        if (scope_pin['manifest_sha256'] != contract.scope_manifest_sha256
                or load_early_scope_policy(root).sha256 != contract.scope_policy_sha256):
            raise ValueError
        release_bytes = read("base-release.json")
        release = load_run_versions(release_bytes, expected_sha256=contract.base_release_sha256)
        stack_bytes = read("executable-stack.json")
        stack = _document(stack_bytes)
        if (type(release) is not InstalledVersions
                or digest(stack_bytes) != contract.executable_stack_sha256 or release.executable_stack_sha256 != contract.executable_stack_sha256
                or release.runtime_image_digest != contract.runtime_image_digest
                or stack["runtime_image"].split("@")[-1] != contract.runtime_image_digest
):
            raise ValueError
        if any(value not in release.resource_manifest_sha256s for value in
                (contract.scope_policy_sha256, contract.scope_manifest_sha256)):
            raise ValueError
        sources = {row["path"]: row for row in stack["files"]}
        for name, row in sources.items():
            destination = ("python/" + name.removeprefix("src/") if name.startswith("src/genosapiens/")
                else name if name.startswith(("config/", "schemas/", "docs/"))
                else "check_bounded_host_capacity.py" if name == "scripts/operator_compatibility_repair/check_bounded_host_capacity.py"
                else name if name in ("track1", "scripts/resolve_track1_resume.py") else None)
            if destination is not None and (files[destination]["sha256"] != row["sha256"]
                    or files[destination]["size_bytes"] != row["size_bytes"]):
                raise ValueError
        interface = digest(canonical({name: sources["src/genosapiens/" + name]["sha256"] for name in EVIDENCE_INTERFACE_FILES}))
        if interface != release.evidence_interface_sha256 or interface != contract.evidence_interface_sha256:
            raise ValueError
        if (digest(read("schemas/track1/mechanism-ready-handoff.schema.json"))
                != release.canonical_handoff_schema_sha256
                or contract.canonical_handoff_schema_sha256 != release.canonical_handoff_schema_sha256
                or digest(read("bounded-wgs-execution.json")) != contract.execution_limits_sha256
                or contract.execution_limits_sha256 != release.policy_bindings.bounded_wgs_execution_sha256
                or digest(read("common-allele-prefilter.json")) != contract.common_allele_policy_sha256
                or contract.common_allele_policy_sha256 != release.policy_bindings.common_allele_policy_sha256):
            raise ValueError
        return contract, release, release_bytes, files, read
    except (OSError, ValueError, TypeError, KeyError, RuntimeError):
        raise WgsOperatorPackageError("installed package binding is invalid") from None

def build_package(*, repository_root, destination, reference_versions, resource_versions, policy_bindings,
                  scope_pin, public_index_pin=None):
    """Copy source files and bind required versions without case execution."""
    from genosapiens.track1_public_validation_coordinator import executable_stack_bytes
    from genosapiens.wgs_analysis_policy import current_interfaces
    from genosapiens.wgs_work_contract import load_focused_work_contract
    from genosapiens.wgs_preparation import PreparationRunOptions
    repository, destination = Path(repository_root), Path(destination)
    if destination.exists():
        raise WgsOperatorPackageError('installation destination already exists')
    stack_bytes = executable_stack_bytes(repository)
    stack = _document(stack_bytes)
    work = load_focused_work_contract(repository)
    interfaces = current_interfaces(repository_root=repository)
    reference_fields = ('reference_identity_id', 'reference_manifest_sha256',
        'reference_bindings_sha256', 'assessment_scope_sha256', 'evaluator_pin_sha256',
        'canonical_handoff_path_sha256', 'path_selection_policy_sha256')
    resource_fields = ('population_manifest_sha256', 'clinical_manifest_sha256',
        'exception_manifest_sha256', 'precision_manifest_sha256',
        'common_allele_public_index_manifest_sha256', 'ontology_manifest_sha256',
        'scope_policy_sha256', 'scope_manifest_sha256')
    if set(reference_versions) != set(reference_fields) or set(resource_versions) != set(resource_fields):
        raise WgsOperatorPackageError('installation version inventory differs')
    if public_index_pin is not None:
        from genosapiens.installed_public_index import PIN_PATH, validate_index_pin
        validate_index_pin(public_index_pin)
        if (public_index_pin['manifest_sha256'] != resource_versions['common_allele_public_index_manifest_sha256']
                or public_index_pin['manifest_sha256'] != policy_bindings.common_allele_public_index_manifest_sha256):
            raise WgsOperatorPackageError('installation index versions differ')
    from genosapiens.track1_early_scope import validate_scope_pin, load_early_scope_policy
    validate_scope_pin(scope_pin)
    if (scope_pin['manifest_sha256'] != resource_versions['scope_manifest_sha256']
            or load_early_scope_policy(repository).sha256 != resource_versions['scope_policy_sha256']):
        raise WgsOperatorPackageError('installation scope versions differ')
    def file_sha(name):
        return digest((repository / name).read_bytes())
    selectors = dict(
        analysis_eligibility_policy_sha256=file_sha('config/track1/analysis-eligibility-v1.json'),
        context_writer_profile_sha256=file_sha('config/track1/context-writer-runtime-v1.json'),
        analysis_profile_policy_sha256=file_sha('config/track1/analysis-profiles.json'),
        population_work_policy_sha256=work.policy.sha256,
        focused_work_contract_sha256=work.sha256)
    resources = tuple(sorted({*resource_versions.values(), *policy_bindings.component_sha256s,
        *work.bindings().values(), *selectors.values(),
        file_sha('config/public-resources/gnomad-joint-faf95-precision-v1.json')}))
    release = InstalledVersions(**reference_versions, policy_bindings=policy_bindings,
        executable_stack_sha256=digest(stack_bytes), runtime_image_digest=stack['runtime_image'].split('@')[-1],
        resource_manifest_sha256s=resources, **interfaces)
    contract = OperatorContract(base_release_sha256=release.sha256,
        executable_stack_sha256=release.executable_stack_sha256, runtime_image_digest=release.runtime_image_digest,
        execution_limits_sha256=policy_bindings.bounded_wgs_execution_sha256,
        common_allele_policy_sha256=policy_bindings.common_allele_policy_sha256,
        default_preparation_run_options_sha256=PreparationRunOptions(False).sha256,
        audited_preparation_run_options_sha256=PreparationRunOptions(True).sha256,
        shared_driver_fields=(), analysis_profile='focused', **selectors, **interfaces, **resource_versions)
    destination.mkdir(mode=0o700)
    _write(destination / 'track1-early-scope-pin.json', canonical(scope_pin), 0o400)
    if public_index_pin is not None:
        _write(destination / PIN_PATH, canonical(public_index_pin), 0o400)
    for row in stack["files"]:
        name = row["path"]
        target = ("python/" + name.removeprefix("src/") if name.startswith("src/genosapiens/")
            else name if name.startswith(("config/", "schemas/", "docs/"))
            else "check_bounded_host_capacity.py" if name == "scripts/operator_compatibility_repair/check_bounded_host_capacity.py"
                else name if name in ("track1", "scripts/resolve_track1_resume.py") else None)
        if target is None:
            continue
        raw = (repository / name).read_bytes()
        if digest(raw) != row["sha256"] or len(raw) != row["size_bytes"]:
            raise WgsOperatorPackageError("installed source changed during build")
        _write(destination / target, raw, 0o500 if target == "track1" or target.endswith(".py") and "/" not in target else 0o400)
    for name, raw in {"operator-contract.json": contract.as_bytes(), "base-release.json": release.as_bytes(),
            "executable-stack.json": stack_bytes,
            "bounded-wgs-execution.json": (repository / "config/track1/bounded-wgs-execution-v5.json").read_bytes(),
            "common-allele-prefilter.json": (repository / "config/track1/common-allele-prefilter-v1.json").read_bytes()}.items():
        _write(destination / name, raw, 0o400)
    files = tuple(dict(relative_path=str(path.relative_to(destination)), sha256=digest(path.read_bytes()),
        size_bytes=path.stat().st_size, mode=stat.S_IMODE(path.stat().st_mode))
        for path in sorted(destination.rglob('*')) if path.is_file())
    manifest = canonical(dict(schema=SCHEMA, version=1, classification='public',
        contract_sha256=digest(contract.as_bytes()), files=files))
    _write(destination / 'package-bindings.json', manifest, 0o400)
    authenticate_package(destination, digest(manifest))
    return digest(manifest)
