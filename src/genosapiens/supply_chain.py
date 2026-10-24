"""Supply-chain policy validation and deterministic image inventories."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import shlex
import sys
from typing import Any, Mapping, Sequence

from genosapiens import __version__


LOCK_SCHEMA = "genosapiens.dependency-lock"
POLICY_SCHEMA = "genosapiens.supply-chain-policy"
PINNED_VERSION = re.compile(r"^[0-9]+(?:\.[0-9A-Za-z]+)+(?:[-+][0-9A-Za-z.]+)?$")
SHA256 = re.compile(r"^[0-9a-f]{64}$")
REQUIREMENT = re.compile(r"^([A-Za-z0-9_.-]+)==([^=<>!~\s]+)$")


class SupplyChainError(ValueError):
    """Raised when a dependency or bundled resource violates policy."""


def _json_object(path: str | Path) -> dict[str, Any]:
    source = Path(path)
    try:
        value = json.loads(source.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise SupplyChainError(f"cannot load {source}: {error}") from error
    if not isinstance(value, dict):
        raise SupplyChainError(f"{source} must contain a JSON object")
    return value


def _object_list(value: Any, field: str) -> list[Mapping[str, Any]]:
    if not isinstance(value, list) or any(not isinstance(item, Mapping) for item in value):
        raise SupplyChainError(f"{field} must be an array of objects")
    return value


def validate_lock(
    lock: Mapping[str, Any],
    policy: Mapping[str, Any],
    expected_role: str,
) -> None:
    if lock.get("schema") != LOCK_SCHEMA:
        raise SupplyChainError(f"lock schema must be {LOCK_SCHEMA!r}")
    if policy.get("schema") != POLICY_SCHEMA:
        raise SupplyChainError(f"policy schema must be {POLICY_SCHEMA!r}")
    if lock.get("image_role") != expected_role:
        raise SupplyChainError(f"lock image_role must be {expected_role!r}")
    if lock.get("python") != "3.12.11":
        raise SupplyChainError("Python runtime must be pinned to 3.12.11")

    allowed_licenses = policy.get("allowed_dependency_licenses")
    if not isinstance(allowed_licenses, list) or any(
        not isinstance(item, str) for item in allowed_licenses
    ):
        raise SupplyChainError("policy allowed_dependency_licenses must be strings")

    dependency_names: set[str] = set()
    for dependency in _object_list(lock.get("dependencies"), "dependencies"):
        name = dependency.get("name")
        version = dependency.get("version")
        hashes = dependency.get("sha256")
        license_id = dependency.get("license")
        if not isinstance(name, str) or not name:
            raise SupplyChainError("every dependency requires a name")
        normalized_name = re.sub(r"[-_.]+", "-", name).lower()
        if normalized_name in dependency_names:
            raise SupplyChainError(f"duplicate dependency: {name}")
        dependency_names.add(normalized_name)
        if not isinstance(version, str) or not PINNED_VERSION.fullmatch(version):
            raise SupplyChainError(f"dependency {name} has an unpinned version")
        if not isinstance(hashes, list) or not hashes or any(
            not isinstance(item, str) or not SHA256.fullmatch(item) for item in hashes
        ):
            raise SupplyChainError(f"dependency {name} requires SHA-256 hashes")
        if license_id not in allowed_licenses:
            raise SupplyChainError(f"dependency {name} has unapproved license {license_id!r}")

    for tool in _object_list(lock.get("tools", []), "tools"):
        name = tool.get("name")
        version = tool.get("version")
        hashes = tool.get("sha256")
        license_id = tool.get("license")
        if not isinstance(name, str) or not name:
            raise SupplyChainError("every tool requires a name")
        normalized_name = re.sub(r"[-_.]+", "-", name).lower()
        if normalized_name in dependency_names:
            raise SupplyChainError(f"duplicate dependency or tool: {name}")
        dependency_names.add(normalized_name)
        if not isinstance(version, str) or not PINNED_VERSION.fullmatch(version):
            raise SupplyChainError(f"tool {name} has an unpinned version")
        if not isinstance(hashes, list) or not hashes or any(
            not isinstance(item, str) or not SHA256.fullmatch(item) for item in hashes
        ):
            raise SupplyChainError(f"tool {name} requires SHA-256 hashes")
        if license_id not in allowed_licenses:
            raise SupplyChainError(f"tool {name} has unapproved license {license_id!r}")

    resource_policy = policy.get("annotation_resources")
    if not isinstance(resource_policy, Mapping):
        raise SupplyChainError("policy annotation_resources must be an object")
    resource_names: set[str] = set()
    for resource in _object_list(
        lock.get("annotation_resources"), "annotation_resources"
    ):
        name = resource.get("name")
        version = resource.get("version")
        distribution = resource.get("distribution")
        if not isinstance(name, str) or not name:
            raise SupplyChainError("every annotation resource requires a name")
        if name in resource_names:
            raise SupplyChainError(f"duplicate annotation resource: {name}")
        resource_names.add(name)
        approved = resource_policy.get(name)
        if not isinstance(approved, Mapping):
            raise SupplyChainError(f"annotation resource {name} is not approved")
        if not isinstance(version, str) or not PINNED_VERSION.fullmatch(version):
            raise SupplyChainError(f"annotation resource {name} has an unpinned version")
        allowed_distributions = approved.get("allowed_distributions")
        if not isinstance(allowed_distributions, list) or distribution not in allowed_distributions:
            raise SupplyChainError(
                f"annotation resource {name} cannot be distributed as {distribution!r}"
            )


def _requirement_records(
    requirements: str,
) -> dict[str, tuple[str, tuple[str, ...]]]:
    installed: dict[str, tuple[str, tuple[str, ...]]] = {}
    for number, raw_line in enumerate(requirements.splitlines(), start=1):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        try:
            tokens = shlex.split(line)
        except ValueError as error:
            raise SupplyChainError(f"malformed requirement line {number}: {error}") from error
        match = REQUIREMENT.fullmatch(tokens[0]) if tokens else None
        hashes = tuple(
            sorted(token.removeprefix("--hash=sha256:") for token in tokens[1:])
        )
        if (
            match is None
            or not hashes
            or any(
                not token.startswith("--hash=sha256:") or not SHA256.fullmatch(value)
                for token, value in zip(tokens[1:], hashes, strict=True)
            )
        ):
            raise SupplyChainError(f"requirement line {number} must pin version and SHA-256")
        name, version = match.groups()
        normalized = re.sub(r"[-_.]+", "-", name).lower()
        if normalized in installed:
            raise SupplyChainError(f"duplicate requirement: {name}")
        installed[normalized] = (version, hashes)
    return installed


def _locked_requirement_records(
    lock: Mapping[str, Any],
) -> dict[str, tuple[str, tuple[str, ...]]]:
    locked: dict[str, tuple[str, tuple[str, ...]]] = {}
    for dependency in _object_list(lock.get("dependencies"), "dependencies"):
        name = dependency.get("name")
        version = dependency.get("version")
        hashes = dependency.get("sha256")
        if not isinstance(name, str) or not isinstance(version, str) or not isinstance(hashes, list):
            raise SupplyChainError("dependency lock must be validated before requirements")
        normalized = re.sub(r"[-_.]+", "-", name).lower()
        locked[normalized] = (version, tuple(sorted(str(value) for value in hashes)))
    return locked


def validate_requirements(requirements: str, lock: Mapping[str, Any]) -> None:
    installed = _requirement_records(requirements)
    locked = _locked_requirement_records(lock)
    if installed != locked:
        raise SupplyChainError("installed requirements differ from the dependency lock")


def validate_requirement_subset(requirements: str, lock: Mapping[str, Any]) -> None:
    installed = _requirement_records(requirements)
    locked = _locked_requirement_records(lock)
    if any(locked.get(name) != record for name, record in installed.items()):
        raise SupplyChainError("build requirements differ from the dependency lock")


def validate_files(
    lock_path: str | Path,
    policy_path: str | Path,
    role: str,
    requirements_path: str | Path | None = None,
    build_requirements_path: str | Path | None = None,
) -> None:
    lock = _json_object(lock_path)
    validate_lock(lock, _json_object(policy_path), role)
    if requirements_path is not None:
        source = Path(requirements_path)
        try:
            requirements = source.read_text(encoding="utf-8")
        except OSError as error:
            raise SupplyChainError(f"cannot load {source}: {error}") from error
        validate_requirements(requirements, lock)
    if build_requirements_path is not None:
        source = Path(build_requirements_path)
        try:
            requirements = source.read_text(encoding="utf-8")
        except OSError as error:
            raise SupplyChainError(f"cannot load {source}: {error}") from error
        validate_requirement_subset(requirements, lock)


def _debian_packages() -> list[dict[str, Any]]:
    status_path = Path("/var/lib/dpkg/status")
    if not status_path.exists():
        return []
    packages: list[dict[str, Any]] = []
    for paragraph in status_path.read_text(encoding="utf-8").split("\n\n"):
        fields: dict[str, str] = {}
        for line in paragraph.splitlines():
            if ": " in line:
                key, value = line.split(": ", 1)
                fields[key] = value
        if fields.get("Status") != "install ok installed":
            continue
        name = fields.get("Package")
        version = fields.get("Version")
        if not name or not version:
            continue
        license_evidence: list[str] = []
        copyright_path = Path("/usr/share/doc") / name / "copyright"
        if copyright_path.exists():
            try:
                for line in copyright_path.read_text(
                    encoding="utf-8", errors="replace"
                ).splitlines():
                    if line.startswith("License: "):
                        license_evidence.append(line.removeprefix("License: ").strip())
            except OSError:
                pass
        packages.append(
            {
                "name": name,
                "version": version,
                "architecture": fields.get("Architecture"),
                "license_evidence": sorted(set(license_evidence)),
            }
        )
    return sorted(packages, key=lambda item: (item["name"], item["version"]))


def build_inventory(image_name: str, lock: Mapping[str, Any]) -> dict[str, Any]:
    application_name = "genosapiens-track1"
    return {
        "schema": "genosapiens.image-inventory",
        "image": image_name,
        "application": {
            "name": application_name,
            "version": __version__,
            "license": "CC-BY-4.0",
        },
        "python": lock["python"],
        "locked_python_dependencies": lock["dependencies"],
        "locked_tools": lock.get("tools", []),
        "annotation_resources": lock["annotation_resources"],
        "os_packages": _debian_packages(),
    }


def _spdx_id(name: str, version: str) -> str:
    value = hashlib.sha256(f"{name}\0{version}".encode()).hexdigest()[:16]
    return f"SPDXRef-Package-{value}"


def build_spdx(image_name: str, inventory: Mapping[str, Any]) -> dict[str, Any]:
    fingerprint = hashlib.sha256(
        json.dumps(inventory, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    application_id = "SPDXRef-Application"
    python_id = "SPDXRef-Python"
    packages: list[dict[str, Any]] = [
        {
            "SPDXID": application_id,
            "name": inventory["application"]["name"],
            "versionInfo": __version__,
            "downloadLocation": "NOASSERTION",
            "filesAnalyzed": False,
            "licenseConcluded": "CC-BY-4.0",
            "licenseDeclared": "CC-BY-4.0",
        },
        {
            "SPDXID": python_id,
            "name": "python",
            "versionInfo": str(inventory["python"]),
            "downloadLocation": "NOASSERTION",
            "filesAnalyzed": False,
            "licenseConcluded": "PSF-2.0",
            "licenseDeclared": "PSF-2.0",
        },
    ]
    relationships: list[dict[str, str]] = [
        {
            "spdxElementId": "SPDXRef-DOCUMENT",
            "relationshipType": "DESCRIBES",
            "relatedSpdxElement": application_id,
        },
        {
            "spdxElementId": application_id,
            "relationshipType": "DEPENDS_ON",
            "relatedSpdxElement": python_id,
        },
    ]
    for component in (
        list(inventory["locked_python_dependencies"])
        + list(inventory.get("locked_tools", []))
    ):
        component_id = _spdx_id(component["name"], component["version"])
        packages.append(
            {
                "SPDXID": component_id,
                "name": component["name"],
                "versionInfo": component["version"],
                "downloadLocation": "NOASSERTION",
                "filesAnalyzed": False,
                "licenseConcluded": component["license"],
                "licenseDeclared": component["license"],
                "checksums": [
                    {"algorithm": "SHA256", "checksumValue": checksum}
                    for checksum in component["sha256"]
                ],
            }
        )
        relationships.append(
            {
                "spdxElementId": application_id,
                "relationshipType": "DEPENDS_ON",
                "relatedSpdxElement": component_id,
            }
        )
    for package in inventory["os_packages"]:
        package_id = _spdx_id(package["name"], package["version"])
        packages.append(
            {
                "SPDXID": package_id,
                "name": package["name"],
                "versionInfo": package["version"],
                "downloadLocation": "NOASSERTION",
                "filesAnalyzed": False,
                "licenseConcluded": "NOASSERTION",
                "licenseDeclared": "NOASSERTION",
                "comment": "Debian license evidence: "
                + (", ".join(package["license_evidence"]) or "not parsed"),
            }
        )
        relationships.append(
            {
                "spdxElementId": python_id,
                "relationshipType": "DEPENDS_ON",
                "relatedSpdxElement": package_id,
            }
        )
    created = datetime.fromtimestamp(0, tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    return {
        "spdxVersion": "SPDX-2.3",
        "dataLicense": "CC0-1.0",
        "SPDXID": "SPDXRef-DOCUMENT",
        "name": f"{image_name}-sbom",
        "documentNamespace": f"https://genosapiens.invalid/spdx/{fingerprint}",
        "creationInfo": {
            "created": created,
            "creators": ["Tool: genosapiens-supply-chain-0.1.0"],
        },
        "packages": packages,
        "relationships": relationships,
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="genosapiens-supply-chain")
    commands = parser.add_subparsers(dest="command", required=True)
    validate = commands.add_parser("validate")
    validate.add_argument("--lock", required=True)
    validate.add_argument("--policy", required=True)
    validate.add_argument("--role", required=True)
    validate.add_argument("--requirements")
    validate.add_argument("--build-requirements")
    for command in ("inventory", "sbom"):
        output = commands.add_parser(command)
        output.add_argument("--lock", required=True)
        output.add_argument("--image", required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        if args.command == "validate":
            validate_files(
                args.lock,
                args.policy,
                args.role,
                args.requirements,
                args.build_requirements,
            )
            return 0
        lock = _json_object(args.lock)
        inventory = build_inventory(args.image, lock)
        value = inventory if args.command == "inventory" else build_spdx(args.image, inventory)
        json.dump(value, sys.stdout, indent=2, sort_keys=True)
        sys.stdout.write("\n")
        return 0
    except SupplyChainError as error:
        print(f"supply-chain error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
