"""Preserve the known failed preparation before a human-authorized retry."""

from pathlib import Path
from pathlib import PurePosixPath
import argparse
import hashlib
import json
import os
import stat

DRIVER = "2f3163441a15dfe1efd06cd84c8cd31430045afed87c606f9bbd7f84aa994e31"
SUPPLEMENT = "511ba6d8b2352334ce55e1c8fb94ea38daf20ed247af2039bbeb55b0baec3e16"
EMPTY_SHA256 = "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"


def _sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as source:
        while chunk := source.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


_INTERRUPTED_FILES = (
    "original-input-bindings.json",
    "case-manifest.json",
    "reference-declaration-clarification.json",
    "complete-concordance.json",
    "operator-logs/preparation.stdout",
    "operator-logs/preparation.stderr",
)
_INTERRUPTED_INPUTS = (
    "HPO_list.txt",
    "WGS_EX2312012_HGWCNDSX7.vcf.gz",
    "WGS_EX2312012_HGWCNDSX7.vcf.gz.tbi",
)


def _private_file(path, *, maximum_bytes=None):
    metadata = path.lstat()
    if (path.is_symlink() or not stat.S_ISREG(metadata.st_mode)
            or metadata.st_uid != os.getuid()
            or stat.S_IMODE(metadata.st_mode) != 0o600
            or metadata.st_nlink != 1
            or maximum_bytes is not None and metadata.st_size > maximum_bytes):
        raise ValueError("Incident file is not private and immutable.")


def _private_directory(path):
    metadata = path.lstat()
    if (path.is_symlink() or not stat.S_ISDIR(metadata.st_mode)
            or metadata.st_uid != os.getuid()
            or stat.S_IMODE(metadata.st_mode) != 0o700):
        raise ValueError("Incident directory is not private and caller-owned.")


def _unique_json(path, *, maximum_bytes):
    _private_file(path, maximum_bytes=maximum_bytes)
    def unique(pairs):
        value = dict(pairs)
        if len(value) != len(pairs):
            raise ValueError("Incident JSON has duplicate fields.")
        return value
    value = json.loads(path.read_bytes(), object_pairs_hook=unique)
    if type(value) is not dict:
        raise ValueError("Incident JSON must be an object.")
    return value


def _safe_relative(value):
    path = PurePosixPath(value)
    return (type(value) is str and not path.is_absolute()
            and path != PurePosixPath(".") and ".." not in path.parts)


def _inventory_closure(state_root, workspace_name):
    inventory_path = state_root / "inventory-preparation.json"
    inventory = _unique_json(inventory_path, maximum_bytes=16 * 1024**2)
    artifacts = inventory.get("artifacts")
    if (inventory.get("schema") != "genosapiens.case-inventory"
            or type(inventory.get("case_id")) is not str
            or type(artifacts) is not list):
        raise ValueError("Interrupted preparation inventory is invalid.")
    # The predecessor driver sealed this inventory immediately after copying
    # inputs, before it created the case manifest and concordance report.
    expected = {
        "original-input-bindings.json",
        "operator-logs/preparation.stdout",
        "operator-logs/preparation.stderr",
        *("inputs/" + name for name in _INTERRUPTED_INPUTS),
    }
    records = {}
    for record in artifacts:
        if type(record) is not dict or not _safe_relative(record.get("relative_path")):
            raise ValueError("Interrupted preparation inventory is invalid.")
        relative = PurePosixPath(record["relative_path"])
        try:
            index = relative.parts.index(workspace_name)
        except ValueError:
            continue
        local = PurePosixPath(*relative.parts[index + 1:]).as_posix()
        if local in expected:
            if (local in records or record.get("classification") != "controlled"
                    or record.get("retention") != "delete"
                    or record.get("status") != "present"
                    or type(record.get("size_bytes")) is not int
                    or record["size_bytes"] < 0
                    or type(record.get("sha256")) is not str
                    or len(record["sha256"]) != 64
                    or any(value not in "0123456789abcdef"
                           for value in record["sha256"])):
                raise ValueError("Interrupted preparation inventory is invalid.")
            records[local] = record
    if set(records) != expected:
        raise ValueError("Interrupted preparation inventory closure is incomplete.")
    for relative, record in records.items():
        path = state_root / relative
        _private_file(path)
        if path.stat().st_size != record["size_bytes"] or _sha(path) != record["sha256"]:
            raise ValueError("Interrupted preparation product hash changed.")
    return _sha(inventory_path)


def _verify_second_interrupted_state(root, release):
    first = root / "preparation-interrupted-qc-memory-01"
    _private_directory(first)
    first_ledger_path = first / "archive.json"
    first_ledger = _unique_json(first_ledger_path, maximum_bytes=65536)
    moved = first_ledger.get("moved")
    if (first_ledger.get("base_release_sha256") != release
            or first_ledger.get("reason")
            != "preparation-terminated-during-full-vcf-qc"
            or type(moved) is not list
            or any(type(item) is not dict for item in moved)
            or {item.get("path") for item in moved}
            != {"inputs", "inventory-preparation.json", *_INTERRUPTED_FILES}):
        raise ValueError("The first QC interruption archive is unrelated.")
    _private_directory(first / "inputs")
    _private_directory(first / "operator-logs")
    _inventory_closure(first, root.name)

    for directory in (root / "inputs", root / "operator-logs"):
        _private_directory(directory)
    if {path.name for path in (root / "inputs").iterdir()} != set(_INTERRUPTED_INPUTS):
        raise ValueError("Interrupted preparation input set is not exact.")
    for name in (*_INTERRUPTED_FILES, *("inputs/" + item for item in _INTERRUPTED_INPUTS)):
        _private_file(root / name)
    for name in ("preparation.stdout", "preparation.stderr"):
        path = root / "operator-logs" / name
        if path.stat().st_size != 0 or _sha(path) != EMPTY_SHA256:
            raise ValueError("Interrupted preparation logs are not the known empty pair.")
    for name in (
        "controlled-preparation-result.json", "qc.json", "normalization.json",
        "supplied-inventory.json", "runs", "controlled-analysis-result.json",
        "inventory-analysis.json", "operator-packet-status.json",
        "preparation-interrupted-qc-memory-02",
    ):
        if (root / name).exists() or (root / name).is_symlink():
            raise ValueError("Unexpected or downstream products exist.")
    clarification = _unique_json(
        root / "reference-declaration-clarification.json", maximum_bytes=65536,
    )
    concordance = _unique_json(
        root / "complete-concordance.json", maximum_bytes=262144,
    )
    if (clarification.get("schema")
            != "genosapiens.track1-reference-declaration-clarification"
            or clarification.get("status") != "applied"
            or clarification.get("original_reference_declaration_count") != 2
            or clarification.get("alternate_reference_value_disclosed") is not False
            or concordance.get("schema") != "genosapiens.reference-concordance"
            or concordance.get("reference_role") != "challenge_alignment"
            or any(concordance.get(field) is not True for field in (
                "header_valid", "reference_basename_match", "contig_names_match",
                "contig_lengths_match", "contig_order_match", "accepted",
            ))
            or type(concordance.get("total_records")) is not int
            or concordance["total_records"] <= 0
            or concordance.get("checked_records") != concordance["total_records"]
            or concordance.get("mismatching_records") != 0
            or concordance.get("unresolved_records") != 0
            or clarification.get("original_input_sha256")
            != concordance.get("validated_input_sha256")
            or clarification.get("record_payload_sha256")
            != concordance.get("record_payload_sha256")):
        raise ValueError("Interrupted preparation receipts do not prove completed concordance.")
    inventory_sha256 = _inventory_closure(root, root.name)
    return first_ledger_path, inventory_sha256


def _directory_tree_sha256(path):
    _private_directory(path)
    records = []
    for item in sorted(path.rglob("*")):
        relative = item.relative_to(path).as_posix()
        if item.is_symlink():
            raise ValueError("Incident tree contains a symlink.")
        if item.is_dir():
            _private_directory(item)
            records.append({"path": relative, "kind": "directory"})
        else:
            _private_file(item)
            records.append({
                "path": relative,
                "kind": "file",
                "sha256": _sha(item),
                "size_bytes": item.stat().st_size,
            })
    content = (json.dumps(records, sort_keys=True, separators=(",", ":")) + "\n").encode()
    return hashlib.sha256(content).hexdigest()


def _archive_second_interrupted_qc(root, release):
    first_ledger_path, inventory_sha256 = _verify_second_interrupted_state(root, release)
    products = [
        ("inputs", root / "inputs", "directory"),
        ("inventory-preparation.json", root / "inventory-preparation.json", "file"),
        *((name, root / name, "file") for name in _INTERRUPTED_FILES),
    ]
    archive = root / "preparation-interrupted-qc-memory-02"
    archive.mkdir(mode=0o700)
    moved = []
    for name, source, kind in products:
        identity = (
            _directory_tree_sha256(source) if kind == "directory" else _sha(source)
        )
        destination = archive / name
        destination.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
        os.rename(source, destination)
        moved.append({"path": name, "kind": kind, "sha256": identity})
    ledger = {
        "schema": "genosapiens.track1-preparation-incident-archive",
        "version": 1,
        "classification": "controlled",
        "llm_policy": "no-llm",
        "retention": "delete",
        "incident_id": "preparation-interrupted-qc-memory-02",
        "base_release_sha256": release,
        "reason": "preparation-terminated-during-full-vcf-qc",
        "predecessor_incident_id": "preparation-interrupted-qc-memory-01",
        "predecessor_archive_sha256": _sha(first_ledger_path),
        "inventory_sha256": inventory_sha256,
        "moved": moved,
    }
    with (archive / "archive.json").open("x") as stream:
        os.chmod(stream.name, 0o600)
        json.dump(ledger, stream, sort_keys=True, separators=(",", ":"))
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    return "archived-interrupted-qc-02"


def _archive_interrupted_qc(root, release):
    if (root / "preparation-interrupted-qc-memory-01").exists():
        return _archive_second_interrupted_qc(root, release)
    first_archive = root / "preparation-before-reference-clarification-01/archive.json"
    if (first_archive.is_symlink() or not first_archive.is_file()
            or first_archive.stat().st_size > 65536):
        raise ValueError("The authenticated preceding preparation failure archive is required.")
    first = json.loads(first_archive.read_bytes())
    if (first.get("base_release_sha256") != release
            or first.get("reason") != "organizer-reference-declaration-clarification"):
        raise ValueError("The preceding preparation archive is unrelated.")
    logs = [root / "operator-logs/preparation.stdout",
            root / "operator-logs/preparation.stderr"]
    for path in logs:
        if (path.is_symlink() or not path.is_file() or path.stat().st_uid != os.getuid()
                or path.stat().st_size != 0 or _sha(path) != EMPTY_SHA256):
            raise ValueError("Interrupted preparation logs are not the known empty pair.")
    required = (
        "inputs", "original-input-bindings.json", "case-manifest.json",
        "inventory-preparation.json", "reference-declaration-clarification.json",
        "complete-concordance.json",
    )
    for name in required:
        path = root / name
        if (path.is_symlink() or not path.exists() or path.stat().st_uid != os.getuid()
                or (path.is_dir() and any(item.is_symlink() for item in path.rglob("*")))):
            raise ValueError("Interrupted preparation products are incomplete or unsafe.")
    inputs = root / "inputs"
    expected_inputs = {
        "HPO_list.txt", "WGS_EX2312012_HGWCNDSX7.vcf.gz",
        "WGS_EX2312012_HGWCNDSX7.vcf.gz.tbi",
    }
    if ({path.name for path in inputs.iterdir()} != expected_inputs
            or any(not path.is_file() for path in inputs.iterdir())):
        raise ValueError("Interrupted preparation input set is not exact.")
    if any((root / name).exists() for name in (
            "controlled-preparation-result.json", "qc.json", "normalization.json",
            "supplied-inventory.json", "runs", "controlled-analysis-result.json",
            "inventory-analysis.json", "operator-packet-status.json")):
        raise ValueError("Unexpected or downstream products exist; preserve and inspect them.")
    clarification_path = root / "reference-declaration-clarification.json"
    concordance_path = root / "complete-concordance.json"
    if clarification_path.stat().st_size > 65536 or concordance_path.stat().st_size > 262144:
        raise ValueError("Interrupted preparation technical receipts are oversized.")
    clarification = json.loads(clarification_path.read_bytes())
    concordance = json.loads(concordance_path.read_bytes())
    if (clarification.get("schema") != "genosapiens.track1-reference-declaration-clarification"
            or clarification.get("status") != "applied"
            or clarification.get("original_reference_declaration_count") != 2
            or clarification.get("alternate_reference_value_disclosed") is not False
            or concordance.get("schema") != "genosapiens.reference-concordance"
            or concordance.get("reference_role") != "challenge_alignment"
            or concordance.get("header_valid") is not True
            or concordance.get("reference_basename_match") is not True
            or concordance.get("contig_names_match") is not True
            or concordance.get("contig_lengths_match") is not True
            or concordance.get("contig_order_match") is not True
            or concordance.get("accepted") is not True
            or concordance.get("total_records") != concordance.get("checked_records")
            or not isinstance(concordance.get("total_records"), int)
            or concordance["total_records"] <= 0
            or concordance.get("mismatching_records") != 0
            or concordance.get("unresolved_records") != 0
            or clarification.get("original_input_sha256") != concordance.get("validated_input_sha256")
            or clarification.get("record_payload_sha256") != concordance.get("record_payload_sha256")):
        raise ValueError("Interrupted preparation receipts do not prove completed concordance.")
    products = []
    for name in required:
        source = root / name
        products.append((name, source, "directory" if source.is_dir() else "file",
                         None if source.is_dir() else _sha(source)))
    products.extend(("operator-logs/" + source.name, source, "file", EMPTY_SHA256)
                    for source in logs)
    archive = root / "preparation-interrupted-qc-memory-01"
    archive.mkdir(mode=0o700)
    moved = []
    for name, source, kind, source_sha256 in products:
        destination = archive / name
        destination.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
        os.rename(source, destination)
        moved.append({"path": name, "kind": kind, "sha256": source_sha256})
    with (archive / "archive.json").open("x") as stream:
        os.chmod(stream.name, 0o600)
        json.dump({"classification": "controlled", "llm_policy": "no-llm",
                   "retention": "delete", "base_release_sha256": release,
                   "reason": "preparation-terminated-during-full-vcf-qc",
                   "moved": moved}, stream, sort_keys=True, indent=2)
        stream.write("\n")
    return "archived-interrupted-qc"


def archive_failure(root, release):
    root = Path(root)
    if (root.name != "track1-base-20260917-" + release[:16]
            or any(path.is_symlink() for path in (root, *root.parents))
            or stat.S_IMODE(root.stat().st_mode) != 0o700
            or root.stat().st_uid != os.getuid()):
        raise ValueError("Preparation workspace is not a private release-bound directory.")
    receipt = root / "controlled-preparation-result.json"
    if not receipt.exists():
        active = (root / "inventory-preparation.json").exists() or any(
            (root / "operator-logs" / name).exists()
            for name in ("preparation.stdout", "preparation.stderr")
        )
        return _archive_interrupted_qc(root, release) if active else "not-needed"
    if receipt.is_symlink() or not receipt.is_file() or receipt.stat().st_size > 65536:
        raise ValueError("Preparation failure receipt is invalid.")
    record = json.loads(receipt.read_bytes())
    required = {
        "schema": "genosapiens.track1-controlled-preparation",
        "classification": "controlled", "llm_policy": "no-llm",
        "retention": "delete", "status": "failed", "rescue_authorized": False,
        "challenge_variant_record_accessed": True,
        "early_snapshot_published": False,
        "failed_stage": "complete-original-reference-concordance",
        "error_class": "ReferenceConcordanceError",
        "matching_diagnosis_supplement_receipt_sha256": SUPPLEMENT,
        "driver_sha256": DRIVER,
    }
    if any(record.get(key) != value for key, value in required.items()):
        raise ValueError("Only the authenticated reference-concordance failure may be archived.")
    if record.get("base_release_sha256") is not None:
        raise ValueError("A release-bound passing preparation result cannot be archived.")
    if any((root / name).exists() for name in (
            "controlled-analysis-result.json", "inventory-analysis.json",
            "operator-packet-status.json")):
        raise ValueError("Downstream results exist; preserve and inspect them before retrying.")
    log_paths = [root / "operator-logs/preparation.stdout",
                 root / "operator-logs/preparation.stderr"]
    for path in [receipt, *log_paths]:
        if (path.is_symlink() or not path.is_file()
                or path.stat().st_uid != os.getuid()):
            raise ValueError("Preparation failure artifacts must be regular caller-owned files.")
    candidates = (
        "controlled-preparation-result.json", "inputs", "original-input-bindings.json",
        "case-manifest.json", "inventory-preparation.json", "complete-concordance.json",
        "reference-declaration-clarification.json", "normalization.json",
        "supplied-inventory.json", "qc.json", "runs",
    )
    products = []
    for name in candidates:
        source = root / name
        if source.is_symlink():
            raise ValueError("Preparation products must be nonsymlinked.")
        if not source.exists():
            continue
        if (source.stat().st_uid != os.getuid()
                or (source.is_dir() and any(path.is_symlink() for path in source.rglob("*")))):
            raise ValueError("Preparation products must be caller-owned and nonsymlinked.")
        products.append((name, source, "directory" if source.is_dir() else "file",
                         None if source.is_dir() else _sha(source)))
    products.extend(("operator-logs/" + source.name, source, "file", _sha(source))
                    for source in log_paths)
    archive = root / "preparation-before-reference-clarification-01"
    archive.mkdir(mode=0o700)
    archive_logs = archive / "operator-logs"
    archive_logs.mkdir(mode=0o700)
    moved = []
    for name, source, kind, source_sha256 in products:
        destination = archive / name
        destination.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
        os.rename(source, destination)
        moved.append({"path": name, "kind": kind, "sha256": source_sha256})
    ledger = archive / "archive.json"
    with ledger.open("x") as stream:
        os.chmod(stream.name, 0o600)
        json.dump({"classification": "controlled", "llm_policy": "no-llm",
                   "retention": "delete", "base_release_sha256": release,
                   "reason": "organizer-reference-declaration-clarification",
                   "moved": moved}, stream, sort_keys=True, indent=2)
        stream.write("\n")
    return "archived"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workspace", type=Path)
    parser.add_argument("release")
    arguments = parser.parse_args()
    status = archive_failure(arguments.workspace, arguments.release)
    print(json.dumps({"schema": "genosapiens.track1-preparation-failure-archive-status",
                      "status": status, "controlled_products": True,
                      "llm_policy": "no-llm"}))


if __name__ == "__main__":
    main()
