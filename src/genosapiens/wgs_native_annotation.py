"""Bounded native annotation publications, reused verbatim by scientific replay."""

from dataclasses import replace
import hashlib
import json
import os
import re
from pathlib import Path

from genosapiens.annotation_admission import validate_annotation_admission
from genosapiens.bounded_json import read_bounded_bytes
from genosapiens.bounded_wgs_config import bounded_wgs_execution_config_as_json
from genosapiens.evidence_model import AlleleEvidence
from genosapiens.evidence_serialization import from_json
from genosapiens.vep_annotation import VepRequest, run_vep_annotation
from genosapiens.wgs_ledger import Artifact, canonical, iter_artifact_lines, publish_bytes, MAX_ROW_BYTES


class NativeAnnotationError(RuntimeError):
    """Static native publication diagnostic."""


def native_index_rows(root, index):
    """Read finite descriptors only after authenticating the entire index."""
    rows = iter_artifact_lines(root, index, maximum_line_bytes=MAX_ROW_BYTES)
    next(rows)  # The scientific replay checks the parent header against its inputs.
    previous = 0
    for content in rows:
        row = json.loads(content)
        if (set(row) != {"ordinal", "identity_id", "directory", "artifacts"}
                or type(row["ordinal"]) is not int or row["ordinal"] <= previous
                or row["directory"] != f'native-{row["ordinal"]:012d}'
                or type(row["artifacts"]) is not list or not row["artifacts"]):
            raise NativeAnnotationError("native annotation index is invalid")
        names = set()
        for item in row["artifacts"]:
            if (set(item) != {"name", "sha256", "size_bytes"}
                    or re.fullmatch(r"[a-zA-Z0-9][a-zA-Z0-9_.-]*", item["name"]) is None
                    or item["name"] in names or re.fullmatch(r"[0-9a-f]{64}", item["sha256"]) is None
                    or type(item["size_bytes"]) is not int or item["size_bytes"] < 0):
                raise NativeAnnotationError("native annotation artifact descriptor is invalid")
            names.add(item["name"])
        previous = row["ordinal"]
        yield row


class NativeAnnotationSession:
    """One native allele at a time; an authenticated index closes exact coverage."""

    def __init__(self, *, root, resources, config, runtime_image, contract_path,
                 bundle_sha256, plan_sha256, stage_sha256, fixed=False, index=None):
        self.root, self.resources, self.config = Path(root), resources, config
        self.runtime_image, self.contract_path = runtime_image, contract_path
        self.bundle_sha256, self.fixed, self.index = bundle_sha256, fixed, index
        self.header = {"schema": "genosapiens.wgs-native-annotation-index", "version": 1,
            "plan_sha256": plan_sha256, "stage_sha256": stage_sha256}
        self.stream, self.rows, self.digest, self.size = None, None, hashlib.sha256(), 0
        self.previous = 0

    def __enter__(self):
        if self.fixed and self.index is not None:
            self.rows = iter_artifact_lines(self.root, self.index, maximum_line_bytes=MAX_ROW_BYTES)
            if json.loads(next(self.rows)) != self.header:
                raise NativeAnnotationError("native annotation parents changed")
        return self

    def _write(self, row):
        content = canonical(row)
        if len(content) > MAX_ROW_BYTES or self.size + len(content) > self.config.annotation.max_shard_bytes:
            raise NativeAnnotationError("native annotation index capacity exceeded")
        self.stream.write(content)
        self.digest.update(content)
        self.size += len(content)

    def enrich(self, value):
        enriched, intake = self.enrich_intake(value)
        return enriched, intake.proofs()

    def enrich_intake(self, value):
        """Validate the original native source artifacts before constructing proofs."""
        if value.transport != "native_vep" or value.retained_ordinal <= self.previous:
            raise NativeAnnotationError("native annotation input order changed")
        self.previous = value.retained_ordinal
        name = f"native-{value.retained_ordinal:012d}"
        destination = self.root / name
        request = VepRequest(value.allele)
        if self.fixed:
            if self.rows is None:
                raise NativeAnnotationError("saved native annotation coverage is missing")
            try:
                row = json.loads(next(self.rows))
            except StopIteration:
                raise NativeAnnotationError("saved native annotation coverage is incomplete") from None
            if (set(row) != {"ordinal", "identity_id", "directory", "artifacts"}
                    or row["ordinal"] != value.retained_ordinal or row["identity_id"] != value.identity_id
                    or row["directory"] != name):
                raise NativeAnnotationError("saved native annotation identity changed")
            expected = row["artifacts"]
        else:
            if self.stream is None:
                descriptor = os.open(self.root / ".native-index.pending",
                    os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600)
                self.stream = os.fdopen(descriptor, "wb")
                self._write(self.header)
            contract = json.loads(read_bounded_bytes(self.contract_path))
            start_name = name + "-start.json"
            if (self.root / start_name).exists() or (self.root / start_name).is_symlink():
                raise NativeAnnotationError("native annotation attempt is already consumed")
            publish_bytes(self.root, start_name, canonical({**self.header,
                "ordinal": value.retained_ordinal, "identity_id": value.identity_id,
                "admission_input_sha256": hashlib.sha256(value.as_bytes()).hexdigest(),
                "attempt_number": 1, "event": "native-annotation-started"}))
            run_vep_annotation((request,), contract_path=self.contract_path,
                reference_root=self.resources.reference_root, bundle_manifest_sha256=self.bundle_sha256,
                runtime_image=self.runtime_image, upstream_image_digest=contract["tool_image_digest"],
                output_directory=destination, constraint_reference=self.resources.vep_constraint,
                reference_compatibility=self.resources.vep_compatibility,
                timeout=self.config.child_process.wall_timeout_seconds,
                buffer_size=self.config.annotation.buffer_size, execution_config=self.config)
            for filename, source in (("reference.fai", self.resources.reference_root /
                    (contract["fasta"]["output_path"] + ".fai")),
                    ("cache.inventory.json", self.resources.reference_root /
                     contract["cache"]["output_root"] / "cache.inventory.json")):
                content = read_bounded_bytes(source, maximum_bytes=self.config.annotation.max_shard_bytes)
                descriptor = os.open(destination / filename, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
                with os.fdopen(descriptor, "wb") as output:
                    output.write(content)
            expected = None
        if destination.is_symlink() or not destination.is_dir():
            raise NativeAnnotationError("native annotation directory is invalid")
        artifacts, descriptors = {}, []
        remaining = self.config.annotation.max_shard_bytes
        for path in sorted(destination.iterdir()):
            content = read_bounded_bytes(path, maximum_bytes=remaining, reference_root=destination)
            if not content and path.suffix in {".stdout", ".stderr", ".warnings"}:
                continue  # Optional empty logs carry no evidence to replay.
            remaining -= len(content)
            artifacts[path.name] = content
            descriptors.append({"name": path.name, "sha256": hashlib.sha256(content).hexdigest(), "size_bytes": len(content)})
        if expected is not None and descriptors != expected:
            raise NativeAnnotationError("saved native annotation artifacts changed")
        binding = json.loads(artifacts["binding.json"])
        if (binding["bundle_manifest_sha256"] != self.bundle_sha256
                or binding["runtime_image"] != self.runtime_image
                or binding["configuration"].get("bounded_execution_sha256") != hashlib.sha256(
                    bounded_wgs_execution_config_as_json(self.config)).hexdigest()):
            raise NativeAnnotationError("native annotation execution binding changed")
        intake = validate_annotation_admission(artifacts,
            expected_binding_sha256=hashlib.sha256(artifacts["binding.json"]).hexdigest(), requests=(request,))
        allele, = intake.alleles
        if not self.fixed:
            self._write({"ordinal": value.retained_ordinal, "identity_id": value.identity_id,
                         "directory": name, "artifacts": descriptors})
        # Annotation replay authenticates unchanged source/call fields and every
        # assessment, including explicit no-match and unresolved native geometry.
        return replace(value, allele=allele, annotation_status="annotated"), intake

    def __exit__(self, kind, value, traceback):
        try:
            if kind is None:
                if self.rows is not None and next(self.rows, None) is not None:
                    raise NativeAnnotationError("saved native annotation coverage contains extras")
                if self.stream is not None:
                    self.stream.flush()
                    os.fsync(self.stream.fileno())
                    os.link(self.root / ".native-index.pending", self.root / "native-index.jsonl")
                    self.index = Artifact("native-index.jsonl", self.digest.hexdigest(), self.size)
        finally:
            if self.stream is not None:
                self.stream.close()
                (self.root / ".native-index.pending").unlink()
            if self.rows is not None:
                self.rows.close()
