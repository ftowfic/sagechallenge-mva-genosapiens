"""Instrumented chromosome pilot and continuation gate for gnomAD reduction."""

from __future__ import annotations

import base64
from contextlib import contextmanager
from dataclasses import asdict, dataclass
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import threading
import time
from typing import Any, Iterator, Sequence
from urllib.error import URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from genosapiens.gnomad_hail import (
    GnomadHailError,
    HAIL_VERSION,
    PUBLIC_GNOMAD_BUCKET,
    PUBLIC_GNOMAD_SOURCE,
    PYSPARK_VERSION,
    RUN_ID,
    _prepare_hail_scratch_parent,
    _validate_runtime_confinement,
    initialize_local_hail,
    prepare_hail_scratch,
)
from genosapiens.gnomad_projection import (
    GnomadProjectionContract,
    gnomad_vcf_metadata,
    load_gnomad_projection_contract,
    project_hail_joint_table_with_layout,
    project_hail_projection_to_vcf,
    resolve_gnomad_projection_layout,
)
from genosapiens.storage import DEFAULT_RESERVE_BYTES


PILOT_SCHEMA = "genosapiens.gnomad-pilot"
PILOT_COMPLETION_SCHEMA = "genosapiens.gnomad-pilot-completion"
PILOT_CHROMOSOMES = ("chr21", "chr22")
OUTPUT_CAP_BYTES = 150_000_000_000
PILOT_MINIMUM_HEADROOM_BYTES = 100_000_000_000
BCFTOOLS_VERSION = "1.13+htslib-1.13+ds"
SAMPLE_INTERVAL_SECONDS = 0.2
COMPLETION_FILE = "gnomad-pilot-complete.json"
SHA256_DIGEST = "sha256:"
SAMPLED_MEMORY_METHOD = "sampled cgroup current at 200 ms intervals"


@dataclass(frozen=True, slots=True)
class PilotChromosomeMeasurement:
    chromosome: str
    rows: int
    remote_bytes: int
    runtime_ns: int
    bcf_bytes: int
    csi_bytes: int
    output_bytes: int
    bcf_path: str
    csi_path: str
    bcf_sha256: str
    csi_sha256: str

    def __post_init__(self) -> None:
        if self.chromosome not in PILOT_CHROMOSOMES:
            raise GnomadHailError(f"unsupported pilot chromosome: {self.chromosome}")
        for name in (
            "rows",
            "remote_bytes",
            "runtime_ns",
            "bcf_bytes",
            "csi_bytes",
            "output_bytes",
        ):
            value = getattr(self, name)
            if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
                raise GnomadHailError(f"pilot {name} must be a positive integer")
        if self.output_bytes != self.bcf_bytes + self.csi_bytes:
            raise GnomadHailError("pilot output bytes must equal BCF plus CSI bytes")
        for name in ("bcf_sha256", "csi_sha256"):
            value = getattr(self, name)
            if len(value) != 64 or any(item not in "0123456789abcdef" for item in value):
                raise GnomadHailError(f"pilot {name} must be a SHA-256 digest")


@dataclass(frozen=True, slots=True)
class GnomadShardMaterialization:
    """Identity and measurements for one atomically staged BCF shard."""

    chromosome: str
    rows: int
    remote_bytes: int
    runtime_ns: int
    bcf_bytes: int
    csi_bytes: int
    output_bytes: int
    bcf_path: str
    csi_path: str
    bcf_sha256: str
    csi_sha256: str

    def __post_init__(self) -> None:
        if not isinstance(self.chromosome, str) or not self.chromosome:
            raise GnomadHailError("shard chromosome must be a non-empty string")
        for name in (
            "rows",
            "remote_bytes",
            "runtime_ns",
            "bcf_bytes",
            "csi_bytes",
            "output_bytes",
        ):
            value = getattr(self, name)
            if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
                raise GnomadHailError(f"shard {name} must be a positive integer")
        if self.output_bytes != self.bcf_bytes + self.csi_bytes:
            raise GnomadHailError("shard output bytes must equal BCF plus CSI bytes")
        for name in ("bcf_sha256", "csi_sha256"):
            value = getattr(self, name)
            if len(value) != 64 or any(item not in "0123456789abcdef" for item in value):
                raise GnomadHailError(f"shard {name} must be a SHA-256 digest")


@dataclass(frozen=True, slots=True)
class GnomadSourceInventory:
    object_count: int
    total_bytes: int
    sha256: str


@dataclass(frozen=True, slots=True)
class GnomadPilotMeasured:
    rows: int
    remote_bytes: int
    output_bytes: int
    output_bytes_per_variant: float
    runtime_ns: int
    peak_memory_bytes: int
    memory_limit_bytes: int | None
    peak_scratch_bytes: int
    available_before_bytes: int
    available_after_bytes: int
    sample_interval_ms: int


@dataclass(frozen=True, slots=True)
class GnomadPilotExtrapolation:
    expected_full_rows: int
    remaining_rows: int
    projected_output_bytes: int
    projected_remote_bytes: int
    projected_runtime_ns: int
    projected_peak_scratch_bytes: int
    projected_peak_storage_bytes: int
    projected_free_from_pre_pilot_bytes: int
    projected_free_from_post_pilot_bytes: int
    projected_free_bytes: int
    method: str


@dataclass(frozen=True, slots=True)
class GnomadPilotGate:
    output_cap_bytes: int
    reserve_bytes: int
    output_within_cap: bool
    free_space_preserved: bool
    allowed: bool
    additional_bytes_required: int
    reasons: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class GnomadPilotReport:
    schema: str
    mode: str
    input_classification: str
    source: str
    projection_id: str
    projection_contract_sha256: str
    source_root_metadata_url: str
    source_root_metadata_object_generation: str
    source_root_metadata_md5_base64: str
    source_root_metadata_verified: bool
    source_object_inventory_sha256: str
    source_object_count: int
    source_object_bytes: int
    source_object_inventory_verified: bool
    source_written_hail_version: str
    source_partitions: int
    chromosomes: tuple[str, ...]
    measurements: tuple[PilotChromosomeMeasurement, ...]
    measured: GnomadPilotMeasured
    extrapolation: GnomadPilotExtrapolation
    gate: GnomadPilotGate
    hail_version: str
    pyspark_version: str
    bcftools_version: str
    image_digest: str
    workers: int
    spark_master: str
    output_format: str
    remote_bytes_method: str
    memory_method: str
    scratch_method: str
    scratch_root: str
    run_root: str
    report_path: str


def _scaled_ceiling(value: int, target_rows: int, measured_rows: int) -> int:
    return (value * target_rows + measured_rows - 1) // measured_rows


def build_gnomad_pilot_report(
    *,
    measurements: Sequence[PilotChromosomeMeasurement],
    projection_id: str,
    projection_contract_sha256: str,
    source_root_metadata_url: str,
    source_root_metadata_object_generation: str,
    source_root_metadata_md5_base64: str,
    source_root_metadata_verified: bool,
    source_object_inventory_sha256: str,
    source_object_count: int,
    source_object_bytes: int,
    source_object_inventory_verified: bool,
    source_written_hail_version: str,
    image_digest: str,
    workers: int,
    memory_method: str,
    expected_full_rows: int,
    source_partitions: int,
    setup_remote_bytes: int,
    total_runtime_ns: int,
    peak_memory_bytes: int,
    memory_limit_bytes: int | None,
    peak_scratch_bytes: int,
    available_before_bytes: int,
    available_after_bytes: int,
    reserve_bytes: int,
    scratch_root: str | Path,
    run_root: str | Path,
) -> GnomadPilotReport:
    values = tuple(measurements)
    if tuple(item.chromosome for item in values) != PILOT_CHROMOSOMES:
        raise GnomadHailError("pilot measurements must contain chr21 then chr22")
    if (
        not projection_id
        or len(projection_contract_sha256) != 64
        or any(item not in "0123456789abcdef" for item in projection_contract_sha256)
        or not source_root_metadata_url.startswith("https://")
        or not source_root_metadata_object_generation.isdigit()
        or not source_root_metadata_md5_base64
        or source_root_metadata_verified is not True
        or len(source_object_inventory_sha256) != 64
        or any(
            item not in "0123456789abcdef"
            for item in source_object_inventory_sha256
        )
        or source_object_inventory_verified is not True
        or not source_written_hail_version
        or len(image_digest) != len(SHA256_DIGEST) + 64
        or not image_digest.startswith(SHA256_DIGEST)
        or any(item not in "0123456789abcdef" for item in image_digest[7:])
        or not isinstance(workers, int)
        or isinstance(workers, bool)
        or not 1 <= workers <= 16
        or not isinstance(memory_method, str)
        or not memory_method
    ):
        raise GnomadHailError("pilot source identity is invalid")
    for name, value in (
        ("expected_full_rows", expected_full_rows),
        ("source_partitions", source_partitions),
        ("source_object_count", source_object_count),
        ("source_object_bytes", source_object_bytes),
        ("total_runtime_ns", total_runtime_ns),
        ("peak_memory_bytes", peak_memory_bytes),
        ("available_before_bytes", available_before_bytes),
        ("available_after_bytes", available_after_bytes),
    ):
        if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
            raise GnomadHailError(f"{name} must be a positive integer")
    for name, value in (
        ("setup_remote_bytes", setup_remote_bytes),
        ("peak_scratch_bytes", peak_scratch_bytes),
    ):
        if not isinstance(value, int) or isinstance(value, bool) or value < 0:
            raise GnomadHailError(f"{name} must be a non-negative integer")
    if memory_limit_bytes is not None and memory_limit_bytes <= 0:
        raise GnomadHailError("memory_limit_bytes must be positive when present")
    if reserve_bytes < DEFAULT_RESERVE_BYTES:
        raise GnomadHailError("pilot must preserve at least 100 GB")

    pilot_rows = sum(item.rows for item in values)
    if expected_full_rows <= pilot_rows:
        raise GnomadHailError("expected full rows must exceed pilot rows")
    remaining_rows = expected_full_rows - pilot_rows
    output_bytes = sum(item.output_bytes for item in values)
    remote_bytes = setup_remote_bytes + sum(item.remote_bytes for item in values)
    chromosome_runtime_ns = sum(item.runtime_ns for item in values)
    if total_runtime_ns < chromosome_runtime_ns:
        raise GnomadHailError("total runtime cannot be shorter than chromosome runtime")
    projected_output = output_bytes + max(
        _scaled_ceiling(item.output_bytes, remaining_rows, item.rows)
        for item in values
    )
    projected_remote = remote_bytes + max(
        _scaled_ceiling(item.remote_bytes, remaining_rows, item.rows)
        for item in values
    )
    projected_runtime = total_runtime_ns + max(
        _scaled_ceiling(item.runtime_ns, remaining_rows, item.rows)
        for item in values
    )
    projected_peak_scratch = _scaled_ceiling(
        peak_scratch_bytes, expected_full_rows, pilot_rows
    )
    projected_peak_storage = projected_output + projected_peak_scratch
    projected_free_from_pre = available_before_bytes - projected_peak_storage
    remaining_output_bytes = projected_output - output_bytes
    projected_free_from_post = (
        available_after_bytes - remaining_output_bytes - projected_peak_scratch
    )
    projected_free = min(projected_free_from_pre, projected_free_from_post)
    output_within_cap = projected_output <= OUTPUT_CAP_BYTES
    free_space_preserved = projected_free >= reserve_bytes
    reasons = []
    if not output_within_cap:
        reasons.append("projected_output_exceeds_150_gb")
    if not free_space_preserved:
        reasons.append("projected_free_space_below_reserve")

    measured = GnomadPilotMeasured(
        rows=pilot_rows,
        remote_bytes=remote_bytes,
        output_bytes=output_bytes,
        output_bytes_per_variant=round(output_bytes / pilot_rows, 6),
        runtime_ns=total_runtime_ns,
        peak_memory_bytes=peak_memory_bytes,
        memory_limit_bytes=memory_limit_bytes,
        peak_scratch_bytes=peak_scratch_bytes,
        available_before_bytes=available_before_bytes,
        available_after_bytes=available_after_bytes,
        sample_interval_ms=round(SAMPLE_INTERVAL_SECONDS * 1000),
    )
    extrapolation = GnomadPilotExtrapolation(
        expected_full_rows=expected_full_rows,
        remaining_rows=remaining_rows,
        projected_output_bytes=projected_output,
        projected_remote_bytes=projected_remote,
        projected_runtime_ns=projected_runtime,
        projected_peak_scratch_bytes=projected_peak_scratch,
        projected_peak_storage_bytes=projected_peak_storage,
        projected_free_from_pre_pilot_bytes=projected_free_from_pre,
        projected_free_from_post_pilot_bytes=projected_free_from_post,
        projected_free_bytes=projected_free,
        method=(
            "pilot totals plus the maximum chr21/chr22 per-row rate applied "
            "to remaining source rows; scratch scales from the measured pilot "
            "peak to all source rows; projections are rounded upward"
        ),
    )
    gate = GnomadPilotGate(
        output_cap_bytes=OUTPUT_CAP_BYTES,
        reserve_bytes=reserve_bytes,
        output_within_cap=output_within_cap,
        free_space_preserved=free_space_preserved,
        allowed=output_within_cap and free_space_preserved,
        additional_bytes_required=max(0, reserve_bytes - projected_free),
        reasons=tuple(reasons),
    )
    run_path = Path(run_root)
    return GnomadPilotReport(
        schema=PILOT_SCHEMA,
        mode="chromosome-21-22-pilot",
        input_classification="public",
        source=PUBLIC_GNOMAD_SOURCE,
        projection_id=projection_id,
        projection_contract_sha256=projection_contract_sha256,
        source_root_metadata_url=source_root_metadata_url,
        source_root_metadata_object_generation=(
            source_root_metadata_object_generation
        ),
        source_root_metadata_md5_base64=source_root_metadata_md5_base64,
        source_root_metadata_verified=source_root_metadata_verified,
        source_object_inventory_sha256=source_object_inventory_sha256,
        source_object_count=source_object_count,
        source_object_bytes=source_object_bytes,
        source_object_inventory_verified=source_object_inventory_verified,
        source_written_hail_version=source_written_hail_version,
        source_partitions=source_partitions,
        chromosomes=PILOT_CHROMOSOMES,
        measurements=values,
        measured=measured,
        extrapolation=extrapolation,
        gate=gate,
        hail_version=HAIL_VERSION,
        pyspark_version=PYSPARK_VERSION,
        bcftools_version=BCFTOOLS_VERSION,
        image_digest=image_digest,
        workers=workers,
        spark_master=f"local[{workers}]",
        output_format="BCFv2.2 plus CSI, one staged shard per chromosome",
        remote_bytes_method="Hadoop FileSystem gs-scheme bytes-read counter delta",
        memory_method=memory_method,
        scratch_method=(
            "sampled peak allocated st_blocks beneath transient run paths"
        ),
        scratch_root=str(Path(scratch_root)),
        run_root=str(run_path),
        report_path=str(run_path / "gnomad-pilot-report.json"),
    )


def pilot_report_as_json(report: GnomadPilotReport) -> str:
    return json.dumps(asdict(report), indent=2, sort_keys=True) + "\n"


def _sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    try:
        with Path(path).open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
    except OSError as error:
        raise GnomadHailError(f"cannot hash file {path}: {error}") from error
    return digest.hexdigest()


def _fsync_file(path: Path) -> None:
    try:
        with path.open("rb") as handle:
            os.fsync(handle.fileno())
    except OSError as error:
        raise GnomadHailError(f"cannot sync pilot file {path}: {error}") from error


def _fsync_directory(path: Path) -> None:
    try:
        directory_fd = os.open(path, os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    except OSError as error:
        raise GnomadHailError(f"cannot sync pilot directory {path}: {error}") from error


def _verify_source_root_metadata(contract: GnomadProjectionContract) -> None:
    digest = hashlib.md5(usedforsecurity=False)
    request = Request(
        contract.source.root_metadata_url,
        headers={"User-Agent": "genosapiens-reference-builder/0.1.0"},
    )
    try:
        with urlopen(request, timeout=60) as response:
            for chunk in iter(lambda: response.read(1024 * 1024), b""):
                digest.update(chunk)
    except (OSError, URLError) as error:
        raise GnomadHailError(f"cannot verify pinned gnomAD root metadata: {error}") from error
    actual = base64.b64encode(digest.digest()).decode("ascii")
    if actual != contract.source.root_metadata_md5_base64:
        raise GnomadHailError("pinned gnomAD root metadata MD5 does not match")


def _source_inventory(contract: GnomadProjectionContract) -> GnomadSourceInventory:
    bucket_root = f"gs://{PUBLIC_GNOMAD_BUCKET}/"
    if not contract.source.hail_path.startswith(bucket_root):
        raise GnomadHailError("cannot inventory a non-public gnomAD source")
    prefix = contract.source.hail_path.removeprefix(bucket_root).rstrip("/") + "/"

    objects: list[dict[str, str | int | None]] = []
    names: set[str] = set()
    seen_page_tokens: set[str] = set()
    page_token: str | None = None
    while True:
        parameters = {
            "fields": "items(name,generation,size,md5Hash),nextPageToken",
            "maxResults": "1000",
            "prefix": prefix,
        }
        if page_token is not None:
            parameters["pageToken"] = page_token
        url = (
            "https://storage.googleapis.com/storage/v1/b/"
            f"{PUBLIC_GNOMAD_BUCKET}/o?{urlencode(parameters)}"
        )
        request = Request(
            url, headers={"User-Agent": "genosapiens-reference-builder/0.1.0"}
        )
        try:
            with urlopen(request, timeout=60) as response:
                document = json.load(response)
        except (OSError, URLError, ValueError) as error:
            raise GnomadHailError(
                f"cannot inventory pinned gnomAD source: {error}"
            ) from error
        if not isinstance(document, dict) or not isinstance(
            document.get("items"), list
        ):
            raise GnomadHailError("gnomAD source inventory response is invalid")
        for value in document["items"]:
            if not isinstance(value, dict):
                raise GnomadHailError("gnomAD source inventory object is invalid")
            name = value.get("name")
            generation = value.get("generation")
            size = value.get("size")
            md5_hash = value.get("md5Hash")
            if (
                not isinstance(name, str)
                or not name.startswith(prefix)
                or name in names
                or not isinstance(generation, str)
                or not generation.isdigit()
                or not isinstance(size, str)
                or not size.isdigit()
                or (md5_hash is not None and not isinstance(md5_hash, str))
            ):
                raise GnomadHailError("gnomAD source inventory object is invalid")
            names.add(name)
            objects.append(
                {
                    "generation": generation,
                    "md5_base64": md5_hash,
                    "name": name,
                    "size_bytes": int(size),
                }
            )
        next_token = document.get("nextPageToken")
        if next_token is None:
            break
        if not isinstance(next_token, str) or not next_token:
            raise GnomadHailError("gnomAD source inventory page token is invalid")
        if next_token in seen_page_tokens:
            raise GnomadHailError("gnomAD source inventory page token repeated")
        seen_page_tokens.add(next_token)
        page_token = next_token

    objects.sort(key=lambda value: str(value["name"]))
    root_name = f"{prefix}metadata.json.gz"
    root = next((value for value in objects if value["name"] == root_name), None)
    if root is None:
        raise GnomadHailError("gnomAD source inventory lacks root metadata")
    if (
        root["generation"] != contract.source.root_metadata_object_generation
        or root["md5_base64"] != contract.source.root_metadata_md5_base64
    ):
        raise GnomadHailError("live gnomAD root metadata identity does not match")
    if len(objects) <= contract.source.expected_partitions:
        raise GnomadHailError("gnomAD source inventory has too few objects")
    encoded = json.dumps(objects, separators=(",", ":"), sort_keys=True).encode()
    return GnomadSourceInventory(
        object_count=len(objects),
        total_bytes=sum(int(value["size_bytes"]) for value in objects),
        sha256=hashlib.sha256(encoded).hexdigest(),
    )


def _allocated_bytes(roots: Sequence[Path]) -> int:
    total = 0
    for root in roots:
        if not root.exists():
            continue
        for current, directories, files in os.walk(root, followlinks=False):
            current_path = Path(current)
            for path in (current_path, *(current_path / name for name in files)):
                try:
                    total += path.lstat().st_blocks * 512
                except FileNotFoundError:
                    continue
            for name in directories:
                path = current_path / name
                if path.is_symlink():
                    try:
                        total += path.lstat().st_blocks * 512
                    except FileNotFoundError:
                        continue
    return total


def _cgroup_memory(
    cgroup_root: Path = Path("/sys/fs/cgroup"),
) -> tuple[int, int | None, int | None, str]:
    candidates = (
        (
            cgroup_root / "memory.current",
            cgroup_root / "memory.max",
            cgroup_root / "memory.peak",
            "cgroup v2 memory.peak",
        ),
        (
            cgroup_root / "memory/memory.usage_in_bytes",
            cgroup_root / "memory/memory.limit_in_bytes",
            cgroup_root / "memory/memory.max_usage_in_bytes",
            "cgroup v1 memory.max_usage_in_bytes",
        ),
    )
    for current_path, limit_path, peak_path, peak_method in candidates:
        if not current_path.is_file():
            continue
        try:
            current = int(current_path.read_text(encoding="ascii").strip())
            raw_limit = limit_path.read_text(encoding="ascii").strip()
            limit = None if raw_limit == "max" else int(raw_limit)
            peak = (
                int(peak_path.read_text(encoding="ascii").strip())
                if peak_path.is_file()
                else None
            )
        except (OSError, ValueError) as error:
            raise GnomadHailError(f"cannot read container memory cgroup: {error}") from error
        if (
            current < 0
            or (limit is not None and limit <= 0)
            or (peak is not None and peak < 0)
        ):
            raise GnomadHailError("container memory cgroup values are invalid")
        method = peak_method if peak is not None else SAMPLED_MEMORY_METHOD
        return current, limit, peak, method
    raise GnomadHailError("container memory cgroup accounting is unavailable")


class _ResourceSampler:
    def __init__(self, scratch_roots: Sequence[Path]) -> None:
        self.scratch_roots = tuple(scratch_roots)
        self.peak_memory_bytes = 0
        self.memory_limit_bytes: int | None = None
        self.memory_method = SAMPLED_MEMORY_METHOD
        self.peak_scratch_bytes = 0
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None
        self._error: BaseException | None = None

    def sample(self) -> None:
        memory, limit, kernel_peak, method = _cgroup_memory()
        self.peak_memory_bytes = max(
            self.peak_memory_bytes,
            kernel_peak if kernel_peak is not None else memory,
        )
        self.memory_limit_bytes = limit
        self.memory_method = method
        self.peak_scratch_bytes = max(
            self.peak_scratch_bytes, _allocated_bytes(self.scratch_roots)
        )

    def _run(self) -> None:
        while not self._stop.wait(SAMPLE_INTERVAL_SECONDS):
            try:
                self.sample()
            except BaseException as error:
                self._error = error
                self._stop.set()

    def start(self) -> None:
        self.sample()
        self._thread = threading.Thread(
            target=self._run,
            name="gnomad-pilot-resource-sampler",
            daemon=True,
        )
        self._thread.start()

    def close(self) -> None:
        self._stop.set()
        if self._thread is not None:
            self._thread.join()
        if self._error is not None:
            raise GnomadHailError(f"resource sampling failed: {self._error}")
        self.sample()


def _gcs_bytes_read(hail_module: Any) -> int:
    try:
        statistics = (
            hail_module.spark_context()
            ._jvm.org.apache.hadoop.fs.FileSystem.getAllStatistics()
        )
        iterator = statistics.iterator()
        total = 0
        while iterator.hasNext():
            item = iterator.next()
            if item.getScheme() == "gs":
                total += int(item.getBytesRead())
        return total
    except Exception as error:
        raise GnomadHailError(f"cannot read Hadoop GCS byte counters: {error}") from error


def _run_bcftools(arguments: Sequence[str]) -> str:
    command = ["bcftools", *arguments]
    try:
        result = subprocess.run(
            command,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except OSError as error:
        raise GnomadHailError(f"cannot execute bcftools: {error}") from error
    if result.returncode != 0:
        detail = result.stderr.strip().splitlines()
        message = detail[-1] if detail else f"exit status {result.returncode}"
        raise GnomadHailError(f"bcftools failed: {message}")
    return result.stdout.strip()


def _require_bcftools() -> None:
    found = _run_bcftools(("--version-only",))
    if found != BCFTOOLS_VERSION:
        raise GnomadHailError(
            f"bcftools runtime is {found!r}; expected {BCFTOOLS_VERSION!r}"
        )


def materialize_gnomad_shard(
    *,
    hail_module: Any,
    table: Any,
    contract: GnomadProjectionContract,
    projection_layout: Any,
    chromosome: str,
    export_root: Path,
    output_root: Path,
    sampler: _ResourceSampler,
) -> GnomadShardMaterialization:
    started = time.perf_counter_ns()
    remote_before = _gcs_bytes_read(hail_module)
    interval = hail_module.parse_locus_interval(
        chromosome, reference_genome=contract.source.assembly
    )
    chromosome_table = hail_module.filter_intervals(table, [interval])
    projected = project_hail_joint_table_with_layout(
        chromosome_table, hail_module, projection_layout
    )
    vcf_table = project_hail_projection_to_vcf(projected, hail_module, contract)
    vcf_path = export_root / f"gnomad.joint.v4.1.{chromosome}.vcf.bgz"
    bcf_path = output_root / f"gnomad.joint.v4.1.{chromosome}.bcf"
    csi_path = Path(f"{bcf_path}.csi")
    # Keep the temporary filename BCF-suffixed so bcftools discovers its CSI.
    bcf_partial = output_root / f".{bcf_path.stem}.partial.bcf"
    csi_partial = Path(f"{bcf_partial}.csi")
    hail_module.export_vcf(
        vcf_table,
        str(vcf_path),
        metadata=gnomad_vcf_metadata(contract),
        tabix=False,
    )
    remote_after = _gcs_bytes_read(hail_module)
    sampler.sample()
    if remote_after <= remote_before:
        raise GnomadHailError(f"{chromosome} Hail export recorded no remote bytes")

    _run_bcftools(
        (
            "view",
            "--no-version",
            "--output-type",
            "b",
            "--output-file",
            str(bcf_partial),
            str(vcf_path),
        )
    )
    _run_bcftools(("index", "--csi", "--force", str(bcf_partial)))
    rows_text = _run_bcftools(("index", "--nrecords", str(bcf_partial)))
    try:
        rows = int(rows_text)
    except ValueError as error:
        raise GnomadHailError(
            f"bcftools returned an invalid row count for {chromosome}: {rows_text!r}"
        ) from error
    _run_bcftools(("view", "--header-only", str(bcf_partial)))
    _fsync_file(bcf_partial)
    _fsync_file(csi_partial)
    sampler.sample()
    vcf_path.unlink()
    (export_root / f".{vcf_path.name}.crc").unlink(missing_ok=True)
    os.replace(bcf_partial, bcf_path)
    os.replace(csi_partial, csi_path)
    _fsync_directory(output_root)
    sampler.sample()
    bcf_bytes = bcf_path.stat().st_size
    csi_bytes = csi_path.stat().st_size
    bcf_sha256 = _sha256_file(bcf_path)
    csi_sha256 = _sha256_file(csi_path)
    runtime_ns = time.perf_counter_ns() - started
    return GnomadShardMaterialization(
        chromosome=chromosome,
        rows=rows,
        remote_bytes=remote_after - remote_before,
        runtime_ns=runtime_ns,
        bcf_bytes=bcf_bytes,
        csi_bytes=csi_bytes,
        output_bytes=bcf_bytes + csi_bytes,
        bcf_path=str(bcf_path),
        csi_path=str(csi_path),
        bcf_sha256=bcf_sha256,
        csi_sha256=csi_sha256,
    )


def _materialize_chromosome(
    *,
    hail_module: Any,
    table: Any,
    contract: GnomadProjectionContract,
    projection_layout: Any,
    chromosome: str,
    export_root: Path,
    output_root: Path,
    sampler: _ResourceSampler,
) -> PilotChromosomeMeasurement:
    """Materialize a fixed pilot chromosome with the shared shard writer."""
    shard = materialize_gnomad_shard(
        hail_module=hail_module,
        table=table,
        contract=contract,
        projection_layout=projection_layout,
        chromosome=chromosome,
        export_root=export_root,
        output_root=output_root,
        sampler=sampler,
    )
    return PilotChromosomeMeasurement(
        chromosome=shard.chromosome,
        rows=shard.rows,
        remote_bytes=shard.remote_bytes,
        runtime_ns=shard.runtime_ns,
        bcf_bytes=shard.bcf_bytes,
        csi_bytes=shard.csi_bytes,
        output_bytes=shard.output_bytes,
        bcf_path=shard.bcf_path,
        csi_path=shard.csi_path,
        bcf_sha256=shard.bcf_sha256,
        csi_sha256=shard.csi_sha256,
    )


def _write_report(report: GnomadPilotReport) -> None:
    path = Path(report.report_path)
    temporary = path.with_suffix(".json.partial")
    completion = path.parent / COMPLETION_FILE
    completion_temporary = completion.with_suffix(".json.partial")
    try:
        with temporary.open("x", encoding="utf-8") as handle:
            handle.write(pilot_report_as_json(report))
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
        _fsync_directory(path.parent)
        completion_document = {
            "report_path": path.name,
            "report_sha256": _sha256_file(path),
            "schema": PILOT_COMPLETION_SCHEMA,
        }
        with completion_temporary.open("x", encoding="utf-8") as handle:
            json.dump(completion_document, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(completion_temporary, completion)
        _fsync_directory(path.parent)
    except OSError as error:
        raise GnomadHailError(f"cannot publish pilot report {path}: {error}") from error


def _remove_failed_run(run_root: Path, scratch_root: str | Path) -> None:
    scratch = Path(scratch_root).resolve(strict=False)
    if not run_root.parent.resolve(strict=False).is_relative_to(scratch):
        raise GnomadHailError("refusing to clean a pilot run outside public scratch")
    try:
        if run_root.is_symlink() or (run_root.exists() and not run_root.is_dir()):
            run_root.unlink()
        elif run_root.exists():
            if not run_root.resolve(strict=True).is_relative_to(scratch):
                raise GnomadHailError(
                    "refusing to clean a pilot run outside public scratch"
                )
            shutil.rmtree(run_root)
    except OSError as error:
        raise GnomadHailError(f"cannot clean failed pilot run {run_root}: {error}") from error


@contextmanager
def _pilot_run_lock(
    reference_root: str | Path,
    scratch_root: str | Path,
    run_id: str,
) -> Iterator[None]:
    _, parent_descriptor = _prepare_hail_scratch_parent(
        reference_root, scratch_root
    )
    try:
        lock_descriptor = os.open(
            f".{run_id}.lock",
            os.O_RDWR | os.O_CREAT | os.O_NOFOLLOW,
            0o600,
            dir_fd=parent_descriptor,
        )
    except OSError as error:
        raise GnomadHailError(f"cannot claim gnomAD pilot run: {error}") from error
    finally:
        os.close(parent_descriptor)

    try:
        try:
            fcntl.flock(lock_descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise GnomadHailError(f"pilot run is already active: {run_id}") from error
        yield
    finally:
        os.close(lock_descriptor)


def run_gnomad_pilot(
    *,
    contract_path: str | Path,
    reference_root: str | Path,
    scratch_root: str | Path,
    run_id: str,
    workers: int,
    reserve_bytes: int,
    image_digest: str,
) -> GnomadPilotReport:
    contract = load_gnomad_projection_contract(contract_path)
    if contract.source.hail_path != PUBLIC_GNOMAD_SOURCE:
        raise GnomadHailError("gnomAD Hail source must match the pinned public table")
    if not isinstance(reserve_bytes, int) or isinstance(reserve_bytes, bool):
        raise GnomadHailError("reserve_bytes must be an integer")
    if reserve_bytes < DEFAULT_RESERVE_BYTES:
        raise GnomadHailError("pilot must preserve at least 100 GB")
    if not isinstance(run_id, str) or not RUN_ID.fullmatch(run_id):
        raise GnomadHailError("run_id must be a lowercase path-safe identifier")
    if (
        not isinstance(image_digest, str)
        or len(image_digest) != len(SHA256_DIGEST) + 64
        or not image_digest.startswith(SHA256_DIGEST)
        or any(item not in "0123456789abcdef" for item in image_digest[7:])
    ):
        raise GnomadHailError("image_digest must be a SHA-256 digest")
    _require_bcftools()
    available_before = shutil.disk_usage(reference_root).free
    required_before = reserve_bytes + PILOT_MINIMUM_HEADROOM_BYTES
    if available_before < required_before:
        raise GnomadHailError(
            "gnomAD pilot requires "
            f"{required_before - available_before} additional bytes before starting"
        )
    with _pilot_run_lock(reference_root, scratch_root, run_id):
        total_started = time.perf_counter_ns()
        projection_contract_sha256 = _sha256_file(contract_path)
        _verify_source_root_metadata(contract)
        source_inventory_before = _source_inventory(contract)

        candidate_run_root = Path(scratch_root) / "gnomad-hail" / run_id
        if candidate_run_root.is_symlink():
            _remove_failed_run(candidate_run_root, scratch_root)
        elif candidate_run_root.exists():
            if (candidate_run_root / COMPLETION_FILE).is_file():
                raise GnomadHailError(
                    f"pilot run path already exists: {candidate_run_root}"
                )
            _remove_failed_run(candidate_run_root, scratch_root)
        layout = prepare_hail_scratch(reference_root, scratch_root, run_id)
        export_root = layout.run_root / "export"
        output_root = layout.run_root / "pilot-output"
        try:
            export_root.mkdir()
            output_root.mkdir()
        except BaseException:
            _remove_failed_run(layout.run_root, scratch_root)
            raise

        transient_roots = (
            layout.hail_tmp,
            layout.spark_local,
            layout.spark_warehouse,
            layout.checkpoints,
            layout.logs,
            layout.home,
            layout.cache,
            layout.ivy,
            layout.java_tmp,
            layout.python_tmp,
            export_root,
        )
        sampler = _ResourceSampler(transient_roots)
        hail_module: Any | None = None
        measurements: list[PilotChromosomeMeasurement] = []
        source_partitions = 0
        setup_remote_bytes = 0
        try:
            sampler.start()
            try:
                try:
                    hail_module = initialize_local_hail(layout, workers)
                    remote_before_setup = _gcs_bytes_read(hail_module)
                    table = hail_module.read_table(contract.source.hail_path)
                    source_partitions = table.n_partitions()
                    if source_partitions != contract.source.expected_partitions:
                        raise GnomadHailError(
                            "gnomAD source partition count differs from the pinned "
                            "contract"
                        )
                    source_globals = hail_module.eval(table.index_globals())
                    projection_layout = resolve_gnomad_projection_layout(
                        contract, source_globals
                    )
                    remote_after_setup = _gcs_bytes_read(hail_module)
                    if remote_after_setup < remote_before_setup:
                        raise GnomadHailError(
                            "Hadoop GCS byte counter decreased during setup"
                        )
                    setup_remote_bytes = remote_after_setup - remote_before_setup
                    for chromosome in PILOT_CHROMOSOMES:
                        measurements.append(
                            _materialize_chromosome(
                                hail_module=hail_module,
                                table=table,
                                contract=contract,
                                projection_layout=projection_layout,
                                chromosome=chromosome,
                                export_root=export_root,
                                output_root=output_root,
                                sampler=sampler,
                            )
                        )
                    _validate_runtime_confinement(hail_module, layout)
                finally:
                    if hail_module is not None:
                        hail_module.stop()
                available_after = shutil.disk_usage(reference_root).free
                source_inventory_after = _source_inventory(contract)
            finally:
                sampler.close()
            if source_inventory_after != source_inventory_before:
                raise GnomadHailError(
                    "gnomAD source object inventory changed during pilot"
                )
            total_runtime_ns = time.perf_counter_ns() - total_started
            report = build_gnomad_pilot_report(
                measurements=measurements,
                projection_id=contract.projection_id,
                projection_contract_sha256=projection_contract_sha256,
                source_root_metadata_url=contract.source.root_metadata_url,
                source_root_metadata_object_generation=(
                    contract.source.root_metadata_object_generation
                ),
                source_root_metadata_md5_base64=(
                    contract.source.root_metadata_md5_base64
                ),
                source_root_metadata_verified=True,
                source_object_inventory_sha256=source_inventory_after.sha256,
                source_object_count=source_inventory_after.object_count,
                source_object_bytes=source_inventory_after.total_bytes,
                source_object_inventory_verified=True,
                source_written_hail_version=contract.source.written_hail_version,
                image_digest=image_digest,
                workers=workers,
                memory_method=sampler.memory_method,
                expected_full_rows=contract.source.expected_rows,
                source_partitions=source_partitions,
                setup_remote_bytes=setup_remote_bytes,
                total_runtime_ns=total_runtime_ns,
                peak_memory_bytes=sampler.peak_memory_bytes,
                memory_limit_bytes=sampler.memory_limit_bytes,
                peak_scratch_bytes=sampler.peak_scratch_bytes,
                available_before_bytes=available_before,
                available_after_bytes=available_after,
                reserve_bytes=reserve_bytes,
                scratch_root=scratch_root,
                run_root=layout.run_root,
            )
            _write_report(report)
        except BaseException:
            _remove_failed_run(layout.run_root, scratch_root)
            raise
        return report
