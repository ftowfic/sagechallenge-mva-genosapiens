"""Confined local Hail/Spark execution for public gnomAD reduction."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from importlib import metadata
import json
import os
from pathlib import Path
import re
import tempfile
from typing import Any, Mapping

from genosapiens.gnomad_projection import (
    GnomadProjectionContract,
    load_gnomad_projection_contract,
    project_hail_joint_table,
)


HAIL_VERSION = "0.2.139"
PYSPARK_VERSION = "3.5.7"
GCS_CONNECTOR_VERSION = "hadoop3-2.2.24"
PUBLIC_GNOMAD_BUCKET = "gcp-public-data--gnomad"
PUBLIC_GNOMAD_SOURCE = (
    f"gs://{PUBLIC_GNOMAD_BUCKET}/release/4.1/ht/joint/"
    "gnomad.joint.v4.1.sites.ht"
)
BOUNDED_PUBLIC_INTERVAL = "chr22:16050000-16050100"
BOUNDED_PUBLIC_MAX_ROWS = 10
PUBLIC_SCRATCH_RELATIVE = Path("scratch/public-build")
RUN_ID = re.compile(r"^[a-z0-9][a-z0-9._-]*$")


class GnomadHailError(RuntimeError):
    """Raised when local Hail execution violates its public-data contract."""


@dataclass(frozen=True, slots=True)
class HailScratchLayout:
    scratch_root: Path
    run_root: Path
    hail_tmp: Path
    spark_local: Path
    spark_warehouse: Path
    checkpoints: Path
    logs: Path
    home: Path
    cache: Path
    ivy: Path
    java_tmp: Path
    python_tmp: Path


@dataclass(frozen=True, slots=True)
class HailScanReport:
    schema: str
    mode: str
    input_classification: str
    source: str | None
    interval: str | None
    row_limit: int
    rows_scanned: int
    source_partitions: int | None
    hail_version: str
    pyspark_version: str
    spark_master: str
    scratch_root: str
    temporary_paths: tuple[str, ...]


def _resolved(path: Path) -> Path:
    try:
        return path.resolve(strict=False)
    except OSError as error:
        raise GnomadHailError(f"cannot resolve scratch path {path}: {error}") from error


def _require_confined(path: Path, root: Path) -> None:
    if not _resolved(path).is_relative_to(_resolved(root)):
        raise GnomadHailError(f"temporary path escapes public scratch: {path}")


def _open_directory_chain(path: Path) -> int:
    flags = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW
    descriptor = os.open(path.anchor, flags)
    try:
        for name in path.parts[1:]:
            try:
                os.mkdir(name, dir_fd=descriptor)
            except FileExistsError:
                pass
            child = os.open(name, flags, dir_fd=descriptor)
            os.close(descriptor)
            descriptor = child
    except BaseException:
        os.close(descriptor)
        raise
    return descriptor


def _prepare_hail_scratch_parent(
    reference_root: str | Path,
    scratch_root: str | Path,
) -> tuple[Path, int]:
    reference = Path(reference_root)
    scratch = Path(scratch_root)
    if not reference.is_absolute() or not scratch.is_absolute():
        raise GnomadHailError("reference and scratch roots must be absolute")
    if ".." in reference.parts or ".." in scratch.parts:
        raise GnomadHailError("reference and scratch roots cannot contain '..'")
    expected = reference / PUBLIC_SCRATCH_RELATIVE
    if scratch != expected:
        raise GnomadHailError(
            f"Hail scratch_root must be the public build path {expected}"
        )

    try:
        reference_descriptor = _open_directory_chain(reference)
        os.close(reference_descriptor)
        scratch_descriptor = _open_directory_chain(scratch)
        try:
            try:
                os.mkdir("gnomad-hail", dir_fd=scratch_descriptor)
            except FileExistsError:
                pass
            parent_descriptor = os.open(
                "gnomad-hail",
                os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
                dir_fd=scratch_descriptor,
            )
        finally:
            os.close(scratch_descriptor)
    except OSError as error:
        raise GnomadHailError(f"cannot initialize Hail scratch: {error}") from error
    return scratch / "gnomad-hail", parent_descriptor


def prepare_hail_scratch(
    reference_root: str | Path,
    scratch_root: str | Path,
    run_id: str,
) -> HailScratchLayout:
    if not RUN_ID.fullmatch(run_id):
        raise GnomadHailError("run_id must be a lowercase path-safe identifier")
    parent, parent_descriptor = _prepare_hail_scratch_parent(
        reference_root, scratch_root
    )
    scratch = Path(scratch_root)
    run_root = parent / run_id
    values = {
        "scratch_root": scratch,
        "run_root": run_root,
        "hail_tmp": run_root / "hail-tmp",
        "spark_local": run_root / "spark-local",
        "spark_warehouse": run_root / "spark-warehouse",
        "checkpoints": run_root / "checkpoints",
        "logs": run_root / "logs",
        "home": run_root / "home",
        "cache": run_root / "cache",
        "ivy": run_root / "ivy",
        "java_tmp": run_root / "java-tmp",
        "python_tmp": run_root / "python-tmp",
    }
    try:
        os.mkdir(run_id, dir_fd=parent_descriptor)
        run_descriptor = os.open(
            run_id,
            os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
            dir_fd=parent_descriptor,
        )
        try:
            for name in (
                "hail-tmp",
                "spark-local",
                "spark-warehouse",
                "checkpoints",
                "logs",
                "home",
                "cache",
                "ivy",
                "java-tmp",
                "python-tmp",
            ):
                os.mkdir(name, dir_fd=run_descriptor)
        finally:
            os.close(run_descriptor)
        _require_confined(run_root, scratch)
    except FileExistsError as error:
        raise GnomadHailError(f"Hail run path already exists: {run_root}") from error
    except OSError as error:
        raise GnomadHailError(f"cannot initialize Hail scratch: {error}") from error
    finally:
        os.close(parent_descriptor)
    return HailScratchLayout(**values)


def hail_environment(layout: HailScratchLayout) -> dict[str, str]:
    java_options = (
        f"-Djava.io.tmpdir={layout.java_tmp} "
        f"-Duser.home={layout.home} "
        "-Duser.name=genosapiens"
    )
    return {
        "HAIL_QUERY_BACKEND": "spark",
        "HADOOP_USER_NAME": "genosapiens",
        "HOME": str(layout.home),
        "IVY_HOME": str(layout.ivy),
        "LOGNAME": "genosapiens",
        "MPLCONFIGDIR": str(layout.cache / "matplotlib"),
        "SPARK_LOCAL_HOSTNAME": "localhost",
        "SPARK_LOCAL_IP": "127.0.0.1",
        "SPARK_LOCAL_DIRS": str(layout.spark_local),
        "TEMP": str(layout.python_tmp),
        "TMP": str(layout.python_tmp),
        "TMPDIR": str(layout.python_tmp),
        "USER": "genosapiens",
        "XDG_CACHE_HOME": str(layout.cache),
        "_JAVA_OPTIONS": java_options,
    }


def hail_spark_conf(layout: HailScratchLayout) -> dict[str, str]:
    java_options = (
        f"-Djava.io.tmpdir={layout.java_tmp} "
        f"-Duser.home={layout.home} "
        "-Duser.name=genosapiens"
    )
    return {
        "spark.driver.bindAddress": "127.0.0.1",
        "spark.driver.extraJavaOptions": java_options,
        "spark.driver.host": "127.0.0.1",
        "spark.eventLog.enabled": "false",
        "spark.executor.extraJavaOptions": java_options,
        "spark.hadoop.fs.AbstractFileSystem.gs.impl": (
            "com.google.cloud.hadoop.fs.gcs.GoogleHadoopFS"
        ),
        "spark.hadoop.fs.gs.auth.null.enable": "true",
        "spark.hadoop.fs.gs.auth.service.account.enable": "false",
        "spark.hadoop.fs.gs.impl": (
            "com.google.cloud.hadoop.fs.gcs.GoogleHadoopFileSystem"
        ),
        "spark.jars.ivy": str(layout.ivy),
        "spark.local.dir": str(layout.spark_local),
        "spark.sql.streaming.checkpointLocation": str(layout.checkpoints),
        "spark.sql.warehouse.dir": str(layout.spark_warehouse),
        "spark.ui.enabled": "false",
        "spark.ui.showConsoleProgress": "false",
    }


def hail_init_kwargs(layout: HailScratchLayout, workers: int) -> dict[str, Any]:
    if not isinstance(workers, int) or isinstance(workers, bool) or not 1 <= workers <= 16:
        raise GnomadHailError("workers must be an integer between 1 and 16")
    return {
        "app_name": "genosapiens-public-gnomad",
        "backend": "spark",
        "copy_log_on_error": False,
        "local_tmpdir": f"file://{layout.spark_local}",
        "log": str(layout.logs / "hail.log"),
        "master": f"local[{workers}]",
        "quiet": True,
        "show_progress": False,
        "spark_conf": hail_spark_conf(layout),
        "tmp_dir": f"file://{layout.hail_tmp}",
    }


def _package_version(name: str) -> str:
    try:
        return metadata.version(name)
    except metadata.PackageNotFoundError as error:
        raise GnomadHailError(f"required runtime package is missing: {name}") from error


def _load_hail() -> Any:
    versions = {"hail": HAIL_VERSION, "pyspark": PYSPARK_VERSION}
    for name, expected in versions.items():
        actual = _package_version(name)
        if actual != expected:
            raise GnomadHailError(
                f"{name} runtime must be {expected}, found {actual}"
            )
    try:
        import hail as hail_module
    except ImportError as error:
        raise GnomadHailError(f"cannot import Hail runtime: {error}") from error
    return hail_module


def initialize_local_hail(layout: HailScratchLayout, workers: int) -> Any:
    for name, value in hail_environment(layout).items():
        os.environ[name] = value
    tempfile.tempdir = str(layout.python_tmp)
    hail_module = _load_hail()
    try:
        hail_module.init(**hail_init_kwargs(layout, workers))
    except Exception as error:
        raise GnomadHailError(f"cannot initialize confined Hail runtime: {error}") from error
    return hail_module


def _temporary_paths(layout: HailScratchLayout) -> tuple[str, ...]:
    return tuple(
        str(path)
        for path in (
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
        )
    )


def _validate_runtime_confinement(
    hail_module: Any,
    layout: HailScratchLayout,
) -> None:
    for path in _temporary_paths(layout):
        _require_confined(Path(path), layout.scratch_root)
    conf = dict(hail_module.spark_context().getConf().getAll())
    expected = hail_spark_conf(layout)
    for name, value in expected.items():
        actual = conf.get(name)
        valid = value in actual if name.endswith("extraJavaOptions") and actual else actual == value
        if not valid:
            raise GnomadHailError(
                f"active Spark configuration {name} is not confined to public scratch"
            )
    if hail_module.tmp_dir() != f"file://{layout.hail_tmp}":
        raise GnomadHailError("active Hail temporary directory is not confined")


def _report(
    *,
    mode: str,
    source: str | None,
    interval: str | None,
    row_limit: int,
    rows_scanned: int,
    source_partitions: int | None,
    workers: int,
    layout: HailScratchLayout,
) -> HailScanReport:
    return HailScanReport(
        schema="genosapiens.gnomad-hail-scan",
        mode=mode,
        input_classification="public",
        source=source,
        interval=interval,
        row_limit=row_limit,
        rows_scanned=rows_scanned,
        source_partitions=source_partitions,
        hail_version=HAIL_VERSION,
        pyspark_version=PYSPARK_VERSION,
        spark_master=f"local[{workers}]",
        scratch_root=str(layout.scratch_root),
        temporary_paths=_temporary_paths(layout),
    )


def run_hail_smoke(
    *,
    reference_root: str | Path,
    scratch_root: str | Path,
    run_id: str,
    workers: int,
    rows: int = 32,
) -> HailScanReport:
    if not isinstance(rows, int) or isinstance(rows, bool) or not 1 <= rows <= 1_000:
        raise GnomadHailError("rows must be an integer between 1 and 1000")
    layout = prepare_hail_scratch(reference_root, scratch_root, run_id)
    hail_module = initialize_local_hail(layout, workers)
    try:
        table = hail_module.utils.range_table(rows, n_partitions=min(workers, rows))
        scanned = table.count()
        _validate_runtime_confinement(hail_module, layout)
        return _report(
            mode="synthetic-runtime-smoke",
            source=None,
            interval=None,
            row_limit=rows,
            rows_scanned=scanned,
            source_partitions=None,
            workers=workers,
            layout=layout,
        )
    finally:
        hail_module.stop()


def run_bounded_public_scan(
    *,
    contract_path: str | Path,
    reference_root: str | Path,
    scratch_root: str | Path,
    run_id: str,
    workers: int,
) -> HailScanReport:
    contract: GnomadProjectionContract = load_gnomad_projection_contract(contract_path)
    if contract.source.hail_path != PUBLIC_GNOMAD_SOURCE:
        raise GnomadHailError("gnomAD Hail source must match the pinned public table")

    layout = prepare_hail_scratch(reference_root, scratch_root, run_id)
    hail_module = initialize_local_hail(layout, workers)
    try:
        table = hail_module.read_table(contract.source.hail_path)
        partitions = table.n_partitions()
        if partitions != contract.source.expected_partitions:
            raise GnomadHailError(
                "gnomAD source partition count differs from the pinned contract"
            )
        locus_interval = hail_module.parse_locus_interval(
            BOUNDED_PUBLIC_INTERVAL, reference_genome=contract.source.assembly
        )
        bounded = hail_module.filter_intervals(table, [locus_interval]).head(
            BOUNDED_PUBLIC_MAX_ROWS
        )
        projected = project_hail_joint_table(bounded, hail_module, contract)
        scanned = projected.count()
        _validate_runtime_confinement(hail_module, layout)
        return _report(
            mode="bounded-public-source",
            source=contract.source.hail_path,
            interval=BOUNDED_PUBLIC_INTERVAL,
            row_limit=BOUNDED_PUBLIC_MAX_ROWS,
            rows_scanned=scanned,
            source_partitions=partitions,
            workers=workers,
            layout=layout,
        )
    finally:
        hail_module.stop()


def hail_report_as_json(report: HailScanReport) -> str:
    return json.dumps(asdict(report), indent=2, sort_keys=True) + "\n"
