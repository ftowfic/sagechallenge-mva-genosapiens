"""Docker-isolated public case execution; evaluator truth never enters a worker."""
from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import hashlib
import http.client
import json
import os
from pathlib import Path
import socket
import sys
import time
from uuid import uuid4

from genosapiens.track1_public_validation import (
    RUNTIME_IMAGE, materialize_runner_view, load_public_policy_family,
    load_public_resource_context, load_public_model_loci,
)
from genosapiens.track1_validation_plan import InputKind, RankerCaseView, RankerInput, RankerResource

# Only public resource products used by this path are exposed, never the whole
# reference root (which can also contain controlled work or Track 2 packages).
PUBLIC_REFERENCE_TREES = (
    "releases/model-loci/grc-ensembl-nhgri/grch38-ensembl-116-nhgri-retrieved-2026-09-16",
    "derived/model-loci/22b69fc15bbbb97dfeac76209e57be0962270ca6a1eaa8011151c9daa2b9fcaf",
    "manifests/model-loci",
    "releases/clinvar/2026-08-22",
    "derived/gnomad/4.1/joint-frequency/gnomad-joint-v4.1-frequency-projection-v1",
    "derived/gnomad/4.1/all-sites/gnomad-v4.1-all-sites-an-v1",
    "manifests/clinvar",
    "manifests/gnomad/4.1",
    "releases/grch38/GCA_000001405.15-masked-v2-plus-hs38d1-current",
    "releases/ensembl/vep/116_GRCh38",
    "releases/open-targets/26.06",
    "releases/orphadata/2026-06-23",
    "releases/gencc/2026-08-30",
    "releases/ontologies/core-ontology-references-v1",
    "releases/repeats/broad-str-analysis/d41c85c97da3cd8c8ac9d7d47883570c95069874",
    "releases/challenge-evaluator/SageBio-rare-disease-real-kid-mva-hackathon-2026/1c761cc23d90aebe6a011fd5b0b99517df42408c",
    "derived/phenotype/hpo-ic/v1/060cf7c49671f9a98f9a8578c35f0fb8eeae62166d79cb40c2bcc81855abcc89",
    "derived/rare-disease-evidence/orphadata-2026-06-23_gencc-2026-08-30_ot-26.06_v6",
    "manifests/grch38",
    "manifests/ensembl/vep",
    "manifests/open-targets",
    "manifests/rare-disease-evidence",
    "manifests/ontologies",
    "manifests/phenotype",
    "manifests/variant-class-resources",
    "manifests/common-allele-public-indexes",
    "derived/common-allele-public-indexes",
    "manifests/track1-common-allele-exceptions",
    "derived/track1-common-allele-exceptions",
)
# The public replay has finite authenticated shards but transient typed rows
# reached 10 GiB across repeated sensitivity verification. Allow 20% headroom
# above the 12 GiB ceiling, with MemorySwap equal to Memory below.
WORKER_MEMORY_BYTES = 12 * 1024**3 * 6 // 5
# The repaired installed-source worker body took 3,373 seconds for both
# development policies and their mandatory exact replays. Keep a finite
# 108-minute limit with headroom for nonempty results and concurrent I/O;
# see docs/track1-execution-investigation-2026-09-16.md. This is an
# operational limit, independent of the scientific experiment budgets.
WORKER_TIMEOUT_SECONDS = 6480
WORKER_STATUS_INTERVAL_SECONDS = 30


def verify_scientific_stack_projection(root, content):
    """Verify every mounted scientific source/config against the full stack.

    Hash metadata for excluded evaluator-side files is safe to cross the worker
    boundary. Their contents, including acceptance truth, remain outside it.
    """
    try:
        if type(content) is not bytes or len(content) > 8 * 1024 * 1024:
            raise ValueError
        value = json.loads(content)
        canonical = (json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True) + "\n").encode("ascii")
        if (canonical != content or set(value) != {"schema", "runtime_image", "files"}
                or value["schema"] != "genosapiens.track1-executable-stack"
                or value["runtime_image"] != RUNTIME_IMAGE):
            raise ValueError
        rows = {}
        for row in value["files"]:
            path = Path(row["path"])
            if (set(row) != {"path", "size_bytes", "sha256"} or not path.parts
                    or path.is_absolute() or path.as_posix() != row["path"] or ".." in path.parts
                    or row["path"] in rows or type(row["size_bytes"]) is not int or row["size_bytes"] < 0
                    or type(row["sha256"]) is not str or len(row["sha256"]) != 64
                    or any(c not in "0123456789abcdef" for c in row["sha256"])):
                raise ValueError
            rows[row["path"]] = row
        if tuple(rows) != tuple(sorted(rows)):
            raise ValueError
        root = Path(root)
        paths = (*tuple((root / "src/genosapiens").glob("*.py")),
            *tuple((root / "config/public-resources").glob("*.json")),
            *tuple((root / "config/reference-identities").glob("*.json")),
            *(root / "config/track1" / name for name in (
                "public-ranking-development-policy-family.json", "public-population-development-derivation.json",
                "challenge-evaluator.json", "epcr-policy.json", "epcr-numerical-suitability.json",
                "bounded-wgs-execution-v5.json", "common-allele-prefilter-v1.json")))
        sources = {path.relative_to(root).as_posix() for path in paths if path.suffix == ".py"}
        if sources != {name for name in rows if name.startswith("src/genosapiens/")}:
            raise ValueError
        for path in paths:
            relative = path.relative_to(root).as_posix()
            row = rows[relative]
            if any(parent.is_symlink() for parent in (path, *path.parents)) or not path.is_file():
                raise ValueError
            digest = hashlib.sha256()
            count = 0
            with path.open("rb") as stream:
                for block in iter(lambda: stream.read(65536), b""):
                    count += len(block)
                    digest.update(block)
            if count != row["size_bytes"] or digest.hexdigest() != row["sha256"]:
                raise ValueError
    except (ValueError, TypeError, KeyError, OSError):
        raise ValueError("scientific stack projection changed") from None


class _WorkerWaitPending(Exception):
    """The bounded Docker wait has not observed termination yet."""


class _WorkerDeadlineExceeded(Exception):
    """The owned worker reached its measured controller deadline."""


class _WorkerOwnedStop(Exception):
    """The coordinator requested termination of this owned worker."""


@dataclass(frozen=True, slots=True)
class RunnerCaseIdentity:
    case_id: str
    split: str
    scope: str = "runner-supplied-call"


class _DockerConnection(http.client.HTTPConnection):
    def __init__(self, *, timeout=60):
        super().__init__("localhost", timeout=timeout)

    def connect(self):
        self.sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        self.sock.settimeout(self.timeout)
        self.sock.connect("/var/run/docker.sock")


def _engine(method, path, value=None, *, timeout=60):
    from genosapiens.track1_public_validation_coordinator import PublicTrack1CoordinatorError
    connection = _DockerConnection(timeout=timeout)
    try:
        body = None if value is None else json.dumps(value, allow_nan=False).encode("ascii")
        connection.request(method, "/v1.44" + path, body=body, headers={"Content-Type": "application/json"})
        response = connection.getresponse()
        content = response.read(1_000_001)
        if response.status >= 300 or len(content) > 1_000_000:
            raise ValueError
        return {} if not content else json.loads(content)
    except TimeoutError:
        if method == "POST" and path.endswith("/wait?condition=not-running"):
            raise _WorkerWaitPending from None
        raise PublicTrack1CoordinatorError("isolated public container operation failed") from None
    except Exception:
        raise PublicTrack1CoordinatorError("isolated public container operation failed") from None
    finally:
        connection.close()


def _workspace_apparent_bytes(scratch):
    """Sample worker-written files by metadata only; no scientific reads."""
    try:
        if not Path(scratch).is_dir() or Path(scratch).is_symlink():
            return None
        total = 0
        def fail(error):
            raise error
        for directory, directories, files in os.walk(scratch, followlinks=False, onerror=fail):
            for name in (*directories, *files):
                path = Path(directory) / name
                if path.is_symlink():
                    return None
            total += sum((Path(directory) / name).stat().st_size for name in files)
        return total
    except OSError:
        return None


def _worker_observation(container, scratch, started):
    """Allowlisted whole-worker measurements; unavailable values stay null."""
    cpu, memory = None, None
    observer_status = "available"
    try:
        stats = _engine("GET", f"/containers/{container}/stats?stream=false", timeout=5)
        # Docker may return a zero-filled stats object after its cgroup is gone.
        # That is unavailable measurement evidence, not measured zero usage.
        if type(stats.get("read")) is not str or stats["read"].startswith("0001-"):
            raise ValueError
        nanoseconds = stats.get("cpu_stats", {}).get("cpu_usage", {}).get("total_usage")
        if type(nanoseconds) is int and nanoseconds >= 0:
            cpu = str(nanoseconds / 1_000_000_000)
        usage = stats.get("memory_stats", {}).get("usage")
        if type(usage) is int and usage >= 0:
            memory = usage
    except Exception:
        observer_status = "unavailable"
    disk = _workspace_apparent_bytes(scratch)
    return {
        "schema": "genosapiens.track1-public-worker-observation",
        "classification": "public", "scope": "whole-isolated-worker",
        "elapsed_seconds": str(max(0, time.monotonic() - started)),
        "elapsed_method": "controller-monotonic-start-to-observation",
        "elapsed_unit": "seconds", "elapsed_quality": "measured",
        "cpu_seconds": cpu, "cpu_method": "container-cgroup-cumulative-nanoseconds",
        "cpu_unit": "seconds", "cpu_quality": "observed" if cpu is not None else "unavailable",
        "memory_bytes": memory, "memory_method": "sampled-cgroup-usage-including-cache",
        "memory_unit": "bytes", "memory_quality": "observed" if memory is not None else "unavailable",
        "workspace_apparent_bytes": disk, "workspace_method": "sampled-scratch-file-stat-size",
        "workspace_scope": "case-scratch-including-persisted-invocations",
        "workspace_unit": "bytes", "workspace_quality": "observed" if disk is not None else "unavailable",
        "peak_memory_bytes": None, "peak_disk_bytes": None,
        "observer_status": observer_status,
        "controlled_use_authorized": False,
    }


def _wait_for_worker(container, scratch, observe, *, started=None, cancellation_event=None):
    """Bound waits to 30 seconds and report liveness independently of success."""
    started = time.monotonic() if started is None else started
    while True:
        if cancellation_event is not None and cancellation_event.is_set():
            raise _WorkerOwnedStop
        remaining = WORKER_TIMEOUT_SECONDS - (time.monotonic() - started)
        if remaining <= 0:
            raise _WorkerDeadlineExceeded
        try:
            result = _engine("POST", f"/containers/{container}/wait?condition=not-running",
                             timeout=min(WORKER_STATUS_INTERVAL_SECONDS, remaining))
        except _WorkerWaitPending:
            observation = _worker_observation(container, scratch, started)
            observation["lifecycle_status"] = "active"
            observe(observation)
            continue
        if (type(result) is not dict or type(result.get("StatusCode")) is not int
                or result.get("Error") is not None):
            from genosapiens.track1_public_validation_coordinator import PublicTrack1CoordinatorError
            raise PublicTrack1CoordinatorError("isolated public container operation failed")
        observation = _worker_observation(container, scratch, started)
        observation["lifecycle_status"] = "terminated"
        observe(observation)
        return result


def worker_container_spec(repository_root, reference_root, runner_root, scratch_root, request_sha256):
    from genosapiens.track1_public_validation_coordinator import PublicTrack1CoordinatorError
    root, references, runner, scratch = map(Path, (repository_root, reference_root, runner_root, scratch_root))
    if os.getuid() == 0 or not scratch.is_relative_to(references / "scratch/public-build") or any(
        parent.is_symlink() for path in (scratch, runner) for parent in (path, *path.parents)
    ):
        raise PublicTrack1CoordinatorError("invalid isolated public worker environment")
    mounts = [
        {"Type": "bind", "Source": str(runner), "Target": str(root), "ReadOnly": True},
        {"Type": "bind", "Source": str(scratch), "Target": str(scratch), "ReadOnly": False},
    ]
    for relative in PUBLIC_REFERENCE_TREES:
        target = references / relative
        source = (
            references / "derived/common-allele-public-indexes-cache-bounded-v1"
            if relative == "derived/common-allele-public-indexes" else target
        )
        if not source.is_dir() or any(parent.is_symlink() for parent in (source, *source.parents)):
            raise PublicTrack1CoordinatorError("isolated public reference is unavailable")
        mounts.append({"Type": "bind", "Source": str(source), "Target": str(target), "ReadOnly": True})
    return {
        # The existing runtime pin is a local OCI content digest, not a registry
        # manifest digest. Engine Image accepts the independently inspected ID.
        "Image": "sha256:" + RUNTIME_IMAGE.rsplit("@sha256:", 1)[1],
        "User": f"{os.getuid()}:{os.getgid()}",
        "Entrypoint": ["python"],
        "WorkingDir": str(root),
        "Env": [f"PYTHONPATH={root / 'src'}", "PYTHONDONTWRITEBYTECODE=1", "PYTHONUNBUFFERED=1"],
        "Cmd": ["-m", "genosapiens.track1_public_worker", str(root), str(references), str(scratch), request_sha256],
        "Labels": {"org.genosapiens.purpose": "public-validation-scientific-case"},
        "HostConfig": {
            "NetworkMode": "none", "ReadonlyRootfs": True, "CapDrop": ["ALL"],
            "SecurityOpt": ["no-new-privileges"], "Memory": WORKER_MEMORY_BYTES, "MemorySwap": WORKER_MEMORY_BYTES,
            "NanoCpus": 4_000_000_000, "PidsLimit": 128,
            "Mounts": mounts, "Tmpfs": {"/tmp": "rw,nosuid,nodev,size=64m"},
        },
    }


def run_isolated_public_case(case, view, resources, policies, journal, *, policy_id=None, cancellation_event=None):
    """Evaluator-side launch and join; only truth-free request bytes cross Docker."""
    from genosapiens.track1_public_validation_coordinator import (
        _bytes, _sha256, _checkpoint_json, _load_attempt, _failed_attempt,
        publish_validation_artifact, publish_validation_file, PublicTrack1CoordinatorError, _persist_replay_started,
        executable_stack_bytes,
    )
    root = resources.repository_root
    references = resources.reference_root
    work = Path(journal) / "isolated-cases" / case.case_id
    runner = work / "runner"
    scratch = work / "scratch"
    scratch.mkdir(parents=True, exist_ok=True, mode=0o700)
    # Code/configuration are exact current bytes. No plan, roster, recipe, seed,
    # truth or whole fixture repository is included in this mounted tree.
    code_paths = tuple(sorted((root / "src/genosapiens").glob("*.py")))
    config_paths = (
        *tuple(sorted((root / "config/public-resources").glob("*.json"))),
        *tuple(sorted((root / "config/reference-identities").glob("*.json"))),
        *(root / "config/track1" / name for name in (
            "public-ranking-development-policy-family.json", "challenge-evaluator.json",
            "public-population-development-derivation.json",
            "epcr-policy.json", "epcr-numerical-suitability.json",
            "bounded-wgs-execution-v5.json", "common-allele-prefilter-v1.json",
        )),
    )
    for path in (*code_paths, *config_paths):
        if path.is_symlink() or not path.is_file():
            raise PublicTrack1CoordinatorError("isolated worker code is invalid")
        publish_validation_artifact(runner / path.relative_to(root), path.read_bytes())
    materialized = materialize_runner_view(view, resources)
    for path, digest, size in zip(materialized.input_paths, materialized.input_sha256s, materialized.input_sizes, strict=True):
        publish_validation_file(runner / path.relative_to(root), path, expected_sha256=digest, expected_size_bytes=size)
    for descriptor, content in zip(view.resources[len(resources.common_resource_bytes):], materialized.resource_bytes[len(resources.common_resource_bytes):], strict=True):
        relative = (
            Path("tests/fixtures/track1/materialized/hg002-public-background.json")
            if descriptor.kind == "runner-public-background" else
            Path("tests/fixtures/track1/roster/source-formats") / f"{descriptor.descriptor_sha256}.json"
        )
        publish_validation_artifact(runner / relative, content)
    request = _bytes({
        "schema": "genosapiens.track1-public-scientific-case-request",
        "case_id": case.case_id, "split": case.split, "policy_id": policy_id,
        "view": asdict(view),
        "executable_stack": json.loads(executable_stack_bytes(root)),
    })
    invocation = uuid4().hex
    request_path = scratch / "requests" / f"{invocation}.json"
    publish_validation_artifact(request_path, request)
    checkpoint = scratch / "checkpoints" / invocation
    container = None
    exit_code = None
    stage = "container-create"
    operation_failed = False
    termination_reason = "unknown"
    observations = 0

    def observe(value):
        nonlocal observations
        observations += 1
        publish_validation_artifact(Path(journal) / "worker-observations" / invocation / f"{observations:06d}.json", _bytes({
            **value, "request_sha256": _sha256(request),
            "last_established_stage": stage,
            "stage_scope": "controller-container-operation",
            "scientific_stage": None,
        }))
        print(f"public worker lifecycle={value['lifecycle_status']} elapsed_seconds={value['elapsed_seconds']} stage={stage} observer={value['observer_status']}", file=sys.stderr, flush=True)
    try:
        spec = worker_container_spec(root, references, runner, scratch, _sha256(request))
        spec["Cmd"].append(str(request_path))
        publish_validation_artifact(Path(journal) / "worker-specifications" / f"{invocation}.json", _bytes(spec))
        container = _engine("POST", "/containers/create", spec)["Id"]
        inspection = _engine("GET", f"/containers/{container}/json")
        publish_validation_artifact(Path(journal) / "worker-inspections" / f"{invocation}.json", _bytes(inspection))
        stage = "container-start"
        worker_started = time.monotonic()
        _engine("POST", f"/containers/{container}/start")
        stage = "container-wait"
        result = _wait_for_worker(container, scratch, observe, started=worker_started,
                                  cancellation_event=cancellation_event)
        exit_code = result["StatusCode"]
        termination_reason = "normal-exit" if exit_code == 0 else "nonzero-exit"
        try:
            terminal = _engine("GET", f"/containers/{container}/json", timeout=5)
            if terminal.get("State", {}).get("OOMKilled") is True:
                termination_reason = "oom-killed"
        except Exception:
            # Lack of an optional observation does not invalidate a completed
            # wait/result or invent an OOM cause.
            pass
    except _WorkerDeadlineExceeded:
        operation_failed = True
        termination_reason = "controller-deadline-exceeded-owned-stop"
    except _WorkerOwnedStop:
        operation_failed = True
        termination_reason = "controller-owned-stop"
    except KeyboardInterrupt:
        operation_failed = True
        termination_reason = "controller-interrupted-owned-stop"
        if cancellation_event is not None:cancellation_event.set()
    except Exception:
        operation_failed = True
    finally:
        if container is not None:
            # No caller/source logs are fetched or forwarded. This destroys only
            # the container started here; durable bind-mounted artifacts remain.
            try:
                _engine("DELETE", f"/containers/{container}?force=true")
            except Exception:
                operation_failed = True
                stage = "container-removal"
                if termination_reason == "controller-deadline-exceeded-owned-stop":
                    termination_reason = "controller-deadline-exceeded-stop-unconfirmed"
                elif termination_reason in {"controller-owned-stop", "controller-interrupted-owned-stop"}:
                    termination_reason = "controller-owned-stop-unconfirmed"
    publish_validation_artifact(Path(journal) / "worker-lifecycle" / f"{invocation}.json", _bytes({
        "schema": "genosapiens.track1-public-worker-termination",
        "classification": "public", "request_sha256": _sha256(request),
        "case_id": case.case_id, "split": case.split,
        "stage": stage, "exit_code": exit_code,
        "container_operation_failed": operation_failed,
        "termination_reason": termination_reason,
        "observation_count": observations,
        "controlled_use_authorized": False,
    }))
    selected = tuple(runtime for runtime in policies.runtimes if policy_id is None or runtime.candidate.policy_id == policy_id)
    output_path = checkpoint / "result.json"
    expected_context = _bytes({"request_sha256": _sha256(request)})
    context_path = checkpoint / "execution-context.json"
    if any(parent.is_symlink() for parent in (checkpoint, *checkpoint.parents)) or checkpoint.exists() and (
        context_path.is_symlink() or not context_path.is_file()
        or context_path.read_bytes() != expected_context
    ):
        raise PublicTrack1CoordinatorError("isolated worker execution context is invalid")
    if exit_code == 0:
        if output_path.is_symlink():
            raise PublicTrack1CoordinatorError("isolated worker result is invalid")
        output = _checkpoint_json(output_path.read_bytes())
        if set(output) != {"attempts", "replays", "request_sha256"} or output["request_sha256"] != _sha256(request):
            raise PublicTrack1CoordinatorError("isolated worker result is invalid")
        attempts = tuple(replace(_load_attempt(value), scope=case.scope) for value in output["attempts"])
        replays = tuple(tuple(value) for value in output["replays"])
    else:
        # Preserve completed worker attempts if it died during replay or later
        # work. An absent completion record is never inferred from console logs.
        attempts = []
        for runtime in selected:
            saved = checkpoint / "attempts" / f"{case.split}.{case.case_id}.{runtime.candidate.policy_id}.1.json"
            if any(parent.is_symlink() for parent in (saved, *saved.parents)):
                raise PublicTrack1CoordinatorError("isolated worker checkpoint is invalid")
            if saved.is_file() and not saved.is_symlink():
                envelope = _checkpoint_json(saved.read_bytes())
                if (
                    set(envelope) != {"schema", "classification", "validation_only", "execution_context_sha256", "attempt"}
                    or envelope["schema"] != "genosapiens.track1-public-validation-attempt-checkpoint"
                    or envelope["classification"] != "public" or envelope["validation_only"] is not True
                    or envelope["execution_context_sha256"] != _sha256(_bytes({"request_sha256": _sha256(request)}))
                ):
                    raise PublicTrack1CoordinatorError("isolated worker checkpoint is invalid")
                value = envelope["attempt"]
                attempts.append(replace(_load_attempt(value), scope=case.scope))
            else:
                failure = _failed_attempt(case, view, runtime.candidate, "isolated-public-worker")
                if termination_reason in {"controller-deadline-exceeded-owned-stop", "controller-deadline-exceeded-stop-unconfirmed"}:
                    failure = replace(failure, status="timed_out")
                attempts.append(failure)
        attempts = tuple(attempts)
        replays = []
        for attempt in attempts:
            saved = checkpoint / 'determinism' / f'{attempt.split}.{attempt.case_id}.{attempt.policy_id}.json'
            if any(parent.is_symlink() for parent in (saved,*saved.parents)):
                raise PublicTrack1CoordinatorError('isolated replay checkpoint is invalid')
            if not saved.exists():continue
            replay=(attempt.split,attempt.case_id,attempt.policy_id,attempt.execution_sha256)
            if _checkpoint_json(saved.read_bytes()) != json.loads(_bytes({
                    'schema':'genosapiens.track1-public-validation-determinism-checkpoint',
                    'classification':'public','execution_context_sha256':_sha256(expected_context),
                    'candidate_sha256':None,'replay':replay})):
                raise PublicTrack1CoordinatorError('isolated replay checkpoint is invalid')
            replays.append(replay)
        replays=tuple(replays)
    if (
        len(attempts) != len(selected)
        or tuple(item.policy_id for item in attempts) != tuple(item.candidate.policy_id for item in selected)
        or any(item.policy_bindings_sha256 != runtime.candidate.bindings.sha256
               for item, runtime in zip(attempts, selected, strict=True))
        or any(item.case_id != case.case_id or item.split != case.split for item in attempts)
        or any(item.input_sha256s != tuple(value.sha256 for value in view.inputs) or item.resource_descriptor_sha256s != tuple(value.descriptor_sha256 for value in view.resources) for item in attempts)
    ):
        raise PublicTrack1CoordinatorError("isolated worker accounting is invalid")
    for attempt in attempts:
        marker = checkpoint / "replay-started" / f"{attempt.split}.{attempt.case_id}.{attempt.policy_id}.json"
        if any(parent.is_symlink() for parent in (marker, *marker.parents)):
            raise PublicTrack1CoordinatorError("isolated replay marker is invalid")
        if marker.exists():
            content = marker.read_bytes()
            value = _checkpoint_json(content)
            expected = {
                "schema": "genosapiens.track1-public-validation-replay-started",
                "classification": "public", "execution_context_sha256": _sha256(expected_context),
                "case_id": attempt.case_id, "split": attempt.split, "policy_id": attempt.policy_id,
                "attempt_number": attempt.attempt_number, "execution_sha256": attempt.execution_sha256,
                "request_sha256": None, "source_marker_sha256": None,
            }
            if value != expected:
                raise PublicTrack1CoordinatorError("isolated replay marker is invalid")
            _persist_replay_started(Path(journal), attempt,
                request_sha256=_sha256(request), source_marker_sha256=_sha256(content))
    objects = checkpoint / "scientific-objects"
    if objects.is_symlink():
        raise PublicTrack1CoordinatorError("isolated scientific artifact is invalid")
    if objects.is_dir():
        for path in objects.iterdir():
            if path.is_symlink() or not path.is_file():
                raise PublicTrack1CoordinatorError("isolated scientific artifact is invalid")
            publish_validation_file(Path(journal) / "scientific-objects" / path.name, path,
                expected_sha256=path.name, expected_size_bytes=path.stat().st_size)
    return attempts, replays


def _main(arguments):
    from genosapiens.track1_public_validation_coordinator import (
        _bytes, _sha256, _checkpoint_json, _run_development_case, _run_held_out_case,
        _persist_case_attempts, publish_validation_artifact, _failed_attempt,
    )
    root, references, scratch = map(Path, arguments[:3])
    request_path = Path(arguments[4])
    if request_path.parent != scratch / "requests" or request_path.is_symlink():
        raise ValueError
    content = request_path.read_bytes()
    if _sha256(content) != arguments[3]:
        raise ValueError
    request = _checkpoint_json(content)
    if set(request) != {"schema", "case_id", "split", "policy_id", "view", "executable_stack"} or request["schema"] != "genosapiens.track1-public-scientific-case-request":
        raise ValueError
    if request["split"] not in {"development", "held_out"} or any(name != "lo" for _, name in socket.if_nameindex()):
        raise ValueError
    view = RankerCaseView(
        tuple(RankerInput(InputKind(value["kind"]), value["sha256"]) for value in request["view"]["inputs"]),
        tuple(RankerResource(**value) for value in request["view"]["resources"]),
    )
    stack_content = _bytes(request["executable_stack"])
    verify_scientific_stack_projection(root, stack_content)
    model_loci = load_public_model_loci(root, references)
    policies = load_public_policy_family(root, model_loci=model_loci)
    if request["split"] == "development" and request["policy_id"] is not None:
        runtimes = tuple(item for item in policies.runtimes if item.candidate.policy_id == request["policy_id"])
        if len(runtimes) != 1:
            raise ValueError
        policies = replace(policies, runtimes=runtimes)
    selected = tuple(item for item in policies.runtimes
                     if request["policy_id"] is None or item.candidate.policy_id == request["policy_id"])
    if not selected:
        raise ValueError
    case = RunnerCaseIdentity(request["case_id"], request["split"])
    journal = scratch / "checkpoints" / request_path.stem
    publish_validation_artifact(journal / "execution-context.json", _bytes({"request_sha256": _sha256(content)}))
    try:
        resources = replace(load_public_resource_context(root, references, scratch), executable_stack_content=stack_content)
    except Exception as error:
        attempts = tuple(_failed_attempt(case, view, item.candidate,
                                         "public-resource-authentication", error=error)
                         for item in selected)
        replays = ()
        publish_validation_artifact(journal / "resource-authentication-failure.json", _bytes({
            "schema": "genosapiens.track1-public-worker-resource-failure",
            "classification": "public", "request_sha256": _sha256(content),
            "error_code": "public-resource-authentication-failed",
            "preparation_stage": attempts[0].preparation_stage,
            "controlled_use_authorized": False,
        }))
    else:
        if request["split"] == "development":
            attempts, replays = _run_development_case(case, view, resources, policies, journal)
        else:
            runtime = next(item for item in policies.runtimes if item.candidate.policy_id == request["policy_id"])
            attempt, replay = _run_held_out_case(case, view, runtime, resources, policies, journal)
            attempts, replays = (attempt,), (() if replay is None else (replay,))
    _persist_case_attempts(journal, attempts, replays, None)
    publish_validation_artifact(journal / "result.json", _bytes({
        "attempts": tuple(asdict(item) for item in attempts),
        "replays": replays, "request_sha256": _sha256(content),
    }))


if __name__ == "__main__":
    import sys
    try:
        _main(sys.argv[1:])
    except Exception:
        print("isolated public scientific case failed", file=sys.stderr)
        raise SystemExit(1)
