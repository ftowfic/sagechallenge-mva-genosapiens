#!/usr/bin/env python3
"""Execute the bounded public Track 1 validation workload.

The script deliberately accepts only repository/public-reference roots and
runner-safe roster case identifiers.  Evaluator truth is joined by the outer
coordinator after scientific execution; it is never passed to preparation or
ranking.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import sys
import time
from uuid import uuid4

from genosapiens.track1_public_validation import (
    load_public_policy_family,
    load_public_resource_context,
    PublicTrack1ValidationError,
    preflight_public_population_policy,
)
from genosapiens.track1_roster import load_submission_acceptance_roster
from genosapiens.track1_public_validation_coordinator import (
    FinalizedPublicValidation,
    PublicValidationStoppedOutcome,
    finalize_public_validation,
    run_public_candidate_validation,
    load_validation_outcome_checkpoint,
    publish_validation_artifact,
    record_finalization_failure,
    PublicTrack1CoordinatorError,
    executable_stack_bytes,
    require_current_public_acceptance_evidence,
)
from genosapiens.track1_validation_plan import (
    load_validation_plan_from_pin, load_validation_artifact_bundle_from_pin,
)
from genosapiens.track1_validation_runner import (
    development_ranker_view,
)


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repository-root", type=Path, required=True)
    parser.add_argument("--reference-root", type=Path, required=True)
    parser.add_argument("--scratch-root", type=Path, required=True)
    parser.add_argument("--case-id")
    parser.add_argument("--policy-id")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--checkpoint-root", type=Path)
    parser.add_argument("--resume-outcome-sha256")
    parser.add_argument("--continuation", type=Path)
    return parser.parse_args()


def _publish(path: Path, content: bytes) -> str:
    return publish_validation_artifact(path, content)


def _report_progress(attempt, required_model_coverage=()) -> None:
    substage = ("" if attempt.preparation_stage is None else
                f" preparation_stage={attempt.preparation_stage}")
    if attempt.execution_stage is not None:
        substage += f" execution_stage={attempt.execution_stage}"
    print(
        f"finished {attempt.split} {attempt.case_id} {attempt.policy_id} status={attempt.status}{substage}",
        file=sys.stderr, flush=True,
    )
    if attempt.population_lookup_counts:
        print("population_lookup " + " ".join(
            f"{name}={count}" for name, count in attempt.population_lookup_counts
        ), file=sys.stderr, flush=True)
    missing = sorted(set(required_model_coverage) - {
        name for name, count in attempt.model_coverage_counts if count > 0
    })
    if missing:
        print("missing_required_model_coverage " + ",".join(missing), file=sys.stderr, flush=True)


def _verify_smoke_execution(journal, attempt, case, view, runtime, resources, policies):
    """Authenticate current products without another scientific invocation."""
    from genosapiens.wgs_public_execution import OUTPUT_ROLES
    from genosapiens.track1_public_validation_coordinator import _verify_saved_attempt
    if (attempt.status != "succeeded" or set(dict(attempt.stage_statuses)) != OUTPUT_ROLES
            or any(status != "succeeded" for _, status in attempt.stage_statuses)):
        raise ValueError
    _verify_saved_attempt(journal, attempt, case, view, runtime, resources, policies)


def _run_preparation_smoke(root, resources, roster, policies, case_id, output):
    """Use the cohort's normal isolated worker and retain scientific bytes locally."""
    from genosapiens import track1_public_worker as worker
    from genosapiens.track1_public_validation_coordinator import (
        _bytes, _execution_context_bytes, _persist_case_attempts, _load_attempt,
    )
    started = time.monotonic()
    case = next(c for c in roster.cases if c.case_id == case_id and c.split == "development")
    plan = load_validation_plan_from_pin(repository_root=root)
    load_validation_artifact_bundle_from_pin(repository_root=root)
    journal = resources.scratch_root / "preparation-smokes" / uuid4().hex
    _publish(journal / "execution-context.json", _execution_context_bytes(root, plan, roster, resources, policies))
    view = development_ranker_view(roster, case.case_id,
                                  runner_resources=resources.common_runner_resources)
    attempts, replays = worker.run_isolated_public_case(case, view, resources, policies, journal)
    _persist_case_attempts(journal, attempts, replays, None)
    for attempt in attempts:
        _report_progress(attempt, case.required_model_coverage)
        path = journal / "attempts" / f"{attempt.split}.{attempt.case_id}.{attempt.policy_id}.1.json"
        if _load_attempt(json.loads(path.read_bytes())["attempt"]) != attempt:
            raise PublicTrack1CoordinatorError("public preparation smoke checkpoint is invalid")
    status = "failed"
    try:
        from genosapiens.track1_roster import load_submission_acceptance_profile, _case_bytes
        profile, _profile_bytes = load_submission_acceptance_profile(root)
        expectation = next(row for row in profile["cases"] if row["case_id"] == case_id)
        if (len(attempts) != len(policies.runtimes) or len(replays) != len(policies.runtimes)
                or any(a.status != "succeeded"
                       or a.required_enumerations != a.completed_enumerations for a in attempts)):
            raise ValueError
        runtimes = {item.candidate.policy_id: item for item in policies.runtimes}
        for attempt in attempts:
            if any(dict(attempt.model_coverage_counts).get(name, 0) < 1
                   for name in case.required_model_coverage):
                raise ValueError
            if expectation["positive_exact_export"]:
                ranks = [i + 1 for i, group in enumerate(attempt.eligible_member_groups)
                         if hashlib.sha256(_case_bytes(group)).hexdigest() == expectation["expected_genotype_sha256"]]
                if not ranks or min(ranks) > profile["positive_export_contract"]["eligible_rank_at_most"]:
                    raise ValueError
            _verify_smoke_execution(journal, attempt, case, view,
                                    runtimes[attempt.policy_id], resources, policies)
        inspections = list((journal / "worker-inspections").iterdir())
        if len(inspections) != 1:
            raise ValueError
        inspection = json.loads(inspections[0].read_bytes())
        host, config = inspection["HostConfig"], inspection["Config"]
        if (inspection["Image"] != "sha256:" + worker.RUNTIME_IMAGE.rsplit("@sha256:", 1)[1]
                or host["NetworkMode"] != "none" or host["ReadonlyRootfs"] is not True
                or config["User"].split(":")[0] in {"", "0", "root"}):
            raise ValueError
        expected_public = {str(resources.reference_root / p) for p in worker.PUBLIC_REFERENCE_TREES}
        mounts = inspection["Mounts"]
        if (not expected_public <= {m["Source"] for m in mounts if m["RW"] is False}
                or any(m["Source"] == "/var/run/docker.sock" or m["Source"] == str(root)
                       or m["Source"] == str(resources.reference_root) or "/data/" in m["Source"] for m in mounts)):
            raise ValueError
        runner = journal / "isolated-cases" / case.case_id / "runner"
        if any("causal-truth" in p.name or "recipe" in p.name or p.name in {"roster.json", "submission-acceptance-profile.json"}
               for p in runner.rglob("*")):
            raise ValueError
        status = "passed"
    except Exception:
        pass
    inventory = []
    for path in sorted(journal.rglob("*")):
        if path.is_file() and not path.is_symlink():
            with path.open("rb") as stream:
                digest = hashlib.file_digest(stream, "sha256").hexdigest()
            inventory.append({"path": str(path.relative_to(journal)),
                "size_bytes": path.stat().st_size, "sha256": digest})
    receipt = _bytes({
        "schema": "genosapiens.track1-public-preparation-smoke", "classification": "public",
        "status": status, "split": "development", "case_id": case.case_id,
        "validation_only": True, "controlled_use_authorized": False,
        "execution_context_sha256": hashlib.sha256((journal / "execution-context.json").read_bytes()).hexdigest(),
        "measurements": {"total_wall_seconds": str(time.monotonic() - started),
                         "attempts": tuple(asdict(a.measurements) for a in attempts)},
        "attempt_statuses": tuple({"policy_id": a.policy_id, "status": a.status,
                                   "preparation_stage": a.preparation_stage,
                                   "execution_stage": a.execution_stage} for a in attempts),
        "artifacts": inventory,
        "measurement_limits": "Attempt CPU includes preparation and subprocess CPU; memory is process high-water and disk is preparation output, not cgroup or filesystem peaks.",
    })
    _publish(journal / "smoke.json", receipt)
    if output is not None:
        _publish(output, receipt)
    print(f"public preparation smoke status={status} policies={len(attempts)}", flush=True)
    return 0 if status == "passed" else 1


def main() -> int:
    arguments = _arguments()
    if arguments.continuation is not None and (arguments.case_id is not None or arguments.resume_outcome_sha256 is not None):
        raise ValueError("continuation requires a fresh campaign")
    root = arguments.repository_root.resolve(strict=True)
    policy_family_sha256 = hashlib.sha256((root / "config/track1/public-ranking-development-policy-family.json").read_bytes()).hexdigest()
    references = arguments.reference_root.resolve(strict=True)
    scratch = arguments.scratch_root.resolve(strict=True)
    if not scratch.is_relative_to(references / "scratch/public-build"):
        raise PublicTrack1CoordinatorError("invalid public validation scratch root")
    startup_id = uuid4().hex
    startup_root = scratch / "startup"
    startup_stack_sha256 = hashlib.sha256(executable_stack_bytes(root)).hexdigest()
    startup_stage = "population-policy-binding"
    def startup_record(status):
        content = (json.dumps({
            "schema": "genosapiens.track1-public-validation-startup",
            "classification": "public", "status": status,
            "stage": startup_stage,
            "executable_stack_sha256": startup_stack_sha256,
            "policy_family_sha256": policy_family_sha256,
            "controlled_use_authorized": False,
        }, sort_keys=True, separators=(",", ":")) + "\n").encode("ascii")
        _publish(startup_root / f"{startup_id}.{status}.json", content)
    startup_record("started")
    try:
        if arguments.case_id is None:
            startup_stage = "current-acceptance-evidence"
            require_current_public_acceptance_evidence(root, executable_stack_sha256=startup_stack_sha256)
        startup_stage = "population-policy-binding"
        preflight_public_population_policy(root)
        startup_stage = "public-resource-authentication"
        resources = load_public_resource_context(root, references, scratch)
        policies = load_public_policy_family(root, model_loci=resources.model_loci)
    except Exception:
        startup_record("failed")
        raise
    startup_record("authenticated")
    roster = load_submission_acceptance_roster(repository_root=root)
    if arguments.case_id is None:
        if arguments.policy_id is not None or arguments.output is None:
            raise ValueError("full execution requires --output and no --policy-id")
        plan = load_validation_plan_from_pin(repository_root=root)
        journal = arguments.checkpoint_root or Path(str(arguments.output.with_suffix("")) + ".checkpoint")

        if arguments.resume_outcome_sha256 is not None:
            outcome = load_validation_outcome_checkpoint(
                journal,
                expected_sha256=arguments.resume_outcome_sha256,
                repository_root=root, plan=plan, roster=roster,
                resources=resources, policies=policies,
            )
        else:
            outcome = run_public_candidate_validation(
                repository_root=root,
                plan=plan,
                roster=roster,
                resources=resources,
                policies=policies,
                progress=lambda attempt: _report_progress(attempt, next(
                    case.required_model_coverage for case in roster.cases if case.case_id == attempt.case_id)),
                workers=arguments.workers,
                checkpoint_root=journal,
                continuation=None if arguments.continuation is None else json.loads(arguments.continuation.read_bytes()),
            )
        result = (
            outcome
            if type(outcome) is PublicValidationStoppedOutcome or outcome.candidate is None
            else finalize_public_validation(
                outcome,
                repository_root=root,
                roster=roster,
                resources=resources,
                checkpoint_root=journal,
            )
        )
        output = arguments.output
        publications = {}
        try:
            if type(result) is FinalizedPublicValidation:
                base = output.with_suffix("")
                components = (
                    ("candidate-validation", result.candidate_validation.as_bytes()),
                    ("candidate", result.candidate_validation.candidate.as_bytes()),
                    ("acceptance", result.acceptance.as_bytes()),
                )
                if result.base_release is not None:
                    components = (*components, ("base-release", result.base_release.as_bytes()))
                for name, content in components:
                    path = Path(str(base) + f".{name}.json")
                    publications[str(path)] = _publish(path, content)
            publications[str(output)] = _publish(output, result.as_bytes())
        except Exception:
            record_finalization_failure(journal, "packet-publication", outcome.sha256)
            raise PublicTrack1CoordinatorError("public validation packet publication failed") from None
        print(json.dumps({
            "status": result.status,
            "publications": publications,
            "sha256": result.sha256,
        }, sort_keys=True, separators=(",", ":")))
        return 0 if result.status == "accepted-and-frozen" else 1

    if arguments.workers != 1 or arguments.resume_outcome_sha256 is not None or arguments.checkpoint_root is not None or arguments.policy_id is not None:
        raise ValueError("development preparation smoke requires both policies and one worker")
    return _run_preparation_smoke(root, resources, roster, policies, arguments.case_id, arguments.output)


if __name__ == "__main__":
    try:
        exit_code = main()
    except (PublicTrack1CoordinatorError, PublicTrack1ValidationError) as error:
        print(str(error), file=sys.stderr)
        exit_code = 1
    except Exception:
        print("public validation execution failed", file=sys.stderr)
        exit_code = 1
    raise SystemExit(exit_code)
