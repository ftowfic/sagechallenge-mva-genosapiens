"""Evaluator-owned finite queries over authenticated completed executions.

One reader owns each complete replay. Only the declared target's recovery and
the fixed population probes survive that reader; rich H-sized results never
enter an attempt, checkpoint, or process-global cache.
"""
from dataclasses import asdict
import hashlib

from genosapiens.evidence_model import _sha
from genosapiens.evidence_serialization import _canonical
from genosapiens.wgs_evidence_store import EvidenceStoreError

MAX_QUERIES = 1024
MAX_PROBES = 256
MAX_CACHE_BYTES = 8 * 1024**2


class AuthenticatedValidationRecovery:
    def __init__(self, *, objects, resources, policies, scratch_root, cases,
                 population_selector, frequency_identities, resources_for_attempt=None):
        self.objects, self.resources, self.policies = objects, resources, policies
        self.scratch_root, self.selector = scratch_root, population_selector
        self.cases = {case.case_id: case for case in cases}
        self.probes = tuple(sorted(set(frequency_identities)))
        if (len(self.cases) != len(cases) or len(cases) > MAX_QUERIES or len(self.probes) > MAX_PROBES
                or len(_canonical(tuple((c.case_id, c.truth_members) for c in cases))) > MAX_CACHE_BYTES):
            raise EvidenceStoreError("validation recovery query capacity exceeded", capacity=True)
        for key in self.probes:
            _sha(key)
        self.resources_for_attempt = resources_for_attempt
        self.cache, self.cache_bytes = {}, 0

    def _get(self, attempt):
        from genosapiens.track1_validation_runner import ValidationAttempt
        from genosapiens.wgs_public_store import open_public_execution
        if type(attempt) is not ValidationAttempt or attempt.recovery_binding is None:
            raise EvidenceStoreError("completed normalized attempt required for recovery")
        case = self.cases.get(attempt.case_id)
        if case is None or (case.split, case.scope) != (attempt.split, attempt.scope):
            raise EvidenceStoreError("validation recovery case differs from evaluator declaration")
        attempt_bytes = _canonical(asdict(attempt)).encode("ascii")
        key = hashlib.sha256(attempt_bytes).hexdigest()
        if key not in self.cache:
            if len(self.cache) >= MAX_QUERIES:
                raise EvidenceStoreError("validation recovery cache exceeds capacity", capacity=True)
            resources = (self.resources if self.resources_for_attempt is None
                else self.resources_for_attempt(attempt))
            with open_public_execution(self.objects, expected_sha256=attempt.execution_sha256,
                    resources=resources, policies=self.policies, scratch_root=self.scratch_root,
                    population_selector=self.selector) as execution:
                manifest = execution.manifest
                if (manifest.case_id != attempt.case_id
                        or manifest.candidate_stack_sha256 != attempt.policy_bindings_sha256
                        or manifest.runner_input_sha256s != attempt.input_sha256s
                        or manifest.runner_resource_descriptor_sha256s != attempt.resource_descriptor_sha256s
                        or tuple((name, "succeeded") for name, _ in manifest.output_bindings) != attempt.stage_statuses
                        or any(getattr(attempt, name) != value for name, value in execution.metrics.items())):
                    raise EvidenceStoreError("attempt differs from authenticated normalized execution")
                target = tuple(sorted(case.truth_members))
                recovery = execution.recovery.target_recovery(target)
                frequencies = tuple(execution.recovery.frequency_observations(self.probes))
                size = len(_canonical((target, recovery, frequencies, manifest)))
                if self.cache_bytes + size > MAX_CACHE_BYTES:
                    raise EvidenceStoreError("validation recovery cache exceeds capacity", capacity=True)
                value = (target, recovery, frequencies, manifest)
            # The outer read scope must close successfully before caching.
            self.cache[key] = value
            self.cache_bytes += size
        return self.cache[key]

    def __call__(self, attempt, target):
        saved, recovery, _, _ = self._get(attempt)
        if target != saved:
            raise EvidenceStoreError("validation target differs from declared evaluator query")
        return recovery

    def frequency_observations(self, attempt, identities):
        requested = set(identities)
        if not requested <= set(self.probes):
            raise EvidenceStoreError("undeclared validation population probe")
        if attempt.recovery_binding is None:
            return ()
        return tuple(value for value in self._get(attempt)[2] if value[0] in requested)

    def manifest(self, attempt):
        return self._get(attempt)[3]
