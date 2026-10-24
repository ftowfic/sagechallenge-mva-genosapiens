"""Streaming logical-policy digest; normalized policy storage is authoritative."""
import hashlib
from genosapiens.evidence_model import _typed
from genosapiens.evidence_serialization import _canonical
from genosapiens.population_policy import PopulationPolicy
from genosapiens.wgs_failures import WgsFailure


class SharedPolicyError(WgsFailure):
    def __init__(self, message, *, category="malformed_evidence"):
        super().__init__("bounded_products", category,
            operation="bounded_product_validate", exception_kind="project_failure")
        RuntimeError.__init__(self, message)



def _bytes(value):
    return (_canonical(value) + "\n").encode("ascii")



class _LogicalDigest:
    def __init__(self, policy_id):
        self.policy_id = policy_id
        self.digest = hashlib.sha256(b'{"policy":{"contexts":[')
        self.contexts = self.exceptions = 0
        self.in_exceptions = False

    def add_context(self, value):
        if self.in_exceptions:
            raise SharedPolicyError("population context order changed")
        if self.contexts:
            self.digest.update(b",")
        self.digest.update(_bytes(value)[:-1])
        self.contexts += 1

    def _exceptions(self):
        if not self.in_exceptions:
            self.digest.update(b'],"exceptions":[')
            self.in_exceptions = True

    def add_exception(self, value):
        self._exceptions()
        if self.exceptions:
            self.digest.update(b",")
        self.digest.update(_bytes(value)[:-1])
        self.exceptions += 1

    def finish(self):
        self._exceptions()
        self.digest.update(b'],"policy_id":' + _bytes(self.policy_id)[:-1]
                           + b'},"schema":"genosapiens.population-policy"}\n')
        return self.digest.hexdigest()



def logical_policy_sha256(policy):
    _typed(policy, PopulationPolicy)
    digest = _LogicalDigest(policy.policy_id)
    for value in policy.contexts:
        digest.add_context(value)
    for value in policy.exceptions:
        digest.add_exception(value)
    return digest.finish()

