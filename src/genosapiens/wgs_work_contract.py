"""Current focused scientific work contract, separate from historical receipts."""
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path

from genosapiens.bounded_json import read_bounded_bytes
from genosapiens.population_filter_candidate import (
    PopulationWorkPolicy, WORK_POLICY_PATH, WORK_POLICY_SHA256, load_population_work_policy,
)

PATH = 'config/track1/focused-work-contract-v1.json'
ROLES = frozenset({'population_work_policy', 'focused_work_contract'})
FRESH_INVOCATION = 'focused-fresh-subsample'


def fresh_sample_binding(sample_name, input_sha256):
    """Bind the original public/synthetic sample without relabeling its source."""
    if (type(sample_name) is not str or not sample_name or len(sample_name) > 1024
            or any(c.isspace() or ord(c) < 32 for c in sample_name)
            or type(input_sha256) is not str or len(input_sha256) != 64
            or any(c not in '0123456789abcdef' for c in input_sha256)):
        raise WorkContractError('fresh selected sample identity is invalid')
    content = (json.dumps(dict(schema='genosapiens.public-wgs-qualification-sample',
        sample_name=sample_name, input_sha256=input_sha256),
        sort_keys=True, separators=(',', ':')) + '\n').encode('ascii')
    return hashlib.sha256(content).hexdigest()


class WorkContractError(ValueError):
    """The scientific work policy, format or domain is not current."""


@dataclass(frozen=True, slots=True)
class FocusedWorkContract:
    policy: PopulationWorkPolicy
    sha256: str
    content: bytes

    def bindings(self):
        return dict(population_work_policy=self.policy.sha256, focused_work_contract=self.sha256)


def load_focused_work_contract(repository_root=None):
    root = Path(__file__).resolve().parents[2] if repository_root is None else Path(repository_root)
    try:
        content = read_bounded_bytes(root / PATH, maximum_bytes=65536)
        expected = dict(schema='genosapiens.focused-scientific-work-contract', version=1,
            analysis_profile='focused', population_work_policy_sha256=WORK_POLICY_SHA256,
            selected_reconstruction=dict(schema='genosapiens.selected-reconstruction-index', version=1,
                context_version=2),
            work_filter=dict(schema='genosapiens.context-population-work-filter', version=1),
            domain='original-closed-context-with-authenticated-retained-member-selection',
            mixed_pairs='retain-original-members-unless-equivalent-omission-is-proven',
            qualification_scope=FRESH_INVOCATION, full_wgs_runtime_validated=False)
        canonical = (json.dumps(expected, sort_keys=True, separators=(',', ':')) + '\n').encode('ascii')
        if content != canonical:
            raise ValueError
        policy = load_population_work_policy(read_bounded_bytes(root / WORK_POLICY_PATH, maximum_bytes=65536),
            expected_sha256=WORK_POLICY_SHA256)
        return FocusedWorkContract(policy, hashlib.sha256(content).hexdigest(), content)
    except (OSError, ValueError, TypeError):
        raise WorkContractError('focused work contract is missing, stale or invalid') from None


def require_work_bindings(values, contract):
    if type(contract) is not FocusedWorkContract or any(values.get(k) != v for k, v in contract.bindings().items()):
        raise WorkContractError('focused work policy or domain binding differs')


def publication_work_bindings(parents, authorized_resources):
    """Retain authenticated D41 parents at the final publication boundary.

    Historical publications without D41 stay readable. A D41 release or artifact
    requires both current bindings, including when the evaluated domain is empty.
    """
    if not ROLES.intersection(parents) and WORK_POLICY_SHA256 not in authorized_resources:
        return {}
    work = load_focused_work_contract()
    require_work_bindings(parents, work)
    if not set(work.bindings().values()) <= set(authorized_resources):
        raise WorkContractError('publication work contract is not authorized')
    return work.bindings()
