"""Invocation-owned checked metadata with immutable cached values."""
import json

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_published_genotype import canonical
from genosapiens.track2_tissue_records import identity


class TissueAdmission:
    """The caller owns the read-only connection and closes it after invocation."""

    def __init__(self, descriptor, policy, connection, *, counters=None):
        require(descriptor['query_policy'] == policy, 'Tissue admission policy changed.')
        self.descriptor_bytes = canonical(descriptor)
        self.policy_bytes = canonical(policy)
        self.descriptor_sha256 = identity(descriptor)
        self.policy_sha256 = identity(policy)
        self.connection = connection
        self.counters = counters if counters is not None else {}
        self._drugs = {}
        self._genetic = {}
        self._features = None
        self._membership = None

    def check(self, descriptor, policy):
        require(identity(descriptor) == self.descriptor_sha256, 'Tissue admission source changed.')
        require(identity(policy) == self.policy_sha256, 'Tissue admission policy changed.')

    def _count(self, name):
        self.counters[name] = self.counters.get(name, 0) + 1

    def features(self):
        if self._features is None:
            self._features = canonical([dict(r) for r in self.connection.execute(
                'SELECT pr_gene_id,symbol,measured_landmark FROM lincs_features WHERE measured_landmark=1 ORDER BY pr_gene_id')])
            self._count('feature_selections')
        return json.loads(self._features)

    def membership(self):
        from genosapiens.track2_tissue_resources import planning_membership
        if self._membership is None:
            self._membership = canonical(planning_membership(json.loads(self.descriptor_bytes), self.connection))
            self._count('membership_selections')
        return json.loads(self._membership)

    def drugs(self, cell):
        from genosapiens.track2_lincs_selection import select_drug_conditions
        key = (self.descriptor_sha256, self.policy_sha256, cell)
        if key not in self._drugs:
            self._drugs[key] = canonical(select_drug_conditions(self.connection, cell, None, json.loads(self.policy_bytes)))
            self._count('drug_metadata_selections')
        else:
            self._count('drug_metadata_reuses')
        return json.loads(self._drugs[key])

    def genetic(self, gene):
        from genosapiens.track2_lincs_selection import select_tissue_genetic_signatures
        key = (self.descriptor_sha256, self.policy_sha256, gene)
        if key not in self._genetic:
            self._genetic[key] = canonical(select_tissue_genetic_signatures(
                self.connection, gene, json.loads(self.policy_bytes), self.membership()))
            self._count('genetic_metadata_selections')
        else:
            self._count('genetic_metadata_reuses')
        return json.loads(self._genetic[key])
