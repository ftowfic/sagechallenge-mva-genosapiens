"""Explicit scientific scope policies, independent of upstream admission."""
from dataclasses import dataclass
from enum import Enum
import hashlib

from genosapiens.wgs_ledger import canonical


class ScienceScopeMode(str, Enum):
    EXHAUSTIVE = 'exhaustive'
    RAPID_PHENOTYPE = 'rapid_phenotype'
    RAPID_STRICT_ASSOCIATION = 'rapid_strict_association'


@dataclass(frozen=True, slots=True)
class ScienceScopePolicy:
    mode: ScienceScopeMode = ScienceScopeMode.EXHAUSTIVE
    target_contexts: int = 1500

    def __post_init__(self):
        if (type(self.mode) is not ScienceScopeMode or type(self.target_contexts) is not int
                or not 1 <= self.target_contexts <= 2**31 - 1):
            raise ValueError('invalid scientific scope policy')

    @property
    def rapid(self):
        return self.mode is not ScienceScopeMode.EXHAUSTIVE

    def as_bytes(self):
        return canonical(dict(schema='genosapiens.track1-scientific-scope', version=1,
            mode=self.mode.value, target_contexts=self.target_contexts if self.rapid else None,
            molecular_protection='applicable-clinical-curated-functional-supported-native',
            uncertain_protection='retain', complete_membership=True,
            missing_phenotype='retain-unless-explicit-strict-association-deferral',
            boundary_ties='retain-all-exact-score-ties', protected_overflow='retain-and-disclose',
            inheritance_uncertainty='not-negative-evidence',
            strict_association_categories={
                'gencc_submissions': ['definitive', 'strong'],
                'evidence_clingen': ['definitive', 'strong'],
                'evidence_gene2phenotype': ['definitive', 'strong'],
                'evidence_genomics_england': ['green'],
                'orphadata_product6': ['validated'],
                'evidence_orphanet': ['validated'],
            }, strict_requires_direct_positive_causal_assertion=True,
            strict_weak_categories=['moderate', 'limited', 'supportive', 'amber'],
            strict_unknown_categories='separate-deferral-reason',
            deferred_interpretation='unevaluated-not-scientifically-invalid'))

    @property
    def sha256(self):
        return hashlib.sha256(self.as_bytes()).hexdigest()


DEFAULT_SCIENCE_SCOPE = ScienceScopePolicy()
