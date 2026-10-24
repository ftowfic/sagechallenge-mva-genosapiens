"""Complete saved hypothesis universe verified through its original owners."""
import sqlite3

from genosapiens.wgs_context_evidence import VerifiedContextEvidence
from genosapiens.wgs_evidence_replay import EvidenceReplay
from genosapiens.wgs_evidence_store import EvidenceStoreError
from genosapiens.wgs_hypothesis_evidence import CompleteHypothesisEvidence
from genosapiens.wgs_pair_verification import VerifiedCompletePairEvidence


class VerifiedHypothesisEvidence(CompleteHypothesisEvidence):
    def __init__(self, reader, reference, members, facts, *, policy, gene_evidence,
                 policy_sha256, observation_ids, outside_context_partner_uncertainty=True, prepared_gene_evidence=None, prepared_policies=None):
        if (reference.kind != "complete-hypothesis-universe-v1"
                or type(members) is not VerifiedContextEvidence
                or members.writer is not reader or facts.writer is not reader
                or facts.published is None):
            raise EvidenceStoreError("foreign complete hypothesis evidence")
        reader._require_reference(reference)
        try:
            super().__init__(EvidenceReplay(reader), members, facts, policy=policy, gene_evidence=gene_evidence,
                policy_sha256=policy_sha256, observation_ids=observation_ids,
                outside_context_partner_uncertainty=outside_context_partner_uncertainty, name=reference.name,
                prepared_gene_evidence=prepared_gene_evidence,prepared_policies=prepared_policies)
            if self.sha256 != reference.sha256:
                self.close()
                raise EvidenceStoreError("complete hypothesis evidence root changed")
        except (ValueError, TypeError, KeyError, sqlite3.Error):
            raise EvidenceStoreError("invalid complete hypothesis evidence") from None

    def _prepare_pairs(self, writer, members, facts, *, policy, gene_evidence, name):
        reader = writer.reader
        return VerifiedCompletePairEvidence(reader, reader.reference(name, kind="complete-pair-evidence-v1"),
            members, facts, policy=policy, gene_evidence=gene_evidence,prepared_gene_evidence=self.prepared_gene_evidence,prepared_policies=self.prepared_policies)
