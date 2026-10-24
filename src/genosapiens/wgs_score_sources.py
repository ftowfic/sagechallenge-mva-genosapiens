"""Source-complete pair-member score facts computed by the existing adapters.

Source proofs are retained once per context/member. Numeric projections omit
only those proofs after exact atomic ownership closure. Pair combination must
still check ownership across members and retain technical-charge order.
"""
from dataclasses import dataclass
import hashlib
import json

from genosapiens.evidence_fusion import (
    EvidenceCharge, FusionPolicy, admission_decision_charges, requirement_member_charges,
    clinical_member_charges, technical_member_charges, disease_context_charges,
    canonical_evidence_charges, policy_as_json, _sum,
    EvidenceFusionError,
)
from genosapiens.evidence_model import _typed
from genosapiens.evidence_serialization import _canonical
from genosapiens.wgs_evidence_replay import EvidenceReplay
from genosapiens.wgs_evidence_store import EvidenceStoreError
from genosapiens.wgs_score_decisions import ScoreChargeValue, ScoreDecisionInputs
from genosapiens.genotype_evaluation import prepare_gene_evidence
from genosapiens.wgs_gene_evidence import CaseGeneEvidence
from genosapiens.wgs_validation_scope import member_validation_scope

MAX_SOURCE_WINDOW_BYTES = 64 * 1024**2


@dataclass(frozen=True, slots=True, repr=False)
class ChargeOwnership:
    lineage_id: str
    evidence_key: str
    charge_sha256: str


@dataclass(frozen=True, slots=True, repr=False)
class SourceScoreProfile:
    # Nontechnical sums are exact, per kind, after ordinary ownership closure.
    # Member slots are assigned only when a complete pair is bound.
    subject_id: str | None
    nontechnical: tuple[ScoreChargeValue, ...]
    technical: tuple[tuple[str, str, str | None], ...]
    ownership: tuple[ChargeOwnership, ...]


def project_member_charges(charges, allele_id):
    charges = canonical_evidence_charges(charges)
    if any(c.member_id != allele_id for c in charges):
        raise EvidenceStoreError("pair member charge has a foreign scope")
    if any(c.evidence_key.startswith(("genotype-model-", "population-policy-")) for c in charges):
        raise EvidenceStoreError("member/shared score charge overlaps a pair-dependent namespace")
    nontechnical = []
    for kind in sorted({c.kind for c in charges if c.kind not in ("technical", "conflict")}):
        raw = _sum(tuple(c for c in charges if c.kind == kind))
        nontechnical.append(ScoreChargeValue(kind, None if raw is None else str(raw), None))
    return SourceScoreProfile(allele_id, tuple(nontechnical),
        tuple((c.lineage_id, c.evidence_key, c.value) for c in charges if c.kind in ("technical", "conflict")),
        tuple(ChargeOwnership(c.lineage_id,c.evidence_key,hashlib.sha256(_canonical(c).encode("ascii")).hexdigest())
              for c in charges))


def pair_score_inputs(left, right, shared, *, completeness, phase, phenotype_units,
                      model_state, population_contribution, policy):
    """Bind actual member slots and every atomic owner before numeric reuse.

    The only collapsed terms are exact commutative nontechnical sums. Each
    technical amount remains in original ownership order for saturation checks.
    Inheritance and population adapters use disjoint fixed key namespaces;
    their source records remain available through the model/population owners.
    """
    if (left.subject_id is None or right.subject_id is None or left.subject_id >= right.subject_id
            or shared.subject_id is not None):
        raise EvidenceStoreError("pair score source slots are invalid")
    owners, charges, technical = {}, [], []
    for slot, profile in (("0",left),("1",right),(None,shared)):
        for owner in profile.ownership:
            key = (owner.lineage_id,owner.evidence_key)
            if key in owners:
                # Different scopes cannot be identical full charges. Never
                # average or silently discard a cross-member ownership clash.
                raise EvidenceFusionError("one atomic evidence item has multiple score owners or values")
            owners[key] = owner.charge_sha256
        charges.extend(ScoreChargeValue(c.kind,c.value,slot) for c in profile.nontechnical)
        technical.extend((lineage,key,ScoreChargeValue("technical",value,slot))
                         for lineage,key,value in profile.technical)
    charges.extend(charge for _,_,charge in sorted(technical,key=lambda row:row[:2]))
    return pair_score_numeric_inputs(tuple(charges), completeness=completeness, phase=phase,
        phenotype_units=phenotype_units, model_state=model_state,
        population_contribution=population_contribution, policy=policy)


def pair_score_numeric_inputs(charges, *, completeness, phase, phenotype_units, model_state,
                              population_contribution, policy):
    """Finish the numeric input after the caller's complete ownership closure."""
    charges = list(charges)
    mapping = next(fact for fact in policy.source_facts if fact.fact == "inheritance.model." + model_state)
    charges.append(ScoreChargeValue("inheritance",mapping.value,None))
    if population_contribution is not None:
        charges.append(ScoreChargeValue("population",population_contribution,None))
    return ScoreDecisionInputs(2,True,completeness,phase,phenotype_units,tuple(charges))


@member_validation_scope()
def member_score_charges(members, key, model_member, qualification, policy):
    """The pair template's actual source dependencies; no qualification labels."""
    allele = members.member(key).allele
    decisions, decision_bytes = [], 0
    for value in members.iter_decisions(subject_id=key):
        decision_bytes += len(_canonical(value))
        if decision_bytes > MAX_SOURCE_WINDOW_BYTES:
            raise EvidenceStoreError("member admission score proofs exceed capacity", capacity=True)
        decisions.append(value)
    alleles = {key: allele}
    charges = canonical_evidence_charges((
        *admission_decision_charges(tuple(decisions), {d.decision_id for d in decisions}, {key}, policy, alleles),
        *requirement_member_charges((qualification,), policy, alleles),
        *clinical_member_charges((qualification,), policy),
        *technical_member_charges((model_member,), policy, alleles),
    ))
    if len(_canonical(charges)) > MAX_SOURCE_WINDOW_BYTES:
        raise EvidenceStoreError("member score source evidence exceeds capacity", capacity=True)
    return charges


class CompletePairScoreSources:
    def __init__(self, writer, universe, *, policy, policy_sha256, gene_evidence, name, gene_resource=None):
        self.writer, self.connection, self.universe = writer, writer.connection, universe
        self.name, self.table = name, "pair_score_sources_" + hashlib.sha256(name.encode()).hexdigest()
        self.closed, self.created, self.member_count = False, False, 0
        self.tables = []
        self.policy = _typed(policy, FusionPolicy)
        self.policy_sha256 = policy_sha256
        if (hashlib.sha256(policy_as_json(self.policy)).hexdigest() != policy_sha256
                or universe.connection is not self.connection or universe.closed or universe.pairs is None):
            raise EvidenceStoreError("foreign pair score source owners")
        try:
            if gene_resource is None:
                prepared = (universe.prepared_gene_evidence if gene_evidence is universe.prepared_gene_evidence.evidence
                            else prepare_gene_evidence(gene_evidence))
                gene_resource = CaseGeneEvidence(writer,prepared,name=name + ".genes")
                gene_evidence = prepared.evidence
            gene_resource.require(writer,gene_evidence)
            self.gene_resource = gene_resource
            self.connection.execute(f"CREATE TABLE {self.table} (id TEXT PRIMARY KEY,charges TEXT,profile TEXT)")
            self.tables.append(self.table)
            self.connection.execute(f"CREATE TABLE {self.table}_owners (lineage TEXT,key TEXT,count INTEGER,PRIMARY KEY(lineage,key))")
            self.tables.append(self.table + "_owners")
            self.created = True
            shared, audits = disease_context_charges(universe.pairs.context, gene_evidence, self.policy,
                                                     prepared_gene_evidence=gene_resource.prepared)
            self.shared = canonical_evidence_charges(shared)
            self.shared_profile = project_member_charges(self.shared,None)
            self._index_owners(self.shared_profile)
            if len(_canonical((self.shared,audits,gene_resource.reference))) > MAX_SOURCE_WINDOW_BYTES:
                raise EvidenceStoreError("shared score source evidence exceeds capacity", capacity=True)
            self.shared_sha256 = writer.add(name + ".shared", "pair-score-shared-sources-v1",
                {"charges":self.shared, "profile":self.shared_profile, "audits":audits, "gene_evidence":gene_resource.reference})
            def rows():
                for key, model, qualification in universe.pairs.member_templates():
                    charges = member_score_charges(universe.members, key, model, qualification, self.policy)
                    profile = project_member_charges(charges, key)
                    self._index_owners(profile)
                    self.connection.execute(f"INSERT INTO {self.table} VALUES (?,?,?)",
                        (key, _canonical(charges), _canonical(profile)))
                    self.member_count += 1
                    yield {"allele_id":key, "charges":charges, "profile":profile}
            self.members_sha256 = writer.add_sequence(name + ".members", "pair-score-member-sources-v1", rows())
            if self.member_count != universe.pairs.summary.pairable_count:
                raise EvidenceStoreError("pair score member source coverage changed")
            # Private proof computed from every actual canonical source charge.
            # A false result requires per-pair closure, never global rejection:
            # some colliding members may represent the same unpairable event.
            self.all_owners_disjoint = self.connection.execute(
                f"SELECT 1 FROM {self.table}_owners WHERE count>1 LIMIT 1").fetchone() is None
            self.sha256 = writer.add(name, "pair-score-source-evidence-v1", {
                "schema":"genosapiens.wgs-pair-score-source-evidence", "version":1,
                "hypotheses":{"name":universe.name,"sha256":universe.sha256},
                "policy":self.policy, "policy_sha256":policy_sha256, "member_count":self.member_count,
                "members":{"name":name + ".members","sha256":self.members_sha256},
                "shared":{"name":name + ".shared","sha256":self.shared_sha256},
            })
        except BaseException:
            writer.failed = True
            self.close()
            raise

    def _index_owners(self, profile):
        self.connection.executemany(f"INSERT INTO {self.table}_owners VALUES (?,?,1) "
            "ON CONFLICT(lineage,key) DO UPDATE SET count=count+1",
            ((owner.lineage_id,owner.evidence_key) for owner in profile.ownership))

    def _require_open(self):
        if self.closed:
            raise EvidenceStoreError("pair score source owner is closed")

    def charges(self, allele_id):
        self._require_open()
        row = self.connection.execute(f"SELECT charges FROM {self.table} WHERE id=?", (allele_id,)).fetchone()
        if row is None:
            raise EvidenceStoreError("missing pair score source member")
        return _typed(json.loads(row[0]), tuple[EvidenceCharge, ...], decode=True)

    def profiles(self):
        self._require_open()
        for key, content in self.connection.execute(f"SELECT id,profile FROM {self.table} ORDER BY id"):
            yield key, _typed(json.loads(content), SourceScoreProfile, decode=True)

    def profile(self, allele_id):
        self._require_open()
        row = self.connection.execute(f"SELECT profile FROM {self.table} WHERE id=?", (allele_id,)).fetchone()
        if row is None:
            raise EvidenceStoreError("missing pair score source member")
        return _typed(json.loads(row[0]), SourceScoreProfile, decode=True)

    def close(self):
        for table in reversed(self.tables):
            self.connection.execute("DROP TABLE " + table)
        self.tables.clear()
        self.created = False
        self.closed = True

    def __enter__(self):
        self._require_open()
        return self

    def __exit__(self, *_):
        self.close()


class VerifiedPairScoreSources(CompletePairScoreSources):
    def __init__(self, reader, reference, universe, *, policy, policy_sha256, gene_evidence, gene_resource=None):
        if reference.kind != "pair-score-source-evidence-v1":
            raise EvidenceStoreError("foreign pair score source root")
        reader._require_reference(reference)
        super().__init__(EvidenceReplay(reader), universe, policy=policy, policy_sha256=policy_sha256,
                         gene_evidence=gene_evidence, name=reference.name,gene_resource=gene_resource)
