"""Complete post-primary diagnosis ordering with bounded context decisions.

The derivative retains a prior proof once per context and a complete locator
tape. All source evidence, primary scores and hypothesis identities remain in
the independently authenticated primary artifact. The external family merge
uses the same exact ordering kernel as production case ordering.
"""
from contextlib import ExitStack, contextmanager
from dataclasses import dataclass
from fractions import Fraction
import hashlib
from itertools import chain
from pathlib import Path

from genosapiens.evidence_model import HypothesisContext, _exact, _sha, _typed
from genosapiens.evidence_serialization import _canonical, _loads
from genosapiens.track1_sensitivity import DisclosedDiagnosis, prepare_diagnosis_match_context, prior_matches_context
from genosapiens.wgs_case_order import CompleteCaseOrder
from genosapiens.wgs_case_resources import bound_reference
from genosapiens.wgs_context_order import _identity_bytes
from genosapiens.wgs_evidence_replay import EvidenceReplay
from genosapiens.wgs_evidence_store import EvidenceManifest, EvidenceReader, EvidenceWriter, EvidenceStoreError
from genosapiens.wgs_scientific_artifact import ScientificArtifactReader

SCHEMA = "genosapiens.wgs-normalized-diagnosis-derivative"
KIND = "complete-diagnosis-derivative-v1"
MAX_DESCRIPTOR_BYTES = 8 * 1024**2
MAX_DECISION_BYTES = 64 * 1024**2


def context_identity(reader, row):
    header = reader.reader.reconstruct(bound_reference(reader.reader, row["context"], "scientific-case-context-v1"))
    members = reader.reader.reconstruct(bound_reference(reader.reader, header["members"], "context-evidence-v1"))
    context = _typed(members["seed"]["context"], HypothesisContext, decode=True)
    if context.context_id != row["context_id"]:
        raise EvidenceStoreError("diagnosis primary context identity changed")
    return context


@dataclass(frozen=True, slots=True, repr=False)
class DiagnosisDerivative:
    evidence: EvidenceManifest
    root_sha256: str
    scientific_manifest_sha256: str
    primary_manifest_sha256: str
    diagnosis: DisclosedDiagnosis
    context_count: int
    hypothesis_count: int
    matched_hypothesis_count: int
    schema: str = SCHEMA
    version: int = 1

    def __post_init__(self):
        for value in (self.root_sha256, self.scientific_manifest_sha256, self.primary_manifest_sha256):
            _sha(value)
        if (type(self.evidence) is not EvidenceManifest or type(self.diagnosis) is not DisclosedDiagnosis
                or self.schema != SCHEMA or type(self.version) is not int or self.version != 1
                or any(type(n) is not int or not 0 <= n < 2**63 for n in
                    (self.context_count, self.hypothesis_count, self.matched_hypothesis_count))
                or self.matched_hypothesis_count > self.hypothesis_count):
            raise EvidenceStoreError("invalid normalized diagnosis derivative")
        self.as_bytes()

    @property
    def status(self):
        return "applied" if self.matched_hypothesis_count else "not_applicable"

    @property
    def reason(self):
        return None if self.matched_hypothesis_count else "no_prior_match"

    def as_bytes(self):
        value = (_canonical({"classification": "controlled", "llm_policy": "no-llm", "derivative": self}) + "\n").encode("ascii")
        if len(value) > MAX_DESCRIPTOR_BYTES:
            raise EvidenceStoreError("diagnosis descriptor exceeds capacity", capacity=True)
        return value

    @property
    def sha256(self):
        return hashlib.sha256(self.as_bytes()).hexdigest()


def diagnosis_derivative_from_bytes(content, *, expected_sha256):
    _sha(expected_sha256)
    if type(content) is not bytes or len(content) > MAX_DESCRIPTOR_BYTES or hashlib.sha256(content).hexdigest() != expected_sha256:
        raise EvidenceStoreError("diagnosis descriptor binding differs")
    value = _loads(content)
    if (type(value) is not dict or set(value) != {"classification", "llm_policy", "derivative"}
            or value["classification"] != "controlled" or value["llm_policy"] != "no-llm"):
        raise EvidenceStoreError("invalid diagnosis descriptor envelope")
    result = _typed(value["derivative"], DiagnosisDerivative, decode=True)
    if result.as_bytes() != content:
        raise EvidenceStoreError("noncanonical diagnosis descriptor")
    return result


class CompleteDiagnosisDerivative(CompleteCaseOrder):
    def __init__(self, writer, primary, *, primary_product, diagnosis):
        if type(primary) is not ScientificArtifactReader:
            raise EvidenceStoreError("verified primary owner required for diagnosis derivative")
        primary._require_open()
        if primary_product.scientific_manifest_sha256 != primary.artifact.sha256:
            raise EvidenceStoreError("foreign diagnosis primary product")
        self.primary, self.primary_product, self.diagnosis = primary, primary_product, _typed(diagnosis, DisclosedDiagnosis)
        self.matched_hypothesis_count = 0
        resources = primary.resources
        super().__init__(writer, resources.fusion_policy, policy_sha256=resources.fusion_policy_sha256, name="diagnosis-order")
        try:
            table = self.prefix + "_priors"
            self.connection.execute(f"CREATE TABLE {table} (ordinal INTEGER PRIMARY KEY,contribution TEXT,content TEXT)")
            self.tables.append(table)
            prepared = prepare_diagnosis_match_context(disease_gene_evidence=resources.genes.prepared.evidence, diagnosis=self.diagnosis)
            strength = str(_exact(resources.prior_policy.strength))
            for row in primary.case.contexts():
                context = context_identity(primary, row)
                proofs = prior_matches_context(context, context.evidence_refs, resources.genes.prepared.evidence,
                    self.diagnosis, prepared=prepared)
                contribution = strength if proofs else "0"
                value = dict(ordinal=row["ordinal"], context_id=row["context_id"], context=row["context"],
                    proofs=proofs, contribution=contribution, prior_eligible=bool(proofs),
                    reason="source_native_or_equivalent" if proofs else "no_native_or_equivalent_match",
                    hypothesis_count=row["hypothesis_count"])
                content = _canonical(value)
                if len(content) > MAX_DECISION_BYTES:
                    raise EvidenceStoreError("diagnosis context proof exceeds capacity", capacity=True)
                self.connection.execute(f"INSERT INTO {table} VALUES (?,?,?)", (row["ordinal"], contribution, content))
                original = primary.case.orders["primary"].context(row["ordinal"])
                self.connection.execute(f"INSERT INTO {self.prefix}_contexts VALUES (?,?,?)",
                    (row["ordinal"], row["context_id"], _canonical(original)))
                self.context_count += 1
                self.matched_hypothesis_count += row["hypothesis_count"] if proofs else 0
            original = primary.case.orders["primary"]
            for family in original.families():
                # Context priors have precisely two possible values. Each
                # filtered run preserves the original family's lexical IDs.
                for contribution in sorted({"0", strength}, key=Fraction):
                    def rows():
                        for locator in original.tape.rows(start=family["start"], count=family["count"]):
                            saved, = self.connection.execute(f"SELECT contribution FROM {table} WHERE ordinal=?", (locator[0],)).fetchone()
                            if saved == contribution:
                                yield _identity_bytes(primary.case.identities.identity(locator)), *locator
                    values = rows()
                    first = next(values, None)
                    if first is None:
                        continue
                    total = str(Fraction(family["total"]) + Fraction(contribution))
                    family_id = self._family(total, tuple(family["ties"]))
                    _, count = self._write_run(family_id, chain((first,), values))
                    self.connection.execute(f"UPDATE {self.prefix}_families SET count=count+? WHERE id=?", (count, family_id))
                    self.hypothesis_count += count
            if self.hypothesis_count != primary.artifact.hypothesis_count or self.context_count != primary.artifact.context_count:
                raise EvidenceStoreError("diagnosis complete primary membership differs")
            self.publish()
            decisions_sha = writer.add_sequence("diagnosis-contexts", "diagnosis-derivative-contexts-v1",
                (self.decision(n) for n in range(1, self.context_count + 1)))
            self.root_sha256 = writer.add("diagnosis-derivative", KIND, dict(schema=SCHEMA, version=1,
                scientific_manifest_sha256=primary.artifact.sha256, primary_manifest_sha256=primary_product.sha256,
                primary_order=primary.case.orders["primary"].reference, resources=resources.reference,
                diagnosis=self.diagnosis, prior_policy=resources.prior_policy, order=self.reference,
                context_count=self.context_count, hypothesis_count=self.hypothesis_count,
                matched_hypothesis_count=self.matched_hypothesis_count,
                decisions={"name": "diagnosis-contexts", "sha256": decisions_sha}))
        except BaseException:
            writer.failed = True
            self.close()
            raise

    def decision(self, ordinal):
        if self.closed:
            raise EvidenceStoreError("diagnosis derivative is closed")
        row = self.connection.execute(f"SELECT content FROM {self.prefix}_priors WHERE ordinal=?", (ordinal,)).fetchone()
        if row is None:
            raise EvidenceStoreError("diagnosis context is unavailable")
        return _loads(row[0])

    def ordered(self):
        count = 0
        for family in self.families():
            for locator in self.tape.rows(start=family["start"], count=family["count"]):
                count += 1
                yield dict(rank=count, locator=locator, hypothesis_id=self.primary.case.identities.identity(locator),
                    total=family["total"], ties=tuple(family["ties"]))
        if count != self.hypothesis_count:
            raise EvidenceStoreError("diagnosis order coverage changed")


def write_diagnosis_derivative(reader, primary_product, diagnosis, *, output_directory):
    reader._require_open()
    root = Path(output_directory)
    root.mkdir(mode=0o700, parents=True, exist_ok=True)
    if any(root.iterdir()):
        raise EvidenceStoreError("diagnosis output directory must be empty")
    with reader.reader.read_scope(), EvidenceWriter(root, scratch_directory=root) as writer:
        with CompleteDiagnosisDerivative(writer, reader, primary_product=primary_product, diagnosis=diagnosis) as derivative:
            values = dict(root_sha256=derivative.root_sha256, scientific_manifest_sha256=reader.artifact.sha256,
                primary_manifest_sha256=primary_product.sha256, diagnosis=diagnosis,
                context_count=derivative.context_count, hypothesis_count=derivative.hypothesis_count,
                matched_hypothesis_count=derivative.matched_hypothesis_count)
        evidence = writer.finish()
    result = DiagnosisDerivative(evidence=evidence, **values)
    from genosapiens.wgs_ledger import publish_bytes
    publish_bytes(root, "diagnosis-derivative.json", result.as_bytes())
    return result


@contextmanager
def open_diagnosis_derivative(value, reader, primary_product, *, directory, scratch_directory):
    if (type(value) is not DiagnosisDerivative or value.scientific_manifest_sha256 != reader.artifact.sha256
            or value.primary_manifest_sha256 != primary_product.sha256):
        raise EvidenceStoreError("foreign normalized diagnosis derivative")
    with reader.reader.read_scope(), EvidenceReader(value.evidence, directory, scratch_directory=scratch_directory) as evidence:
        with evidence.read_scope(), evidence.complete_root_scope(), ExitStack() as stack:
            root = evidence.reference("diagnosis-derivative", kind=KIND)
            if root.sha256 != value.root_sha256:
                raise EvidenceStoreError("diagnosis derivative root differs")
            derivative = stack.enter_context(CompleteDiagnosisDerivative(EvidenceReplay(evidence), reader,
                primary_product=primary_product, diagnosis=value.diagnosis))
            if (derivative.root_sha256 != value.root_sha256 or any(getattr(derivative, name) != getattr(value, name)
                    for name in ("context_count", "hypothesis_count", "matched_hypothesis_count"))):
                raise EvidenceStoreError("diagnosis derivative metadata differs from complete replay")
            yield derivative
