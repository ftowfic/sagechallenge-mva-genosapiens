"""Bounded operational evidence for the finite public acceptance experiment.

The scientific manifests retain all scored hypotheses and complete contexts.
This companion retains admission outcomes outside those contexts and actual
annotation work, so acceptance metrics do not require a whole-callset replay.
It binds completed products but cannot authorize a scientific invocation.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from genosapiens.candidate_admission import AdmissionDecision, ShardedCandidateAdmission
from genosapiens.evidence_model import _sha
from genosapiens.wgs_failures import WgsFailure
from genosapiens.wgs_admission import LedgerAdmissionInput
from genosapiens.wgs_provenance import execution_sha256
from genosapiens.wgs_scientific import _bytes, _load_manifest, _read_rows, _ShardWriter, _write_manifest
from genosapiens.wgs_shards import ScientificShard
from genosapiens.wgs_scientific_artifact import ScientificContextAnnotation, ScientificArtifactReader


@dataclass(frozen=True)
class ValidationEligibilityWindow:
    first_ordinal: int
    last_ordinal: int
    policy_sha256: str
    decisions_sha256: str
    input_dispositions: tuple[tuple[str,int], ...]
    occurrence_dispositions: tuple[tuple[str,int], ...]
    schema: str = 'genosapiens.wgs-validation-eligibility-window'
    version: int = 1

    def __post_init__(self):
        from genosapiens.wgs_analysis_eligibility import DISPOSITIONS
        if (self.schema!='genosapiens.wgs-validation-eligibility-window' or type(self.version) is not int or self.version!=1
                or any(type(v) is not int or v<1 for v in (self.first_ordinal,self.last_ordinal))
                or self.last_ordinal<self.first_ordinal
                or any(tuple(k for k,v in rows)!=tuple(sorted(DISPOSITIONS))
                    or any(type(v) is not int or v<0 for k,v in rows)
                    for rows in (self.input_dispositions,self.occurrence_dispositions))
                or sum(v for k,v in self.input_dispositions)!=self.last_ordinal-self.first_ordinal+1):
            raise ValueError('invalid eligibility audit window')
        _sha(self.policy_sha256); _sha(self.decisions_sha256)

    @classmethod
    def from_batch(cls,batch):
        import hashlib
        from genosapiens.wgs_analysis_eligibility import DISPOSITIONS
        rows=batch.decisions
        return cls(rows[0].retained_ordinal,rows[-1].retained_ordinal,rows[0].eligibility_policy_sha256,
            hashlib.sha256(b''.join(r.as_bytes() for r in rows)).hexdigest(),
            tuple((key,sum(r.disposition==key for r in rows)) for key in sorted(DISPOSITIONS)),
            tuple((key,sum(len(r.source_allele_ids) for r in rows if r.disposition==key)) for key in sorted(DISPOSITIONS)))


@dataclass(frozen=True, slots=True, repr=False)
class ValidationAdmission:
    """The admission facts consumed by public validation readback."""

    inputs: tuple[LedgerAdmissionInput, ...]
    decisions: tuple[AdmissionDecision, ...]
    eligibility: ValidationEligibilityWindow | None = None

    def __post_init__(self):
        if (type(self.inputs) is not tuple or type(self.decisions) is not tuple
                or any(type(value) is not LedgerAdmissionInput for value in self.inputs)
                or any(type(value) is not AdmissionDecision for value in self.decisions)
                or self.eligibility is not None and (type(self.eligibility) is not ValidationEligibilityWindow
                    or not self.inputs or self.eligibility.first_ordinal!=self.inputs[0].retained_ordinal
                    or self.eligibility.last_ordinal!=self.inputs[-1].retained_ordinal)):
            raise WgsFailure("bounded_products", "malformed_evidence", operation="bounded_product_validate", exception_kind="project_failure")


@dataclass(frozen=True, slots=True, repr=False)
class ValidationAuditRow:
    ordinal: int
    kind: Literal["admission", "annotation"]
    admission: ValidationAdmission | None = None
    annotation: ScientificContextAnnotation | None = None

    def __post_init__(self):
        if (type(self.ordinal) is not int or self.ordinal < 1
                or self.kind not in {"admission", "annotation"}
                or self.kind == "admission" and (
                    type(self.admission) is not ValidationAdmission
                    or self.annotation is not None)
                or self.kind == "annotation" and (
                    self.admission is not None
                    or type(self.annotation) is not ScientificContextAnnotation)):
            raise WgsFailure("bounded_products", "malformed_evidence", operation="bounded_product_validate", exception_kind="project_failure")

    @property
    def key(self):
        return f"{self.ordinal:012d}:{self.kind}"


@dataclass(frozen=True, slots=True, repr=False)
class ValidationAuditManifest:
    parent_bindings: tuple[tuple[str, str], ...]
    shards: tuple[ScientificShard, ...]
    record_count: int
    input_count: int
    decision_count: int
    context_count: int
    stream_sha256: str
    schema: Literal["genosapiens.wgs-public-validation-audit"] = "genosapiens.wgs-public-validation-audit"
    analysis_eligibility: dict | None = None
    version: Literal[4] = 4
    classification: Literal["controlled"] = "controlled"
    llm_policy: Literal["no-llm"] = "no-llm"

    def __post_init__(self):
        if (self.schema != "genosapiens.wgs-public-validation-audit" or type(self.version) is not int
                or self.version != 4 or self.classification != "controlled" or self.llm_policy != "no-llm"
                or self.parent_bindings != tuple(sorted(self.parent_bindings))
                or set(dict(self.parent_bindings)) != ({"prepared_case", "policy_bindings", "execution_config", "scientific"}
                    | (set() if self.analysis_eligibility is None else {'analysis_eligibility_policy','analysis_eligibility_manifest','analysis_eligibility_rescue_manifest'}))
                or any(type(value) is not int or value < 0 for value in (
                    self.record_count, self.input_count, self.decision_count, self.context_count))
                or self.record_count != sum(row.row_count for row in self.shards)
                or any(row.kind != "validation_audit" for row in self.shards)):
            raise WgsFailure("bounded_products", "malformed_evidence", operation="bounded_product_validate", exception_kind="project_failure")
        for value in (self.stream_sha256, *(value for _, value in self.parent_bindings)):
            _sha(value)

    def as_bytes(self):
        content = _bytes({"classification": self.classification, "llm_policy": self.llm_policy,
                          "manifest": self})
        if len(content) > 8 * 1024 * 1024:
            raise WgsFailure("bounded_products", "response_capacity", operation="bounded_product_validate", exception_kind="project_failure")
        return content

    @property
    def sha256(self):
        import hashlib
        return hashlib.sha256(self.as_bytes()).hexdigest()


def load_validation_audit(content, *, expected_sha256):
    return _load_manifest(content, expected_sha256=expected_sha256,
                          result_type=ValidationAuditManifest, result_key="manifest")


class ValidationAuditWriter:
    def __init__(self, directory, *, execution_config, prepared_sha256, policy_sha256):
        self.root, self.config = Path(directory), execution_config
        self.parents = {"prepared_case": prepared_sha256, "policy_bindings": policy_sha256,
                        "execution_config": execution_sha256(execution_config)}
        self.writer = _ShardWriter(self.root, "validation_audit",
            max_rows=execution_config.source_shards.max_uncompressed_rows,
            max_bytes=execution_config.source_shards.max_uncompressed_bytes)
        self.count = self.inputs = self.decisions = self.contexts = 0
        self.finished = False
        self.eligibility_inputs=0
        self.eligibility_mode=None
        self.analysis_eligibility=None

    def _put(self, row):
        if self.finished:
            raise WgsFailure("bounded_products", "malformed_evidence", operation="bounded_product_validate", exception_kind="project_failure")
        self.writer.write(row, key=row.key)
        self.count += 1

    def admission(self, batch):
        if type(batch) is not ShardedCandidateAdmission:
            raise WgsFailure("bounded_products", "malformed_evidence", operation="bounded_product_validate", exception_kind="project_failure")
        self._put(ValidationAuditRow(self.count + 1, "admission",
            admission=ValidationAdmission(batch.inputs, batch.decisions)))
        self.inputs += len(batch.inputs)
        self.decisions += len(batch.decisions)

    def admission_projected(self, subjects, decisions, *, admission_reader):
        """Explicitly requested audits resolve original inputs with read counts."""
        import hashlib
        from genosapiens.wgs_analysis_eligibility import EligibilityBatch
        from genosapiens.wgs_ledger import canonical
        refs=admission_reader.refs_for_subjects(subjects,expected_rows_sha256=hashlib.sha256(
            canonical([value.admission_ref[4] for value in subjects])).hexdigest())
        records=admission_reader.get_records(refs,purpose='audit')
        eligibility=EligibilityBatch(tuple(admission_reader._compact(value.retained_ordinal).decision for value in subjects))
        self.admission_encoded(tuple(row.value for row in records),decisions,eligibility=eligibility)

    def admission_encoded(self, inputs, decisions, *, eligibility=None):
        """Keep the existing audit format from canonical worker-owned rows."""
        import json
        from genosapiens.evidence_model import _typed
        values = tuple(_typed(json.loads(content), AdmissionDecision, decode=True) for content in decisions)
        current=eligibility is not None
        if self.eligibility_mode is None: self.eligibility_mode=current
        if self.eligibility_mode!=current: raise ValueError('mixed eligibility audit domains')
        window=None if eligibility is None else ValidationEligibilityWindow.from_batch(eligibility)
        self._put(ValidationAuditRow(self.count + 1, 'admission',
            admission=ValidationAdmission(inputs, values,window)))
        self.inputs += len(inputs)
        if current: self.eligibility_inputs+=len(inputs)
        self.decisions += len(values)

    def annotations(self, scientific_reader):
        if type(scientific_reader) is not ScientificArtifactReader or self.contexts:
            raise WgsFailure("bounded_products", "malformed_evidence", operation="bounded_product_validate", exception_kind="project_failure")
        self.analysis_eligibility=scientific_reader.resources.analysis_eligibility
        if self.analysis_eligibility is not None:
            from genosapiens.wgs_analysis_eligibility import eligibility_component_parents
            if self.inputs!=self.eligibility_inputs or self.inputs!=self.analysis_eligibility['input_rows']:
                raise ValueError('incomplete eligibility audit coverage')
            self.parents.update(eligibility_component_parents(self.analysis_eligibility))
        elif self.eligibility_inputs:
            raise ValueError('missing sealed eligibility audit domain')
        for annotation in scientific_reader.annotation_contexts():
            self._put(ValidationAuditRow(self.count + 1, "annotation", annotation=annotation))
            self.contexts += 1

    def finish(self, scientific_sha256):
        if self.finished:
            raise WgsFailure("bounded_products", "malformed_evidence", operation="bounded_product_validate", exception_kind="project_failure")
        shards = self.writer.close()
        manifest = ValidationAuditManifest(tuple(sorted({**self.parents,
            "scientific": scientific_sha256}.items())), shards, self.count, self.inputs,
            self.decisions, self.contexts, self.writer.stream.hexdigest(),analysis_eligibility=self.analysis_eligibility)
        _write_manifest(self.root, "validation-audit-manifest.json", manifest.as_bytes())
        self.finished = True
        return manifest

    def close(self):
        self.writer._close()


def iter_validation_audit(manifest, *, directory, execution_config, expected_parents, scientific_reader):
    if (type(manifest) is not ValidationAuditManifest
            or type(scientific_reader) is not ScientificArtifactReader
            or dict(expected_parents).get("scientific") != scientific_reader.artifact.sha256
            or manifest.context_count != scientific_reader.artifact.positive_context_count
            or manifest.parent_bindings != tuple(sorted(expected_parents.items()))
            or dict(manifest.parent_bindings)["execution_config"] != execution_sha256(execution_config)
            or any(row.size_bytes > execution_config.source_shards.max_uncompressed_bytes
                   or row.row_count > execution_config.source_shards.max_uncompressed_rows for row in manifest.shards)):
        raise WgsFailure("bounded_products", "malformed_evidence", operation="bounded_product_validate", exception_kind="project_failure")
    count = inputs = decisions = contexts = 0
    if _bytes(manifest.analysis_eligibility)!=_bytes(scientific_reader.resources.analysis_eligibility):
        raise ValueError('validation eligibility component differs from science')
    eligibility_rows=None
    if manifest.analysis_eligibility is not None:
        from genosapiens.wgs_analysis_eligibility import verify_analysis_eligibility_component,iter_eligibility_decisions,EligibilityBatch
        bound=verify_analysis_eligibility_component(scientific_reader.directory.parent,manifest.analysis_eligibility,full=False,
            eligibility_root=scientific_reader.eligibility_root)
        eligibility_rows=iter_eligibility_decisions(scientific_reader.eligibility_root,bound,
            scratch_directory=scientific_reader.reader.scratch_directory)
    def rows():
        return _read_rows(Path(directory), manifest.shards, "validation_audit", ValidationAuditRow,
            expected_stream_sha256=manifest.stream_sha256, row_key=lambda value: value.key)

    from contextlib import closing
    from contextlib import ExitStack
    with ExitStack() as stack:
        annotations=stack.enter_context(closing(scientific_reader.annotation_contexts()))
        if eligibility_rows is not None: stack.callback(eligibility_rows.close)
        for _, _, row in rows():
            count += 1
            if row.ordinal != count:
                raise WgsFailure("bounded_products", "malformed_evidence", operation="bounded_product_validate", exception_kind="project_failure")
            if row.kind == "admission":
                if eligibility_rows is not None:
                    batch=EligibilityBatch(tuple(next(eligibility_rows) for _ in row.admission.inputs))
                    assessed=batch.require_inputs(row.admission.inputs)
                    if (row.admission.eligibility!=ValidationEligibilityWindow.from_batch(batch)
                            or any(d.subject_id not in {v.subject_id for v in assessed} for d in row.admission.decisions)):
                        raise ValueError('validation assessed decision coverage changed')
                elif row.admission.eligibility is not None:
                    raise ValueError('unbound validation eligibility window')
                inputs += len(row.admission.inputs)
                decisions += len(row.admission.decisions)
            else:
                if row.annotation != next(annotations, None):
                    raise WgsFailure("bounded_products", "malformed_evidence", operation="bounded_product_validate", exception_kind="project_failure")
                contexts += 1
        if next(annotations, None) is not None:
            raise WgsFailure("bounded_products", "malformed_evidence", operation="bounded_product_validate", exception_kind="project_failure")
        if eligibility_rows is not None and next(eligibility_rows,None) is not None:
            raise ValueError('validation eligibility input coverage incomplete')
    if (count, inputs, decisions, contexts) != (manifest.record_count, manifest.input_count,
                                               manifest.decision_count, manifest.context_count):
        raise WgsFailure("bounded_products", "malformed_evidence", operation="bounded_product_validate", exception_kind="project_failure")
    for _, _, row in rows():
        yield row
