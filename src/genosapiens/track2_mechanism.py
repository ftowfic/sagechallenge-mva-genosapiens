"""Controlled mechanism-readiness assessment from explicit evidence dimensions."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping


INPUT_FIELDS = {
    "allele_evidence",
    "conflicts",
    "direction",
    "exact_disease_validity",
    "genotype_completeness",
    "scenario_id",
    "inheritance_consistency",
    "joint_evidence",
    "phase",
    "technical_support",
    "transcript_relevance",
}
DIRECTIONS = {
    "dominant_negative",
    "gain_of_function",
    "hypomorphic",
    "loss_of_function",
    "other",
    "unknown",
}
MODEL_EFFECTS = DIRECTIONS
TRANSCRIPT_FATES = {
    "altered_splicing",
    "nonsense_mediated_decay",
    "stable_altered_transcript",
    "stable_transcript",
    "unknown",
}
NON_VEP_SUPPORT = {
    "cellular_phenotype_assay",
    "curated_functional_evidence",
    "natural_history_evidence",
    "protein_function_assay",
    "rescue_experiment",
    "transcript_assay",
}


class MechanismReadinessError(ValueError):
    """Raised when mechanism confidence is based on incomplete or forbidden inputs."""


@dataclass(frozen=True, slots=True)
class AlleleMechanismEvidence:
    allele_id: str
    predicted_consequence: str | None
    curated_consequence: str | None
    curation_status: str
    evidence_ids: tuple[str, ...]

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any]) -> "AlleleMechanismEvidence":
        required = {
            "allele_id",
            "predicted_consequence",
            "curated_consequence",
            "curation_status",
            "evidence_ids",
        }
        if set(value) != required:
            raise MechanismReadinessError("allele mechanism evidence fields are incomplete")
        if not isinstance(value.get("evidence_ids"), list):
            raise MechanismReadinessError("allele mechanism evidence IDs must be a list")
        record = cls(
            allele_id=value["allele_id"],
            predicted_consequence=value["predicted_consequence"],
            curated_consequence=value["curated_consequence"],
            curation_status=value["curation_status"],
            evidence_ids=tuple(value["evidence_ids"]),
        )
        if (
            not isinstance(record.allele_id, str)
            or not record.allele_id
            or record.curation_status not in {"curated", "predicted", "unknown"}
            or not record.evidence_ids
            or any(not isinstance(item, str) or not item for item in record.evidence_ids)
            or (
                record.curation_status == "curated"
                and not record.curated_consequence
            )
            or (
                record.curation_status == "predicted"
                and not record.predicted_consequence
            )
        ):
            raise MechanismReadinessError("allele mechanism evidence is invalid")
        return record


@dataclass(frozen=True, slots=True)
class MechanismReadiness:
    scenario_id: str
    genotype_completeness: str
    phase: str
    technical_support: str
    exact_disease_validity: str
    inheritance_consistency: str
    transcript_relevance: str
    allele_evidence: tuple[AlleleMechanismEvidence, ...]
    joint_evidence: str
    conflicts: tuple[str, ...]
    direction: str
    missing_direction: bool
    readiness_state: str
    mechanism_confidence: str
    confidence_inputs: tuple[str, ...]
    prediction_acceptance: dict | None = None

    @property
    def permits_planning(self):
        return self.readiness_state in {'ready', 'prediction_accepted'}

    def as_track2_payload(self) -> dict[str, Any]:
        return {
            "conflicts": list(self.conflicts),
            "genotype_completeness": self.genotype_completeness,
            "missing_direction": self.missing_direction,
            "phase": self.phase,
            "readiness_state": self.readiness_state,
            "subtype_applicability": self.exact_disease_validity,
            "technical_support": self.technical_support,
            "transcript_relevance": self.transcript_relevance,
        }


@dataclass(frozen=True, slots=True)
class AlleleMechanismModel:
    allele_id: str
    transcript_id: str
    vep_consequence: str
    vep_severity: str
    transcript_fate: str
    residual_product: str
    molecular_function: str
    evidence_support_types: tuple[str, ...]
    evidence_ids: tuple[str, ...]

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any]) -> "AlleleMechanismModel":
        required = {
            "allele_id",
            "evidence_ids",
            "evidence_support_types",
            "molecular_function",
            "residual_product",
            "transcript_fate",
            "transcript_id",
            "vep_consequence",
            "vep_severity",
        }
        if set(value) != required:
            raise MechanismReadinessError("allele mechanism model fields are incomplete")
        record = cls(
            allele_id=value["allele_id"],
            transcript_id=value["transcript_id"],
            vep_consequence=value["vep_consequence"],
            vep_severity=value["vep_severity"],
            transcript_fate=value["transcript_fate"],
            residual_product=value["residual_product"],
            molecular_function=value["molecular_function"],
            evidence_support_types=tuple(value["evidence_support_types"]),
            evidence_ids=tuple(value["evidence_ids"]),
        )
        if (
            not record.allele_id
            or not record.transcript_id
            or not record.vep_consequence
            or not record.vep_severity
            or record.transcript_fate not in TRANSCRIPT_FATES
            or not record.residual_product
            or record.molecular_function not in MODEL_EFFECTS
            or not record.evidence_support_types
            or not record.evidence_ids
        ):
            raise MechanismReadinessError("allele mechanism model is incomplete")
        support = set(record.evidence_support_types)
        non_vep = support & NON_VEP_SUPPORT
        predicted = support == {'policy_prediction'} and record.molecular_function == 'loss_of_function'
        if record.molecular_function != "unknown" and not non_vep and not predicted:
            raise MechanismReadinessError(
                "VEP severity alone cannot infer molecular function"
            )
        if record.transcript_fate != "unknown" and not (
            support & {"transcript_assay", "curated_functional_evidence"}
        ):
            raise MechanismReadinessError(
                "VEP severity alone cannot infer transcript fate or NMD"
            )
        if record.residual_product != "unknown" and not (
            support
            & {
                "protein_function_assay",
                "transcript_assay",
                "curated_functional_evidence",
            }
        ):
            raise MechanismReadinessError(
                "VEP severity alone cannot infer residual product"
            )
        return record


@dataclass(frozen=True, slots=True)
class JointMechanismModel:
    allele_ids: tuple[str, ...]
    joint_effect: str
    cellular_lesion_class: str
    discriminating_baseline_lesion: str
    cell_process: str
    tissue: str
    developmental_stage: str
    disease_stage: str
    natural_history: str
    temporal_phase: str
    response_shape: str
    feedback: str
    evidence_support_types: tuple[str, ...]
    evidence_ids: tuple[str, ...]

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any]) -> "JointMechanismModel":
        required = {
            "allele_ids",
            "cell_process",
            "cellular_lesion_class",
            "developmental_stage",
            "discriminating_baseline_lesion",
            "disease_stage",
            "evidence_ids",
            "evidence_support_types",
            "feedback",
            "joint_effect",
            "natural_history",
            "response_shape",
            "temporal_phase",
            "tissue",
        }
        if set(value) != required:
            raise MechanismReadinessError("joint mechanism model fields are incomplete")
        record = cls(
            allele_ids=tuple(value["allele_ids"]),
            joint_effect=value["joint_effect"],
            cellular_lesion_class=value["cellular_lesion_class"],
            discriminating_baseline_lesion=value["discriminating_baseline_lesion"],
            cell_process=value["cell_process"],
            tissue=value["tissue"],
            developmental_stage=value["developmental_stage"],
            disease_stage=value["disease_stage"],
            natural_history=value["natural_history"],
            temporal_phase=value["temporal_phase"],
            response_shape=value["response_shape"],
            feedback=value["feedback"],
            evidence_support_types=tuple(value["evidence_support_types"]),
            evidence_ids=tuple(value["evidence_ids"]),
        )
        text_fields = (
            record.cellular_lesion_class,
            record.discriminating_baseline_lesion,
            record.cell_process,
            record.tissue,
            record.developmental_stage,
            record.disease_stage,
            record.natural_history,
            record.temporal_phase,
            record.response_shape,
            record.feedback,
        )
        if (
            not record.allele_ids
            or len(record.allele_ids) != len(set(record.allele_ids))
            or record.joint_effect not in MODEL_EFFECTS
            or any(not item for item in text_fields)
            or not record.evidence_support_types
            or not record.evidence_ids
        ):
            raise MechanismReadinessError("joint mechanism model is incomplete")
        support = set(record.evidence_support_types)
        non_vep = support & NON_VEP_SUPPORT
        predicted = support == {'policy_prediction'} and record.joint_effect == 'loss_of_function'
        if record.joint_effect != "unknown" and not non_vep and not predicted:
            raise MechanismReadinessError("VEP severity alone cannot infer a joint effect")
        if (
            record.cellular_lesion_class != "unknown"
            or record.discriminating_baseline_lesion != "unknown"
        ) and not (support & {"cellular_phenotype_assay", "curated_functional_evidence"}):
            raise MechanismReadinessError(
                "VEP severity alone cannot infer the model lesion"
            )
        return record


@dataclass(frozen=True, slots=True)
class MechanismModel:
    scenario_id: str
    allele_models: tuple[AlleleMechanismModel, ...]
    joint_model: JointMechanismModel
    explicit_unknowns: tuple[str, ...]

    def as_track2_payload(self) -> dict[str, Any]:
        return {
            "allele_models": [
                {
                    "allele_id": item.allele_id,
                    "transcript_fate": item.transcript_fate,
                    "residual_product": item.residual_product,
                    "molecular_function": item.molecular_function,
                }
                for item in self.allele_models
            ],
            "explicit_unknowns": list(self.explicit_unknowns),
            "joint_model": self.joint_model.joint_effect,
            "lesion": self.joint_model.cellular_lesion_class,
            "natural_history": self.joint_model.natural_history,
            "stage": {
                "developmental": self.joint_model.developmental_stage,
                "disease": self.joint_model.disease_stage,
                "temporal_phase": self.joint_model.temporal_phase,
            },
            "tissue": self.joint_model.tissue,
        }


def build_mechanism_model(data: Mapping[str, Any]) -> MechanismModel:
    if set(data) != {"allele_models", "explicit_unknowns", "scenario_id", "joint_model"}:
        raise MechanismReadinessError("mechanism model fields are incomplete or unknown")
    if not isinstance(data.get("allele_models"), list) or not all(
        isinstance(item, Mapping) for item in data["allele_models"]
    ):
        raise MechanismReadinessError("allele models must be an explicit list")
    alleles = tuple(AlleleMechanismModel.from_mapping(item) for item in data["allele_models"])
    if not isinstance(data.get("joint_model"), Mapping):
        raise MechanismReadinessError("joint mechanism model is required")
    joint = JointMechanismModel.from_mapping(data["joint_model"])
    if set(joint.allele_ids) != {item.allele_id for item in alleles}:
        raise MechanismReadinessError("joint model does not resolve every allele")
    unknowns = data.get("explicit_unknowns")
    if not isinstance(unknowns, list) or any(
        not isinstance(item, str) or not item for item in unknowns
    ):
        raise MechanismReadinessError("mechanism model unknowns must be explicit")
    unknown_values = {
        "joint_model.joint_effect": joint.joint_effect,
        "joint_model.cellular_lesion_class": joint.cellular_lesion_class,
        "joint_model.discriminating_baseline_lesion": joint.discriminating_baseline_lesion,
        "joint_model.cell_process": joint.cell_process,
        "joint_model.tissue": joint.tissue,
        "joint_model.developmental_stage": joint.developmental_stage,
        "joint_model.disease_stage": joint.disease_stage,
        "joint_model.natural_history": joint.natural_history,
        "joint_model.temporal_phase": joint.temporal_phase,
        "joint_model.response_shape": joint.response_shape,
        "joint_model.feedback": joint.feedback,
    }
    required_unknowns = {name for name, value in unknown_values.items() if value == "unknown"}
    for index, allele in enumerate(alleles):
        for field in ("transcript_fate", "residual_product", "molecular_function"):
            if getattr(allele, field) == "unknown":
                required_unknowns.add(f"allele_models[{index}].{field}")
    if not required_unknowns.issubset(set(unknowns)):
        raise MechanismReadinessError("mechanism model omits explicit unknowns")
    if not isinstance(data.get("scenario_id"), str) or not data["scenario_id"]:
        raise MechanismReadinessError("mechanism model hypothesis identity is required")
    return MechanismModel(
        scenario_id=data["scenario_id"],
        allele_models=alleles,
        joint_model=joint,
        explicit_unknowns=tuple(sorted(set(unknowns))),
    )


def assess_mechanism_readiness(data: Mapping[str, Any]) -> MechanismReadiness:
    if set(data) - {'prediction_acceptance'} != INPUT_FIELDS:
        forbidden = sorted(set(data) & {"epcr", "track1_score", "ranking_score"})
        if forbidden:
            raise MechanismReadinessError(
                f"{forbidden[0]} cannot populate mechanism confidence"
            )
        raise MechanismReadinessError("mechanism-readiness inputs are incomplete or unknown")
    allele_data = data.get("allele_evidence")
    if not isinstance(allele_data, list) or not allele_data or not all(
        isinstance(item, Mapping) for item in allele_data
    ):
        raise MechanismReadinessError("allele-level mechanism evidence is required")
    alleles = tuple(AlleleMechanismEvidence.from_mapping(item) for item in allele_data)
    if len({item.allele_id for item in alleles}) != len(alleles):
        raise MechanismReadinessError("allele mechanism evidence IDs must be unique")
    conflicts = data.get("conflicts")
    if not isinstance(conflicts, list) or any(
        not isinstance(item, str) or not item for item in conflicts
    ):
        raise MechanismReadinessError("mechanism conflicts must be an explicit list")
    if not isinstance(data.get("scenario_id"), str) or not data["scenario_id"]:
        raise MechanismReadinessError("mechanism hypothesis identity is required")
    allowed = {
        "genotype_completeness": {"complete", "incomplete", "unknown"},
        "phase": {"confirmed_cis", "confirmed_trans", "not_applicable", "unknown"},
        "technical_support": {"supported", "uncertain", "unsupported"},
        "exact_disease_validity": {"broad_only", "conflicting", "exact_supported", "unknown"},
        "inheritance_consistency": {"consistent", "inconsistent", "unknown"},
        "transcript_relevance": {"canonical_relevant", "irrelevant", "relevant", "unknown"},
        "joint_evidence": {"conflicting", "supported", "unknown", "prediction_accepted"},
        "direction": DIRECTIONS,
    }
    for field, values in allowed.items():
        if data.get(field) not in values:
            raise MechanismReadinessError(f"invalid mechanism-readiness field: {field}")
    blocking = (
        data["technical_support"] == "unsupported"
        or data["exact_disease_validity"] == "conflicting"
        or data["inheritance_consistency"] == "inconsistent"
        or data["transcript_relevance"] == "irrelevant"
        or data["joint_evidence"] == "conflicting"
        or bool(conflicts)
    )
    complete = (
        data["genotype_completeness"] == "complete"
        and data["phase"] != "unknown"
        and data["technical_support"] == "supported"
        and data["exact_disease_validity"] == "exact_supported"
        and data["inheritance_consistency"] == "consistent"
        and data["transcript_relevance"] in {"canonical_relevant", "relevant"}
        and all(item.curation_status == "curated" for item in alleles)
        and data["joint_evidence"] == "supported"
        and data["direction"] != "unknown"
    )
    prediction=data.get('prediction_acceptance')
    prediction_complete=False
    if prediction is not None:
        from genosapiens.track2_adjudication import validate_automatic_decision, AdjudicationError
        if not isinstance(prediction,dict) or set(prediction)!={'decision','evidence'}:
            raise MechanismReadinessError('Prediction acceptance has incomplete provenance.')
        decision=prediction['decision'];evidence=prediction['evidence']
        try:
            validate_automatic_decision(decision,domain='mechanism',
                subject={k:v for k,v in data.items() if k!='prediction_acceptance'},
                evidence=evidence,checks=evidence['checks'],execution=decision['execution'])
        except (AdjudicationError, KeyError, TypeError) as error:
            raise MechanismReadinessError('Prediction acceptance cannot replay its policy.') from error
        prediction_complete=(decision['outcome']=='accepted'
            and data['joint_evidence']=='prediction_accepted' and data['direction']=='loss_of_function'
            and data['genotype_completeness']=='complete' and data['phase']!='unknown'
            and data['technical_support']=='supported' and data['exact_disease_validity']=='exact_supported'
            and data['inheritance_consistency']=='consistent'
            and data['transcript_relevance'] in {'canonical_relevant','relevant'}
            and all(a.curation_status=='predicted' and a.curated_consequence is None for a in alleles))
        if decision['outcome']=='accepted' and not prediction_complete:
            raise MechanismReadinessError('Prediction acceptance conflicts with its mechanism facts.')
    if data['joint_evidence']=='prediction_accepted' and not prediction_complete:
        raise MechanismReadinessError('Predicted joint evidence requires automatic acceptance.')
    if blocking:
        readiness_state = "not_ready"
        confidence = "low"
    elif complete:
        readiness_state = "ready"
        confidence = "high"
    elif prediction_complete:
        readiness_state = 'prediction_accepted'
        confidence = 'policy_prediction'
    else:
        readiness_state = "not_evaluable"
        confidence = "not_evaluable"
    return MechanismReadiness(
        scenario_id=data["scenario_id"],
        genotype_completeness=data["genotype_completeness"],
        phase=data["phase"],
        technical_support=data["technical_support"],
        exact_disease_validity=data["exact_disease_validity"],
        inheritance_consistency=data["inheritance_consistency"],
        transcript_relevance=data["transcript_relevance"],
        allele_evidence=alleles,
        joint_evidence=data["joint_evidence"],
        conflicts=tuple(conflicts),
        direction=data["direction"],
        missing_direction=data["direction"] == "unknown",
        readiness_state=readiness_state,
        mechanism_confidence=confidence,
        confidence_inputs=(
            "genotype_completeness",
            "phase",
            "technical_support",
            "exact_disease_validity",
            "inheritance_consistency",
            "transcript_relevance",
            "allele_evidence",
            "joint_evidence",
            "conflicts",
            "direction",
        ),
        prediction_acceptance=prediction,
    )


def readiness_inputs(readiness):
    from dataclasses import asdict
    data={k:v for k,v in asdict(readiness).items() if k in INPUT_FIELDS}
    if readiness.prediction_acceptance is not None:
        data['prediction_acceptance']=readiness.prediction_acceptance
    return data
