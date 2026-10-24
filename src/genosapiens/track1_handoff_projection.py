"""Public declarations projected from authenticated bounded evidence records."""

from collections.abc import Mapping

from genosapiens.track1_handoff import public_evidence_reference_id, HandoffError
from genosapiens.wgs_provenance import PublicEvidenceRecord


def project_public_evidence_reference(record: PublicEvidenceRecord):
    """Preserve publisher-native provenance without inventing support lineage."""
    if type(record) is not PublicEvidenceRecord:
        raise HandoffError("persisted public evidence record is required")
    raw = record.source_fields
    ref = record.reference
    if record.owner_kind == "public":
        row = raw.get("record", raw) if isinstance(raw, Mapping) else {}
        source = record.record
        dataset, release, locator = (source["source_dataset"], source["source_release"],
                                     source["source_record_id"])
        kind = ("clinical" if dataset in {"ClinVar VCF", "Open Targets variant", "Open Targets evidence_eva"}
                else "population" if dataset in {"gnomAD joint sites", "gnomAD all-sites AN"}
                else "molecular")
    else:
        row = raw
        table = record.record.get("table")
        if record.owner_kind == "disease_gene" and table in {"native_repeat_locus", "native_repeat_disease"}:
            dataset, release, locator = (record.source_binding["source_dataset"],
                                        record.source_binding["source_release"], ref.record_id)
        else:
            dataset, release, locator = row["source_dataset"], row["source_release"], row["source_record_id"]
        kind = ("phenotype" if record.owner_kind == "disease_profile" else
                {"disease_alias": "disease", "disease_mapping": "disease",
                 "disease_prevalence": "prevalence", "native_repeat_locus": "molecular"}.get(table, "gene"))
    return {
        "evidence_id": public_evidence_reference_id(evidence_id=ref.evidence_id,
            source_record_id=ref.record_id, manifest_sha256=ref.public_manifest_sha256),
        "record_type": kind,
        "source_disease_id": row.get("source_disease_id"),
        "mapped_disease_id": row.get("mapped_disease_id"),
        "disease_specificity": row.get("disease_specificity"),
        "source_record_id": locator, "lineage_key": row.get("lineage_key"),
        "dependency_key": row.get("dependency_key"), "assertion_key": row.get("assertion_key"),
        "source_dataset": dataset, "source_release": release,
        "manifest_sha256": ref.public_manifest_sha256,
    }
