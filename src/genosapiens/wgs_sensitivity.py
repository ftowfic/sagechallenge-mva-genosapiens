"""Source-resource validation shared by normalized scientific consumers.

All persisted primary/sensitivity products use wgs_scientific_products.
"""
from genosapiens.disease_gene_evidence import disease_gene_evidence_as_bytes, rebuild_disease_gene_evidence
from genosapiens.disease_profiles import build_disease_profiles, disease_profiles_as_bytes
from genosapiens.evidence_model import DiseaseGeneEvidence, DiseaseProfile
from genosapiens.hpo_ic import HpoICStore
from genosapiens.track1_sensitivity import _authenticated_ic_hash, _ic_hash


class WgsSensitivityError(RuntimeError):
    """Static bounded primary/sensitivity failure."""



def _checked_ic_hash(store, authenticated_resources):
    """Reuse the loaded immutable store's proof when the caller owns one."""
    if authenticated_resources is not None:
        return _authenticated_ic_hash(authenticated_resources, store)
    return _ic_hash(store)



def _profiles(
    profile: DiseaseProfile | None,
    *,
    profile_evidence_id: str | None,
    ic_store: HpoICStore,
):
    if profile is None:
        if profile_evidence_id is not None:
            raise WgsSensitivityError("profile evidence binding is invalid")
        return None
    if (
        type(profile) is not DiseaseProfile
        or profile.evidence_id != profile_evidence_id
    ):
        raise WgsSensitivityError("profile evidence binding is invalid")
    rebuilt = build_disease_profiles(
        ic_store.graph,
        (value.source_fields for value in profile.store.annotations),
        source_binding=profile.store.source_binding,
    )
    if disease_profiles_as_bytes(rebuilt) != disease_profiles_as_bytes(
        profile.store
    ):
        raise WgsSensitivityError("phenotype profile is not source-bound")
    return profile.store



def _validate_gene_sources(values: tuple[DiseaseGeneEvidence, ...]) -> None:
    if type(values) is not tuple:
        raise WgsSensitivityError("disease-gene evidence is invalid")
    for value in values:
        if (
            type(value) is not DiseaseGeneEvidence
            or disease_gene_evidence_as_bytes(
                rebuild_disease_gene_evidence(value.store)
            ) != disease_gene_evidence_as_bytes(value.store)
        ):
            raise WgsSensitivityError("disease-gene evidence is not source-bound")
