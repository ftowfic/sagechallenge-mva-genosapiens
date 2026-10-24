"""Pure Broad str-analysis disease-locus admission, not a repeat caller or ranker.

Public API:
* load_pinned_repeat_catalog: authenticate supplied catalog/contract/manifest bytes.
* evaluate_repeat_catalog: replay that authentication and match one AlleleEvidence.
* RepeatReferenceInput: original reference contract/manifest/alias bytes and pins.
* load_synthetic_repeat_catalog: separate hash-bound validation fixture profile,
  never evidence from the Broad publisher or a controlled-use authorization.

No downloads, files, phenotype inputs, callbacks, thresholds, or gene-symbol joins.
Returned proofs retain native sizes/uncertainty and HypothesisContext alternatives.
Consumers must recompute proofs, not accept serialized outcomes as authorization.
The input allele's caller provenance is retained, not independently re-called here.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping
from dataclasses import dataclass
import hashlib
import re
from types import MappingProxyType
from typing import Literal
from urllib.parse import quote

from genosapiens.allele_identity import GenomicInterval, RepeatAllele
from genosapiens import longtr_admission as longtr
from genosapiens.disease_gene_evidence import (
    SCHEMA as DISEASE_SCHEMA, DiseaseGeneContext, DiseaseGeneEvidenceStore,
    GeneAssertion, GeneRequirements, SourceRecord, _bytes, _freeze, _hash, _moi,
    rebuild_disease_gene_evidence,
)
from genosapiens.evidence_model import (
    AlleleEvidence, AlleleReference, DiseaseGeneEvidence, EvidenceReference,
    HypothesisContext, LocusBinding, SourcedRecord, _Controlled, _sha, _typed,
)
from genosapiens.evidence_serialization import _loads
from genosapiens.manifest import PublicResourceManifest, SourceIdentity, Transformation
from genosapiens.reference_concordance import ReferenceIdentityContract, verify_reference_manifest_contract
from genosapiens.variant_class_resources import VariantClassBundleContract, _verify_manifest_contract


ADAPTER_VERSION = "broad-str-analysis-disease-locus"
CATALOG_COMMIT = "d41c85c97da3cd8c8ac9d7d47883570c95069874"
CATALOG_RESOURCE_ID = "broad-str-analysis-disease-catalog-grch38-d41c85c"
CATALOG_SHA256 = "6ad5c05b361899f9c6229aa841827cc00dedc7a9d5a9a35d8fac9e23710bd9cd"
CATALOG_SIZE_BYTES = 63431
CATALOG_RECORD_COUNT = 83
CATALOG_URL = (f"https://raw.githubusercontent.com/broadinstitute/str-analysis/{CATALOG_COMMIT}/"
               "str_analysis/variant_catalogs/variant_catalog_without_offtargets.GRCh38.json")
SOURCE_FAMILY = "broad-str-analysis"
SYNTHETIC_RESOURCE_ID = "genosapiens-synthetic-repeat-catalog"
SYNTHETIC_ADAPTER_VERSION = "genosapiens-synthetic-repeat-catalog"
SYNTHETIC_PROFILE_SCHEMA = "genosapiens.synthetic-repeat-catalog-profile"
IUPAC = MappingProxyType({
    "A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG",
    "W": "AT", "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT",
    "V": "ACG", "N": "ACGT",
})
POLICY = _freeze({
    "adapter": ADAPTER_VERSION,
    "locus": "exact-source-LocusId; no-symbol-or-prefix-join",
    "coordinates": "catalog-zero-based-half-open-to-native-one-based-closed: start+1,end",
    "region": "single-region-or-explicit-MainReferenceRegion-member; ordered-structure-token-selection",
    "structure": "IUPAC-literals-and-parenthesized-IUPAC-motifs-with-star-or-plus",
    "motif": "exact-or-positionwise-IUPAC-subset; partial-overlap-unresolved; no-rotation-or-reverse-complement",
    "motif_lists": "selected-structure-or-explicit-PathogenicMotifs-or-BenignMotifs-alternative-at-validated-main-region; association-only",
    "motif_uncertainty": "no-full-match-and-any-unresolved-comparison-is-unresolved; no-pooling-partial-matches; no-overriding-region-structure-RepeatUnit-conflicts",
    "motif_labels": "both-lists-support-context; overlapping-labels-retained-without-clinical-classification",
    "gene": "explicit-GeneId-ENSG-only; preserve-version; ENST-or-other-retained-unmapped",
    "disease": "explicit-Diseases-entry: OMIM-six-digits-or-source-namespaced-Symbol; no-Name-join",
    "inheritance": "disease-entry-only; AD/AR/XD/XR-to-shared-HPO; unknown-not-count",
    "size": "retain-all-values-bounds-units-and-read-support; no-thresholds-or-conversions",
})


class RepeatCatalogError(ValueError):
    """Static diagnostic, never source or case content."""


@dataclass(frozen=True, slots=True, repr=False)
class RepeatReferenceInput:
    """Hash-bound existing reference artifacts, not a caller compatibility flag.

    Alias bytes use genosapiens.contig-aliases. Only coordinate compatibility
    is established here, never sequence identity, accepted case REF, or rescue.
    """

    contract_bytes: bytes
    manifest_bytes: bytes
    aliases_bytes: bytes
    expected_contract_sha256: str
    expected_manifest_sha256: str


@dataclass(frozen=True, slots=True, repr=False)
class RepeatAssociationProof(_Controlled):
    """Separate structure/list comparisons and the source of association support.

    motif_interpretation remains the selected LocusStructure comparison, even
    when a listed alternative supplies admission. Reasons identify every matched
    source field; list labels never establish a clinical classification.
    """
    allele_ref: AlleleReference
    native_identity: RepeatAllele
    input_source: SourcedRecord
    reference_manifest_sha256: str
    adapter_version: str
    configuration_sha256: str
    evidence_refs: tuple[EvidenceReference, ...]
    association_status: Literal["match", "no_match", "unresolved"]
    outcome: Literal["admitted", "evaluated_no_match", "unresolved"]
    context: HypothesisContext | None
    gene_id: str | None
    locus_id: str | None
    selected_region: GenomicInterval | None
    selected_motif: str | None
    motif_interpretation: str
    pathogenic_motif_interpretation: str
    benign_motif_interpretation: str
    reasons: tuple[str, ...]
    reference_binding_sha256: str | None
    catalog_manifest_sha256: str
    catalog_contract_sha256: str
    catalog_evidence_sha256: str

    def as_bytes(self) -> bytes:
        """Deterministic controlled audit only; replay, never trust this outcome."""
        return _bytes({"schema": self.adapter_version, "classification": "controlled", "llm_policy": "no-llm", "proof": self})


def _json_bytes(content, expected_sha256=None):
    if type(content) is not bytes or not content or len(content) > 8 * 1024 * 1024:
        raise RepeatCatalogError("invalid bounded source bytes")
    if expected_sha256 is not None:
        _sha(expected_sha256)
        if hashlib.sha256(content).hexdigest() != expected_sha256:
            raise RepeatCatalogError("source checksum mismatch")
    return _loads(content.decode("utf-8", errors="strict"))


def _optional_text(row, key):
    value = row.get(key)
    if value is not None and type(value) is not str:
        raise RepeatCatalogError("invalid catalog field type")
    return value if value and value.strip() else None


def _catalog_evidence(raw, binding, source_family, evidence_id):
    """Shared source-shape interpretation after either closed loader authenticates bytes."""
    records, assertions, unassigned = {}, [], []
    views = defaultdict(list)
    for index, locus in enumerate(raw):
        if type(locus) is not dict or not _optional_text(locus, "LocusId"):
            raise RepeatCatalogError("invalid catalog locus")
        root_id = f"#/{index}"
        records[root_id] = SourceRecord(root_id, "native_repeat_locus", _freeze(locus))
        source_gene = _optional_text(locus, "GeneId")
        gene = source_gene if source_gene and re.fullmatch(r"ENSG[0-9]{11}(?:\.[1-9][0-9]*)?", source_gene) else None
        diagnostics = () if gene else ("unknown_gene_identity_no_display_join",)
        if source_family == SYNTHETIC_RESOURCE_ID:
            diagnostics += ("independently_synthetic_catalog_not_publisher_evidence",)
        if source_gene and gene is None:
            diagnostics += ("source_GeneId_not_an_Ensembl_gene_identifier",)
        diseases = locus.get("Diseases")
        if diseases is not None and type(diseases) is not list:
            raise RepeatCatalogError("invalid catalog disease list")
        for key in ("LocusStructure", "RepeatUnit", "MainReferenceRegion", "Inheritance"):
            _optional_text(locus, key)
        regions = locus.get("ReferenceRegion")
        if regions is not None and not (type(regions) is str or type(regions) is list
                                       and all(type(region) is str for region in regions)):
            raise RepeatCatalogError("invalid catalog reference regions")
        variant_types = locus.get("VariantType")
        if variant_types is not None and not (type(variant_types) is str or type(variant_types) is list
                                             and all(type(value) is str for value in variant_types)):
            raise RepeatCatalogError("invalid catalog variant types")
        for key in ("PathogenicMotifs", "BenignMotifs"):
            motifs = locus.get(key)
            if motifs is not None and not (type(motifs) is list and all(type(motif) is str for motif in motifs)):
                raise RepeatCatalogError("invalid catalog motif list")
        for number, disease in enumerate(diseases or ()):
            if type(disease) is not dict:
                raise RepeatCatalogError("invalid catalog disease assertion")
            record_id = f"{root_id}/Diseases/{number}"
            records[record_id] = SourceRecord(record_id, "native_repeat_disease", _freeze(disease))
            omim, symbol = _optional_text(disease, "OMIM"), _optional_text(disease, "Symbol")
            native = "OMIM:" + omim if omim and re.fullmatch(r"[1-9][0-9]{5}", omim) else (
                source_family + ":Symbol:" + quote(symbol, safe="") if symbol else None)
            if native is None:
                unassigned.append(record_id)
                continue
            flags = (*diagnostics, "catalog_association_not_curated_validity", "catalog_support_numerically_unscored")
            if omim and not native.startswith("OMIM:"):
                flags += ("unresolved_OMIM_mapping_source_symbol_retained",)
            inherited = _optional_text(disease, "Inheritance")
            moi, moi_flags = _moi(None, {"XR": "XLR", "XD": "XLD"}.get(inherited, inherited))
            requirement = GeneRequirements(
                moi, inherited, "mode_of_inheritance" if moi else "unknown", None,
                None, None, None, None, (root_id, record_id),
                (*moi_flags, "unknown_required_allele_count", "unknown_mechanism_direction"),
            )
            assertions.append(GeneAssertion(
                record_id, source_family, gene, "positive", None, "native_disease_locus",
                (requirement,), flags, locus["LocusId"], native, root_id,
            ))
            views[native].append(record_id)
        if not diseases:
            unassigned.append(root_id)
    assertions = tuple(sorted(assertions, key=lambda item: item.record_id))
    views = _freeze({key: tuple(sorted(ids)) for key, ids in sorted(views.items())})
    digest = _hash((DISEASE_SCHEMA, records, assertions, views, (), sorted(unassigned), binding))
    store = DiseaseGeneEvidenceStore(MappingProxyType(records), assertions, views, (), tuple(sorted(unassigned)), binding, digest)
    return DiseaseGeneEvidence(evidence_id, binding["bundle_manifest_sha256"], store)


def load_pinned_repeat_catalog(
    catalog_bytes: bytes, *, contract_bytes: bytes, manifest_bytes: bytes,
    expected_contract_sha256: str, expected_manifest_sha256: str,
    evidence_id: str = "repeat-catalog",
) -> DiseaseGeneEvidence:
    """Authenticate the selected immutable Broad catalog and interpret its real schema.

    manifest_bytes is its existing public-resource manifest, NOT a fabricated
    validation flag or the multi-resource bundle manifest. Its hash is the exact
    EvidenceReference.public_manifest_sha256. The contract hash is checked both
    against supplied bytes and the manifest transformation. No publication occurs.
    Original UTF-8 bytes are retained for source replay; raw JSON fields are also
    addressable by array JSON pointers, including every empty/unknown source value.
    """
    try:
        contract = VariantClassBundleContract.from_mapping(_json_bytes(contract_bytes, expected_contract_sha256))
        resource = contract.resource(CATALOG_RESOURCE_ID)
        manifest = PublicResourceManifest.from_mapping(_json_bytes(manifest_bytes, expected_manifest_sha256))
        _verify_manifest_contract(manifest, resource, expected_contract_sha256)
        artifact = resource.artifact("json")
        if (resource.release != "commit." + CATALOG_COMMIT or resource.lane != "repeat_expansion"
                or resource.validation["expected_record_count"] != CATALOG_RECORD_COUNT
                or artifact.url != CATALOG_URL or artifact.sha256 != CATALOG_SHA256
                or artifact.size_bytes != CATALOG_SIZE_BYTES or len(catalog_bytes) != CATALOG_SIZE_BYTES):
            raise RepeatCatalogError("catalog resource pin mismatch")
        raw = _json_bytes(catalog_bytes, CATALOG_SHA256)
        if type(raw) is not list or len(raw) != CATALOG_RECORD_COUNT:
            raise RepeatCatalogError("catalog record count mismatch")
        binding = _freeze({
            "classification": "public", "adapter_version": ADAPTER_VERSION,
            "bundle_manifest_sha256": expected_manifest_sha256,
            "contract_sha256": expected_contract_sha256, "source_sha256": CATALOG_SHA256,
            "source_dataset": CATALOG_RESOURCE_ID, "source_release": resource.release,
            "catalog_utf8": catalog_bytes.decode("utf-8"), "contract_utf8": contract_bytes.decode("utf-8"),
            "manifest_utf8": manifest_bytes.decode("utf-8"),
        })
        return _catalog_evidence(raw, binding, SOURCE_FAMILY, evidence_id)
    except (ValueError, RuntimeError, TypeError, KeyError, AttributeError, OverflowError, RecursionError):
        raise RepeatCatalogError("invalid or mismatched pinned repeat catalog") from None


def load_synthetic_repeat_catalog(
    catalog_bytes: bytes, *, profile_bytes: bytes, manifest_bytes: bytes,
    expected_profile_sha256: str, expected_manifest_sha256: str,
    evidence_id: str = "synthetic-repeat-catalog",
) -> DiseaseGeneEvidence:
    """Load an independently synthetic validation fixture without production overrides.

    The closed JSON profile contains schema=SYNTHETIC_PROFILE_SCHEMA, profile_id
    (synthetic-<identifier>), assembly=GRCh38, catalog_sha256, size_bytes, and
    record_count (zero allowed). Supply independently selected profile/manifest
    hashes, not hashes recovered from a purported evidence store.

    The existing PublicResourceManifest must use resource_id=SYNTHETIC_RESOURCE_ID,
    release=profile_id, a placeholder source URL
    https://example.invalid/genosapiens/synthetic-repeat-catalog/<profile_id>.json,
    and exactly releases/synthetic-repeat-catalog/<profile_id>/catalog.json as its
    output. Source/output hashes and size must equal the profile and actual bytes;
    transformation is SYNTHETIC_PROFILE_SCHEMA with the profile hash as parameters.

    Hash checks establish reproducibility, not independent origin, publisher
    evidence, empirical acceptance, or controlled-use readiness. Source family,
    Symbol namespace and adapter version stay explicitly synthetic in outputs.
    """
    try:
        profile = _json_bytes(profile_bytes, expected_profile_sha256)
        if (type(profile) is not dict or set(profile) != {
                "schema", "profile_id", "assembly", "catalog_sha256", "size_bytes", "record_count"}
                or profile["schema"] != SYNTHETIC_PROFILE_SCHEMA or profile["assembly"] != "GRCh38"
                or type(profile["profile_id"]) is not str
                or not re.fullmatch(r"synthetic-[a-z0-9][a-z0-9._-]{0,63}", profile["profile_id"])
                or type(profile["size_bytes"]) is not int or profile["size_bytes"] <= 0
                or type(profile["record_count"]) is not int or profile["record_count"] < 0):
            raise RepeatCatalogError("invalid synthetic catalog profile")
        raw = _json_bytes(catalog_bytes, profile["catalog_sha256"])
        manifest = PublicResourceManifest.from_mapping(_json_bytes(manifest_bytes, expected_manifest_sha256))
        profile_id = profile["profile_id"]
        if (len(catalog_bytes) != profile["size_bytes"] or type(raw) is not list or len(raw) != profile["record_count"]
                or manifest.resource_id != SYNTHETIC_RESOURCE_ID or manifest.release != profile_id
                or manifest.assembly != "GRCh38"
                or manifest.source != SourceIdentity(
                    f"https://example.invalid/genosapiens/synthetic-repeat-catalog/{profile_id}.json", profile["catalog_sha256"])
                or manifest.transformation != Transformation(SYNTHETIC_PROFILE_SCHEMA, expected_profile_sha256)
                or len(manifest.outputs) != 1
                or manifest.outputs[0].path != f"releases/synthetic-repeat-catalog/{profile_id}/catalog.json"
                or manifest.outputs[0].sha256 != profile["catalog_sha256"]
                or manifest.outputs[0].size_bytes != profile["size_bytes"]):
            raise RepeatCatalogError("synthetic catalog manifest or bytes mismatch")
        binding = _freeze({
            "classification": "public", "adapter_version": SYNTHETIC_ADAPTER_VERSION,
            "bundle_manifest_sha256": expected_manifest_sha256, "contract_sha256": expected_profile_sha256,
            "source_sha256": profile["catalog_sha256"], "source_dataset": SYNTHETIC_RESOURCE_ID,
            "source_release": profile_id, "catalog_utf8": catalog_bytes.decode("utf-8"),
            "contract_utf8": profile_bytes.decode("utf-8"), "manifest_utf8": manifest_bytes.decode("utf-8"),
        })
        return _catalog_evidence(raw, binding, SYNTHETIC_RESOURCE_ID, evidence_id)
    except (ValueError, RuntimeError, TypeError, KeyError, AttributeError, OverflowError, RecursionError):
        raise RepeatCatalogError("invalid or mismatched synthetic repeat catalog") from None


def _reference_interval(native, source):
    contract = ReferenceIdentityContract.from_mapping(_json_bytes(source.contract_bytes, source.expected_contract_sha256))
    manifest = PublicResourceManifest.from_mapping(_json_bytes(source.manifest_bytes, source.expected_manifest_sha256))
    verify_reference_manifest_contract(contract, manifest, contract_sha256=source.expected_contract_sha256)
    if (native.reference_manifest_sha256 != source.expected_manifest_sha256
            or contract.primary_assembly_accession != "GCA_000001405.15"):
        raise RepeatCatalogError("repeat reference binding mismatch")
    identity = native.identity
    role = contract.reference_role(identity.reference_role.value)
    if role.identity_id != identity.reference_identity_id:
        raise RepeatCatalogError("repeat reference role mismatch")
    output = next(output for output in manifest.outputs if output.path == contract.expected_output("contig_aliases").path)
    aliases = _json_bytes(source.aliases_bytes, output.sha256)
    if (len(source.aliases_bytes) != output.size_bytes or aliases["schema"] != "genosapiens.contig-aliases"
            or aliases["contract_id"] != contract.contract_id
            or aliases["sequence_count"] != contract.sequence_count
            or aliases["reference_roles"] != {role.role: {
                "identity_id": role.identity_id, "naming": role.naming,
                "sequence_relationship": role.sequence_relationship,
            } for role in contract.reference_roles}):
        raise RepeatCatalogError("repeat reference alias binding mismatch")
    forward, reverse = aliases["alignment_to_tool"], aliases["tool_to_alignment"]
    if (not isinstance(forward, dict) or not isinstance(reverse, dict)
            or len(set(forward.values())) != len(forward)
            or reverse != {value: key for key, value in forward.items()}):
        raise RepeatCatalogError("invalid reference coordinate aliases")
    interval = identity.interval
    alignment = interval.contig if identity.reference_role.value == "challenge_alignment" else reverse.get(interval.contig)
    if alignment not in forward:
        return None
    length = aliases["sequence_lengths"][alignment]
    if type(length) is not int or interval.end > length:
        raise RepeatCatalogError("repeat interval outside reference bounds")
    return GenomicInterval(forward[alignment], interval.start, interval.end)


def _selected_region(locus):
    regions = locus.get("ReferenceRegion")
    regions = (regions,) if type(regions) is str else regions
    if not regions:
        return None, None, "missing_reference_region"
    main = locus.get("MainReferenceRegion")
    if len(regions) > 1 and not main:
        return None, None, "missing_multi_region_selection"
    main = main or regions[0]
    if regions.count(main) != 1:
        return None, None, "ambiguous_or_nonmember_main_reference_region"
    selected = regions.index(main)
    region = re.fullmatch(r"([^:]+):([0-9]+)-([0-9]+)", main)
    if not region or int(region[3]) <= int(region[2]):
        return None, None, "unsupported_reference_region"
    converted = GenomicInterval(region[1], int(region[2]) + 1, int(region[3]))
    structure = locus.get("LocusStructure")
    if not structure or not re.fullmatch(r"(?:[ACGTRYSWKMBDHVN]+|\([ACGTRYSWKMBDHVN]+\)[*+])+", structure):
        return converted, None, "unsupported_locus_structure"
    motifs = re.findall(r"\(([ACGTRYSWKMBDHVN]+)\)[*+]", structure)
    types = locus.get("VariantType")
    types = (types,) if type(types) is str else types
    if len(motifs) != len(regions) or not types or len(types) != len(regions) or types[selected] != "Repeat":
        return converted, None, "unresolved_structure_region_correspondence"
    motif = motifs[selected]
    repeat_unit = locus.get("RepeatUnit")
    if not repeat_unit or _motif_match(repeat_unit, motif) not in ("exact", "iupac_subset"):
        return converted, motif, "unresolved_catalog_repeat_unit_structure"
    return converted, motif, None


def _motif_match(observed, pattern):
    if not observed or not pattern or any(char not in IUPAC for char in observed + pattern):
        return "unresolved"
    if observed == pattern:
        return "exact"
    if len(observed) != len(pattern) or any(not set(IUPAC[a]) & set(IUPAC[b]) for a, b in zip(observed, pattern)):
        return "different"
    if all(set(IUPAC[a]) <= set(IUPAC[b]) for a, b in zip(observed, pattern)):
        return "iupac_subset"
    return "unresolved_iupac_overlap"


def _motif_list(observed, patterns):
    if patterns is None:
        return "not_supplied"
    comparisons = tuple(_motif_match(observed, pattern) for pattern in patterns)
    if "exact" in comparisons:
        return "exact"
    if "iupac_subset" in comparisons:
        return "iupac_subset"
    return "unresolved" if any(value.startswith("unresolved") for value in comparisons) else "no_match"


def evaluate_repeat_catalog(
    allele: AlleleEvidence, catalog: DiseaseGeneEvidence, *, reference: RepeatReferenceInput | None,
) -> tuple[RepeatAssociationProof, ...]:
    """Replay verified catalog evidence and return every native disease alternative.

    Add catalog to EvidenceDocument.disease_gene_evidence; proofs' context and
    evidence_refs then resolve in the shared typed model. Empty Diseases with a
    genuine GeneId returns association_status=match and outcome=unresolved so the
    main admission evaluator can seek independent positive gene-disease support.
    Missing reference evidence is unresolved. Bad hashes or stale interpretations
    raise RepeatCatalogError, never a silently smaller candidate set.

    After binding the selected main region and its reference structure, accept
    the selected motif OR an explicitly listed PathogenicMotifs/BenignMotifs
    alternative. All use exact/positionwise IUPAC rules without rotations or
    reverse complements. The lists establish locus association, not diagnosis;
    conflicting region/structure bindings cannot be rescued by a motif-list hit.
    """
    try:
        allele = _typed(allele, AlleleEvidence)
        catalog = _typed(catalog, DiseaseGeneEvidence)
        binding = catalog.store.source_binding
        if (binding.get("adapter_version") not in (ADAPTER_VERSION, SYNTHETIC_ADAPTER_VERSION)
                or binding.get("bundle_manifest_sha256") != catalog.public_manifest_sha256
                or rebuild_disease_gene_evidence(catalog.store) != catalog.store):
            raise RepeatCatalogError("stale repeat catalog evidence")
        native = next(item for item in allele.identities if item.identity_id == allele.native_identity_id)
        identity = native.identity
        if type(identity) is not RepeatAllele:
            raise RepeatCatalogError("native repeat allele required")
        interval = _reference_interval(native, reference) if reference is not None else None
        reference_hash = _hash((reference.expected_contract_sha256, reference.expected_manifest_sha256,
                                hashlib.sha256(reference.aliases_bytes).hexdigest())) if reference is not None else None
        synthetic = binding["adapter_version"] == SYNTHETIC_ADAPTER_VERSION
        configuration_hash = _hash((POLICY, SYNTHETIC_ADAPTER_VERSION, binding["contract_sha256"], binding["source_sha256"])) if synthetic else _hash((POLICY, CATALOG_COMMIT, CATALOG_SHA256))
        matches = [record for record in catalog.store.records.values()
                   if record.table == "native_repeat_locus" and record.source_fields["LocusId"] == identity.locus_id]

        def proof(record, status, outcome, reasons, *, context=None, selected=None, motif=None, interpretation="not_evaluated"):
            locus = record.source_fields if record else {}
            source_gene = locus.get("GeneId")
            gene = source_gene if source_gene and re.fullmatch(r"ENSG[0-9]{11}(?:\.[1-9][0-9]*)?", source_gene) else None
            refs = context.evidence_refs if context else ((EvidenceReference(
                catalog.evidence_id, record.record_id, catalog.public_manifest_sha256),) if record else ())
            if synthetic and "independently_synthetic_catalog_not_publisher_evidence" not in reasons:
                reasons = (*reasons, "independently_synthetic_catalog_not_publisher_evidence")
            return RepeatAssociationProof(
                AlleleReference(allele.allele_id, native.identity_id), identity, native.source,
                native.reference_manifest_sha256, binding["adapter_version"], configuration_hash, refs, status, outcome,
                context, gene, locus.get("LocusId"), selected, motif, interpretation,
                _motif_list(identity.repeat_unit, locus.get("PathogenicMotifs")),
                _motif_list(identity.repeat_unit, locus.get("BenignMotifs")), tuple(reasons), reference_hash,
                catalog.public_manifest_sha256, binding["contract_sha256"], catalog.store.evidence_sha256,
            )

        if interval is None:
            return (proof(matches[0] if len(matches) == 1 else None, "unresolved", "unresolved",
                          ("missing_reference_evidence" if reference is None else "unresolved_reference_alias",)),)
        if not matches:
            if longtr.metadata(allele).get("locus_id_source") == "derived-source-region-key":
                return (proof(None, "unresolved", "unresolved", ("source-locus-name-unavailable; no-catalog-ID-inference",)),)
            return (proof(None, "no_match", "evaluated_no_match", ("no_exact_catalog_locus",)),)
        if len(matches) != 1:
            return tuple(proof(record, "unresolved", "unresolved", ("ambiguous_catalog_locus",)) for record in matches)
        record = matches[0]
        selected, motif, reason = _selected_region(record.source_fields)
        if selected is not None and selected != interval:
            return (proof(record, "no_match", "evaluated_no_match", ("different_selected_reference_region", *((reason,) if reason else ())),
                          selected=selected, motif=motif),)
        if reason:
            return (proof(record, "unresolved", "unresolved", (reason,), selected=selected, motif=motif),)
        interpretation = _motif_match(identity.repeat_unit, motif)
        motif_checks = (
            ("LocusStructure", interpretation),
            ("PathogenicMotifs", _motif_list(identity.repeat_unit, record.source_fields.get("PathogenicMotifs"))),
            ("BenignMotifs", _motif_list(identity.repeat_unit, record.source_fields.get("BenignMotifs"))),
        )
        matched_fields = tuple(field for field, result in motif_checks if result in ("exact", "iupac_subset"))
        if not matched_fields:
            unresolved = any(result.startswith("unresolved") for _, result in motif_checks)
            outcome = "unresolved" if unresolved else "evaluated_no_match"
            return (proof(record, "unresolved" if unresolved else "no_match", outcome,
                          ("unresolved_repeat_motif" if unresolved else "different_repeat_motif",),
                          selected=selected, motif=motif, interpretation=interpretation),)
        motif_reasons = tuple("motif_association_via_" + field for field in matched_fields)
        if "PathogenicMotifs" in matched_fields and "BenignMotifs" in matched_fields:
            motif_reasons += ("overlapping_catalog_motif_labels_not_clinical_classification",)
        result = []
        for assertion in catalog.store.assertions:
            if assertion.locus_record_id != record.record_id:
                continue
            req = assertion.requirements[0]
            target = DiseaseGeneContext(assertion.native_disease_id, assertion.gene_id, req.inheritance,
                                        assertion.native_disease_id, locus_id=assertion.locus_id)
            refs = tuple(EvidenceReference(catalog.evidence_id, key, catalog.public_manifest_sha256) for key in req.record_ids)
            context = HypothesisContext(
                "repeat-context:" + _hash((target, req, catalog.public_manifest_sha256)),
                target.disease_id, target.native_disease_id, target.gene_id, target.locus_id, target.inheritance,
                "source_native_repeat_locus", assertion.requirements, refs, (LocusBinding(refs[0], ("LocusId",)),),
            )
            result.append(proof(record, "match", "admitted", (
                "explicit_catalog_disease_locus_assertion", "association_only_not_pathogenicity",
                *motif_reasons, *assertion.diagnostics,
            ), context=context, selected=selected, motif=motif, interpretation=interpretation))
        unassigned = [key for key in catalog.store.unassigned_record_ids if key.startswith(record.record_id + "/Diseases/")]
        if not result or unassigned:
            result.append(proof(record, "match", "unresolved", (
                "unresolved_native_disease_assertion" if unassigned else "independent_disease_gene_support_required"
                if (record.source_fields.get("GeneId") or "").startswith("ENSG") else "no_explicit_disease_or_gene_support",
                *motif_reasons,
            ), selected=selected, motif=motif, interpretation=interpretation))
        return tuple(result)
    except (ValueError, RuntimeError, TypeError, KeyError, AttributeError, StopIteration, OverflowError, RecursionError):
        raise RepeatCatalogError("invalid repeat catalog admission evidence") from None
