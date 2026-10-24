"""Lossless v6 gene evidence plus task 7.6 source-native disease-locus support.

No genotype completeness, VEP interpretation, candidate generation or release
freeze is performed here. Numerical values are explicitly development policy.
Mapped disease views are indexes, not equivalence assertions. Requirements are
source-specific alternatives, never a union of dominant/recessive mechanisms.

Vocabulary references (consulted 2026-09-06; not runtime dependencies):
* https://thegencc.org/faq/ (validity definitions)
* https://github.com/GeneCurationCoalition/gencc-search/blob/main/app/Classification.php
* https://github.com/GeneCurationCoalition/gencc-sub/blob/main/app/Models/Inheritance.php
* https://sciences.orphadata.com/docs/OrphadataFreeAccessProductsDescription.pdf
  v2, June 2025, pp. 12-15, 26, 30-32 (association types, inheritance, mappings)
* https://www.orphacode.org/data/docs/OrphanetNomenclaturePackDescription.pdf
  v6, July 2026, pp. 3-7 (active subtype AggregationLevel, NOT replacement links)
* https://github.com/opentargets/evidence_datasource_parsers/tree/master/modules
  Orphanet.py, Gene2Phenotype.py, ClinGen.py, PanelApp.py (native field semantics)
* https://github.com/EBI-G2P/gene2phenotype_api/blob/main/gene2phenotype_project/
  gene2phenotype_app/fixtures/attribs.json (literal allelic vocabulary)

An Orphanet association type can explicitly assert a mechanism. G2P product
level/structure consequences do NOT establish mechanism direction or NMD fate.
GenCC and ClinGen MOI values, including those stored in allelic_requirements,
are NOT literal allele counts. Unsupported vocabularies remain unknown.
"""

from __future__ import annotations

from collections import defaultdict, deque
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, fields, is_dataclass, replace
from datetime import date
from graphlib import CycleError, TopologicalSorter
import hashlib
import json
import math
from pathlib import Path
import re
from types import MappingProxyType
from urllib.parse import quote, unquote

from genosapiens.bounded_json import BoundedJSONError, load_bounded_json
from genosapiens.ontology_bundle import OntologyBundleError, _file_identity, _file_signature, _require_within_root
from genosapiens.rare_disease import (
    BUNDLE_MANIFEST_SCHEMA, OT_EVIDENCE_DATASETS, RareDiseaseContract,
    RareDiseaseError, TABLE_COLUMNS, TABLE_SCHEMAS,
    _artifact_contract_sha256, _bundle_manifest_path, _source_manifest_path,
)


SCHEMA = "genosapiens.disease-gene-evidence"
# v6 publication tables are unchanged. Native catalog records are not v6 rows.
NATIVE_LOCUS_TABLES = ("native_repeat_locus", "native_repeat_disease")
PUBLIC_TABLES = ("gene_assertion", "disease_inheritance", "disease", "disease_alias", "disease_relation")
SOURCE_COLUMNS = TABLE_COLUMNS["gene_assertion"]
SOURCE_NATIVE_DISEASE_NAMESPACES = (
    "broad-str-analysis:Symbol:", "genosapiens-synthetic-repeat-catalog:Symbol:",
)
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_GENCC_VALIDITY = {
    "GENCC:100001": "definitive", "GENCC:100002": "strong", "GENCC:100003": "moderate",
    "GENCC:100004": "limited", "GENCC:100005": "disputed", "GENCC:100006": "refuted",
    "GENCC:100007": "animal model only", "GENCC:100008": "no known disease relationship",
    "GENCC:100009": "supportive",
}
_MOI = {
    "HP:0000006": ("autosomal dominant", "AD"),
    "HP:0000007": ("autosomal recessive", "AR"),
    "HP:0001417": ("x-linked", "XL"),
    "HP:0001419": ("x-linked recessive", "XLR"),
    "HP:0001423": ("x-linked dominant", "XLD"),
    "HP:0001427": ("mitochondrial", "MIT"),
    "HP:0001450": ("y-linked", "YLD"),
    "HP:0032113": ("semidominant", "SD", "semi-dominant"),
    "HP:0010985": ("gonosomal", "GS"),
    "HP:0001426": ("multifactorial", "MF", "multigenic/multifactorial"),
    "HP:0032382": ("uniparental disomy", "UNP"),
    "HP:0001428": ("somatic mutation", "SOM"),
    "HP:0001466": ("contiguous gene syndrome", "CGS"),
    "HP:0003743": ("genetic anticipation", "GA"),
    "HP:0003745": ("sporadic", "SPR"),
    "HP:0001442": ("somatic mosaicism", "SM"),
    "HP:0012274": ("autosomal dominant inheritance with paternal imprinting", "ADIPI"),
    "HP:0012275": ("autosomal dominant inheritance with maternal imprinting", "ADIMI"),
    "HP:0010984": ("digenic", "DI"),
}
_MOI_TITLES = {
    title.casefold(): curie for curie, titles in _MOI.items()
    for title in (*titles, *(name + " inheritance" for name in titles[:1]))
}
_G2P_ALLELES = {
    "monoallelic_autosomal": ("HP:0000006", 1),
    "biallelic_autosomal": ("HP:0000007", 2),
    "monoallelic_PAR": ("g2p:monoallelic_PAR", 1),
    "biallelic_PAR": ("g2p:biallelic_PAR", 2),
    "monoallelic_X_hemizygous": ("g2p:monoallelic_X_hemizygous", 1),
    "monoallelic_X_heterozygous": ("g2p:monoallelic_X_heterozygous", 1),
    "monoallelic_Y_hemizygous": ("g2p:monoallelic_Y_hemizygous", 1),
    "mitochondrial": ("HP:0001427", None),
}
_ORPHA_TYPES = {
    "Disease-causing germline mutation(s) in": ("causal", "germline", None),
    "Disease-causing germline mutation(s) (loss of function) in": ("causal", "germline", "loss_of_function"),
    "Disease-causing germline mutation(s) (gain of function) in": ("causal", "germline", "gain_of_function"),
    "Disease-causing somatic mutation(s) in": ("causal", "somatic", None),
    "Modifying germline mutation in": ("modifier", "germline", None),
    "Major susceptibility factor in": ("susceptibility", None, None),
    "Role in the phenotype of": ("phenotype_role", None, None),
    "Part of a fusion gene in": ("fusion", None, None),
    "Candidate gene tested in": ("candidate", None, None),
    "Biomarker tested in": ("biomarker", None, None),
}
_FAMILIES = {
    "orphanet": "orphanet", "orphadata": "orphanet", "clingen": "clingen",
    "g2p": "gene2phenotype", "gene2phenotype": "gene2phenotype",
    "genomics_england": "genomics_england", "genomics england panelapp": "genomics_england",
}


class DiseaseGeneEvidenceError(RuntimeError):
    """Malformed evidence or a public snapshot that is not bound to its pins."""


def is_disease_identifier(value: object, *, allow_native_token: bool = False) -> bool:
    """Non-scoring syntax validation; never normalize or equate disease anchors.

    Catalog Symbol identifiers require an exact registered namespace and canonical
    percent encoding of the bounded source symbol. Arbitrary Symbol namespaces,
    raw narrative, malformed/double escapes and control characters do not qualify.
    A syntactically valid identifier still needs separately verified source proof.

    allow_native_token preserves the phenotype scorer's existing unqualified/raw
    identifier vocabulary, not observation validation or broader source trust.
    It cannot bypass the registered source-domain or Symbol-namespace checks.
    """
    if not isinstance(value, str):
        return False
    for namespace in SOURCE_NATIVE_DISEASE_NAMESPACES:
        if value.startswith(namespace):
            token = value[len(namespace):]
            if len(token) > 384:
                return False
            try:
                symbol = unquote(token, encoding="utf-8", errors="strict")
            except UnicodeError:
                return False
            return bool(re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.:/()+ -]{0,127}", symbol)
                        and symbol == symbol.strip() and quote(symbol, safe="") == token)
    if ":Symbol:" in value or any(value.startswith(namespace.split(":", 1)[0] + ":")
                                  for namespace in SOURCE_NATIVE_DISEASE_NAMESPACES):
        return False
    pattern = r"[A-Za-z0-9][A-Za-z0-9_.:-]{0,127}" if allow_native_token else r"[A-Za-z][A-Za-z0-9_]*:[A-Za-z0-9][A-Za-z0-9_.-]{0,127}"
    return re.fullmatch(pattern, value) is not None


@dataclass(frozen=True, slots=True)
class AggregationPolicy:
    version: str = "disease-gene-development-not-release-frozen"
    direct_cap: float = 1.0
    contextual_cap: float = 0.1
    conflict_penalty: float = 0.2
    shared_source_discount: float = 1.0
    weights: tuple[tuple[str, float], ...] = (
        ("definitive", 1.0), ("strong", 0.9), ("moderate", 0.7),
        ("supportive", 0.6), ("limited", 0.4), ("validated", 0.6),
        ("green", 0.7), ("amber", 0.4), ("aggregate", 0.1),
    )

    def __post_init__(self):
        if not isinstance(self.weights, tuple) or any(
            not isinstance(item, tuple) or len(item) != 2 or not isinstance(item[0], str) for item in self.weights
        ):
            raise DiseaseGeneEvidenceError("policy weights must be immutable named pairs")
        numbers = (self.direct_cap, self.contextual_cap, self.conflict_penalty, self.shared_source_discount,
                   *(value for _, value in self.weights))
        if any(type(value) not in (int, float) or not math.isfinite(value) or value < 0 for value in numbers):
            raise DiseaseGeneEvidenceError("policy values must be finite and nonnegative")
        if (
            not isinstance(self.version, str) or "development" not in self.version
            or self.shared_source_discount != 1.0
            or {key for key, _ in self.weights} != {key for key, _ in AggregationPolicy.__dataclass_fields__["weights"].default}
            or len(dict(self.weights)) != len(self.weights)
            or not 0 < self.contextual_cap < min(value for key, value in self.weights if key != "aggregate") - self.conflict_penalty
            or self.direct_cap < max(value for _, value in self.weights)
        ):
            raise DiseaseGeneEvidenceError("invalid development caps, weights or shared-source policy")


DEFAULT_POLICY = AggregationPolicy()


def _plain(value: object) -> object:
    if is_dataclass(value):
        return {field.name: _plain(getattr(value, field.name)) for field in fields(value)}
    if isinstance(value, Mapping):
        return {key: _plain(child) for key, child in value.items()}
    if isinstance(value, (tuple, list)):
        return [_plain(child) for child in value]
    return value


def _bytes(value: object) -> bytes:
    try:
        return (json.dumps(_plain(value), sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode()
    except (TypeError, ValueError, RecursionError):
        raise DiseaseGeneEvidenceError("evidence must be finite JSON") from None


def _freeze(value: object) -> object:
    if isinstance(value, Mapping):
        if any(not isinstance(key, str) for key in value):
            raise DiseaseGeneEvidenceError("metadata keys must be strings")
        return MappingProxyType({key: _freeze(child) for key, child in value.items()})
    if isinstance(value, (list, tuple)):
        return tuple(_freeze(child) for child in value)
    if value is None or type(value) in (str, bool, int, float):
        return value
    raise DiseaseGeneEvidenceError("evidence must be JSON-compatible")


def _hash(value: object) -> str:
    return hashlib.sha256(_bytes(value)).hexdigest()


def _components(links: Iterable[tuple[str, str]]) -> dict[str, str]:
    neighbours = defaultdict(set)
    for left, right in links:
        neighbours[left].add(right)
        neighbours[right].add(left)
    result = {}
    for node in sorted(neighbours):
        if node in result:
            continue
        pending, members = [node], set()
        while pending:
            current = pending.pop()
            if current not in members:
                members.add(current)
                pending.extend(neighbours[current] - members)
        result.update(dict.fromkeys(members, node))
    return result


def _moi(curie: str | None, title: str | None) -> tuple[str | None, tuple[str, ...]]:
    by_id = curie if curie in _MOI else None
    by_title = _MOI_TITLES.get((title or "").strip().casefold())
    if curie and title and by_id != by_title:
        return None, ("inheritance_id_title_disagreement_or_unknown",)
    result = by_id or by_title
    return result, () if result else ("unknown_inheritance",)


@dataclass(frozen=True, slots=True)
class SourceRecord:
    record_id: str
    table: str
    source_fields: Mapping[str, object]


@dataclass(frozen=True, slots=True)
class GeneRequirements:
    inheritance: str | None
    allelic_requirement: str | None
    allelic_requirement_kind: str
    required_allele_count: int | None
    association_requirement: str | None
    consequence_requirement: str | None
    mechanism_direction: str | None
    allele_origin: str | None
    record_ids: tuple[str, ...]
    diagnostics: tuple[str, ...]


def _requirements_disagree(left: GeneRequirements, right: GeneRequirements) -> bool:
    """Compare explicit source requirements; missing values are not matches."""
    for name in ("inheritance", "allelic_requirement", "required_allele_count", "consequence_requirement",
                 "mechanism_direction", "allele_origin"):
        if name == "allelic_requirement" and (
            left.allelic_requirement_kind != "literal_allelic_model"
            or right.allelic_requirement_kind != "literal_allelic_model"
        ):
            continue
        first, second = getattr(left, name), getattr(right, name)
        if first is not None and second is not None and first != second:
            return True
    return False


@dataclass(frozen=True, slots=True)
class GeneAssertion:
    """Gene evidence, or an explicitly scoped native disease-locus assertion.

    A catalog gene mapping does not turn locus evidence into unrestricted gene
    curation. Native locus assertions require the locus in the requested context.
    """
    record_id: str
    source_family: str
    gene_id: str | None
    polarity: str
    validity: str | None
    evidence_kind: str
    requirements: tuple[GeneRequirements, ...]
    diagnostics: tuple[str, ...]
    locus_id: str | None = None
    native_disease_id: str | None = None
    locus_record_id: str | None = None


@dataclass(frozen=True, slots=True)
class DiseaseEvidenceRelation:
    subject_id: str | None
    object_id: str | None
    kind: str
    record_ids: tuple[str, ...]
    diagnostics: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class DiseaseGeneEvidenceStore:
    records: Mapping[str, SourceRecord]
    assertions: tuple[GeneAssertion, ...]
    views: Mapping[str, tuple[str, ...]]
    relations: tuple[DiseaseEvidenceRelation, ...]
    unassigned_record_ids: tuple[str, ...]
    source_binding: Mapping[str, object]
    evidence_sha256: str


TRACK2_ALLELE_POLICY = 'track2-affected-copy-derivation/1.0.0'


def track2_requirement_projection(store):
    """Derive Track 2 copy requirements without changing shared Track 1 rules.

    Raw source records and their identities remain intact. The projection has
    a distinct digest. Each requirement records its derivation in diagnostics.
    Native gene/subtype groups retain specific conflicts. Broad disease rows
    retain their lower specificity under the existing applicability owner.
    """
    from dataclasses import replace
    projected=[]
    for assertion in store.assertions:
        requirements=[]
        for requirement in assertion.requirements:
            inheritance=requirement.inheritance
            count=requirement.required_allele_count
            literal=requirement.allelic_requirement
            kind=requirement.allelic_requirement_kind
            basis='source_stated' if count is not None else 'unresolved'
            literal_conflict=False
            if count is None and literal == 'BIALLELIC, autosomal or pseudoautosomal':
                count=2
                literal_conflict=inheritance not in (None,'HP:0000007')
                inheritance=inheritance or 'HP:0000007'
                kind='mode_of_inheritance'
                basis='vocabulary_interpreted_genomics_england_biallelic'
            elif count is None and inheritance in {'HP:0000006','HP:0000007'}:
                count=1 if inheritance=='HP:0000006' else 2
                basis='inheritance_derived'
            diagnostics=tuple(d for d in requirement.diagnostics
                if not (count is not None and d in {'unknown_required_allele_count','unknown_inheritance'}))
            requirements.append(replace(requirement,inheritance=inheritance,required_allele_count=count,
                allelic_requirement_kind=kind,
                diagnostics=(*diagnostics,TRACK2_ALLELE_POLICY,'affected_copies:'+basis,
                    *(('track2_specific_requirement_conflict',) if literal_conflict else ()))))
        projected.append(replace(assertion,requirements=tuple(requirements)))
    specific={}
    for assertion in projected:
        if assertion.polarity!='positive': continue
        for requirement in assertion.requirements:
            if 'disease_level_inheritance_not_gene_specific' not in requirement.diagnostics:
                key=(assertion.gene_id,assertion_native_disease_id(store,assertion),
                     requirement.mechanism_direction,requirement.consequence_requirement)
                specific.setdefault(key,set()).add((requirement.inheritance,requirement.required_allele_count))
    conflicts={key for key,values in specific.items()
        if len({v[0] for v in values if v[0] is not None})>1
        or len({v[1] for v in values if v[1] is not None})>1}
    projected=[replace(a,requirements=tuple(replace(r,diagnostics=(*r.diagnostics,
        'track2_specific_requirement_conflict')) if (a.gene_id,assertion_native_disease_id(store,a),
        r.mechanism_direction,r.consequence_requirement) in conflicts else r for r in a.requirements))
        for a in projected]
    return replace(store,assertions=tuple(projected),evidence_sha256=_hash({
        'source_evidence_sha256':store.evidence_sha256,'policy':TRACK2_ALLELE_POLICY,'assertions':projected}))


def build_disease_gene_evidence(
    rows: Iterable[Mapping[str, object]], *, source_binding: Mapping[str, object],
    inheritance_rows: Iterable[Mapping[str, object]] = (),
    disease_rows: Iterable[Mapping[str, object]] = (),
    alias_rows: Iterable[Mapping[str, object]] = (),
    relation_rows: Iterable[Mapping[str, object]] = (),
) -> DiseaseGeneEvidenceStore:
    """Pure public/synthetic v6 builder; every row and extra source field survives.

    Only native identifiers and explicit row mappings form views. No display
    joins, preference filtering or transitive equivalence through mapped views.
    Orphadata disease-level inheritance is linked only to its native disease and
    remains labeled disease-level, not an invented gene-specific curation.
    """
    if not isinstance(source_binding, Mapping):
        raise DiseaseGeneEvidenceError("source binding must be an object")
    binding = _freeze(source_binding)
    _bytes(binding)
    records = {}
    by_table = defaultdict(list)
    occurrences = defaultdict(int)
    for table, items in zip(PUBLIC_TABLES, (rows, inheritance_rows, disease_rows, alias_rows, relation_rows), strict=True):
        schema = dict(part.strip().replace('"', '').split() for part in TABLE_SCHEMAS[table].split(","))
        for row in items:
            if not isinstance(row, Mapping) or not set(TABLE_COLUMNS[table]).issubset(row):
                raise DiseaseGeneEvidenceError(f"{table} is missing v6 source fields")
            for name, dtype in schema.items():
                value = row[name]
                valid = value is None
                if dtype == "VARCHAR[]":
                    valid |= isinstance(value, (tuple, list)) and all(item is None or isinstance(item, str) for item in value)
                elif dtype == "BOOLEAN":
                    valid |= type(value) is bool
                elif dtype in ("BIGINT", "INTEGER"):
                    valid |= type(value) is int
                elif dtype == "DOUBLE":
                    valid |= type(value) in (int, float) and math.isfinite(value)
                else:
                    valid |= isinstance(value, str)
                if not valid:
                    raise DiseaseGeneEvidenceError(f"invalid {table} field type: {name}")
            provenance = ["source_release"]
            if table != "disease_relation":
                provenance.append("source_dataset")
            if table not in ("disease", "disease_relation"):
                provenance.append("source_record_id")
            if table == "gene_assertion":
                provenance.append("source_family")
            if any(not row[name] or not row[name].strip() for name in provenance):
                raise DiseaseGeneEvidenceError(f"{table} lacks source record provenance")
            frozen = _freeze(row)
            digest = _hash((table, frozen))
            occurrences[digest] += 1
            record = SourceRecord(f"{digest}:{occurrences[digest]}", table, frozen)
            records[record.record_id] = record
            by_table[table].append(record)
    catalog = defaultdict(list)
    inheritance = defaultdict(list)
    for record in by_table["disease"]:
        catalog[record.source_fields["disease_id"]].append(record)
    for record in by_table["disease_inheritance"]:
        inheritance[record.source_fields["disease_id"]].append(record)
    relations = []
    for record in by_table["disease"]:
        row = record.source_fields
        target_id = row["aggregation_target_id"]
        if not target_id or target_id == row["disease_id"]:
            continue
        targets = catalog[target_id]
        target = targets[0].source_fields if len(targets) == 1 else {}
        valid = (
            row["source_dataset"] == target.get("source_dataset") == "orpha_nomenclature"
            and row["source_release"] == target.get("source_release")
            and row["status"] in ("Active", "Active: Historic")
            and target.get("status") in ("Active", "Active: Historic")
            and row["classification_level_id"] == "36554"
            and row["classification_level"] == "Subtype of disorder"
            and target.get("classification_level_id") == "36547"
            and target.get("classification_level") == "Disorder"
            and row["aggregation_status"] == "Applicable"
        )
        relations.append(DiseaseEvidenceRelation(
            row["disease_id"], target_id, "narrower" if valid else "unknown",
            tuple(sorted({record.record_id, *(item.record_id for item in targets)})),
            ("active_subtype_aggregation",) if valid else ("unresolved_aggregation_relation",),
        ))
    for record in by_table["disease_alias"]:
        row = record.source_fields
        # Compound NTBT/E and BTNT/E are not exact concept equivalence.
        label = (row["mapping_relation"] or "").split(" (", 1)[0].strip()
        kind = {"E": "equivalent", "NTBT": "narrower", "BTNT": "broader"}.get(label, "unknown")
        if row["source_dataset"] != "orphadata_product1" or row["mapping_validation_status"] != "Validated":
            kind = "unknown"
        if not row["disease_id"] or not row["normalised_reference"]:
            kind = "unknown"
        relations.append(DiseaseEvidenceRelation(
            row["disease_id"], row["normalised_reference"], kind, (record.record_id,),
            ("unresolved_mapping_relation",) if kind == "unknown" else (),
        ))
    for record in by_table["disease_relation"]:
        row = record.source_fields
        relations.append(DiseaseEvidenceRelation(
            row["root_disease_id"], row["target_disease_id"], "unknown", (record.record_id,),
            ("replacement_or_unsupported_relation_not_hierarchy",),
        ))
    equivalences = _components((item.subject_id, item.object_id) for item in relations if item.kind == "equivalent")
    hierarchy = defaultdict(set)
    for item in relations:
        if item.kind in ("narrower", "broader"):
            child, parent = (item.subject_id, item.object_id) if item.kind == "narrower" else (item.object_id, item.subject_id)
            hierarchy[equivalences.get(child, child)].add(equivalences.get(parent, parent))
    invalid = set()
    while True:
        try:
            TopologicalSorter({node: parents - invalid for node, parents in hierarchy.items() if node not in invalid}).prepare()
            break
        except CycleError as error:
            invalid.update(error.args[1])
    if invalid:
        relations = [replace(item, kind="unknown", diagnostics=(*item.diagnostics, "cyclic_or_conflicting_disease_relation"))
                     if equivalences.get(item.subject_id, item.subject_id) in invalid
                     or equivalences.get(item.object_id, item.object_id) in invalid else item for item in relations]
    assertions, unassigned = [], []
    views = defaultdict(list)
    for record in by_table["gene_assertion"]:
        row = record.source_fields
        dataset = row["source_dataset"]
        family = _FAMILIES.get(row["source_family"].strip().casefold(), row["source_family"].strip().casefold())
        gene_id = row["ensembl_gene_id"] or row["gene_id"]
        diagnostics = []
        if not gene_id:
            diagnostics.append("unknown_gene_identity_no_display_join")
        if not row["ensembl_gene_id"]:
            diagnostics.append("missing_ensembl_mapping_native_identity_retained")
        if row["disease_specificity"] not in ("subtype", "disorder", "broad"):
            diagnostics.append("unknown_disease_specificity")
        native = row["source_disease_id"]
        if not native:
            diagnostics.append("unknown_native_disease")
        polarity, validity, kind = "unknown", None, "unknown"
        association, origin, direction = None, None, None
        consequence = None
        if dataset == "gencc_submissions":
            by_id = _GENCC_VALIDITY.get(row["validity_id"])
            title = (row["validity"] or "").strip().casefold()
            title = {"disputed evidence": "disputed", "refuted evidence": "refuted"}.get(title, title)
            if by_id and title and by_id != title:
                diagnostics.append("validity_id_title_disagreement")
            elif row["validity_id"] and not by_id:
                diagnostics.append("unknown_validity_identifier")
            else:
                validity = by_id or (title if title in _GENCC_VALIDITY.values() else None)
            kind = "curated"
        elif dataset in ("evidence_clingen", "evidence_gene2phenotype", "evidence_genomics_england", "evidence_orphanet"):
            native_family = _FAMILIES[OT_EVIDENCE_DATASETS[dataset]]
            if family != native_family:
                diagnostics.append("source_family_disagreement")
            family = native_family
            titles = {(row[name] or "").strip().casefold() for name in ("validity", "confidence")} - {""}
            titles = {{"disputed evidence": "disputed", "refuted evidence": "refuted"}.get(title, title) for title in titles}
            if len(titles) == 1:
                validity = next(iter(titles))
            elif titles:
                diagnostics.append("validity_confidence_disagreement")
            kind = "curated"
            if dataset in ("evidence_gene2phenotype", "evidence_orphanet"):
                value = row["consequence_requirement"]
                if value and re.fullmatch(r"SO[:_][0-9]{7}", value):
                    consequence = value.replace("_", ":", 1)
                elif value:
                    diagnostics.append("unsupported_consequence_requirement")
            if dataset == "evidence_orphanet":
                # Its ETL also admits modifiers, so no unspecified consequence
                # is promoted to a proven causative mechanism.
                direction = {"SO:0002054": "loss_of_function", "SO:0002053": "gain_of_function"}.get(consequence)
                origin = "germline"
                if direction is None:
                    kind = "curated_role_unknown"
        elif dataset == "orphadata_product6":
            family = "orphanet"
            association = row["association_type"]
            role, origin, direction = _ORPHA_TYPES.get(association, ("unknown", None, None))
            kind = "curated" if role == "causal" else role
            validity = (row["validity"] or "").strip().casefold() or None
            if role == "unknown":
                diagnostics.append("unknown_association_type")
        elif dataset == "association_by_datasource_direct":
            kind = "aggregate"
            if family in set(_FAMILIES.values()) and row["score"] is not None and row["score"] > 0:
                polarity, validity = "contextual", "aggregate"
            diagnostics.append("aggregate_not_independent_curated_assertion")
        else:
            diagnostics.append("unsupported_source_dataset")
        if kind != "aggregate":
            positive_terms = {"definitive", "strong", "moderate", "supportive", "limited"}
            if dataset in ("orphadata_product6", "evidence_orphanet"):
                positive_terms = {"validated"}
            elif dataset == "evidence_genomics_england":
                positive_terms = {"green", "amber"}
            if validity in ("disputed", "refuted"):
                polarity = "negative"
            elif validity in positive_terms and kind not in ("unknown", "candidate", "biomarker"):
                polarity = "positive"
            elif validity == "animal model only":
                polarity, kind = "contextual", "animal_model_only"
            elif validity in ("no known disease relationship", "not validated"):
                polarity = "no_positive_support"
            else:
                diagnostics.append("unknown_or_unsupported_validity")
        requirements = []
        moi, moi_diagnostics = _moi(row["inheritance_id"], row["inheritance"])
        alleles = row["allelic_requirements"] or ()
        alternatives = [(moi, None, None, (record.record_id,), moi_diagnostics)]
        if dataset == "gencc_submissions":
            # v6 repeats moi_title in this array, not a count vocabulary.
            diagnostics.append("gencc_moi_is_not_literal_allele_count")
            if alleles:
                alternatives = []
                for allele in dict.fromkeys(alleles):
                    token, flags = _moi(None, allele)
                    if token != moi:
                        token, flags = None, ("inheritance_requirement_disagreement",)
                    alternatives.append((token, allele, None, (record.record_id,), (*moi_diagnostics, *flags)))
        elif dataset in ("evidence_clingen", "evidence_genomics_england", "evidence_gene2phenotype") and alleles:
            alternatives = []
            for allele in dict.fromkeys(alleles):
                token, count = _G2P_ALLELES.get(allele, (None, None)) if dataset == "evidence_gene2phenotype" else (None, None)
                inferred, flags = _moi(None, allele)
                token = token or inferred
                flags = () if token else flags
                if moi and token and moi != token:
                    token, count, flags = None, None, ("inheritance_requirement_disagreement",)
                alternatives.append((token, allele, count, (record.record_id,), flags))
        elif dataset in ("orphadata_product6", "evidence_orphanet") and native:
            linked = [item for item in inheritance[native] if item.source_fields["source_dataset"] == "orphadata_product9_ages"]
            if linked:
                alternatives = []
                for item in linked:
                    # Orpha numeric IDs are source-local; interpret only its
                    # documented controlled inheritance titles, never HP IDs.
                    token, flags = _moi(None, item.source_fields["inheritance"])
                    alternatives.append((token, None, None, (record.record_id, item.record_id),
                                         (*flags, "disease_level_inheritance_not_gene_specific")))
        for token, allele, count, refs, flags in alternatives:
            allele_kind = "unknown"
            if dataset == "evidence_gene2phenotype" and allele in _G2P_ALLELES:
                allele_kind = "literal_allelic_model"
            elif _moi(None, allele)[0]:
                allele_kind = "mode_of_inheritance"
            requirements.append(GeneRequirements(
                token, allele, allele_kind, count, association, consequence, direction, origin, tuple(sorted(refs)),
                (*flags, *(("unknown_required_allele_count",) if count is None else ()),
                 *(("unknown_mechanism_direction",) if direction is None else ())),
            ))
        assertion = GeneAssertion(record.record_id, family, gene_id, polarity, validity, kind,
                                  tuple(sorted(requirements, key=_bytes)), tuple(diagnostics))
        assertions.append(assertion)
        ids = {value for value in (native, row["mapped_disease_id"]) if value}
        if not ids or not gene_id:
            unassigned.append(record.record_id)
        for disease_id in ids:
            views[disease_id].append(record.record_id)
    assertions = tuple(sorted(assertions, key=lambda item: item.record_id))
    views = MappingProxyType({key: tuple(sorted(value)) for key, value in sorted(views.items())})
    relations = tuple(sorted(relations, key=_bytes))
    # Record IDs bind the complete raw rows. Hash interpretations once rather
    # than serializing all raw payloads on every context evaluation.
    digest = _hash((SCHEMA, sorted(records), assertions, views, relations, sorted(unassigned), binding))
    return DiseaseGeneEvidenceStore(MappingProxyType(dict(sorted(records.items()))), assertions, views,
                                    relations, tuple(sorted(unassigned)), binding, digest)


@dataclass(frozen=True, slots=True)
class DiseaseGeneContext:
    disease_id: str
    gene_id: str | None
    inheritance: str | None
    native_disease_id: str | None = None
    allelic_requirement: str | None = None
    required_allele_count: int | None = None
    consequence_requirement: str | None = None
    mechanism_direction: str | None = None
    allele_origin: str | None = None
    locus_id: str | None = None

    def __post_init__(self):
        for field in fields(self):
            value = getattr(self, field.name)
            if field.name == "required_allele_count":
                if value is not None and (type(value) is not int or value <= 0):
                    raise DiseaseGeneEvidenceError("context allele count must be a positive integer or unknown")
            elif value is not None and (not isinstance(value, str) or not value.strip()):
                raise DiseaseGeneEvidenceError("context identifiers/requirements must be nonempty strings or unknown")
        if not self.disease_id or not (self.gene_id or self.locus_id):
            raise DiseaseGeneEvidenceError("context requires disease and native gene or locus identifiers")
        if self.inheritance:
            normalized, _ = _moi(self.inheritance if self.inheritance.startswith("HP:") else None, None
                                 if self.inheritance.startswith("HP:") else self.inheritance)
            if normalized:
                object.__setattr__(self, "inheritance", normalized)
            elif self.inheritance not in {value[0] for value in _G2P_ALLELES.values()}:
                raise DiseaseGeneEvidenceError("unsupported context inheritance; use unknown, not a wildcard")
        if self.consequence_requirement and re.fullmatch(r"SO_[0-9]{7}", self.consequence_requirement):
            object.__setattr__(self, "consequence_requirement", self.consequence_requirement.replace("_", ":", 1))


@dataclass(frozen=True, slots=True)
class RequirementEvaluation:
    requirements: GeneRequirements
    checks: Mapping[str, str]


@dataclass(frozen=True, slots=True)
class GeneEvidenceContribution:
    record_id: str
    source_family: str
    disease_applicability: str
    relation_record_ids: tuple[str, ...]
    requirement_evaluations: tuple[RequirementEvaluation, ...]
    raw_support: float | None
    contribution: float
    diagnostics: tuple[str, ...]
    discounted_by_record_ids: tuple[str, ...] = ()
    support_channel: str = "unscored"


@dataclass(frozen=True, slots=True)
class GeneEvidencePenalty:
    penalty_id: str
    reason: str
    amount: float
    positive_record_ids: tuple[str, ...]
    negative_record_ids: tuple[str, ...]
    component: str = "technical_and_conflict_penalties"


@dataclass(frozen=True, slots=True)
class DiseaseGeneAggregation:
    context: DiseaseGeneContext
    policy: AggregationPolicy
    source_binding: Mapping[str, object]
    evidence_sha256: str
    ledger: tuple[GeneEvidenceContribution, ...]
    penalties: tuple[GeneEvidencePenalty, ...]
    direct_support: float
    contextual_support: float
    positive_support: float
    conflict_penalty: float
    diagnostics: tuple[str, ...]


def assertion_native_disease_id(store: DiseaseGeneEvidenceStore, assertion: GeneAssertion) -> str | None:
    """Use the source adapter's explicit native anchor, never a mapped view."""
    if assertion.locus_id is not None:
        return assertion.native_disease_id
    return store.records[assertion.record_id].source_fields["source_disease_id"]


def assertion_matches_target(assertion: GeneAssertion, context: DiseaseGeneContext) -> bool:
    """Locus and gene namespaces are distinct; neither is a wildcard for the other."""
    if assertion.locus_id is not None or context.locus_id is not None:
        return (assertion.locus_id == context.locus_id and context.locus_id is not None
                and (context.gene_id is None or assertion.gene_id == context.gene_id))
    return assertion.gene_id == context.gene_id


def rebuild_disease_gene_evidence(store: DiseaseGeneEvidenceStore) -> DiseaseGeneEvidenceStore:
    """Replay a known source adapter, not arbitrary stored interpretations.

    Native catalogs additionally rehash retained original bytes against their
    resource/manifest pins. Existing v6 stores keep their original gene semantics.
    Callers must compare the returned store to the supplied one.
    """
    if 'track2_concept_supplement' in store.source_binding:
        from genosapiens.track2_concept_supplement import supplement_store
        rows = {table: [record.source_fields for record in store.records.values() if record.table == table]
                for table in PUBLIC_TABLES}
        if any(record.table not in {*PUBLIC_TABLES, 'track2_concept'} for record in store.records.values()):
            raise DiseaseGeneEvidenceError('unsupported supplemented disease evidence table')
        base = build_disease_gene_evidence(rows['gene_assertion'], source_binding={
            k:v for k,v in store.source_binding.items() if k != 'track2_concept_supplement'},
            inheritance_rows=rows['disease_inheritance'], disease_rows=rows['disease'],
            alias_rows=rows['disease_alias'], relation_rows=rows['disease_relation'])
        rebuilt = supplement_store(base)
        if rebuilt.source_binding['track2_concept_supplement'] != store.source_binding['track2_concept_supplement']:
            raise DiseaseGeneEvidenceError('Track 2 concept source changed')
        return rebuilt
    if "catalog_utf8" in store.source_binding or any(record.table in NATIVE_LOCUS_TABLES for record in store.records.values()):
        from genosapiens.repeat_catalog_admission import (
            SYNTHETIC_ADAPTER_VERSION, load_pinned_repeat_catalog, load_synthetic_repeat_catalog,
        )

        binding = store.source_binding
        try:
            if binding.get("adapter_version") == SYNTHETIC_ADAPTER_VERSION:
                return load_synthetic_repeat_catalog(
                    binding["catalog_utf8"].encode("utf-8"),
                    profile_bytes=binding["contract_utf8"].encode("utf-8"),
                    manifest_bytes=binding["manifest_utf8"].encode("utf-8"),
                    expected_profile_sha256=binding["contract_sha256"],
                    expected_manifest_sha256=binding["bundle_manifest_sha256"],
                ).store
            return load_pinned_repeat_catalog(
                binding["catalog_utf8"].encode("utf-8"),
                contract_bytes=binding["contract_utf8"].encode("utf-8"),
                manifest_bytes=binding["manifest_utf8"].encode("utf-8"),
                expected_contract_sha256=binding["contract_sha256"],
                expected_manifest_sha256=binding["bundle_manifest_sha256"],
            ).store
        except (KeyError, TypeError, AttributeError, ValueError):
            raise DiseaseGeneEvidenceError("invalid native locus source binding") from None
    if any(record.table not in PUBLIC_TABLES for record in store.records.values()):
        raise DiseaseGeneEvidenceError("unsupported disease evidence source table")
    rows = {table: [record.source_fields for record in store.records.values() if record.table == table]
            for table in PUBLIC_TABLES}
    return build_disease_gene_evidence(
        rows["gene_assertion"], source_binding=store.source_binding,
        inheritance_rows=rows["disease_inheritance"], disease_rows=rows["disease"],
        alias_rows=rows["disease_alias"], relation_rows=rows["disease_relation"],
    )


class PreparedDiseaseGeneStore:
    """Private exact-owner indexes over an already validated immutable store.

    Candidate indexes preserve original assertion/relation order. All target,
    applicability and requirement decisions remain in the original kernels.
    """
    def __init__(self, store):
        self.store = store
        genes, loci, relations = defaultdict(list), defaultdict(list), defaultdict(list)
        edges, source_edges = defaultdict(list), defaultdict(set)
        for assertion in store.assertions:
            genes[assertion.gene_id].append(assertion)
            loci[assertion.locus_id].append(assertion)
        for ordinal, relation in enumerate(store.relations):
            for key in relation.record_ids:
                relations[key].append(ordinal)
            if not relation.subject_id or not relation.object_id:
                continue
            if relation.kind in ("equivalent", "narrower"):
                edges[relation.subject_id].append((relation.object_id, relation.kind, relation.record_ids))
                source_edges[relation.subject_id].add(relation.object_id)
            if relation.kind in ("equivalent", "broader"):
                edges[relation.object_id].append((relation.subject_id,
                    "equivalent" if relation.kind == "equivalent" else "narrower", relation.record_ids))
                source_edges[relation.object_id].add(relation.subject_id)
        self.genes = MappingProxyType({key:tuple(rows) for key,rows in genes.items()})
        self.loci = MappingProxyType({key:tuple(rows) for key,rows in loci.items()})
        self.relations = MappingProxyType({key:tuple(rows) for key,rows in relations.items()})
        self.edges = MappingProxyType({key:tuple(rows) for key,rows in edges.items()})
        self.source_edges = MappingProxyType({key:frozenset(rows) for key,rows in source_edges.items()})
        self.equivalences = MappingProxyType(_components((item.subject_id,item.object_id)
            for item in store.relations if item.kind == "equivalent"))

    def require(self, store):
        if store is not self.store:
            raise DiseaseGeneEvidenceError("prepared disease gene store owner changed")

    def referenced_relations(self, refs):
        for ordinal in sorted({ordinal for key in refs for ordinal in self.relations.get(key,())}):
            relation = self.store.relations[ordinal]
            if relation.record_ids and set(relation.record_ids) <= refs:
                yield relation

    def source_assertions(self, context):
        return self.store.assertions if context.gene_id is None else self.genes.get(context.gene_id,())

    def aggregation_assertions(self, context):
        # With a locus and no gene, the original target predicate accepts every
        # matching locus gene. Without a locus, even a None gene is exact.
        return (self.loci.get(context.locus_id,()) if context.locus_id is not None
                else self.genes.get(context.gene_id,()))


def aggregate_disease_gene_evidence(
    store: DiseaseGeneEvidenceStore, context: DiseaseGeneContext, *, policy: AggregationPolicy = DEFAULT_POLICY,
    _prepared=None,
) -> DiseaseGeneAggregation:
    """Evaluate one disease/native-gene/inheritance requirement context.

    native_disease_id optionally anchors a source-specific interpretation inside
    an explicit mapped view. This is not equivalence to that view. Otherwise a
    mapped-only row remains disease-applicability unknown. Missing inheritance
    is never compatible by default; unresolved context support shares the small
    contextual cap with aggregates, broad and indirect evidence. Unknown count
    or mechanism does not manufacture a requirement or invalidate known MOI.

    Positive support is NOT net of conflict_penalty. The single penalty has an
    owner component/id for later fusion to charge once, not in both components.
    Native catalog assertions participate in locus applicability/requirements,
    but have no numerical validity weight. Their null raw support is not negative
    evidence and does not undo source-backed admission.
    """
    if not isinstance(context, DiseaseGeneContext) or not isinstance(policy, AggregationPolicy):
        raise DiseaseGeneEvidenceError("typed context and aggregation policy are required")
    if _prepared is not None:
        if type(_prepared) is not PreparedDiseaseGeneStore:
            raise DiseaseGeneEvidenceError("invalid prepared disease gene store")
        _prepared.require(store)
    assertions = store.assertions if _prepared is None else _prepared.aggregation_assertions(context)
    anchor = context.native_disease_id or context.disease_id
    view_ids = set(store.views.get(context.disease_id, ()))
    if anchor != context.disease_id and not any(
        assertion_matches_target(item, context) and item.record_id in view_ids
        and assertion_native_disease_id(store, item) == anchor
        for item in assertions
    ):
        raise DiseaseGeneEvidenceError("native context is not explicitly represented in this disease view")
    edges = defaultdict(list) if _prepared is None else _prepared.edges
    equivalences = (_components((item.subject_id, item.object_id) for item in store.relations if item.kind == "equivalent")
                    if _prepared is None else _prepared.equivalences)
    for relation in store.relations if _prepared is None else ():
        if not relation.subject_id or not relation.object_id:
            continue
        if relation.kind in ("equivalent", "narrower"):
            edges[relation.subject_id].append((relation.object_id, relation.kind, relation.record_ids))
        if relation.kind in ("equivalent", "broader"):
            edges[relation.object_id].append((relation.subject_id, "equivalent" if relation.kind == "equivalent" else "narrower",
                                              relation.record_ids))
    # Typed reachability: retain proof references, never follow unknown links.
    reached = {anchor: (0, ())}
    queue = deque([anchor])
    while queue:
        current = queue.popleft()
        distance, refs = reached[current]
        for target, kind, proof in sorted(edges.get(current, ())):
            candidate = (distance + (kind == "narrower"), tuple(sorted(set(refs) | set(proof))))
            if target not in reached or candidate[0] < reached[target][0]:
                reached[target] = candidate
                queue.append(target)
    ledger, details = [], {}
    weights = dict(policy.weights)
    for assertion in assertions:
        if not assertion_matches_target(assertion, context):
            continue
        row = store.records[assertion.record_id].source_fields
        native = assertion_native_disease_id(store, assertion)
        distance, proof = reached.get(native, (None, ()))
        in_view = assertion.record_id in view_ids
        if distance is None and not in_view:
            continue
        applicability = "exact" if native == anchor else "equivalent" if distance == 0 else "broader" if distance else "unknown"
        flags = list(assertion.diagnostics)
        if applicability == "unknown":
            flags.append("mapped_view_not_equivalence_applicability_unknown")
        if anchor != context.disease_id:
            flags.append("source_native_interpretation_in_mapped_view")
        evaluations = []
        for requirement in assertion.requirements:
            checks = {}
            for name in ("inheritance", "allelic_requirement", "required_allele_count", "consequence_requirement",
                         "mechanism_direction", "allele_origin"):
                source, requested = getattr(requirement, name), getattr(context, name)
                if name == "allelic_requirement" and requirement.allelic_requirement_kind != "literal_allelic_model":
                    source = None
                checks[name] = "unknown" if source is None or requested is None else "match" if source == requested else "different"
            evaluations.append(RequirementEvaluation(requirement, MappingProxyType(checks)))
        compatible = [item for item in evaluations if "different" not in item.checks.values()]
        moi_match = any(item.checks["inheritance"] == "match" for item in compatible)
        if not compatible:
            flags.append("different_inheritance_or_requirements")
        elif not moi_match:
            flags.append("inheritance_applicability_unknown")
        indirect = row.get("is_ontology_propagated") is True
        scope_known = (row.get("association_scope") in ("source_direct", "ot_mapped_evidence", "ot_direct_aggregate")
                       and row.get("disease_specificity") in ("subtype", "disorder", "broad"))
        if indirect:
            flags.append("ontology_propagated")
        if not scope_known or row.get("is_ontology_propagated") is None:
            flags.append("unknown_association_scope")
        broad = applicability == "broader" or row.get("disease_specificity") == "broad"
        if broad:
            flags.append("broad_evidence")
        curated = assertion.evidence_kind == "curated"
        direct = distance == 0 and curated and not broad and not indirect and scope_known and row["is_ontology_propagated"] is False
        eligible = bool(compatible) and distance is not None
        # Unknown disease mappings remain visible but cannot imply applicability.
        raw = weights.get(assertion.validity) if assertion.polarity in ("positive", "contextual") else None
        if assertion.locus_id is not None:
            raw = None
        if assertion.polarity == "negative":
            flags.append("disputed_or_refuted_no_positive_support")
            if not scope_known or row["is_ontology_propagated"] is not False:
                flags.append("validity_conflict_scope_unknown_or_indirect")
        elif assertion.polarity == "contextual" and assertion.evidence_kind == "animal_model_only":
            raw = policy.contextual_cap
        entry = GeneEvidenceContribution(assertion.record_id, assertion.source_family, applicability, proof,
                                         tuple(evaluations), raw, 0.0, tuple(flags))
        ledger.append(entry)
        details[assertion.record_id] = (assertion, direct, eligible, moi_match, distance)
    # A compatible, more specific curated validity assertion takes precedence
    # over umbrella evidence, including an umbrella's incompatible inheritance.
    specific = [entry for entry in ledger if details[entry.record_id][1:4] == (True, True, True)
                and details[entry.record_id][0].polarity in ("positive", "negative")]
    specific_polarities = {details[entry.record_id][0].polarity for entry in specific}
    matching_requirements = {
        entry.record_id: tuple(evaluation.requirements for evaluation in entry.requirement_evaluations
                               if evaluation.checks["inheritance"] == "match" and "different" not in evaluation.checks.values())
        for entry in ledger
    }
    represented_assertions = defaultdict(list)
    for entry in ledger:
        assertion, _, eligible, moi_match, _ = details[entry.record_id]
        row = store.records[entry.record_id].source_fields
        native = assertion_native_disease_id(store, assertion)
        if (eligible and moi_match and native
                and row.get("source_dataset") in ("gencc_submissions", "orphadata_product6", *OT_EVIDENCE_DATASETS)
                and row["association_scope"] in ("source_direct", "ot_mapped_evidence")
                and row["is_ontology_propagated"] is False):
            represented_assertions[(assertion.source_family, equivalences.get(native, native))].append(entry.record_id)
    updated = []
    for entry in ledger:
        assertion, direct, eligible, moi_match, distance = details[entry.record_id]
        flags = list(entry.diagnostics)
        represented = []
        if assertion.evidence_kind == "aggregate" and eligible:
            native = store.records[entry.record_id].source_fields["source_disease_id"]
            for peer_id in represented_assertions[(assertion.source_family, equivalences.get(native, native))]:
                if any(not _requirements_disagree(evaluation.requirements, peer)
                       for evaluation in entry.requirement_evaluations if "different" not in evaluation.checks.values()
                       for peer in matching_requirements[peer_id]):
                    represented.append(peer_id)
        if represented:
            flags.append("aggregate_source_already_represented")
            eligible = False
        if specific and distance is not None and distance > 0:
            flags.append("umbrella_evidence_superseded")
            if assertion.polarity in ("positive", "negative") and specific_polarities - {assertion.polarity}:
                flags.append("umbrella_validity_conflict")
            if any(item.checks["inheritance"] == "different" for item in entry.requirement_evaluations):
                flags.append("umbrella_inheritance_conflict")
            eligible = False
        details[entry.record_id] = (assertion, direct, eligible, moi_match, distance)
        updated.append(replace(entry, diagnostics=tuple(flags), discounted_by_record_ids=tuple(sorted(represented))))
    ledger = updated
    positive_by_disease, negative_by_disease = defaultdict(list), defaultdict(list)
    for entry in ledger:
        assertion, direct, eligible, moi_match, distance = details[entry.record_id]
        row = store.records[entry.record_id].source_fields
        if (eligible and moi_match and assertion.evidence_kind == "curated"
                and row["is_ontology_propagated"] is False
                and row["disease_specificity"] in ("subtype", "disorder", "broad")
                and row["association_scope"] in ("source_direct", "ot_mapped_evidence")):
            native = row["source_disease_id"]
            specificity_key = equivalences.get(native, native)
            if assertion.polarity == "positive":
                positive_by_disease[specificity_key].append(entry.record_id)
            elif assertion.polarity == "negative":
                negative_by_disease[specificity_key].append(entry.record_id)
    positive_ids, negative_ids, incompatible_models = set(), set(), set()
    # Compare opposing source models directly, not their shared compatibility
    # with a partially specified request or a transitive unknown-source group.
    for disease_id in positive_by_disease.keys() & negative_by_disease.keys():
        for positive in positive_by_disease[disease_id]:
            for negative in negative_by_disease[disease_id]:
                if any(not _requirements_disagree(left, right)
                       for left in matching_requirements[positive] for right in matching_requirements[negative]):
                    positive_ids.add(positive)
                    negative_ids.add(negative)
                else:
                    incompatible_models.update((positive, negative))
    ledger = [replace(entry, diagnostics=(*entry.diagnostics, "opposing_assertion_has_incompatible_requirements"))
              if entry.record_id in incompatible_models else entry for entry in ledger]
    penalties = ()
    if positive_ids and negative_ids:
        penalties = (GeneEvidencePenalty(
            _hash((SCHEMA, context, "equal_specificity_validity_conflict")),
            "equal_specificity_validity_conflict", policy.conflict_penalty,
            tuple(sorted(positive_ids)), tuple(sorted(negative_ids)),
        ),)
    unresolved_requirements = set()
    for name in ("required_allele_count", "consequence_requirement", "mechanism_direction", "allele_origin"):
        if getattr(context, name) is not None:
            continue
        values = defaultdict(set)
        for entry in ledger:
            assertion, direct, eligible, moi_match, distance = details[entry.record_id]
            if eligible and moi_match and assertion.polarity == "positive":
                for evaluation in entry.requirement_evaluations:
                    value = getattr(evaluation.requirements, name)
                    if value is not None and "different" not in evaluation.checks.values():
                        values[value].add(entry.record_id)
        if len(values) > 1:
            unresolved_requirements.update(record for ids in values.values() for record in ids)
    # Connected dependence includes transitive shared publications; greedily
    # discounting only publications already credited would be order-dependent.
    links = []
    for entry in ledger:
        if not details[entry.record_id][2] or not entry.raw_support:
            continue
        links.append((entry.record_id, "family:" + entry.source_family))
        for value in store.records[entry.record_id].source_fields["publications"] or ():
            if value and (match := re.fullmatch(r"(?:PMID:)?([1-9][0-9]*)", value.strip(), re.IGNORECASE)):
                links.append((entry.record_id, "publication:" + match.group(1)))
    components = _components(links)
    credited_groups = {}
    totals = {"direct": 0.0, "contextual": 0.0}
    output = []
    for entry in sorted(ledger, key=lambda item: (
        not (details[item.record_id][1] and details[item.record_id][3] and item.record_id not in unresolved_requirements),
        -(item.raw_support or 0), store.records[item.record_id].source_fields.get("representation_kind") == "integrated_copy",
        item.record_id,
    )):
        assertion, direct, eligible, moi_match, distance = details[entry.record_id]
        bucket = "direct" if direct and moi_match and entry.record_id not in unresolved_requirements else "contextual"
        flags = list(entry.diagnostics)
        if entry.record_id in unresolved_requirements:
            flags.append("unresolved_requirement_alternatives_no_hybrid_support")
        group = components.get(entry.record_id)
        blocker = credited_groups.get(group)
        contribution = 0.0
        if eligible and entry.raw_support:
            if blocker:
                flags.append("shared_source_discount")
            else:
                cap = policy.direct_cap if bucket == "direct" else policy.contextual_cap
                contribution = min(entry.raw_support, max(0.0, cap - totals[bucket]))
                if contribution < entry.raw_support:
                    flags.append(f"{bucket}_cap")
                if contribution:
                    credited_groups[group] = entry.record_id
                    totals[bucket] += contribution
        output.append(replace(entry, contribution=contribution, diagnostics=tuple(flags), support_channel=bucket,
                              discounted_by_record_ids=(blocker,) if eligible and entry.raw_support and blocker else entry.discounted_by_record_ids))
    totals = {channel: math.fsum(entry.contribution for entry in output if entry.support_channel == channel)
              for channel in totals}
    diagnostics = []
    if not ledger:
        diagnostics.append("no_evidence_unknown_not_negative")
    if any(entry.disease_applicability == "unknown" for entry in ledger):
        diagnostics.append("unknown_disease_applicability")
    if any("inheritance_applicability_unknown" in entry.diagnostics for entry in ledger):
        diagnostics.append("unknown_inheritance_applicability")
    if penalties:
        diagnostics.append("unresolved_equal_specificity_validity_conflict")
    if unresolved_requirements:
        diagnostics.append("unresolved_requirement_alternatives")
    return DiseaseGeneAggregation(
        context, policy, store.source_binding, store.evidence_sha256, tuple(sorted(output, key=lambda item: item.record_id)), penalties,
        totals["direct"], totals["contextual"], totals["direct"] + totals["contextual"],
        sum(item.amount for item in penalties), tuple(diagnostics),
    )


def disease_gene_evidence_as_bytes(value: DiseaseGeneEvidenceStore | DiseaseGeneAggregation) -> bytes:
    """Deterministic JSON; source rows live once in the store, ledgers reference them."""
    if not isinstance(value, (DiseaseGeneEvidenceStore, DiseaseGeneAggregation)):
        raise DiseaseGeneEvidenceError("expected a disease-gene store or aggregation")
    return _bytes({"schema": SCHEMA, "kind": "store" if isinstance(value, DiseaseGeneEvidenceStore) else "aggregation",
                   "data": value})


def read_pinned_disease_gene_evidence(
    *, contract_path: str | Path, reference_root: str | Path, expected_bundle_sha256: str,
    diagnosis: str | None = None,
    track2_concepts: bool = False,
) -> DiseaseGeneEvidenceStore:
    """Read only governing JSON and the five necessary public v6 Parquets.

    The completed pinned bundle binds upstream releases/source manifests without
    reopening original XML/CSV, OT partitions, phenotype data or ontology files.
    Only gene assertions, inheritance and typed disease relations are consumed.
    There are no writes, repairs, downloads or persistent DuckDB databases.
    """
    try:
        return _read_pinned_disease_gene_evidence(
            contract_path=contract_path, reference_root=reference_root, expected_bundle_sha256=expected_bundle_sha256,
            diagnosis=diagnosis, track2_concepts=track2_concepts,
        )
    except (BoundedJSONError, OntologyBundleError, RareDiseaseError, OSError, KeyError, TypeError, ValueError):
        raise DiseaseGeneEvidenceError("invalid or mismatched public disease-gene input provenance") from None


def _read_pinned_disease_gene_evidence(
    *, contract_path: str | Path, reference_root: str | Path, expected_bundle_sha256: str,
    diagnosis: str | None = None,
    track2_concepts: bool = False,
) -> DiseaseGeneEvidenceStore:
    control_paths = (Path(contract_path),)
    signatures = {path: _file_signature(path) for path in control_paths}
    root = Path(reference_root).resolve(strict=True)
    if not isinstance(expected_bundle_sha256, str) or not _SHA256.fullmatch(expected_bundle_sha256):
        raise DiseaseGeneEvidenceError("a pinned bundle SHA-256 is required")
    contract_data, contract_bytes = load_bounded_json(contract_path)
    if not isinstance(contract_data, Mapping):
        raise DiseaseGeneEvidenceError("v6 source contract must be an object")
    contract = RareDiseaseContract.from_mapping(contract_data)
    if contract.projection.schema_version != "v6":
        raise DiseaseGeneEvidenceError("disease-gene evidence requires the immutable v6 source contract")
    contract_sha256 = hashlib.sha256(contract_bytes).hexdigest()
    manifest_path = _bundle_manifest_path(root, contract, contract_sha256)
    legacy_manifest_path = _bundle_manifest_path(root, contract)
    # A restricted reader can mount only the selected versioned manifest.
    # Producer naming must not make that reader depend on an unmounted legacy file.
    if not manifest_path.exists():
        versioned = legacy_manifest_path.with_name(
            f"{contract.bundle_id}.{contract_sha256[:12]}.bundle.json")
        if versioned.exists():
            manifest_path = versioned
    if (
        not manifest_path.exists()
        and legacy_manifest_path.exists()
        and _file_identity(legacy_manifest_path)[0] == expected_bundle_sha256
    ):
        manifest_path = legacy_manifest_path
    _require_within_root(root, manifest_path, "disease-gene bundle manifest")
    signatures[manifest_path] = _file_signature(manifest_path)
    manifest, manifest_bytes = load_bounded_json(manifest_path)
    if hashlib.sha256(manifest_bytes).hexdigest() != expected_bundle_sha256:
        raise DiseaseGeneEvidenceError("disease-gene bundle checksum mismatch")
    if not isinstance(manifest, Mapping):
        raise DiseaseGeneEvidenceError("public v6 bundle manifest must be an object")
    identity = {
        "schema": BUNDLE_MANIFEST_SCHEMA, "classification": "public", "bundle_id": contract.bundle_id,
        "projection_id": contract.projection_id, "projection_schema": "v6",
        "contract_sha256": contract_sha256, "assembly": contract.assembly,
    }
    if any(manifest.get(key) != value for key, value in identity.items()):
        raise DiseaseGeneEvidenceError("bundle manifest is not bound to the public v6 source contract")
    if not re.fullmatch(r"sha256:[0-9a-f]{64}", manifest["image_digest"]):
        raise DiseaseGeneEvidenceError("bundle manifest image identity is invalid")
    date.fromisoformat(manifest["retrieved_on"])
    upstream = manifest["upstream"]
    expected_upstream = {
        "hpo_release": contract.hpo.release, "hpo_ontology_sha256": contract.hpo.ontology_sha256,
        "hpo_annotations_sha256": contract.hpo.annotations_sha256,
        "open_targets_release": contract.open_targets.release,
        "open_targets_croissant_sha256": contract.open_targets.croissant_sha256,
        "open_targets_manifest_path": contract.open_targets.manifest_path,
        "open_targets_manifest_sha256": upstream["open_targets_manifest_sha256"],
    }
    if upstream != expected_upstream or not _SHA256.fullmatch(upstream["open_targets_manifest_sha256"]):
        raise DiseaseGeneEvidenceError("bundle upstream identities differ from the source pins")
    if len(manifest["sources"]) != len(contract.artifacts):
        raise DiseaseGeneEvidenceError("bundle source manifest index is incomplete")
    for source, artifact in zip(manifest["sources"], contract.artifacts, strict=True):
        # The completed bundle binds source provenance. Upstream manifests are
        # not mounted by this reader, so validate their two producer naming
        # forms from the contract rather than their current filesystem presence.
        legacy_source = _source_manifest_path(root, artifact)
        versioned_source = legacy_source.with_name(
            f"{artifact.resource_id}.{_artifact_contract_sha256(artifact)[:12]}.json")
        if source["manifest_path"] not in {
            str(path.relative_to(root)) for path in (legacy_source, versioned_source)
        }:
            raise DiseaseGeneEvidenceError("bundle source manifest path differs from the v6 source contract")
        expected = {
            "role": artifact.role, "resource_id": artifact.resource_id, "release": artifact.release,
            "extraction_date": artifact.extraction_date, "source_revision": artifact.source_revision,
            "source_url": artifact.url, "source_sha256": artifact.sha256,
            "license": artifact.license, "license_url": artifact.license_url,
            "manifest_path": source["manifest_path"],
            "manifest_sha256": source["manifest_sha256"],
        }
        if source != expected or not _SHA256.fullmatch(source["manifest_sha256"]):
            raise DiseaseGeneEvidenceError("bundle source provenance differs from the v6 source contract")
    outputs = {}
    for output in manifest["outputs"]:
        table = output["table"]
        if table not in contract.projection.expected_row_counts or table in outputs:
            raise DiseaseGeneEvidenceError("bundle contains unknown or duplicate tables")
        if (
            set(output) != {"table", "path", "size_bytes", "sha256", "row_count", "retention"}
            or output["path"] != f"{contract.projection.output_prefix}/{table}.parquet"
            or type(output["row_count"]) is not int or output["row_count"] != contract.projection.expected_row_counts[table]
            or type(output["size_bytes"]) is not int or output["size_bytes"] <= 0
            or not _SHA256.fullmatch(output["sha256"]) or output["retention"] != "durable"
        ):
            raise DiseaseGeneEvidenceError("bundle output checksum table differs from its contract")
        outputs[table] = output
    if set(outputs) != set(contract.projection.expected_row_counts):
        raise DiseaseGeneEvidenceError("bundle output checksum table is incomplete")
    validation_counts = manifest["validation"]["row_counts"]
    if (not isinstance(validation_counts, Mapping) or any(type(value) is not int for value in validation_counts.values())
            or validation_counts != dict(contract.projection.expected_row_counts)):
        raise DiseaseGeneEvidenceError("bundle validation counts differ from the source contract")
    files = [(root / outputs[table]["path"], outputs[table]["sha256"], outputs[table]["size_bytes"])
             for table in PUBLIC_TABLES]
    for path, digest, size in files:
        _require_within_root(root, path, "public disease-gene input")
        if any(character in str(path) for character in "*?[]{}"):
            raise DiseaseGeneEvidenceError("public input paths cannot contain glob metacharacters")
        signatures[path] = _file_signature(path)
        if signatures[path][2] != size:
            raise DiseaseGeneEvidenceError("public disease-gene input size differs from completed checksum metadata")
    binding = {
        **identity, "bundle_manifest_path": str(manifest_path.relative_to(root)),
        "bundle_manifest_sha256": expected_bundle_sha256, "sources": manifest["sources"], "upstream": upstream,
        "image_digest": manifest["image_digest"], "retrieved_on": manifest["retrieved_on"],
        "outputs": {table: outputs[table] for table in PUBLIC_TABLES},
    }
    releases = {
        "orphadata_product6": contract.artifact("orphadata_genes").release,
        "gencc_submissions": contract.artifact("gencc_submissions").release,
        **dict.fromkeys(OT_EVIDENCE_DATASETS, contract.open_targets.release),
        "association_by_datasource_direct": contract.open_targets.release,
        "orphadata_product9_ages": contract.artifact("orphadata_natural_history").release,
        "orpha_nomenclature": contract.artifact("orpha_nomenclature").release,
        "orphadata_product1": contract.artifact("orphadata_diseases").release,
    }
    datasets = {
        "gene_assertion": {"orphadata_product6", "gencc_submissions", "association_by_datasource_direct", *OT_EVIDENCE_DATASETS},
        "disease_inheritance": {"orphadata_product9_ages"}, "disease": {"orpha_nomenclature"},
        "disease_alias": {"orphadata_product1"}, "disease_relation": {"orpha_nomenclature"},
    }
    import duckdb

    counts = defaultdict(int)
    tables = {}
    try:
        with duckdb.connect(":memory:", config={
            "threads": "1", "memory_limit": "512MB", "temp_directory": "",
            "autoinstall_known_extensions": "false", "autoload_known_extensions": "false",
        }) as connection:
            ordered=list(zip(PUBLIC_TABLES,files,strict=True))
            if diagnosis is not None:
                # Resolve the existing exact/subtype rules from the small
                # metadata tables before materializing gene assertions.
                ordered.sort(key=lambda item:item[0]=='gene_assertion')
            for table, (path, _, _) in ordered:
                sql="SELECT * FROM read_parquet(?, hive_partitioning = false)"
                parameters=[str(path)]
                if diagnosis is not None and table=='gene_assertion':
                    from genosapiens.track2_diagnosis_scope import resolve_diagnosis
                    metadata=build_disease_gene_evidence([],source_binding=binding,
                        inheritance_rows=tables['disease_inheritance'],disease_rows=tables['disease'],
                        alias_rows=tables['disease_alias'],relation_rows=tables['disease_relation'])
                    if track2_concepts:
                        from genosapiens.track2_concept_supplement import supplement_store
                        metadata=supplement_store(metadata)
                    _,concepts,_=resolve_diagnosis(metadata,diagnosis)
                    concepts=sorted(concepts)
                    # Source and mapped views are retained. The existing scope
                    # predicate still decides scientific applicability.
                    placeholders=','.join('?' for _ in concepts)
                    sql+=f' WHERE source_disease_id IN ({placeholders}) OR mapped_disease_id IN ({placeholders})'
                    parameters.extend([*concepts,*concepts])
                    binding['diagnosis_selection']={'classification':'controlled','diagnosis':diagnosis,'concept_ids':concepts,
                        'coverage':'all_source_or_mapped_assertions_for_resolved_concepts'}
                    del metadata
                    # The completed producer owns whole-table qualification.
                    # This count uses Parquet metadata, not full row readback.
                    count=connection.execute('SELECT count(*) FROM read_parquet(?)',[str(path)]).fetchone()[0]
                    if count!=outputs[table]['row_count']:
                        raise DiseaseGeneEvidenceError('Parquet count differs from the completed source contract')
                cursor = connection.execute(sql, parameters)
                expected_schema = tuple(tuple(part.strip().replace('"', '').split())
                                        for part in TABLE_SCHEMAS[table].split(","))
                if tuple((item[0], str(item[1])) for item in cursor.description) != expected_schema:
                    raise DiseaseGeneEvidenceError(f"{table} Parquet schema differs from v6")
                tables[table] = []
                while batch := cursor.fetchmany(2048):
                    for values in batch:
                        row = dict(zip(TABLE_COLUMNS[table], values, strict=True))
                        dataset = row.get("source_dataset", "orpha_nomenclature")
                        if dataset not in datasets[table] or row["source_release"] != releases[dataset]:
                            raise DiseaseGeneEvidenceError("Parquet row source/release is not bound to its table")
                        counts[(table, dataset)] += 1
                        tables[table].append(row)
                        if diagnosis is not None and table=='gene_assertion' and len(tables[table])>100000:
                            raise DiseaseGeneEvidenceError('diagnosis gene evidence exceeds its bounded input scope')
                if not (diagnosis is not None and table=='gene_assertion') and len(tables[table]) != outputs[table]["row_count"]:
                    raise DiseaseGeneEvidenceError("Parquet row counts differ from the source contract")
    except duckdb.Error:
        raise DiseaseGeneEvidenceError("cannot read the pinned v6 disease-gene Parquets") from None
    expected_counts = {
        "gencc_submissions": contract.artifact("gencc_submissions").expected_records,
        **{name: contract.open_targets.dataset_rows[name] for name in OT_EVIDENCE_DATASETS},
        "association_by_datasource_direct": contract.validation_expectations["open_targets_gene_aggregates"],
    }
    if diagnosis is None and any(counts[("gene_assertion", name)] != count for name, count in expected_counts.items()):
        raise DiseaseGeneEvidenceError("gene source row counts differ from the source contract")
    if diagnosis is None and sum(counts[("gene_assertion", name)] for name in OT_EVIDENCE_DATASETS) != contract.validation_expectations["open_targets_gene_evidence"]:
        raise DiseaseGeneEvidenceError("Open Targets gene evidence count differs from the source contract")
    store = build_disease_gene_evidence(
        tables["gene_assertion"], source_binding=binding, inheritance_rows=tables["disease_inheritance"],
        disease_rows=tables["disease"], alias_rows=tables["disease_alias"], relation_rows=tables["disease_relation"],
    )
    if any(_file_signature(path) != before for path, before in signatures.items()):
        raise DiseaseGeneEvidenceError("disease-gene input changed during the read-only build")
    if track2_concepts:
        from genosapiens.track2_concept_supplement import supplement_store
        store=supplement_store(store)
    return store
