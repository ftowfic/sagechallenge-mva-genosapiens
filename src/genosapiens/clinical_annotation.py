"""Offline, exact-allele clinical annotation; no scoring or frequency decisions.

Only source replay authenticates an assertion. Serialized reports and the model's
explicit_condition_link flag are not authentication tokens. Task 7.6 owns the
source readers; this module adds REF/normalization gates and a lossless projection.
"""

from __future__ import annotations

from collections.abc import Mapping
from contextlib import ExitStack
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Literal

from genosapiens import admission_links as links, clinvar as cv, open_targets as ot
from genosapiens.admission_links import PinnedPublicSource
from genosapiens.evidence_model import (
    AlleleEvidence, ClinicalEvidence, EvidenceReference, IdentityStage,
    PublicSourcedEvidence, SourcedRecord, _Controlled, _typed,
)
from genosapiens.evidence_serialization import _canonical, _loads, _plain
from genosapiens.reference_concordance import (
    _OpenIndexedFasta, _parse_fai_bytes, load_reference_identity,
)
from genosapiens.variant_normalization import (
    SourceAllele, SourceSmallVariantRecord, _normalize_sequence_allele,
    source_allele_id, source_record_id,
)
from genosapiens.vep_annotation import _sequence


SCHEMA = "genosapiens.clinical-annotation"
_DATASETS = {"ClinVar VCF", "Open Targets variant", "Open Targets evidence_eva"}
_CONFLICT_MARKERS = {
    "ClinVar VCF": {
        "CLNSIG": ("conflicting_classifications_of_pathogenicity", "conflicting_interpretations_of_pathogenicity"),
        "CLNREVSTAT": ("criteria_provided,_conflicting_classifications", "criteria_provided,_conflicting_interpretations"),
    },
    "Open Targets evidence_eva": {
        "clinicalSignificances": ("conflicting classifications of pathogenicity", "conflicting interpretations of pathogenicity"),
        "confidence": ("criteria provided, conflicting classifications", "criteria provided, conflicting interpretations"),
        "qualityControls": ("conflicting_clinical_significance",),
    },
}


def affirmative_source_conflict(dataset, fields):
    """Match supported affirmative fields exactly, preserving comma-bearing labels."""
    affirmative = dataset == "ClinVar VCF" and fields.get("CLNSIGCONF") not in (None, ".", "")
    for field, vocabulary in _CONFLICT_MARKERS.get(dataset, {}).items():
        values = fields.get(field)
        if type(values) is str:
            values = values.split("|") if field == "CLNSIG" else (values,)
        affirmative |= any(value is not None and value.lower() in vocabulary for value in values or ())
    return affirmative

_RULES = {
    "identity": "native-sequence-must-be-normalized-on-its-own-reference; no-proxy-or-repair",
    "lookup": "declared-coordinate-alias; public-native-REF-check; exact-REF-and-ALT",
    "clinvar": "2026-08-22-VCF-aggregates-only; no-condition-significance-product",
    "open_targets": "26.06-variant-explicit-fields-to-variantId-to-individual-evidence_eva-row",
    "condition": "diseaseFromSourceId; mapped-disease-is-not-equivalence",
    "date": "EVA-evidenceDate-verbatim; no-release-date-as-evaluation-date",
    "conflict": "affirmative-source-field-vocabulary-or-pathogenic/benign-at-identical-native-condition-and-target",
    "conflict_markers": _CONFLICT_MARKERS,
    "unknown": "no-assertion-is-unknown-not-benign; missing-source-is-not-negative",
    "scope": "raw-annotation-only; no-admission-score-frequency-rescue-or-eligibility-decision",
}


class ClinicalAnnotationError(ValueError):
    """Static diagnostics without source rows, paths, SQL, or chained exceptions."""


@dataclass(frozen=True, slots=True, repr=False)
class ClinicalCoverage(_Controlled):
    allele_id: str
    native_identity_id: str
    source: Literal["clinvar", "open_targets"]
    status: Literal["source_unavailable", "unknown", "aggregate_only", "assertions"]
    record_refs: tuple[EvidenceReference, ...]
    assertion_refs: tuple[EvidenceReference, ...]
    conflict_refs: tuple[EvidenceReference, ...]
    gaps: tuple[str, ...]


@dataclass(frozen=True, slots=True, repr=False)
class ClinicalAnnotationResult(_Controlled):
    alleles: tuple[AlleleEvidence, ...]
    coverage: tuple[ClinicalCoverage, ...]
    bindings: SourcedRecord
    public_sourced_evidence: tuple[PublicSourcedEvidence, ...]


def _reference_checks(files, alleles, pin):
    maps, _, provenance = links._reference(files, pin)
    _, manifest = files.source(pin)
    contract = load_reference_identity(pin.contract_path)
    outputs = {o.path: o for o in manifest.outputs}
    roles = {"public_annotation_benchmark"}
    for allele in alleles:
        seq = _sequence(links._native(allele).identity)
        if seq is not None:
            roles.add(seq.reference_role.value)
    references = {}
    for name in sorted(roles):
        role = contract.reference_role(name)
        fasta = outputs[contract.expected_output(role.fasta_output_role).path]
        fai = outputs[contract.expected_output(role.fai_output_role).path]
        stream = files.open(pin.root / fasta.path, fasta.sha256, fasta.size_bytes)
        index = files.open(pin.root / fai.path, fai.sha256, fai.size_bytes)
        references[name] = _OpenIndexedFasta(stream, _parse_fai_bytes(index.read()))
    public = references["public_annotation_benchmark"]
    gaps = {}
    for allele in alleles:
        native = links._native(allele)
        identity = native.identity
        seq = _sequence(identity)
        reference_identity = seq if seq is not None else identity
        # Even an unsupported native class must bind the declared reference.
        if (native.reference_manifest_sha256 != pin.manifest_sha256
                or contract.reference_role(reference_identity.reference_role.value).identity_id != reference_identity.reference_identity_id):
            raise ValueError
        if seq is None:
            gaps[allele.allele_id] = "native-class-without-exact-clinical-adapter; no-proxy"
            continue
        for item in allele.identities:
            if item.stage is IdentityStage.NORMALIZED and (
                _sequence(item.identity) != seq or item.reference_manifest_sha256 != native.reference_manifest_sha256
            ):
                raise ValueError
        own = references[seq.reference_role.value]
        own_names = (
            maps[(seq.reference_role.value, seq.reference_identity_id)][1]
            if seq.reference_role.value == "public_annotation_benchmark" else {}
        )

        def fetch_own(contig, position, length):
            return own.fetch(own_names.get(contig, contig), position, length)

        if fetch_own(seq.contig, seq.position, len(seq.reference)) != seq.reference:
            raise ValueError
        # Reuse the normalization primitive as a check, without manufacturing a
        # case concordance report or publishing a replacement identity.
        record = SourceSmallVariantRecord(1, seq.contig, seq.position, seq.reference, (seq.alternate,))
        record_id = source_record_id(native.source.source_sha256, record)
        checked = _normalize_sequence_allele(
            SourceAllele(source_allele_id(record_id, 1), record_id, native.source.source_sha256,
                         record, 1, seq.alternate), reference_identity_id=seq.reference_identity_id,
            reference_role=seq.reference_role.value, fetch_reference=fetch_own,
        )
        if (checked.position, checked.reference, checked.alternate) != (seq.position, seq.reference, seq.alternate):
            raise ValueError
        key, reason = links._sequence(allele, maps, public)
        if key is None:
            gaps[allele.allele_id] = reason
        elif not links._DNA.fullmatch(seq.reference) or not links._DNA.fullmatch(seq.alternate):
            gaps[allele.allele_id] = "ambiguous-sequence-exact-link-unresolved"
        elif public.fetch(key[0], key[1], len(key[2])) != key[2]:
            raise ValueError
    return provenance, gaps


def _ref(record):
    manifest = record.public_manifest_sha256
    return EvidenceReference("clinical-source-" + manifest, record.source_record_id, manifest)


def _projection(record):
    row = record.source_fields["record"]
    if record.source_dataset == "ClinVar VCF":
        info = row["info"]
        # Retain delimiter spelling and group structure; do not infer assertions.
        return ClinicalEvidence(
            record, condition_ids=None if info.get("CLNDISDB") is None else (info["CLNDISDB"],),
            significance=None if info.get("CLNSIG") is None else (info["CLNSIG"],),
            review_status=info.get("CLNREVSTAT"), explicit_condition_link=False,
        )
    if record.source_dataset == "Open Targets variant":
        return ClinicalEvidence(record)
    target, condition = links.parse_open_targets_eva_row(row)
    return ClinicalEvidence(
        record, assertion_id=row["id"], condition_ids=None if condition is None else (condition,),
        significance=row["clinicalSignificances"], review_status=row["confidence"],
        explicit_condition_link=condition is not None, target_id=target,
        confidence=row["confidence"], evaluated_on=row["evidenceDate"],
    )


def _polarities(values):
    result = set()
    for value in values or ():
        if value is None:
            continue
        value = value.lower().replace("_", " ")
        if value in ("pathogenic", "likely pathogenic", "pathogenic/likely pathogenic"):
            result.add("pathogenic")
        if value in ("benign", "likely benign", "benign/likely benign"):
            result.add("benign")
    return result


def annotate_clinical(
    alleles: tuple[AlleleEvidence, ...], *, clinvar: PinnedPublicSource | None,
    open_targets: PinnedPublicSource | None, reference: PinnedPublicSource,
) -> ClinicalAnnotationResult:
    """Attach authenticated raw evidence plus per-source, per-allele coverage.

    Source selection is explicit; None produces source_unavailable coverage.
    Existing clinical entries must be exactly reproducible from these sources;
    they are never trusted, dropped, or silently overwritten. Other evidence is
    unchanged. Both source iterators are exhausted before any result is returned.
    """
    try:
        alleles = links._inputs(alleles)
        with ExitStack() as stack:
            files = links._Files(stack)
            ref_binding, identity_gaps = _reference_checks(files, alleles, reference)
            records = {a.allele_id: {} for a in alleles}
            gaps = {(a.allele_id, name): [] for a in alleles for name in ("clinvar", "open_targets")}
            provenance = {}
            pins = {}
            for name, pin, reader in (
                ("clinvar", clinvar, links.clinvar_clinical_report),
                ("open_targets", open_targets, links.open_targets_clinical_report),
            ):
                if pin is None:
                    pins[name] = None
                    continue
                data, manifest = files.source(pin)
                pins[name] = {"contract_sha256": pin.contract_sha256, "manifest_sha256": pin.manifest_sha256}
                # Holding every allowed artifact also protects the projection
                # after each reused iterator's own stable-file scope has closed.
                if name == "clinvar":
                    contract = cv.ClinvarContract.from_mapping(data)
                    cv._verify_manifest_contract(manifest, contract, pin.contract_sha256)
                    allowed = {a.path for a in contract.artifacts}
                else:
                    contract = ot.OpenTargetsContract.from_mapping(data)
                    prefix = Path(contract.output_prefix)
                    allowed = {str(prefix / m.name) for m in contract.metadata}
                    allowed.update(o.path for o in manifest.outputs if Path(o.path).parent in (
                        prefix / "output/variant", prefix / "output/evidence_eva",
                    ))
                for output in manifest.outputs:
                    if output.path in allowed:
                        files.open(pin.root / output.path, output.sha256, output.size_bytes)
                for proof in reader(alleles, source=pin, reference=reference):
                    provenance[name] = proof.provenance
                    gaps[proof.allele_id, name].extend(proof.reasons)
                    for record in proof.source_refs:
                        if record.source_dataset not in _DATASETS:
                            raise ValueError
                        key = (record.public_manifest_sha256, record.source_record_id)
                        previous = records[proof.allele_id].setdefault(key, record)
                        if previous != record:
                            raise ValueError
            attached, coverage, stores = [], [], {}
            for allele in alleles:
                sources = sorted(records[allele.allele_id].values(), key=lambda r: (r.source_dataset, r.source_record_id))
                clinical = [_projection(record) for record in sources]
                conflicts = {_ref(c.source): set() for c in clinical}
                for item in clinical:
                    row = item.source.source_fields["record"]
                    dataset = item.source.source_dataset
                    fields = row["info"] if dataset == "ClinVar VCF" else row
                    affirmative = affirmative_source_conflict(dataset, fields)
                    if affirmative or len(_polarities(item.significance)) == 2:
                        conflicts[_ref(item.source)].add(_ref(item.source))
                linked = [c for c in clinical if c.explicit_condition_link]
                for index, left in enumerate(linked):
                    for right in linked[index + 1:]:
                        if (left.condition_ids, left.target_id) != (right.condition_ids, right.target_id):
                            continue
                        left_signs, right_signs = _polarities(left.significance), _polarities(right.significance)
                        if ("pathogenic" in left_signs and "benign" in right_signs
                                or "benign" in left_signs and "pathogenic" in right_signs):
                            conflicts[_ref(left.source)].add(_ref(right.source))
                            conflicts[_ref(right.source)].add(_ref(left.source))
                clinical = tuple(replace(c, conflict_refs=tuple(sorted(
                    conflicts[_ref(c.source)], key=lambda r: (r.evidence_id, r.record_id),
                ))) for c in clinical)
                if any(c not in clinical for c in allele.clinical or ()):
                    raise ValueError
                attached.append(replace(allele, clinical=clinical or (None if allele.clinical is None else ())))
                for name, pin in (("clinvar", clinvar), ("open_targets", open_targets)):
                    selected = tuple(c for c in clinical if c.source.source_dataset.startswith(
                        "ClinVar" if name == "clinvar" else "Open Targets",
                    ))
                    assertions = tuple(_ref(c.source) for c in selected if c.assertion_id is not None)
                    reasons = []
                    if pin is None:
                        status = "source_unavailable"
                        reasons.append("source-not-selected; clinical-significance-unknown")
                    elif assertions:
                        status = "assertions"
                    elif selected and name == "clinvar":
                        status = "aggregate_only"
                        reasons.append("aggregate-condition-assertion-link-unresolved")
                    else:
                        status = "unknown"
                        reasons.append("no-matched-source-assertion; clinical-significance-unknown")
                    if allele.allele_id in identity_gaps:
                        reasons.append(identity_gaps[allele.allele_id])
                    for c in selected:
                        if c.assertion_id is not None:
                            if not c.explicit_condition_link:
                                reasons.append("native-condition-unresolved; mapped-disease-not-equivalence")
                            if c.target_id is None:
                                reasons.append("source-target-unresolved; native-condition-retained")
                            if not c.significance or all(v is None for v in c.significance):
                                reasons.append("source-significance-unknown")
                            if c.review_status is None:
                                reasons.append("source-review-status-unknown")
                            if c.evaluated_on is None:
                                reasons.append("source-evaluation-date-unknown")
                    if not selected:
                        reasons.extend(gaps[allele.allele_id, name])
                    conflict_refs = {_ref(c.source) for c in selected if c.conflict_refs}
                    coverage.append(ClinicalCoverage(
                        allele.allele_id, allele.native_identity_id, name, status,
                        tuple(_ref(c.source) for c in selected), assertions,
                        tuple(sorted(conflict_refs, key=lambda r: (r.evidence_id, r.record_id))),
                        tuple(sorted(set(reasons))),
                    ))
                for record in sources:
                    manifest = record.public_manifest_sha256
                    previous = stores.setdefault(manifest, {}).setdefault(record.source_record_id, record)
                    if previous != record:
                        raise ValueError
            binding = {
                "schema": SCHEMA, "rules": _RULES, "sources": pins, "reference": ref_binding,
                "source_replay": provenance,
                "implementations": {Path(path).name: links._open_file_sha256(files.open(Path(path)))
                                    for path in (__file__, links.__file__)},
                "inputs": [{"allele_id": a.allele_id, "native": links._native(a),
                            "normalized": tuple(i for i in a.identities if i.stage is IdentityStage.NORMALIZED)}
                           for a in alleles],
            }
            digest = links._digest(binding)
            result = ClinicalAnnotationResult(
                tuple(attached), tuple(coverage), SourcedRecord(SCHEMA, "current", digest, digest, None, _plain(binding)),
                tuple(PublicSourcedEvidence("clinical-source-" + pin, pin, rows) for pin, rows in sorted(stores.items())),
            )
        return result
    except Exception:
        raise ClinicalAnnotationError("offline clinical annotation failed") from None


def verify_clinical_annotation(
    result: ClinicalAnnotationResult, alleles: tuple[AlleleEvidence, ...], *,
    clinvar: PinnedPublicSource | None, open_targets: PinnedPublicSource | None,
    reference: PinnedPublicSource,
) -> None:
    """Consumer replay of actual bound bytes, not a check of self-declared hashes."""
    try:
        if type(result) is not ClinicalAnnotationResult or result != annotate_clinical(
            alleles, clinvar=clinvar, open_targets=open_targets, reference=reference,
        ):
            raise ValueError
    except Exception:
        raise ClinicalAnnotationError("clinical annotation replay failed") from None


def clinical_annotation_as_json(result: ClinicalAnnotationResult) -> str:
    """Explicit controlled serialization; does not authenticate the result."""
    try:
        value = _typed(result, ClinicalAnnotationResult)
        return _canonical({"schema": SCHEMA, "classification": "controlled", "llm_policy": "no-llm", "result": value}) + "\n"
    except Exception:
        raise ClinicalAnnotationError("cannot encode clinical annotation") from None


def clinical_annotation_from_json(content: str | bytes) -> ClinicalAnnotationResult:
    """Decode structure only; verify_clinical_annotation is required before use."""
    try:
        data = _loads(content)
        if not isinstance(data, Mapping) or set(data) != {"schema", "classification", "llm_policy", "result"} or (
            data["schema"], data["classification"], data["llm_policy"]
        ) != (SCHEMA, "controlled", "no-llm"):
            raise ValueError
        return _typed(data["result"], ClinicalAnnotationResult, decode=True)
    except Exception:
        raise ClinicalAnnotationError("invalid clinical annotation JSON") from None
