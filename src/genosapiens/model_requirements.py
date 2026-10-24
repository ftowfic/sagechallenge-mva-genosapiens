"""Decision 7 property kernel, NOT source authentication or partner qualification.

The consumer must replay context requirements, original allele/call evidence,
VEP/native publications, clinical sources and the independently pinned policy.
Constructors and hashes establish structure/identity only. No I/O, callbacks,
annotation execution, source precedence, genotype enumeration or scoring occurs.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction
from graphlib import TopologicalSorter
import hashlib
import re
from typing import Literal, get_args

from genosapiens.allele_identity import RepeatAllele, StructuralAllele, identity_to_dict
from genosapiens import longtr_admission as longtr
from genosapiens.annotation_admission import INTRAGENIC_TERMS, NON_OVERLAP_TERMS
from genosapiens.clinical_annotation import _projection
from genosapiens.disease_gene_evidence import GeneRequirements
from genosapiens.evidence_model import (
    AlleleEvidence, EvidenceReference, HypothesisContext, PublicSourcedEvidence,
    SourcedRecord, _Controlled, _exact, _ids, _typed,
)
from genosapiens.evidence_serialization import _canonical, _plain
from genosapiens.ontology_bundle import OntologyTerm
from genosapiens.vep_annotation import _value
from genosapiens.vep_native_annotation import native_interval_precision


SCHEMA = "genosapiens.model-requirements"
CURATED_SCHEMA = "genosapiens.curated-allele-property"
State = Literal["supported", "contradicted", "unresolved"]
RequirementField = Literal[
    "inheritance", "allelic_requirement", "required_allele_count",
    "association_requirement", "consequence_requirement", "mechanism_direction", "allele_origin",
]
_CURATED_VALUES = {
    "mechanism_direction": ("loss_of_function", "gain_of_function", "dominant_negative"),
    "transcript_fate": ("nonsense_mediated_decay", "escape_nonsense_mediated_decay"),
    "allele_origin": ("germline", "somatic"),
}
_NATIVE_MEANINGS = {
    "transcript_overlap": ("transcript_overlap", "source-required-transcript-overlap"),
    "transcript_covered": ("transcript_covered", "source-required-full-transcript-coverage"),
    "breakpoint_in_transcript": ("breakpoint_in_transcript", "source-required-breakpoint-in-transcript"),
    "reported_dosage_loss": ("reported_dosage_loss", "source-required-reported-dosage-loss"),
    "reported_dosage_gain": ("reported_dosage_gain", "source-required-reported-dosage-gain"),
    "repeat_size_at_least": ("repeat_size_at_least", "source-required-repeat-size"),
}


class ModelRequirementsError(ValueError):
    """Static diagnostics; malformed bindings are not biological unknowns."""


def _hash(value):
    return hashlib.sha256((_canonical(value) + "\n").encode("ascii")).hexdigest()


@dataclass(frozen=True, slots=True, repr=False)
class RequirementSelector(_Controlled):
    """Match an existing complete source alternative, never add a requirement."""

    requirement: GeneRequirements
    field: RequirementField
    source_refs: tuple[EvidenceReference, ...]

    def __post_init__(self):
        _Controlled.__post_init__(self)
        if (getattr(self.requirement, self.field) is None or not self.requirement.record_ids
                or set(self.requirement.record_ids) != {r.record_id for r in self.source_refs}
                or len(set(self.source_refs)) != len(self.source_refs)):
            raise ModelRequirementsError("invalid requirement selector")


@dataclass(frozen=True, slots=True, repr=False)
class TranscriptScope(_Controlled):
    """Empty IDs means every retained source transcript of the context gene.

    Specified IDs are all required; missing versions never satisfy exact versioned
    selectors. Stable matching is explicit and does not discard raw versions.
    Disagreement between applicable transcripts is unresolved, not worst-effect.
    """

    transcript_ids: tuple[str, ...] = ()
    matching: Literal["exact", "ensembl_stable"] = "exact"

    def __post_init__(self):
        _Controlled.__post_init__(self)
        _ids(self.transcript_ids)
        if any(re.fullmatch(r"ENST[0-9]{11}(?:\.[1-9][0-9]*)?", t) is None for t in self.transcript_ids):
            raise ModelRequirementsError("invalid transcript scope")


@dataclass(frozen=True, slots=True, repr=False)
class SequenceConsequenceRule(_Controlled):
    selector: RequirementSelector
    transcripts: TranscriptScope = TranscriptScope()
    relationship: Literal["exact", "is_a"] = "exact"
    contradicting_terms: tuple[str, ...] = ()

    def __post_init__(self):
        _Controlled.__post_init__(self)
        if self.selector.field != "consequence_requirement":
            raise ModelRequirementsError("predictions cannot evaluate curated requirements")
        _ids(self.contradicting_terms)


@dataclass(frozen=True, slots=True, repr=False)
class NativePropertyRule(_Controlled):
    selector: RequirementSelector
    property: Literal["transcript_overlap", "transcript_covered", "breakpoint_in_transcript", "reported_dosage_loss",
                      "reported_dosage_gain", "repeat_size_at_least"]
    transcripts: TranscriptScope = TranscriptScope()
    threshold: str | None = None
    repeat_unit: Literal["repeat-units", "base-pairs"] | None = None
    repeat_quantifier: Literal["any_reported_size", "all_reported_sizes"] = "any_reported_size"

    def __post_init__(self):
        _Controlled.__post_init__(self)
        if self.selector.field not in ("association_requirement", "consequence_requirement"):
            raise ModelRequirementsError("native properties cannot evaluate curated requirements")
        if self.property == "repeat_size_at_least":
            if self.threshold is None or _exact(self.threshold) < 0 or self.repeat_unit is None:
                raise ModelRequirementsError("explicit repeat threshold and unit required")
        elif self.threshold is not None or self.repeat_unit is not None:
            raise ModelRequirementsError("unexpected native rule parameters")
        if self.property not in ("transcript_overlap", "transcript_covered", "breakpoint_in_transcript") and self.transcripts.transcript_ids:
            raise ModelRequirementsError("measurement is not transcript-specific evidence")


@dataclass(frozen=True, slots=True, repr=False)
class CuratedPropertyRule(_Controlled):
    selector: RequirementSelector
    property: Literal["mechanism_direction", "transcript_fate", "allele_origin"]
    transcripts: TranscriptScope = TranscriptScope()

    def __post_init__(self):
        _Controlled.__post_init__(self)
        field = "consequence_requirement" if self.property == "transcript_fate" else self.property
        if self.selector.field != field:
            raise ModelRequirementsError("curated property and requirement field disagree")


@dataclass(frozen=True, slots=True, repr=False)
class ClinicalApplicabilityRule(_Controlled):
    """Assertion applicability only, never clinical strength or qualification."""

    selector: RequirementSelector

    def __post_init__(self):
        _Controlled.__post_init__(self)
        if self.selector.field != "association_requirement":
            raise ModelRequirementsError("clinical applicability is not molecular mechanism")


RequirementRule = SequenceConsequenceRule | NativePropertyRule | CuratedPropertyRule | ClinicalApplicabilityRule


@dataclass(frozen=True, slots=True, repr=False)
class ModelRequirementPolicy(_Controlled):
    """Closed public policy; public records resolve allele.evidence_refs only.

    source binds independently authored rule/transcript/threshold semantics.
    ontology_terms are the existing parser's records, replayed by the consumer
    against ontology_source; no inferred synonyms, replacements or other edges.
    """

    source: SourcedRecord
    rules: tuple[RequirementRule, ...] = ()
    ontology_terms: tuple[OntologyTerm, ...] = ()
    ontology_source: SourcedRecord | None = None
    public_evidence: tuple[PublicSourcedEvidence, ...] = ()
    authentication: Literal["consumer_must_replay_sources"] = "consumer_must_replay_sources"

    def __post_init__(self):
        _Controlled.__post_init__(self)
        if self.source.public_manifest_sha256 is None or self.source.source_fields is None:
            raise ModelRequirementsError("source-bound requirement policy required")
        if len({_canonical(r.selector) for r in self.rules}) != len(self.rules):
            raise ModelRequirementsError("ambiguous requirement rules")
        if bool(self.ontology_terms) != (self.ontology_source is not None):
            raise ModelRequirementsError("ontology records require a source binding")
        if self.ontology_source is not None and self.ontology_source.public_manifest_sha256 is None:
            raise ModelRequirementsError("ontology manifest required")
        terms = {t.term_id: t for t in self.ontology_terms}
        if len(terms) != len(self.ontology_terms) or any(
            not set(t.parents) <= terms.keys() for t in self.ontology_terms
        ):
            raise ModelRequirementsError("invalid ontology graph")
        try:
            TopologicalSorter({k: t.parents for k, t in terms.items()}).prepare()
        except ValueError:
            raise ModelRequirementsError("invalid ontology graph") from None
        if any(type(r) is SequenceConsequenceRule and r.relationship == "is_a" for r in self.rules) and not terms:
            raise ModelRequirementsError("is_a requires a pinned ontology")
        _ids(tuple(p.evidence_id for p in self.public_evidence))

    @property
    def policy_sha256(self) -> str:
        return _hash({"schema": SCHEMA, "policy": self})


@dataclass(frozen=True, slots=True, repr=False)
class RequirementBasis(_Controlled):
    state: State
    reason: str
    source_refs: tuple[SourcedRecord, ...] = ()
    evidence_refs: tuple[EvidenceReference, ...] = ()
    transcript_id: str | None = None
    retained_values: Mapping[str, object] | None = None


@dataclass(frozen=True, slots=True, repr=False)
class RequirementAssessment(_Controlled):
    context_id: str
    allele_id: str
    native_identity_id: str
    requirement_index: int
    requirement: GeneRequirements
    field: RequirementField
    state: State
    reasons: tuple[str, ...]
    requirement_refs: tuple[EvidenceReference, ...]
    bases: tuple[RequirementBasis, ...]
    policy_sha256: str
    input_sha256: str
    rule: RequirementRule | None
    authentication: Literal["consumer_must_replay_sources"] = "consumer_must_replay_sources"
    qualification: Literal["not_assessed"] = "not_assessed"


def _consensus(states):
    values = set(states)
    return next(iter(values)) if len(values) == 1 else "unresolved"


def _identifier_match(actual, expected, *, stable=False):
    if actual is None or expected is None:
        return None
    if actual == expected:
        return True
    if stable and re.fullmatch(r"ENS[GT][0-9]{11}(?:\.[1-9][0-9]*)?", actual) and re.fullmatch(
        r"ENS[GT][0-9]{11}(?:\.[1-9][0-9]*)?", expected,
    ):
        # Gene requirements specifying a version are not silently weakened.
        return actual.split(".")[0] == expected.split(".")[0]
    return False


def _term(value, policy):
    if not policy.ontology_terms:
        return value if value in INTRAGENIC_TERMS | NON_OVERLAP_TERMS else None
    matches = [t.term_id for t in policy.ontology_terms if not t.obsolete and value in (t.term_id, t.name)]
    return matches[0] if len(matches) == 1 else None


def _matches_term(value, wanted, relationship, policy):
    value, wanted = _term(value, policy), _term(wanted, policy)
    if value is None or wanted is None:
        return None
    parents = {t.term_id: t.parents for t in policy.ontology_terms if not t.obsolete}
    pending, seen = [value], set()
    while pending:
        current = pending.pop()
        if current == wanted:
            return True
        if current not in seen and relationship == "is_a":
            pending.extend(parents.get(current, ()))
        seen.add(current)
    return False


def _transcripts(context, allele, scope, *, native_roles=None):
    """Read original annotation, not enrichment. Return gaps alongside matches."""
    selected, gaps, found = [], [], set()
    for tx in allele.transcript_consequences or ():
        source, fields = tx.source, tx.source.source_fields or {}
        if native_roles is not None and fields.get("role") not in native_roles:
            continue
        if source.source_dataset == "Ensembl VEP" and source.source_release == "116" and "csq" in fields:
            raw = {k: _value(v) for k, v in fields["csq"].items()}
            if raw.get("Feature_type") != "Transcript":
                continue
            gene, transcript = raw.get("Gene"), raw.get("Feature")
            terms = tuple((raw.get("Consequence") or "").split("&"))
        elif source.source_dataset in ("Ensembl VEP native", "Ensembl VEP cache context") and source.source_release == "116":
            raw = fields.get("raw_record", {})
            gene, transcript = raw.get("gene_id"), raw.get("transcript_id")
            terms = raw.get("consequence_terms") if fields.get("annotation_scope") == "structural_effect" else None
        else:
            gaps.append(RequirementBasis("unresolved", "unrecognized-annotation-source", (source,)))
            continue
        gene_match = _identifier_match(gene, context.gene_id, stable=context.gene_id is not None and "." not in context.gene_id)
        if gene_match is False:
            continue
        matches = {t for t in scope.transcript_ids if _identifier_match(transcript, t, stable=scope.matching == "ensembl_stable")}
        if scope.transcript_ids and not matches:
            continue
        if gene_match is None or transcript is None:
            gaps.append(RequirementBasis("unresolved", "source-transcript-or-gene-unavailable", (source,)))
            continue
        found.update(matches)
        selected.append((tx, raw, transcript, terms))
    if set(scope.transcript_ids) - found or not selected:
        gaps.append(RequirementBasis("unresolved", "required-source-transcript-unavailable"))
    return selected, gaps


def _sequence_bases(context, allele, rule, policy):
    selected, bases = _transcripts(context, allele, rule.transcripts)
    wanted = rule.selector.requirement.consequence_requirement
    for tx, _, transcript, terms in selected:
        positive = [_matches_term(t, wanted, rule.relationship, policy) for t in terms or ()]
        negative = [_matches_term(t, n, rule.relationship, policy) for t in terms or () for n in rule.contradicting_terms]
        state = "unresolved"
        if True in positive and True not in negative:
            state = "supported"
        elif True in negative and True not in positive and None not in positive:
            state = "contradicted"
        sources = (tx.source,) + ((policy.ontology_source,) if policy.ontology_source else ())
        bases.append(RequirementBasis(state, "source-prediction-only; no-curated-or-clinical-inference", sources,
            transcript_id=transcript, retained_values={"terms": terms, "relationship": rule.relationship}))
    return bases


def _native_bases(context, allele, rule):
    mapping = next(i for i in allele.identities if i.identity_id == allele.native_identity_id)
    identity = mapping.identity
    sources = (mapping.source,) + tuple(c.source for c in allele.calls or ())
    if rule.property.startswith("reported_dosage"):
        bases = []
        for original in allele.identities:
            if type(original.identity) is not StructuralAllele:
                continue
            dosage, state = original.identity.dosage, "unresolved"
            if dosage is not None and dosage.copy_number is not None and dosage.baseline_copy_number is not None:
                match = (dosage.copy_number < dosage.baseline_copy_number if rule.property == "reported_dosage_loss"
                         else dosage.copy_number > dosage.baseline_copy_number)
                state = "supported" if match else "contradicted"
            bases.append(RequirementBasis(state, "reported-dosage-comparison-only; missing-baseline-unresolved",
                (original.source,), retained_values={"identity_id": original.identity_id, "dosage": _plain(dosage)}))
        return bases or [RequirementBasis("unresolved", "native-dosage-unavailable", sources)]
    if rule.property == "repeat_size_at_least":
        if type(identity) is not RepeatAllele or context.locus_id != identity.locus_id:
            return [RequirementBasis("unresolved", "native-repeat-locus-not-established", sources)]
        bases, threshold = [], _exact(rule.threshold)
        for original in ((mapping,) if longtr.is_projection(allele) else allele.identities):
            if type(original.identity) is not RepeatAllele:
                continue
            sizes = original.identity.allele_sizes or ()
            values = {"identity_id": original.identity_id, "sizes": _plain(sizes), "threshold": rule.threshold,
                      "unit": rule.repeat_unit, "quantifier": rule.repeat_quantifier}
            if (original.identity.locus_id != context.locus_id or not sizes
                    or any(s.lower_bound is None or s.upper_bound is None or s.unit.value != rule.repeat_unit for s in sizes)):
                state = "unresolved"
            else:
                operation = any if rule.repeat_quantifier == "any_reported_size" else all
                lower = operation(Fraction(str(s.lower_bound)) >= threshold for s in sizes)
                upper = operation(Fraction(str(s.upper_bound)) >= threshold for s in sizes)
                state = "supported" if lower else "contradicted" if not upper else "unresolved"
            bases.append(RequirementBasis(state, "all-source-bounded-size-interpretations; no-nominal-bound-imputation",
                                          (original.source,), retained_values=values))
        return bases
    roles = (("breakpoint_uncertainty",) if rule.property == "breakpoint_in_transcript"
             else ("affected_interval", "structural_interval"))
    selected, bases = _transcripts(context, allele, rule.transcripts, native_roles=roles)
    precision = native_interval_precision(allele) if type(identity) is StructuralAllele else None
    seen_points = set()
    for tx, raw, transcript, _ in selected:
        fields = tx.source.source_fields or {}
        if "reference_context" in fields:
            raw = fields["reference_context"]["source_fields"]["raw_record"]
        geometry = raw.get("transcript_geometry")
        bounds, reason = None, "source-native-bounds-or-transcript-geometry-unavailable"
        if geometry is not None and precision is not None:
            if rule.property == "breakpoint_in_transcript" and fields.get("role") == "breakpoint_uncertainty":
                original = fields.get("source_geometry", {})
                for ordinal, point in enumerate(identity.breakpoints):
                    if _canonical(original) == _canonical(point) and point.confidence_interval is not None:
                        lo, hi = (point.position + v for v in point.confidence_interval)
                        bounds = (lo, hi, lo, hi)
                        seen_points.add(ordinal)
            elif rule.property in ("transcript_overlap", "transcript_covered") and fields.get("role") in ("affected_interval", "structural_interval"):
                if _canonical(fields.get("interval_evidence")) != _canonical(precision):
                    raise ModelRequirementsError("native precision binding mismatch")
                contributors = precision.source_fields["contributors"]
                allowed = (None, "source-endpoints-imprecise")
                if contributors and all(c["reason"] in allowed and c["recipe"] is not None for c in contributors):
                    recipes = [c["recipe"] for c in contributors]
                    bounds = (max(r["start"]["min"] for r in recipes), min(r["start"]["max"] for r in recipes),
                              max(r["end"]["min"] for r in recipes), min(r["end"]["max"] for r in recipes))
                    if not (bounds[0] <= identity.interval.start <= bounds[1]
                            and bounds[2] <= identity.interval.end <= bounds[3]):
                        bounds, reason = None, "source-bounds-and-retained-nominal-interval-conflict"
        state = "unresolved"
        if bounds is not None:
            smin, smax, emin, emax = bounds
            if smin <= smax and emin <= emax and smin <= emax and geometry["contig"] == fields["query_geometry"]["contig"]:
                smax, emin = min(smax, emax), max(emin, smin)
                left, right = ((geometry["start"], geometry["end"]) if rule.property == "transcript_covered"
                               else (geometry["end"], geometry["start"]))
                state = ("supported" if smax <= left and emin >= right else
                         "contradicted" if smin > left or emax < right else "unresolved")
                reason = "all-source-bounded-geometry-interpretations; overlap-not-molecular-effect"
        bases.append(RequirementBasis(state, reason, (tx.source, mapping.source), transcript_id=transcript,
            retained_values={"bounds": bounds, "transcript_geometry": geometry, "geometry_role": fields.get("role")}))
    if (rule.property == "breakpoint_in_transcript" and type(identity) is StructuralAllele
            and seen_points != set(range(len(identity.breakpoints)))):
        bases.append(RequirementBasis("unresolved", "required-breakpoint-bounds-or-context-unavailable", sources))
    return bases


def _validate_curated_record(record):
    row=record.source_fields
    keys={'schema', 'record_id', 'native_identity', 'reference_manifest_sha256', 'native_disease_id',
          'gene_id', 'locus_id', 'transcript_id', 'property', 'value', 'negated', 'citations'}
    if (row is None or set(row)!=keys or row['schema']!=CURATED_SCHEMA
            or row['record_id']!=record.source_record_id or _hash(row)!=record.source_sha256
            or type(row['negated']) is not bool or type(row['citations']) is not tuple or not row['citations']
            or any(type(row[k]) is not str for k in ('property','value','native_disease_id'))
            or any(row[k] is not None and type(row[k]) is not str for k in ('gene_id','locus_id','transcript_id'))
            or any(type(v) is not str or not v.strip() for v in row['citations'])):
        raise ModelRequirementsError('invalid curated allele source record')
    return row


def _curated_applicable(context, mapping, row, property_name, transcripts):
    applicable=(_canonical(row['native_identity'])==_canonical(identity_to_dict(mapping.identity))
        and row['reference_manifest_sha256']==mapping.reference_manifest_sha256
        and row['native_disease_id']==(context.native_disease_id or context.disease_id)
        and (row['gene_id'],row['locus_id'])==(context.gene_id,context.locus_id)
        and row['property']==property_name)
    matches={t for t in transcripts.transcript_ids if _identifier_match(
        row['transcript_id'],t,stable=transcripts.matching=='ensembl_stable')}
    return applicable and (not transcripts.transcript_ids or bool(matches)),matches


def frozen_functional_check(context, allele, transcripts, policy):
    """Check every retained functional record, including a valid empty snapshot.

    This query uses the existing native identity and transcript applicability
    contract. Its scope is the supplied frozen sources, not all literature.
    """
    mapping=next(i for i in allele.identities if i.identity_id==allele.native_identity_id)
    sources=[];matches=[];invalid=[]
    for public in policy.public_evidence:
        sources.append({'evidence_id':public.evidence_id,'manifest_sha256':public.public_manifest_sha256,
                        'record_count':len(public.records)})
        for record in public.records.values():
            if record.source_dataset!=CURATED_SCHEMA or record.source_release!='current':
                invalid.append({'record_id':record.source_record_id,'reason':'unsupported_functional_source'})
                continue
            try:
                row=_validate_curated_record(record)
            except ModelRequirementsError:
                invalid.append({'record_id':record.source_record_id,'reason':'invalid_functional_source'})
                continue
            applicable,_=_curated_applicable(context,mapping,row,'mechanism_direction',transcripts)
            if applicable:
                matches.append({'record_id':record.source_record_id,'source_sha256':record.source_sha256,
                    'value':row['value'],'negated':row['negated'],'transcript_id':row['transcript_id'],
                    'citations':list(row['citations']),'evidence_id':public.evidence_id})
    contrary=any((r['value']=='loss_of_function' and r['negated']) or
                 (r['value']!='loss_of_function' and not r['negated']) for r in matches)
    return {'state':'invalid' if invalid else 'complete' if sources else 'unavailable',
            'sources':sources,'matches':matches,'invalid_records':invalid,
            'contrary':contrary,'rule':'exact-native-allele-gene-disease-transcript-and-negation/v1'}


def _curated_bases(context, allele, rule, policy):
    mapping = next(i for i in allele.identities if i.identity_id == allele.native_identity_id)
    wanted = getattr(rule.selector.requirement, rule.selector.field)
    bases, found = [], set()
    for public in policy.public_evidence:
        for ref in allele.evidence_refs:
            if (ref.evidence_id, ref.public_manifest_sha256) != (public.evidence_id, public.public_manifest_sha256):
                continue
            record = public.records[ref.record_id]
            if record.source_dataset != CURATED_SCHEMA or record.source_release != "current":
                continue
            row = _validate_curated_record(record)
            applicable,matches = _curated_applicable(context,mapping,row,rule.property,rule.transcripts)
            if not applicable:
                continue
            found.update(matches)
            state = "unresolved"
            if wanted in _CURATED_VALUES[rule.property] and row["value"] == wanted:
                state = "contradicted" if row["negated"] else "supported"
            bases.append(RequirementBasis(state, "exact-allele-curated-assertion; not-gene-model-or-prediction", (record,),
                                         (ref,), row["transcript_id"], {"property": row["property"], "value": row["value"], "negated": row["negated"]}))
    if not bases or set(rule.transcripts.transcript_ids) - found:
        bases.append(RequirementBasis("unresolved", "applicable-curated-allele-evidence-unavailable"))
    return bases


def _clinical_bases(context, allele):
    bases = []
    for item in allele.clinical or ():
        source = item.source
        releases = {"ClinVar VCF": "2026-08-22", "Open Targets variant": "26.06", "Open Targets evidence_eva": "26.06"}
        if releases.get(source.source_dataset) != source.source_release:
            bases.append(RequirementBasis("unresolved", "unrecognized-clinical-source", (source,)))
            continue
        # Reuse the original projection; optional enriched fields grant no trust.
        original = _projection(source)
        state = "unresolved"
        if original.explicit_condition_link:
            anchor = context.native_disease_id or context.disease_id
            gene = _identifier_match(original.target_id, context.gene_id,
                                     stable=context.gene_id is not None and "." not in context.gene_id)
            row = source.source_fields["record"]
            restricted = any(row.get(field) for field in ("allelicRequirements", "alleleOrigins"))
            if anchor in original.condition_ids and gene is True and context.locus_id is None and not restricted:
                state = "supported"
        bases.append(RequirementBasis(state, "exact-native-condition-target-applicability-only; not-clinical-qualification", (source,),
            retained_values={"condition_ids": original.condition_ids, "target_id": original.target_id,
                             "significance": original.significance, "review_status": original.review_status}))
    return bases or [RequirementBasis("unresolved", "clinical-assertion-unavailable; not-benign")]


def assess_model_requirements(
    context: HypothesisContext, allele: AlleleEvidence, *, policy: ModelRequirementPolicy,
    prepared=None,
) -> tuple[RequirementAssessment, ...]:
    """Assess each existing requirement field, with independent source alternatives.

    Consumer authentication is a precondition, NOT a caller-selectable flag.
    Attach original bound transcript/native results before calling. Unrecognized
    fields/values/sources stay unresolved. Required count and inheritance are
    delegated to genotype/model construction and never inferred from one allele.
    No result, including clinical applicability, is a qualification certificate.
    """
    try:
        from genosapiens.genotype_evaluation import prepare_genotype_evaluation
        if type(policy) is not ModelRequirementPolicy:
            raise ValueError
        if prepared is None:
            prepared = prepare_genotype_evaluation(context=context, alleles=(allele,),
                genotypes=None, gene_evidence=(), policies={'requirements': policy}, pairs=())
        prepared.require(context=context, alleles=(allele,), policies={'requirements': policy})
        key = (prepared.context_key(context), allele.allele_id)
        if key in prepared.requirements:
            return prepared.requirements[key]
        inputs = prepared.context_allele_hash(context, allele)
        policy_hash, result = policy.policy_sha256, []
        for index, requirement in enumerate(context.requirements):
            refs = tuple(r for r in context.evidence_refs if r.record_id in requirement.record_ids)
            for field in get_args(RequirementField):
                rule = next((r for r in policy.rules if r.selector.requirement == requirement and r.selector.field == field
                             and set(r.selector.source_refs) == set(refs)), None)
                if getattr(requirement, field) is None:
                    bases = [RequirementBasis("unresolved", "source-requirement-unspecified")]
                elif field in ("inheritance", "allelic_requirement", "required_allele_count"):
                    bases = [RequirementBasis("unresolved", "genotype-model-property-outside-allele-kernel")]
                elif rule is None:
                    bases = [RequirementBasis("unresolved", "unrecognized-or-unbound-source-requirement")]
                elif type(rule) is SequenceConsequenceRule:
                    bases = _sequence_bases(context, allele, rule, policy)
                elif type(rule) is NativePropertyRule:
                    if getattr(requirement, field) not in _NATIVE_MEANINGS[rule.property]:
                        bases = [RequirementBasis("unresolved", "native-rule-cannot-reinterpret-source-requirement")]
                    else:
                        bases = _native_bases(context, allele, rule)
                elif type(rule) is CuratedPropertyRule:
                    bases = _curated_bases(context, allele, rule, policy)
                else:
                    if getattr(requirement, field) not in ("clinical_assertion_applicability", "source-required-clinical-applicability"):
                        bases = [RequirementBasis("unresolved", "clinical-rule-cannot-reinterpret-source-requirement")]
                    else:
                        bases = _clinical_bases(context, allele)
                state = _consensus(b.state for b in bases)
                if type(rule) is ClinicalApplicabilityRule and any(b.state == "supported" for b in bases):
                    # Applicability asks whether an assertion exists, not whether
                    # every identity/coverage record is itself an assertion.
                    state = "supported"
                reasons = tuple(sorted({b.reason for b in bases}))
                if len({b.state for b in bases}) > 1:
                    reasons += ("mixed-or-unavailable-property-evidence",)
                result.append(RequirementAssessment(context.context_id, allele.allele_id, allele.native_identity_id,
                    index, requirement, field, state, reasons, refs, tuple(bases), policy_hash, inputs, rule))
        result = tuple(result)
        prepared.requirements[key] = result
        prepared.requirement_hashes[id(result)] = (result, _hash(result))
        return result
    except (ValueError, TypeError, KeyError, AttributeError, RuntimeError, OverflowError, RecursionError):
        raise ModelRequirementsError("invalid model requirement inputs or source bindings") from None
