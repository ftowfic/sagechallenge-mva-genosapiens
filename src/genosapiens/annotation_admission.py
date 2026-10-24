"""Pure Decision 6a replay of a pinned, already published VEP annotation run.

This is association evidence, not pathogenicity, disease support, or permission
to enumerate only admitted alleles. The caller owns acquisition and independently
selects the binding pin and complete original ordered requests. No VEP execution,
FASTA lookup, cache deserialization, I/O, or caller callback occurs here.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, replace
import hashlib
import json
import re
from types import SimpleNamespace
from typing import Literal

from genosapiens.allele_identity import ReferenceRole, RepeatAllele, StructuralAllele
from genosapiens.evidence_model import (
    AlleleEvidence, IdentityMapping, SourcedRecord, TranscriptConsequence, _Controlled,
    _validation_session, _construct_owned,
)
from genosapiens.evidence_serialization import _plain, from_json, to_json
from genosapiens.reference_concordance import _parse_fai_bytes
from genosapiens import vep_annotation as vep, vep_native_annotation as native
from genosapiens.vep_bundle import INVENTORY_SCHEMA, VepVersions


SCHEMA = "genosapiens.annotation-admission"
SEQUENCE_RULE = "sequence-transcript-overlap"
NATIVE_RULE = "native-structural-overlap"
# SO terms from the pinned release-116 ensembl-variation vocabulary. These are
# location claims only, not a severity ordering or a functional-effect filter.
INTRAGENIC_TERMS = frozenset({
    "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
    "stop_gained", "frameshift_variant", "stop_lost", "start_lost",
    "transcript_amplification", "feature_elongation", "feature_truncation",
    "inframe_insertion", "inframe_deletion", "missense_variant", "protein_altering_variant",
    "splice_donor_5th_base_variant", "splice_region_variant", "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant", "incomplete_terminal_codon_variant",
    "start_retained_variant", "stop_retained_variant", "synonymous_variant",
    "coding_sequence_variant", "mature_miRNA_variant", "5_prime_UTR_variant",
    "3_prime_UTR_variant", "non_coding_transcript_exon_variant", "intron_variant",
    "NMD_transcript_variant", "non_coding_transcript_variant",
})
NON_OVERLAP_TERMS = frozenset({
    "upstream_gene_variant", "downstream_gene_variant", "intergenic_variant",
    "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant",
    "regulatory_region_ablation", "regulatory_region_amplification", "regulatory_region_variant",
})
POLICY_BYTES = vep._json({
    "schema": SCHEMA, "rules": [SEQUENCE_RULE, NATIVE_RULE],
    "vep_version": "116.2", "cache_version": "116",
    "variation_revision": native.SOURCE_REVISIONS["ensembl-variation"],
    "intragenic_terms": sorted(INTRAGENIC_TERMS), "non_overlap_terms": sorted(NON_OVERLAP_TERMS),
    "unknown_terms": "unresolved; retain alongside known-term outcomes",
    "sequence": "original-transcript-CSQ-Gene; no-canonical-MANE-biotype-filter",
    "native": "replay-native-queries-and-parsers; each-role-and-result-separately",
    "certainty": "affected-interval=definite; other-bounded-structural-geometry=possible",
    "unbounded": "unresolved; unknown-source-semantics-also-unresolved; never-invent-an-interval",
    "scope": "context-only-association; no-mechanism-pathogenicity-or-population-inference",
    "repeat": "replay-assessment-only; separate-repeat-locus-rule-owned-by-consumer",
    "comparison": "source,consequence_id,consequences,variant_class; optional-typed-crossrefs-ignored",
}).encode("ascii")
POLICY_SHA256 = hashlib.sha256(POLICY_BYTES).hexdigest()
_SHA = re.compile(r"[0-9a-f]{64}")
_GENE = re.compile(r"ENSG[0-9]{11}(?:\.[1-9][0-9]*)?")
_SOURCE_KEYS = (
    "runtime_image", "upstream_image_digest", "executable_sha256", "contract_sha256",
    "cache_inventory_sha256", "fasta_sha256", "fai_sha256", "resource_reference_identity_id",
    "configuration_sha256", "input_evidence_sha256", "tool_version",
)


class AnnotationAdmissionError(ValueError):
    """Static diagnostics; malformed/stale evidence is not an unresolved verdict."""


@dataclass(frozen=True, slots=True, repr=False)
class AnnotationAssociationProof(_Controlled):
    """Scientific ID excludes full-run bindings and optional typed enrichment.

    Source records and provenance still retain their exact byte/manifest pins.
    Multiple rows/terms/geometries, including contradictory outcomes, survive.
    Construction alone is not authentication: consumers must replay the batch.
    """

    proof_id: str
    allele_id: str
    native: IdentityMapping
    rule: Literal["sequence-transcript-overlap", "native-structural-overlap"]
    outcome: Literal["admitted", "evaluated_no_match", "unresolved"]
    target_namespace: str | None
    target_id: str | None
    assessment_status: str
    geometry_role: str | None
    overlap_certainty: Literal["definite", "possible", "unknown", "none"]
    reasons: tuple[str, ...]
    source_refs: tuple[SourcedRecord, ...]
    provenance: SourcedRecord
    policy_sha256: str = POLICY_SHA256
    context_only: Literal[True] = True

    @property
    def rule_key(self) -> tuple[str, str]:
        return self.allele_id, self.rule


def _equal(a, b):
    # JSON comparison also distinguishes bool/int and list/null, unlike ==.
    return vep._json(_plain(a)) == vep._json(_plain(b))


def _gene(value):
    if value in (None, "", ".", "-"):
        return None
    if type(value) is not str or _GENE.fullmatch(value) is None:
        raise ValueError
    return value


def _transcript_claims(transcripts):
    if transcripts is None:
        return None
    return [{key: _plain(getattr(t, key)) for key in (
        "consequence_id", "source", "consequences", "variant_class",
    )} for t in transcripts]


@_validation_session()
def sequence_annotation_associations(
    allele: AlleleEvidence,
    *,
    assessment_status: str,
    transcripts: tuple[TranscriptConsequence, ...],
    non_transcripts: tuple[SourcedRecord, ...],
    provenance: SourcedRecord,
) -> tuple[AnnotationAssociationProof, ...]:
    """Apply the existing sequence-overlap rule to one closed VEP response."""

    try:
        if (
            type(allele) is not AlleleEvidence
            or type(assessment_status) is not str
            or not assessment_status
            or type(transcripts) is not tuple
            or any(type(value) is not TranscriptConsequence for value in transcripts)
            or type(non_transcripts) is not tuple
            or any(type(value) is not SourcedRecord for value in non_transcripts)
            or type(provenance) is not SourcedRecord
        ):
            raise ValueError
        selected = tuple(
            value for value in allele.identities
            if value.identity_id == allele.native_identity_id
        )
        if len(selected) != 1 or type(selected[0].identity) in {
            StructuralAllele, RepeatAllele,
        }:
            raise ValueError
        native_identity = selected[0]
        identity = [
            allele.allele_id,
            native_identity.identity_id,
            _plain(native_identity.identity),
            native_identity.reference_manifest_sha256,
        ]
        proofs = []
        for source in (
            *(value.source for value in transcripts),
            *non_transcripts,
        ):
            raw = source.source_fields["csq"]
            values = vep._value(raw["Consequence"]).split("&")
            target = _gene(vep._value(raw["Gene"]))
            transcript = (
                vep._value(raw["Feature"])
                if raw["Feature_type"] == "Transcript" else None
            )
            for ordinal, term in enumerate(values):
                if term in NON_OVERLAP_TERMS:
                    outcome = "evaluated_no_match"
                    reason = "source-consequence-does-not-establish-transcript-overlap"
                elif term not in INTRAGENIC_TERMS:
                    outcome = "unresolved"
                    reason = "consequence-outside-pinned-overlap-vocabulary"
                elif not transcript:
                    outcome = "evaluated_no_match"
                    reason = "source-feature-is-not-a-transcript"
                elif target is None:
                    outcome = "unresolved"
                    reason = "original-source-gene-unavailable"
                else:
                    outcome = "admitted"
                    reason = "source-transcript-intragenic-consequence"
                locator = [
                    source.source_fields.get("vcf_record_id"),
                    source.source_fields.get("csq_ordinal"),
                    None,
                    transcript,
                    ordinal,
                    term,
                ]
                certainty = (
                    "definite" if outcome == "admitted" else
                    "none" if outcome == "evaluated_no_match" else "unknown"
                )
                scientific = [
                    SCHEMA, POLICY_SHA256, identity, SEQUENCE_RULE, outcome,
                    target, reason, None, certainty, locator,
                ]
                proofs.append(_construct_owned(AnnotationAssociationProof,
                    "annotation-association-" + vep._digest(
                        vep._json(_plain(scientific))
                    ),
                    allele.allele_id,
                    native_identity,
                    SEQUENCE_RULE,
                    outcome,
                    "Ensembl" if target else None,
                    target,
                    assessment_status,
                    None,
                    certainty,
                    (reason, "context-only"),
                    (source,),
                    provenance,
                ))
        return tuple(sorted(proofs, key=lambda value: value.proof_id))
    except (ValueError, TypeError, KeyError, AttributeError, RuntimeError):
        raise AnnotationAdmissionError(
            "invalid or incomplete sequence annotation admission inputs"
        ) from None


@dataclass(frozen=True, slots=True, repr=False)
class ValidatedAnnotationAdmission:
    """Complete source replay with association construction deferred."""

    alleles: tuple[AlleleEvidence, ...]
    molecular_sources: tuple[tuple[SourcedRecord, ...], ...]
    _sequence_plans: tuple
    _native_plans: tuple

    @_validation_session()
    def proofs(self) -> tuple[AnnotationAssociationProof, ...]:
        proofs = [_construct_owned(AnnotationAssociationProof, *plan)
                  for plan in self._native_plans]
        for allele, status, transcripts, non_transcripts, provenance in self._sequence_plans:
            proofs.extend(sequence_annotation_associations(allele,
                assessment_status=status, transcripts=transcripts,
                non_transcripts=non_transcripts, provenance=provenance))
        return tuple(sorted(proofs, key=lambda p: p.proof_id))


def replay_annotation_admission(artifacts: Mapping[str, bytes], *, expected_binding_sha256: str,
        requests: tuple[vep.VepRequest, ...]) -> tuple[AnnotationAssociationProof, ...]:
    return validate_annotation_admission(artifacts,
        expected_binding_sha256=expected_binding_sha256, requests=requests).proofs()


def validate_annotation_admission(
    artifacts: Mapping[str, bytes], *, expected_binding_sha256: str,
    requests: tuple[vep.VepRequest, ...],
) -> ValidatedAnnotationAdmission:
    return _validate_annotation_admission(artifacts,
        expected_binding_sha256=expected_binding_sha256, requests=requests)


def validate_annotation_admission_files(artifacts, *, expected_binding_sha256, requests):
    """Track 2 file input. The byte-only admission interface remains pure."""
    from genosapiens.track2_annotation_artifacts import AlleleArtifact
    if not isinstance(artifacts, Mapping) or type(artifacts.get('alleles.json')) is not AlleleArtifact:
        raise AnnotationAdmissionError('A file annotation input requires an allele artifact.')
    return _validate_annotation_admission(artifacts,
        expected_binding_sha256=expected_binding_sha256, requests=requests, file_input=True)


def _validate_annotation_admission(artifacts, *, expected_binding_sha256, requests, file_input=False):
    """Replay relevant bytes from run_vep_annotation, never arbitrary assessments.

    Required keys: binding.json, alleles.json, transport.vcf, reference.fai (the
    VEP FASTA's actual index); annotated.vcf when output_sha256 is non-null.
    Native batches additionally require cache.inventory.json and every executed
    native/context input/output. Logs and native-cache-context.pl are optional;
    if supplied they must match their publication pins. Unknown keys fail.

    Reference naming is taken from the independently pinned publication's
    verified compatibility proof, with lengths and cache entries replayed from
    hash-bound bytes. This does not reverify the FASTA/cache installation or
    assert case-reference acceptance. Original request order and full input
    evidence must match input_evidence_sha256. Published optional typed mappings
    may differ; original bound source fields and scientific claims may not.

    Repeat assessments are also replayed for closed-batch membership, but do not
    acquire a structural or sequence rule: their separate locus adapter owns
    admission. Every applicable request emits a rule key, even with no target.
    """
    try:
        if (not isinstance(artifacts, Mapping) or type(requests) is not tuple
                or any(type(r) is not vep.VepRequest for r in requests)
                or len({r.allele.allele_id for r in requests}) != len(requests)
                or any(type(k) is not str or type(v) is not bytes
                    for k, v in artifacts.items() if not (file_input and k == 'alleles.json'))):
            raise ValueError
        artifacts = dict(artifacts)
        used = set()

        def content(name, digest, *, decode=True):
            if type(digest) is not str or _SHA.fullmatch(digest) is None:
                raise ValueError
            value = artifacts[name]
            if hashlib.sha256(value).hexdigest() != digest:
                raise ValueError
            used.add(name)
            return value.decode("utf-8", errors="strict") if decode else value

        binding = native._strict_json(content("binding.json", expected_binding_sha256))
        required = {*_SOURCE_KEYS, "schema", "classification", "llm_policy", "configuration", "bundle_manifest_sha256",
                    "assessments", "transport_sha256", "runtime_version_probe_sha256", "output_sha256", "csq_fields",
                    "constraint_source", "alleles_sha256"}
        if (not required <= set(binding) or set(binding) - required - {"native_annotation", "vep.stdout_sha256", "vep.stderr_sha256"}
                or binding["schema"] != vep.SCHEMA or binding["classification"] != "controlled"
                or binding["llm_policy"] != "no-llm" or binding["tool_version"] != "116.2"
                or re.fullmatch(r"[a-z0-9][a-z0-9._-]*", binding["resource_reference_identity_id"]) is None
                or re.fullmatch(r"[A-Za-z0-9_./:-]+@sha256:[0-9a-f]{64}", binding["runtime_image"]) is None
                or re.fullmatch(r"sha256:[0-9a-f]{64}", binding["upstream_image_digest"]) is None):
            raise ValueError
        for key in (*_SOURCE_KEYS, "bundle_manifest_sha256", "runtime_version_probe_sha256"):
            if key.endswith("sha256") and _SHA.fullmatch(binding[key]) is None:
                raise ValueError
        configuration = binding["configuration"]
        if (set(configuration) - {"bounded_execution_sha256"} != {"schema", "flags", "buffer_size", "timeout_seconds", "reference_compatibility",
                                   "constraint_adapter", "constraint_manifest_sha256", "native_backend", "resource_reference_role"}
                or binding["configuration_sha256"] != vep._digest(vep._json(configuration))
                or configuration["schema"] != vep.SCHEMA or configuration["flags"] != list(vep.FLAGS)
                or not _equal(configuration["native_backend"], native.configuration())
                or configuration["resource_reference_role"] != ReferenceRole.PUBLIC_ANNOTATION_BENCHMARK.value
                or configuration["constraint_adapter"] != vep.CONSTRAINT_SCHEMA
                or _SHA.fullmatch(configuration["constraint_manifest_sha256"]) is None):
            raise ValueError
        if ("bounded_execution_sha256" in configuration
                and _SHA.fullmatch(configuration["bounded_execution_sha256"]) is None):
            raise ValueError
        for key in ("buffer_size", "timeout_seconds"):
            native._integer(configuration[key])
        for request in requests:
            request.__post_init__()
        expected_input = [{"identity_id": r.identity_id, "allele": json.loads(to_json(r.allele))} for r in requests]
        if binding["input_evidence_sha256"] != vep._digest(vep._json(expected_input)):
            raise ValueError
        if file_input:
            published = artifacts['alleles.json'].records(binding['alleles_sha256'])
            used.add('alleles.json')
        else:
            published = native._strict_json(content("alleles.json", binding["alleles_sha256"]))
            if type(published) is not list:
                raise ValueError
        alleles = tuple(from_json(vep._json(a), AlleleEvidence) for a in published)
        del published
        if len(alleles) != len(requests) or len(binding["assessments"]) != len(requests):
            raise ValueError
        index = _parse_fai_bytes(content("reference.fai", binding["fai_sha256"]).encode("utf-8"))
        reference = SimpleNamespace(index=index)
        compatibility = configuration["reference_compatibility"]
        if (set(compatibility) != {"identity_contract_sha256", "identity_contract_id", "primary_assembly_accession",
                                   "vep_fasta_manifest_sha256", "source_manifest_sha256", "scope", "source_naming"}
                or compatibility["primary_assembly_accession"] != "GCA_000001405.15"
                or _SHA.fullmatch(compatibility["identity_contract_sha256"]) is None
                or _SHA.fullmatch(compatibility["vep_fasta_manifest_sha256"]) is None):
            raise ValueError
        maps = {(ReferenceRole.PUBLIC_ANNOTATION_BENCHMARK.value, binding["resource_reference_identity_id"]):
                (compatibility["vep_fasta_manifest_sha256"], dict(vep.resource_native_aliases(reference)))}
        for naming in compatibility["source_naming"]:
            key = (ReferenceRole(naming["role"]).value, naming["identity_id"])
            if (set(naming) != {"role", "identity_id", "aliases_sha256", "aliases"}
                    or key in maps or _SHA.fullmatch(naming["aliases_sha256"]) is None
                    or _SHA.fullmatch(compatibility["source_manifest_sha256"]) is None
                    or type(naming["aliases"]) is not dict
                    or any(type(a) is not str or not a or b not in index for a, b in naming["aliases"].items())):
                raise ValueError
            maps[key] = (compatibility["source_manifest_sha256"], naming["aliases"])
        natives = tuple(next(i for i in r.allele.identities if i.identity_id == r.identity_id) for r in requests)
        native_indexes = [i for i, n in enumerate(natives) if type(n.identity) in (StructuralAllele, RepeatAllele)]
        cache_entries, cache_contigs = {}, set()
        if native_indexes:
            inventory = native._strict_json(content("cache.inventory.json", binding["cache_inventory_sha256"]))
            if (set(inventory) != {"schema", "cache_root", "entries"}
                    or inventory["schema"] not in {INVENTORY_SCHEMA, f"{INVENTORY_SCHEMA}/v1"}
                    or inventory["cache_root"] != "homo_sapiens/116_GRCh38" or type(inventory["entries"]) is not list):
                raise ValueError
            for entry in inventory["entries"]:
                path = entry["path"]
                prefix = inventory["cache_root"] + "/"
                if (set(entry) != {"path", "size_bytes", "sha256"} or not path.startswith(prefix)
                        or any(p in ("", ".", "..") for p in path.split("/"))
                        or _SHA.fullmatch(entry["sha256"]) is None):
                    raise ValueError
                native._integer(entry["size_bytes"], 0)
                path = path[len(prefix):]
                if path in cache_entries:
                    raise ValueError
                cache_entries[path] = entry["sha256"]
                if "/" in path:
                    cache_contigs.add(path.split("/", 1)[0])

        records, queries, new_transcripts = {}, {}, {}
        for i, (request, allele, selected, assessment) in enumerate(zip(requests, alleles, natives, binding["assessments"], strict=True)):
            keys = {"allele_id", "native_identity_id", "requested_identity_id", "source_reference_manifest_sha256",
                    "transport_record_id", "status", "reason"}
            if i in native_indexes:
                keys.add("native_geometry_assessments")
            if (set(assessment) != keys
                    or not _equal(replace(allele, transcript_consequences=None), replace(request.allele, transcript_consequences=None))
                    or any(assessment[k] != v for k, v in {
                        "allele_id": request.allele.allele_id, "native_identity_id": request.identity_id,
                        "requested_identity_id": request.identity_id,
                        "source_reference_manifest_sha256": selected.reference_manifest_sha256,
                    }.items())):
                raise ValueError
            coordinate = maps.get((selected.identity.reference_role.value, selected.identity.reference_identity_id)) if i in native_indexes else None
            sequence = vep._sequence(selected.identity)
            if sequence is not None:
                coordinate = maps.get((sequence.reference_role.value, sequence.reference_identity_id))
            if coordinate is not None and selected.reference_manifest_sha256 != coordinate[0]:
                raise ValueError
            if i in native_indexes:
                if assessment["transport_record_id"] is not None:
                    raise ValueError
                queries[i] = native.native_queries(request, maps, reference, cache_contigs)
                for q in queries[i]:
                    if q.interval_evidence:
                        for contributor in q.interval_evidence.source_fields["contributors"]:
                            reason = contributor["reason"] or ""
                            recipe = contributor["recipe"]
                            if (recipe is not None and recipe["event_type"] != selected.identity.sv_type.value
                                    or reason.endswith("-conflict") or reason == "source-interval-contract-unrecognized-or-invalid"
                                    and contributor["source"]["source_dataset"] == native.INTERVAL_RECIPE_DATASET):
                                raise ValueError
                continue
            reason = None
            if sequence is None:
                reason = "native-identity-has-no-established-sequence-equivalence"
            elif coordinate is None:
                reason = "source-reference-not-coordinate-bound"
            elif any(re.fullmatch(r"[ACGT]+", v) is None for v in (sequence.reference, sequence.alternate)):
                reason = "ambiguous-sequence-bases"
            elif sequence.contig not in coordinate[1]:
                reason = "contig-outside-resource-alias-coverage"
            else:
                contig = coordinate[1][sequence.contig]
                if sequence.position + len(sequence.reference) - 1 > index[contig].length:
                    raise ValueError
            if (assessment["status"] != ("unassessed" if reason else "annotated") or assessment["reason"] != reason
                    or assessment["transport_record_id"] != (None if reason else f"vep-{i + 1:09d}")):
                raise ValueError
            if not reason:
                records[i] = vep.VepVcfRecord(assessment["transport_record_id"], contig, sequence.position,
                                               sequence.reference, (sequence.alternate,))
        order = {name: i for i, name in enumerate(index)}
        ordered = tuple(sorted(records.values(), key=lambda r: (order[r.contig], r.position, r.record_id)))
        transport = "##fileformat=VCFv4.2\n" + "".join(
            f"##contig=<ID={name},length={index[name].length}>\n" for name in dict.fromkeys(r.contig for r in ordered)
        ) + vep._COLUMNS + "\n" + "".join(r.line() for r in ordered)
        if content("transport.vcf", binding["transport_sha256"]) != transport:
            raise ValueError
        provenance = SourcedRecord(SCHEMA, "current", "annotation-replay-" + expected_binding_sha256,
            expected_binding_sha256, binding["bundle_manifest_sha256"], {
                "binding": {k: binding[k] for k in (*_SOURCE_KEYS, "schema", "bundle_manifest_sha256",
                    "configuration", "transport_sha256", "output_sha256", "alleles_sha256")},
                "binding_scope": "relevant-run-fields; source_sha256-pins-complete-binding-bytes",
                "policy_sha256": POLICY_SHA256,
                "policy": native._strict_json(POLICY_BYTES.decode("ascii")),
                "reference_validation": "pinned-publication-and-metadata; not-FASTA-reexecution",
            })
        native_plans, sequence_plans = [], []
        molecular_sources = [() for _ in requests]

        def emit(i, outcome, reason, sources=(), target=None, query=None, locator=None, certainty=None):
            if type(natives[i].identity) is RepeatAllele:
                return
            rule = NATIVE_RULE if i in native_indexes else SEQUENCE_RULE
            certainty = certainty or ("none" if outcome == "evaluated_no_match" else "unknown")
            identity = [requests[i].allele.allele_id, natives[i].identity_id, _plain(natives[i].identity),
                        natives[i].reference_manifest_sha256]
            role = query.role if query else None
            scientific = [SCHEMA, POLICY_SHA256, identity, rule, outcome, target, reason,
                          role, certainty, locator]
            native_plans.append((
                "annotation-association-" + vep._digest(vep._json(_plain(scientific))), requests[i].allele.allele_id,
                natives[i], rule, outcome, "Ensembl" if target else None, target,
                binding["assessments"][i]["status"], role, certainty, (reason, "context-only"),
                tuple(sources), provenance,
            ))

        def terms(i, source, values, target, transcript, query=None, geometry=None):
            for ordinal, term in enumerate(values):
                if term in NON_OVERLAP_TERMS:
                    outcome, reason = "evaluated_no_match", "source-consequence-does-not-establish-transcript-overlap"
                elif term not in INTRAGENIC_TERMS:
                    outcome, reason = "unresolved", "consequence-outside-pinned-overlap-vocabulary"
                elif not transcript:
                    outcome, reason = "evaluated_no_match", "source-feature-is-not-a-transcript"
                elif target is None:
                    outcome, reason = "unresolved", "original-source-gene-unavailable"
                else:
                    outcome, reason = "admitted", "source-transcript-intragenic-consequence"
                fields = source.source_fields
                locator = [fields.get("vcf_record_id"), fields.get("csq_ordinal"),
                           query.query_id if query else None, transcript, ordinal, term]
                emit(i, outcome, reason, (source,) if geometry is None else (geometry, source), target, query,
                     locator, "definite" if outcome == "admitted" else None)

        if ordered:
            source = SourcedRecord("Ensembl VEP", "116", "vep-output", binding["output_sha256"],
                binding["bundle_manifest_sha256"], {k: binding[k] for k in ("schema", *_SOURCE_KEYS, "transport_sha256")})
            parsed = vep.parse_vep_vcf(content("annotated.vcf", binding["output_sha256"]),
                                     expected=ordered, source=source, required_fields=vep.RICH_FIELDS)
            if list(parsed.fields) != binding["csq_fields"]:
                raise ValueError
            for i, record in records.items():
                key = (record.record_id, 1)
                new_transcripts[i] = parsed.transcripts[key]
                molecular_sources[i] = (*(t.source for t in parsed.transcripts[key]),
                                        *parsed.non_transcripts[key])
                sequence_plans.append((requests[i].allele,
                    binding["assessments"][i]["status"], parsed.transcripts[key],
                    parsed.non_transcripts[key], provenance))
        elif binding["output_sha256"] is not None or binding["csq_fields"] is not None:
            raise ValueError

        if native_indexes:
            _replay_native(binding, queries, requests, cache_entries, content, new_transcripts, emit, terms)
        elif "native_annotation" in binding:
            raise ValueError
        for i, request in enumerate(requests):
            previous = request.allele.transcript_consequences
            added = new_transcripts.get(i, ())
            expected = (previous or ()) + added if added or i in records else previous
            if not _equal(_transcript_claims(alleles[i].transcript_consequences), _transcript_claims(expected)):
                raise ValueError
            if i not in native_indexes and i not in records:
                emit(i, "unresolved", binding["assessments"][i]["reason"])
        optional = {"vep.stdout": binding.get("vep.stdout_sha256"), "vep.stderr": binding.get("vep.stderr_sha256")}
        optional.update(binding.get("native_annotation", {}).get("artifacts", {}))
        for name in set(artifacts) - used:
            content(name, optional[name], decode=False)
        return ValidatedAnnotationAdmission(alleles, tuple(molecular_sources),
                                           tuple(sequence_plans), tuple(native_plans))
    except (ValueError, TypeError, KeyError, IndexError, AttributeError, RuntimeError, OverflowError, RecursionError):
        raise AnnotationAdmissionError("invalid or incomplete pinned annotation admission inputs") from None


def _replay_native(binding, queries, requests, cache_entries, content, new_transcripts, emit, terms):
    """Reconstruct publication claims using the delivered query and source parsers."""
    proof = binding["native_annotation"]
    versions = proof["versions"]
    VepVersions.from_mapping(versions)
    if (set(proof) != {"configuration", "versions", "artifacts"}
            or not _equal(proof["configuration"], native.configuration())
            or any(versions[k] != v for k, v in {"cache": "116", "assembly": "GRCh38",
                "source_assembly": "GRCh38.p14", "genebuild": "GENCODE50", "gencode": "GENCODE 50"}.items())):
        raise ValueError
    hashes = proof["artifacts"]
    allowed = {"native.input.tsv", "context.input.jsonl", "native-cache-context.pl"}
    allowed.update(name + "." + suffix for name in ("native", "context") for suffix in ("output.jsonl", "stdout", "stderr", "warnings"))
    if type(hashes) is not dict or set(hashes) - allowed or any(_SHA.fullmatch(v) is None for v in hashes.values()):
        raise ValueError
    valid = tuple(q for group in queries.values() for q in group if q.reason is None)
    contexts = {q.query_id: replace(q, query_id="native-query-" + vep._digest(vep._json([q.query_id, "reference-context"])),
                annotation_scope="transcript_context", sv_type=None) for q in valid if q.annotation_scope == "structural_effect"}
    results = {}
    for context in (False, True):
        expected = tuple(sorted((q for q in (*valid, *contexts.values())
                                 if (q.annotation_scope == "transcript_context") == context),
                                key=lambda q: (q.contig, q.start, q.end, q.query_id)))
        name = "context" if context else "native"
        if not expected:
            if any(k.startswith(name + ".") for k in hashes):
                raise ValueError
            continue
        input_name, output_name = name + (".input.jsonl" if context else ".input.tsv"), name + ".output.jsonl"
        transport = "".join(vep._json(q.transport()) if context else q.line() for q in expected)
        if content(input_name, hashes[input_name]) != transport:
            raise ValueError
        snapshot = hashes if context else {k: v for k, v in hashes.items() if k.startswith("native.")}
        if context and hashes["native-cache-context.pl"] != vep._digest(native.CACHE_HELPER):
            raise ValueError
        source = SourcedRecord("Ensembl VEP cache context" if context else "Ensembl VEP native", "116",
            name + "-output", hashes[output_name], binding["bundle_manifest_sha256"], {
                **{k: binding[k] for k in _SOURCE_KEYS}, "schema": native.SCHEMA,
                "source_revisions": dict(native.SOURCE_REVISIONS), "source_sha256": dict(native.SOURCE_HASHES),
                "versions": versions, "helper_sha256": vep._digest(native.CACHE_HELPER),
                "schema_sha256": vep._digest(vep._json(native.PROTOCOL)),
                "transport_sha256": hashes[input_name], "artifacts": snapshot,
            })
        parsed = native.parse_native_output(content(output_name, hashes[output_name]), expected=expected,
                                            source=source, context=context, cache_entries=cache_entries)
        results.update({k: (*v, source) for k, v in parsed.items()})
    for i, group in queries.items():
        assessment = binding["assessments"][i]
        completed = sum(q.reason is None for q in group)
        status = "unassessed" if not completed else "partially_annotated" if completed < len(group) else (
            "annotated" if any(q.annotation_scope == "structural_effect" for q in group) else "annotated_context")
        gaps = {q.reason for q in group if q.reason}
        reason = None if not gaps else next(iter(gaps)) if len(gaps) == 1 else "multiple-geometry-assessment-gaps"
        if assessment["status"] != status or assessment["reason"] != reason or len(assessment["native_geometry_assessments"]) != len(group):
            raise ValueError
        all_transcripts = []
        for query, stored in zip(group, assessment["native_geometry_assessments"], strict=True):
            reference_context = None
            context_tx = ()
            if query.reason:
                source = SourcedRecord("Ensembl VEP native assessment", "116", query.query_id,
                    vep._digest(vep._json([query.transport(), _plain(query.source_geometry), query.reason])),
                    binding["bundle_manifest_sha256"], {"schema": native.SCHEMA, "versions": versions,
                        "configuration_sha256": binding["configuration_sha256"], "cache_inventory_sha256": binding["cache_inventory_sha256"]})
                tx, raw, count = (), None, None
            else:
                tx, _, raw, source = results[query.query_id]
                if query.query_id in contexts:
                    cq = contexts[query.query_id]
                    context_tx, _, context_raw, context_source = results[cq.query_id]
                    by_transcript = {t.ensembl_transcript_id: t for t in context_tx if t.ensembl_transcript_id is not None}
                    joined, reconstructed = set(), []
                    for effect in tx:
                        cached = by_transcript.get(effect.ensembl_transcript_id)
                        if cached is not None:
                            native._validate_joined_claims(effect.source.source_fields["raw_record"], cached.source.source_fields["raw_record"])
                            joined.add(cached.ensembl_transcript_id)
                            identifier = "native-enriched-" + vep._digest(vep._json([effect.consequence_id, cached.consequence_id]))
                            effect = replace(effect, consequence_id=identifier, source=replace(effect.source,
                                source_record_id=identifier, source_fields={**effect.source.source_fields,
                                    "effect_source_record_id": effect.source.source_record_id, "reference_context": _plain(cached.source)}))
                        reconstructed.append(effect)
                    tx = tuple(reconstructed)
                    context_record = native._bound_source(context_source, cq, context_raw, "geometry-context")
                    reference_context = _plain(replace(context_record, source_fields={**context_record.source_fields,
                        "transcript_count": len(context_tx), "joined_transcripts": sorted(joined),
                        "unjoined_context_transcripts": [t.source.source_fields["raw_record"]["transcript_id"]
                                                         for t in context_tx if t.ensembl_transcript_id not in joined]}))
                count = len(tx)
            fields = {"role": query.role, "annotation_scope": query.annotation_scope,
                "status": "unassessed" if query.reason else "completed", "reason": query.reason,
                "query_id": query.query_id, "query_geometry": {k: v for k, v in query.transport().items() if k != "query_id"},
                "source_geometry": _plain(query.source_geometry), "transcript_count": count,
                "lookup_result": None if query.reason else "matched" if count else "completed-no-overlap",
                "raw_result": raw, "annotation_binding": source.source_fields,
                "interval_evidence": _plain(query.interval_evidence), "reference_context": reference_context,
                "native_identity_id": requests[i].identity_id, "allele_id": requests[i].allele.allele_id}
            geometry = replace(source, source_record_id=query.query_id, source_fields=fields)
            if not _equal(stored, geometry):
                raise ValueError
            all_transcripts.extend(tx)
            geometry_gap = None
            if query.role == "structural_interval":
                recipes = [c["recipe"] for c in query.interval_evidence.source_fields["contributors"] if c["recipe"] is not None]
                if any(r[k][v] is None for r in recipes for k in ("start", "end") for v in ("min", "max")):
                    geometry_gap = "unbounded-source-geometry"
                elif any(r["interval_semantics"] == "unknown" for r in recipes):
                    geometry_gap = "unknown-source-geometry"
            if query.reason or geometry_gap:
                if geometry_gap and tx:
                    for t in tx:
                        raw_tx = t.source.source_fields["raw_record"]
                        emit(i, "unresolved", geometry_gap, (geometry, t.source),
                             _gene(raw_tx.get("gene_id")), query, [query.query_id, raw_tx["transcript_id"]])
                else:
                    emit(i, "unresolved", query.reason or geometry_gap, (geometry,), query=query, locator=query.query_id)
                continue
            certainty = "definite" if query.role == "affected_interval" else "possible"
            if query.annotation_scope == "structural_effect":
                for t in tx:
                    raw_tx = t.source.source_fields["raw_record"]
                    terms(i, t.source, raw_tx["consequence_terms"], _gene(raw_tx.get("gene_id")),
                          raw_tx["transcript_id"], query, geometry)
                overlaps = context_tx
            else:
                overlaps = tx
            for t in overlaps:
                raw_tx = t.source.source_fields["raw_record"]
                target = _gene(raw_tx.get("gene_id"))
                emit(i, "admitted" if target else "unresolved",
                     "validated-native-transcript-overlap" if target else "original-source-gene-unavailable",
                     (geometry, t.source), target, query,
                     [query.query_id, raw_tx["transcript_id"], raw_tx["transcript_geometry"]],
                     certainty if target else "unknown")
            if not overlaps:
                emit(i, "evaluated_no_match", "completed-empty-transcript-geometry", (geometry,), query=query, locator=query.query_id)
        new_transcripts[i] = tuple(all_transcripts)
