"""Offline source-specific association proofs, not candidate or clinical scores.

The three report iterators must be exhausted successfully before accepting a
batch. They verify actual files against independently selected contract/manifest
pins, hold stable file descriptors through evaluation, and never query a remote
service. Only supplied-query geometry and matched records are retained, not a
genome-wide regulatory-link or clinical-variant table. Callers may stream proofs
to their controlled evidence registry; this module does not own that registry.

An ``admitted`` proof establishes ONLY the association rule. The consumer still
has to resolve applicable positive native disease support, reconcile its entire
call inventory, and recompute this report rather than trust serialized outcomes.
Pure row parsers validate syntax/semantics but do NOT authenticate their inputs.
"""

from __future__ import annotations

from collections.abc import Iterator, Mapping
from contextlib import ExitStack
from dataclasses import dataclass, replace
from decimal import Decimal
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import socket
from typing import Literal

from genosapiens import clinvar as cv, open_targets as ot, variant_class_resources as vc
from genosapiens.allele_identity import (
    MosaicAllele, MitochondrialAllele, SequenceAllele, StructuralAllele, RepeatAllele,
)
from genosapiens.evidence_model import AlleleEvidence, IdentityMapping, SourcedRecord, _Controlled
from genosapiens.evidence_serialization import _loads, _plain
from genosapiens.grch38_reference import IndexedFasta
from genosapiens.governance import ArtifactClassification
from genosapiens.llm_guard import LLMPolicy
from genosapiens.manifest import PublicResourceManifest
from genosapiens.reference_concordance import (
    ReferenceRoleBinding, _open_file_sha256, _open_stable_binary, load_reference_identity,
)
from genosapiens.submission_rendering import load_submission_contig_renderer
from genosapiens.vep_annotation import VepRequest, _local
from genosapiens.vep_native_annotation import INTERVAL_RECIPE_DATASET, native_queries


SCHEMA = "genosapiens.admission-links"
ENCODE_CCRE = "encode-registry-v4-ccre-grch38"
ENCODE_RE2G = "encode-re2g-registry-v4-grch38"
_GENE = re.compile(r"ENSG[0-9]{11}(?:\.[1-9][0-9]*)?")
_CCRE_ACCESSION = re.compile(r"EH38E[0-9]+")
_DNA = re.compile(r"[ACGT]+")
_SHA = re.compile(r"[0-9a-f]{64}")
_CCRE_COLUMNS = ("chromosome", "start", "end", "current_id", "legacy_id", "class")
_RE2G_COLUMNS = ("legacy_id", "gene_id", "gene_symbol", "gene_biotype", "method",
                 "context", "accession", "biosample", "score")
_EVA_LISTS = ("alleleOrigins", "allelicRequirements", "clinicalSignificances",
              "cohortPhenotypes", "literature", "qualityControls")
_EVA_TEXT = ("targetId", "id", "targetFromSourceId", "diseaseFromSourceMappedId",
             "confidence", "datasourceId", "datatypeId", "diseaseFromSource",
             "diseaseFromSourceId", "releaseDate", "studyId", "variantFromSourceId",
             "variantFunctionalConsequenceId", "variantHgvsId", "variantId", "variantRsId",
             "diseaseId", "publicationDate", "evidenceDate", "directionOnTrait", "directionOnTarget")
_POLICY = {
    "schema": SCHEMA,
    "encode_join": "Registry-V4-BED6[4]=rE2G-TSV9[0]; target=TSV[1]",
    "encode_target": "raw-target-preserved; supported-Ensembl-or-missing-or-unmapped; no-symbol-or-alias-mapping",
    "coordinates": "zero-based-half-open; native-reference-span-start=position-1",
    "native_geometry": "reuse-vep-native-query-roles; no-new-event-span-or-proxy",
    "native_source_gate": "invalid-binding-fails; unbounded-or-unknown-interval-unresolved; independent-bounded-anchors-retained",
    "clinical": "exact-variant-row-fields-to-source-variantId-to-single-EVA-record",
    "condition": "diseaseFromSourceId-only; mapped-identifiers-retained-not-equated",
    "clinvar": "VCF-aggregates-retained; condition-link-unresolved",
    "unknown": "missing-assertion-is-unknown-not-benign",
}


class AdmissionLinkError(ValueError):
    """Static diagnostics only; no patient identifiers, paths, rows, or SQL."""


@dataclass(frozen=True, slots=True, repr=False)
class PinnedPublicSource:
    """Pins supplied by the selected policy, never recovered from a SourcedRecord.

    root is the common reference namespace containing manifests and releases.
    For Open Targets, only metadata, output/variant and output/evidence_eva must
    be mounted below the contract's release prefix; its independently pinned
    full manifest is still required. For reference
    coordinates, use the reference identity contract and reconstructed manifest.
    """

    contract_path: Path
    contract_sha256: str
    root: Path
    manifest_path: Path
    manifest_sha256: str


@dataclass(frozen=True, slots=True, repr=False)
class AssociationProof(_Controlled):
    """Immutable controlled result using native identity and existing source records.

    target_namespace is Ensembl for these adapters; condition_id is the original
    source identifier, never a disease label or an inferred cross-source mapping.
    Source locators retain separate physical rows even for repeated assertions.
    """

    allele_id: str
    native: IdentityMapping
    rule: Literal["regulatory_link", "explicit_clinical_link"]
    outcome: Literal["admitted", "evaluated_no_match", "unresolved"]
    target_namespace: str | None
    target_id: str | None
    condition_id: str | None
    reasons: tuple[str, ...]
    source_refs: tuple[SourcedRecord, ...]
    provenance: SourcedRecord

    classification = ArtifactClassification.CONTROLLED
    llm_policy = LLMPolicy.NO_LLM


def _json(value) -> str:
    return json.dumps(_plain(value), sort_keys=True, separators=(",", ":"), allow_nan=False, ensure_ascii=True)


def _digest(value) -> str:
    return hashlib.sha256(_json(value).encode()).hexdigest()


def _text(value):
    if type(value) is not str or not value or value.strip() != value or any(c in value for c in "\r\n\t"):
        raise ValueError
    return value


def _integer(value, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError
    return value


def parse_encode_ccre_row(row: str) -> Mapping[str, object]:
    """Registry V4 BED6: column 3 is a D/F anchor; column 4 is the E cCRE.

    Internal current_id/legacy_id keys retain the source-contract labels;
    both accessions are current biological identifiers, with only E joining rE2G.
    """
    try:
        columns = tuple(row.rstrip("\r\n").split("\t"))
        if len(columns) != 6:
            raise ValueError
        for value in columns:
            _text(value)
        if not re.fullmatch(r"[0-9]+", columns[1]) or not re.fullmatch(r"[0-9]+", columns[2]):
            raise ValueError
        start, end = int(columns[1]), int(columns[2])
        if start >= end or not re.fullmatch(r"EH38[DF][0-9]+", columns[3]):
            raise ValueError
        # A missing cCRE accession is retained but cannot manufacture a join.
        if columns[4] != "." and not re.fullmatch(r"EH38E[0-9]+", columns[4]):
            raise ValueError
        return dict(zip(_CCRE_COLUMNS, columns, strict=True), start=start, end=end, columns=columns)
    except (ValueError, TypeError, AttributeError):
        raise AdmissionLinkError("invalid ENCODE cCRE row") from None


def classify_encode_target(value: str) -> str:
    """Classify source identity syntax only, without resolving a display label."""
    _text(value)
    return 'missing' if value == '.' else 'ensembl' if _GENE.fullmatch(value) else 'unmapped'


def parse_encode_re2g_row(row: str) -> Mapping[str, object]:
    """Nine source columns; gene_id is a raw target label, not always Ensembl."""
    try:
        text = row.rstrip("\r\n")
        columns = tuple(text.split("\t"))
        # Tabs delimit the nine fields; embedded CR/LF and untrimmed/empty
        # fields remain invalid. Check line structure once, rather than
        # recreating nine control-character generators on every tissue row.
        if (len(columns) != 9 or '\r' in text or '\n' in text
                or any(not value or value.strip() != value for value in columns)):
            raise ValueError
        if not _CCRE_ACCESSION.fullmatch(columns[0]):
            raise ValueError
        target_kind = classify_encode_target(columns[1])
        # Validate a number, but never round, threshold, or replace its spelling.
        if columns[8] != "." and not Decimal(columns[8]).is_finite():
            raise ValueError
        return dict(zip(_RE2G_COLUMNS, columns, strict=True), columns=columns,
                    source_target_id=columns[1], target_kind=target_kind,
                    ensembl_gene_id=columns[1] if target_kind == 'ensembl' else None)
    except (ValueError, TypeError, AttributeError, ArithmeticError):
        raise AdmissionLinkError("invalid ENCODE rE2G row") from None


def parse_open_targets_variant_row(row: Mapping[str, object]) -> tuple[str, int, str, str]:
    """Bind a publisher key to its explicit fields, including opaque OTVAR keys.

    No rsID/HGVS lookup, indel shifting, or reconstruction of a hashed ID occurs.
    The exact row's fields, not a guessed decoding of its key, identify the allele.
    """
    try:
        key = _text(row["variantId"])
        chrom, pos = _text(row["chromosome"]), _integer(row["position"], 1)
        ref, alt = _text(row["referenceAllele"]), _text(row["alternateAllele"])
        if not _DNA.fullmatch(ref) or not _DNA.fullmatch(alt) or ref == alt:
            raise ValueError
        if key != f"{chrom}_{pos}_{ref}_{alt}":
            opaque = re.fullmatch(r"OTVAR_(?:([^_]+)_([1-9][0-9]*)_)?[0-9a-f]{32}", key)
            if opaque is None or opaque[1] is not None and (opaque[1], int(opaque[2])) != (chrom, pos):
                raise ValueError
        return chrom, pos, ref, alt
    except (ValueError, TypeError, KeyError, AttributeError):
        raise AdmissionLinkError("invalid Open Targets exact variant row") from None


def parse_open_targets_eva_row(row: Mapping[str, object]) -> tuple[str | None, str | None]:
    """One 26.06 germline assertion, never a significance/condition product.

    All fields are retained by the loader. confidence is publisher-native review
    context, not a fabricated ClinVar star count. The source's diseaseId mapping
    remains in the row; it does not replace diseaseFromSourceId or prove equality.
    """
    try:
        if set(row) != {*_EVA_TEXT, *_EVA_LISTS, "score"}:
            raise ValueError
        for name in _EVA_TEXT:
            if row[name] is not None:
                _text(row[name])
        for name in _EVA_LISTS:
            values = row[name]
            if values is not None and (type(values) not in (list, tuple) or any(
                value is not None and type(value) is not str for value in values
            )):
                raise ValueError
        if row["datasourceId"] != "eva" or row["datatypeId"] != "genetic_association":
            raise ValueError
        _text(row["id"])
        _text(row["variantId"])
        if row["targetId"] is not None and not _GENE.fullmatch(row["targetId"]):
            raise ValueError
        # JSON validation preserves nulls and rejects NaN without scoring anything.
        if row["score"] is not None and type(row["score"]) not in (int, float):
            raise ValueError
        _json(row)
        return row["targetId"], row["diseaseFromSourceId"]
    except (ValueError, TypeError, KeyError, AttributeError):
        raise AdmissionLinkError("invalid Open Targets EVA assertion row") from None


def parse_clinvar_vcf_row(row: str) -> Mapping[str, object]:
    """Retain site-level INFO strings verbatim, including conflicting aggregates."""
    try:
        fields = tuple(row.rstrip("\r\n").split("\t"))
        if len(fields) != 8 or not re.fullmatch(r"[1-9][0-9]*", fields[1]):
            raise ValueError
        if not fields[2].isdecimal() or not fields[0] or not fields[3] or not fields[4]:
            raise ValueError
        info = cv._parse_info(fields[7]) if fields[7] != "." else {}
        return {"columns": fields, "info": info}
    except (ValueError, TypeError, AttributeError, RuntimeError):
        raise AdmissionLinkError("invalid ClinVar VCF row") from None


class _Files:
    """One report lifetime: authenticated inputs and stable descriptor checks."""

    def __init__(self, stack, *, reuse_published_outputs=False):
        if type(reuse_published_outputs) is not bool:
            raise ValueError
        self.stack = stack
        self.opened = {}
        self.reuse_published_outputs = reuse_published_outputs
        self.published_outputs = {}

    def open(self, path, sha256=None, size=None):
        path = _local(path)
        if not path.is_file():
            raise ValueError
        if path not in self.opened:
            self.opened[path] = self.stack.enter_context(_open_stable_binary(path, "admission source"))
        stream = self.opened[path]
        published = self.published_outputs.get(path)
        trusted = self.reuse_published_outputs and published is not None
        if trusted and (sha256, size) != published:
            raise ValueError
        if sha256 is not None and (not _SHA.fullmatch(sha256)
                or not trusted and _open_file_sha256(stream) != sha256):
            raise ValueError
        if size is not None and path.stat().st_size != size:
            raise ValueError
        stream.seek(0)
        return stream

    def json(self, path, sha256):
        stream = self.open(path, sha256)
        content = stream.read(8 * 1024 * 1024 + 1)
        if len(content) > 8 * 1024 * 1024:
            raise ValueError
        return _loads(content.decode("utf-8"))

    def source(self, source):
        if type(source) is not PinnedPublicSource:
            raise ValueError
        root = _local(source.root)
        if not _local(source.manifest_path).is_relative_to(root):
            raise ValueError
        contract = self.json(source.contract_path, source.contract_sha256)
        manifest = PublicResourceManifest.from_mapping(self.json(source.manifest_path, source.manifest_sha256))
        if self.reuse_published_outputs:
            for output in manifest.outputs:
                path = _local(root / output.path)
                if not path.is_relative_to(root):
                    raise ValueError
                binding = (output.sha256, output.size_bytes)
                if path in self.published_outputs and self.published_outputs[path] != binding:
                    raise ValueError
                self.published_outputs[path] = binding
        return contract, manifest

    def lines(self, path):
        stream = self.open(path)
        binary = gzip.GzipFile(fileobj=stream) if Path(path).suffix in (".gz", ".bgz") else stream
        text = io.TextIOWrapper(binary, encoding="utf-8", newline="")
        try:
            number = 0
            while row := text.readline(8 * 1024 * 1024 + 1):
                number += 1
                if len(row) > 8 * 1024 * 1024:
                    raise ValueError
                yield number, row
        finally:
            text.detach()
            if binary is not stream:
                binary.close()


def _reference(files, pin):
    data, manifest = files.source(pin)
    contract = load_reference_identity(pin.contract_path)
    if (data["primary_assembly_accession"] != "GCA_000001405.15"
            or contract.assembly != "GRCh38-masked-v2-no-alt-plus-hs38d1"):
        raise ValueError
    role = contract.reference_role("public_annotation_benchmark")
    outputs = {o.path: o for o in manifest.outputs}
    maps = {}
    for output_role in ("contig_aliases", *(r.fai_output_role for r in contract.reference_roles)):
        path = contract.expected_output(output_role).path
        output = outputs[path]
        files.open(pin.root / path, output.sha256, output.size_bytes)
    renderers = {}
    for source_role in contract.reference_roles:
        renderers[source_role.role] = load_submission_contig_renderer(
            identity_path=pin.contract_path, reference_manifest_path=pin.manifest_path,
            reference_root=pin.root, reference_binding=ReferenceRoleBinding(source_role.role, source_role.identity_id),
        )
    public_names = renderers[role.role]
    for source_role in contract.reference_roles:
        names = renderers[source_role.role]
        maps[(source_role.role, source_role.identity_id)] = (
            pin.manifest_sha256,
            {
                name: public_names.reference_aliases[canonical]
                for name, (canonical, _) in names.accepted_names.items()
                if canonical in public_names.reference_aliases
            },
        )
    reference = IndexedFasta(pin.root / contract.expected_output(role.fasta_output_role).path,
                             pin.root / contract.expected_output(role.fai_output_role).path)
    provenance = {"reference_contract_sha256": pin.contract_sha256, "reference_manifest_sha256": pin.manifest_sha256,
                  "reference_roles": data["reference_roles"], "aliases_sha256": public_names.contig_aliases_sha256,
                  "relationship": "coordinate-compatible-only; no sequence-equivalence or case-acceptance claim",
                  "source_reference_role": role.role, "source_reference_identity_id": role.identity_id}
    return maps, reference, provenance


_PUBLIC_QUERY_CHROMOSOMES = {
    **{"chr" + str(number): str(number) for number in range(1, 23)},
    "chrX": "X",
    "chrY": "Y",
    "chrM": "MT",
}


def _public_query_chromosome(contig):
    try:
        return _PUBLIC_QUERY_CHROMOSOMES[contig]
    except (KeyError, TypeError):
        raise ValueError("public clinical query contig is unsupported") from None


def _inputs(alleles):
    if type(alleles) is not tuple or any(type(a) is not AlleleEvidence for a in alleles):
        raise ValueError
    if len({a.allele_id for a in alleles}) != len(alleles):
        raise ValueError
    if any(name != "lo" for _, name in socket.if_nameindex()):
        raise AdmissionLinkError("admission source reports require offline execution")
    return tuple(sorted(alleles, key=lambda a: a.allele_id))


def _native(allele):
    return next(i for i in allele.identities if i.identity_id == allele.native_identity_id)


def _sequence(allele, maps, reference):
    native = _native(allele)
    identity = native.identity
    sequence = identity.sequence if type(identity) in (MosaicAllele, MitochondrialAllele) else identity
    binding = maps.get((sequence.reference_role.value, sequence.reference_identity_id))
    if binding is None:
        return None, "source-reference-not-coordinate-bound"
    if native.reference_manifest_sha256 != binding[0]:
        raise ValueError
    if type(sequence) is not SequenceAllele:
        return None, "native-identity-has-no-exact-sequence-rule; no-proxy"
    chrom = binding[1].get(sequence.contig)
    if chrom is None or chrom not in reference.index:
        return None, "contig-outside-source-alias-coverage"
    if sequence.position + len(sequence.reference) - 1 > reference.index[chrom].length:
        raise ValueError
    # This returns a lookup tuple, never a replacement/proxy SequenceAllele.
    return (chrom, sequence.position, sequence.reference, sequence.alternate), None


def _proof(allele, rule, outcome, reasons, sources, provenance, target=None, condition=None):
    return AssociationProof(allele.allele_id, _native(allele), rule, outcome,
                            "Ensembl" if target else None, target, condition,
                            tuple(reasons), tuple(sources), provenance)


def _provenance(adapter, sources, reference_proof):
    value = {"schema": SCHEMA, "adapter": adapter,
             "configuration_sha256": _digest(_POLICY), "policy": _POLICY,
             "sources": sources, "reference": reference_proof,
             "scope": "association-only; disease-support-and-inventory-accounting-required-by-consumer",
             "interpretation": "no-score-cutoff; no-mechanism-or-tissue-inference; no-frequency-rescue"}
    digest = _digest(value)
    return SourcedRecord(SCHEMA, "current", adapter + "-" + digest, digest, None, value)


def _source_record(dataset, release, path, ordinal, sha, manifest_sha, row):
    return SourcedRecord(dataset, release, f"{path}#row={ordinal}", sha, manifest_sha,
                         {"artifact_path": path, "row_ordinal": ordinal, "record": row})


def _query_tree(rows):
    """Balanced interval index of supplied queries only, with subtree end bounds."""
    if not rows:
        return None
    mid = len(rows) // 2
    left, right = _query_tree(rows[:mid]), _query_tree(rows[mid + 1:])
    return rows[mid], left, right, max(rows[mid][1], left[3] if left else 0, right[3] if right else 0)


def _overlaps(tree, start, end):
    if tree is None or tree[3] <= start:
        return
    row, left, right, _ = tree
    yield from _overlaps(left, start, end)
    if row[0] < end:
        if start < row[1]:
            yield row[2]
        yield from _overlaps(right, start, end)


def encode_query_geometries(alleles, maps, fasta):
    """Reuse the existing source-aware regulatory geometry in scan and index paths."""
    geometries, by_contig, unresolved = {}, {}, {}
    for allele in alleles:
        sequence, reason = _sequence(allele, maps, fasta)
        if sequence is not None:
            chrom, pos, ref, _ = sequence
            rows = ((chrom, pos - 1, pos + len(ref) - 1, {"role": "sequence-reference-span", "certainty": "context-only"}),)
        elif type(_native(allele).identity) in (StructuralAllele, RepeatAllele):
            selected = _native(allele)
            queries = native_queries(VepRequest(allele), maps, fasta, set(fasta.index))
            rows, gaps, gap_sources = [], [], []
            for q in queries:
                recipes = []
                if q.interval_evidence is not None:
                    for contributor in q.interval_evidence.source_fields["contributors"]:
                        source_reason = contributor["reason"] or ""
                        recipe = contributor["recipe"]
                        # Same contributor gates as annotation admission;
                        # a nominal query span does not authenticate it.
                        if (source_reason.endswith("-conflict")
                                or source_reason == "source-interval-contract-unrecognized-or-invalid"
                                and contributor["source"]["source_dataset"] == INTERVAL_RECIPE_DATASET):
                            raise ValueError
                        if recipe is not None:
                            if (recipe["event_type"] != selected.identity.sv_type.value
                                    or recipe["reference"] != {
                                        "role": selected.identity.reference_role.value,
                                        "identity_id": selected.identity.reference_identity_id,
                                        "manifest_sha256": selected.reference_manifest_sha256,
                                    }
                                    or selected.identity.interval is not None
                                    and recipe["contig"] != selected.identity.interval.contig):
                                raise ValueError
                            recipes.append(recipe)
                geometry_gap = None
                if q.role == "structural_interval":
                    if any(r[k][v] is None for r in recipes for k in ("start", "end") for v in ("min", "max")):
                        geometry_gap = "unbounded-source-geometry"
                    elif any(r["interval_semantics"] == "unknown" for r in recipes):
                        geometry_gap = "unknown-source-geometry"
                if q.reason or geometry_gap:
                    gaps.append(q.reason or geometry_gap)
                    if q.interval_evidence is not None:
                        gap_sources.append(q.interval_evidence)
                else:
                    rows.append((q.contig, q.start - 1, q.end, _plain(q)))
            if gaps:
                unresolved[allele.allele_id] = (tuple(gaps), tuple(gap_sources))
        else:
            rows = ()
        if not rows:
            unresolved.setdefault(allele.allele_id, ((reason or "geometry-unresolved",), ()))
        for index, (chrom, start, end, geometry) in enumerate(rows):
            key = (allele.allele_id, index)
            geometries[key] = (allele, geometry)
            by_contig.setdefault(chrom, []).append((start, end, key))
    trees = {chrom: _query_tree(sorted(rows)) for chrom, rows in by_contig.items()}
    return geometries, trees, unresolved, by_contig


def encode_regulatory_report(
    alleles: tuple[AlleleEvidence, ...], *, resources: tuple[PinnedPublicSource, PinnedPublicSource],
    reference: PinnedPublicSource,
) -> Iterator[AssociationProof]:
    """Scan BED and gzip TSV once per batch; join legacy column 4 to column 0.

    Existing native query semantics retain interval/anchor/uncertainty roles.
    Bounded native geometry without established affected-interval semantics
    remains possible/context-only. Invalid contributors fail validation;
    unbounded/unknown intervals remain unresolved without erasing independent
    bounded anchors. No structural allele becomes a sequence proxy.
    Missing links/targets remain unresolved. A complete BED scan with no overlap
    is evaluated_no_match in this resource and query scope only.
    """
    try:
        alleles = _inputs(alleles)
        with ExitStack() as stack:
            files = _Files(stack)
            maps, fasta, ref_proof = _reference(files, reference)
            if len(resources) != 2:
                raise ValueError
            installed, bindings = {}, []
            for pin in resources:
                data, manifest = files.source(pin)
                contract = vc.VariantClassBundleContract.from_mapping(data)
                if (contract.regulatory_join.registry_resource_id != ENCODE_CCRE
                        or contract.regulatory_join.links_resource_id != ENCODE_RE2G
                        or (contract.regulatory_join.registry_legacy_id_column, contract.regulatory_join.links_legacy_id_column) != (4, 0)):
                    raise ValueError
                resource = contract.resource(manifest.resource_id)
                if resource.resource_id not in (ENCODE_CCRE, ENCODE_RE2G) or resource.release != "Registry-V4":
                    raise ValueError
                if resource.resource_id in installed:
                    raise ValueError
                vc._verify_manifest_contract(manifest, resource, pin.contract_sha256)
                paths = {}
                for artifact in resource.artifacts:
                    output = next(o for o in manifest.outputs if o.path == artifact.path)
                    files.open(pin.root / artifact.path, output.sha256, output.size_bytes)
                    paths[artifact.role] = pin.root / artifact.path
                vc.validate_class_resource(resource, paths)
                width = 6 if resource.resource_id == ENCODE_CCRE else 9
                if resource.validation["column_count"] != width:
                    raise ValueError
                artifact = resource.artifact(resource.primary_artifact_role)
                installed[resource.resource_id] = (pin, resource, artifact, paths[artifact.role])
                bindings.append({"resource_id": resource.resource_id, "contract_sha256": pin.contract_sha256,
                                 "manifest_sha256": pin.manifest_sha256, "artifact": _plain(artifact)})
            if resources[0].contract_sha256 != resources[1].contract_sha256:
                raise ValueError
            provenance = _provenance("encode-registry-v4-re2g", sorted(bindings, key=lambda b: b["resource_id"]), ref_proof)
            geometries, trees, unresolved, _ = encode_query_geometries(alleles, maps, fasta)
            hits, overlapping = {}, set()
            pin, resource, artifact, path = installed[ENCODE_CCRE]
            for ordinal, text in files.lines(path):
                row = parse_encode_ccre_row(text)
                chrom, start, end = row["chromosome"], row["start"], row["end"]
                if chrom in fasta.index and end > fasta.index[chrom].length:
                    raise ValueError
                keys = tuple(_overlaps(trees.get(chrom), start, end))
                if not keys:
                    continue
                record = _source_record(resource.resource_id, resource.release, artifact.path, ordinal,
                                        artifact.sha256, pin.manifest_sha256, row)
                for key in keys:
                    overlapping.add(key[0])
                hits.setdefault(row["legacy_id"], []).append((row, record, keys))
            linked = set()
            pin, resource, artifact, path = installed[ENCODE_RE2G]
            for ordinal, text in files.lines(path):
                row = parse_encode_re2g_row(text)
                matches = hits.get(row["legacy_id"], ())
                if not matches:
                    continue
                record = _source_record(resource.resource_id, resource.release, artifact.path, ordinal,
                                        artifact.sha256, pin.manifest_sha256, row)
                for element, element_record, keys in matches:
                    # The contract joins IDs, not coordinate-looking context
                    # strings. Preserve that separate source column verbatim.
                    linked.add(element_record.source_record_id)
                    for key in keys:
                        allele, geometry = geometries[key]
                        geometry_digest = _digest(geometry)
                        geometry_record = SourcedRecord("genosapiens.admission-geometry", "current", geometry_digest,
                                                       geometry_digest, None, geometry)
                        target = row["ensembl_gene_id"]
                        yield _proof(allele, "regulatory_link", "admitted" if target else "unresolved",
                                     ("legacy-ccre-explicit-target-overlap" if target else "source-target-unresolved", "context-only"),
                                     (element_record, record, geometry_record), provenance, target)
            for elements in hits.values():
                for _, record, keys in elements:
                    if record.source_record_id not in linked:
                        for key in keys:
                            yield _proof(geometries[key][0], "regulatory_link", "unresolved",
                                         ("overlapping-element-without-matched-target-link",), (record,), provenance)
            for allele in alleles:
                if allele.allele_id in unresolved:
                    reasons, sources = unresolved[allele.allele_id]
                    yield _proof(allele, "regulatory_link", "unresolved", reasons, sources, provenance)
                elif allele.allele_id not in overlapping:
                    yield _proof(allele, "regulatory_link", "evaluated_no_match",
                                 ("no-element-overlap-in-verified-registry-scope",), (), provenance)
    except Exception:
        raise AdmissionLinkError("offline ENCODE association report failed") from None


def _open_targets_metadata(files, pin):
    data, manifest = files.source(pin)
    contract = ot.OpenTargetsContract.from_mapping(data)
    if (contract.resource, contract.release, contract.assembly) != ("open-targets-platform", "26.06", "GRCh38"):
        raise ValueError
    release = pin.root / contract.output_prefix
    metadata = tuple(ot.ReleaseArtifact(None, m.name, m.size_bytes, m.sha1, m.sha256) for m in contract.metadata)
    for m in metadata:
        files.open(release / m.relative_path, m.sha256, m.size_bytes)
    ot._verify_artifacts(release, metadata)
    artifacts, models, _ = ot._release_model(release, contract)
    ot._verify_manifest_contract(manifest, contract, pin.contract_sha256, artifacts)
    return release, contract, manifest, artifacts, models


def preflight_clinical_source(pin, kind):
    """Qualify topology and current metadata before hashing scientific outputs."""
    with ExitStack() as stack:
        files = _Files(stack)
        if kind == "open_targets":
            release, contract, manifest, artifacts, _models = _open_targets_metadata(files, pin)
            paths = tuple(f"{contract.output_prefix}/{a.relative_path}" for a in artifacts
                          if a.dataset in (None, "variant", "evidence_eva"))
        elif kind == "clinvar":
            data, manifest = files.source(pin)
            contract = cv.ClinvarContract.from_mapping(data)
            cv._verify_manifest_contract(manifest, contract, pin.contract_sha256)
            paths = tuple(a.path for a in contract.artifacts)
        else:
            raise AdmissionLinkError("unknown clinical preflight source")
        outputs = {o.path: o for o in manifest.outputs}
        for relative in paths:
            path = _local(pin.root / relative)
            output = outputs[relative]
            if (not path.is_relative_to(_local(pin.root)) or not path.is_file()
                    or path.stat().st_size != output.size_bytes):
                raise AdmissionLinkError("clinical resource topology mismatch")


def _open_targets(files, pin, *, reuse_published_reconciliation=False):
    release, contract, manifest, artifacts, models = _open_targets_metadata(files, pin)
    outputs = {o.path: o for o in manifest.outputs}
    selected = []
    for artifact in artifacts:
        if artifact.dataset not in (None, "variant", "evidence_eva"):
            continue
        output = outputs[f"{contract.output_prefix}/{artifact.relative_path}"]
        if output.size_bytes != artifact.size_bytes or artifact.sha256 not in (None, output.sha256):
            raise ValueError
        files.open(release / artifact.relative_path, output.sha256, output.size_bytes)
        selected.append(replace(artifact, sha256=output.sha256))
    datasets = tuple(d for d in contract.datasets if d.name in ("variant", "evidence_eva"))
    selected_models = tuple(m for m in models if m.name in ("variant", "evidence_eva"))
    if {d.name for d in datasets} != {"variant", "evidence_eva"}:
        raise ValueError
    for dataset in datasets:
        partitions = tuple(a for a in selected if a.dataset == dataset.name)
        prefix = "output/" + dataset.name
        ot._verify_release_tree(release / prefix, tuple(replace(a, relative_path=Path(a.relative_path).name) for a in partitions))
    if type(reuse_published_reconciliation) is not bool:
        raise ValueError
    verified = (
        tuple(selected) if reuse_published_reconciliation
        else ot._verify_artifacts(release, selected)
    )
    if not reuse_published_reconciliation:
        ot._validate_parquet(
            release, replace(contract, datasets=datasets, semantic_joins=()),
            verified, selected_models,
        )
    for model in selected_models:
        expected_key = "variantId" if model.name == "variant" else "id"
        if model.key_fields != (model.name + "/" + expected_key,):
            raise ValueError
        if model.name == "evidence_eva" and set(model.top_level_fields) != {*_EVA_TEXT, *_EVA_LISTS, "score"}:
            raise ValueError
    return release, contract, verified


def _rows(connection, path, where, parameters):
    cursor = connection.execute("SELECT * FROM read_parquet(?, file_row_number=true) WHERE " + where,
                                [str(path), *parameters])
    columns = tuple(c[0] for c in cursor.description)
    while rows := cursor.fetchmany(256):
        for values in rows:
            row = dict(zip(columns, values, strict=True))
            ordinal = row.pop("file_row_number")
            yield ordinal, row


def open_targets_clinical_report(
    alleles: tuple[AlleleEvidence, ...], *, source: PinnedPublicSource, reference: PinnedPublicSource,
) -> Iterator[AssociationProof]:
    """Exact variant fields -> verified variantId -> individual EVA assertions.

    All selected partitions must exist and match the release integrity index,
    full public manifest, and Croissant schema. Unrelated datasets are not read.
    Missing exact rows/assertions are unknown, even if other clinical evidence or
    an ``explicit_condition_link`` flag was supplied on the allele.

    An assertion proof keeps its bound variant and EVA row first, followed by
    other exact physical variant rows for that allele as raw context, including
    keys without EVA assertions. This adds no proofs or admission outcomes.
    """
    try:
        import duckdb
        alleles = _inputs(alleles)
        if duckdb.__version__ != "1.4.2":
            raise ValueError
        with ExitStack() as stack:
            files = _Files(stack)
            maps, fasta, ref_proof = _reference(files, reference)
            release, contract, artifacts = _open_targets(files, source)
            provenance = _provenance("open-targets-26.06-eva", [{"contract_sha256": source.contract_sha256,
                "manifest_sha256": source.manifest_sha256, "artifacts": [_plain(a) for a in artifacts],
                "schema_source": "pinned-croissant", "partition_scope": ["variant", "evidence_eva"]}], ref_proof)
            lookups, unresolved = {}, {}
            for allele in alleles:
                sequence, reason = _sequence(allele, maps, fasta)
                if sequence is None or not _DNA.fullmatch(sequence[2]) or not _DNA.fullmatch(sequence[3]):
                    unresolved[allele.allele_id] = reason or "ambiguous-sequence-exact-link-unresolved"
                    continue
                chrom, pos, ref, alt = sequence
                native_chrom = _public_query_chromosome(chrom)
                lookups.setdefault((native_chrom, pos, ref, alt), []).append(allele)
            connection = stack.enter_context(duckdb.connect(database=":memory:", config={
                "threads": "1", "memory_limit": "512MB", "preserve_insertion_order": "true",
                "autoinstall_known_extensions": "false", "autoload_known_extensions": "false",
            }))
            connection.execute("CREATE TEMP TABLE queries (chromosome VARCHAR, position INTEGER, referenceAllele VARCHAR, alternateAllele VARCHAR, variantId VARCHAR)")
            if lookups:
                connection.executemany("INSERT INTO queries VALUES (?, ?, ?, ?, ?)",
                                       [(*key, "_".join(map(str, key))) for key in sorted(lookups)])
            variants, variant_records_by_allele = {}, {}
            for artifact in artifacts:
                if artifact.dataset != "variant":
                    continue
                for ordinal, row in _rows(connection, release / artifact.relative_path,
                    "(chromosome, position, referenceAllele, alternateAllele) IN "
                    "(SELECT chromosome, position, referenceAllele, alternateAllele FROM queries) "
                    "OR variantId IN (SELECT variantId FROM queries)", []):
                    key = parse_open_targets_variant_row(row)
                    if key not in lookups:
                        raise ValueError
                    record = _source_record("Open Targets variant", "26.06", artifact.relative_path, ordinal,
                                            artifact.sha256, source.manifest_sha256, row)
                    variant_id = row["variantId"]
                    if variant_id in variants:
                        raise ValueError
                    variants[variant_id] = (key, record)
                    for allele in lookups[key]:
                        variant_records_by_allele.setdefault(allele.allele_id, []).append(record)
            connection.execute("CREATE TEMP TABLE matched_variants (variantId VARCHAR)")
            if variants:
                connection.executemany("INSERT INTO matched_variants VALUES (?)", [(v,) for v in sorted(variants)])
            asserted = set()
            for artifact in artifacts:
                if artifact.dataset != "evidence_eva":
                    continue
                for ordinal, row in _rows(connection, release / artifact.relative_path,
                    "variantId IN (SELECT variantId FROM matched_variants)", []):
                    target, condition = parse_open_targets_eva_row(row)
                    key, variant_record = variants[row["variantId"]]
                    record = _source_record("Open Targets evidence_eva", "26.06", artifact.relative_path, ordinal,
                                            artifact.sha256, source.manifest_sha256, row)
                    for allele in lookups[key]:
                        asserted.add(allele.allele_id)
                        outcome = "admitted" if target or condition else "unresolved"
                        reasons = ("exact-variant-key-to-source-assertion" if target or condition else "source-target-and-condition-unresolved",)
                        if condition is None:
                            reasons += ("native-condition-unresolved; mapped-disease-not-equivalence",)
                        rows = (variant_record, record, *(r for r in variant_records_by_allele[allele.allele_id]
                                if r.source_record_id != variant_record.source_record_id))
                        yield _proof(allele, "explicit_clinical_link", outcome, reasons,
                                     rows, provenance, target, condition)
            for allele in alleles:
                if allele.allele_id not in asserted:
                    rows = tuple(variant_records_by_allele.get(allele.allele_id, ()))
                    yield _proof(allele, "explicit_clinical_link", "unresolved",
                                 (unresolved.get(allele.allele_id, "no-matched-source-assertion; clinical-significance-unknown"),), rows, provenance)
    except Exception:
        raise AdmissionLinkError("offline Open Targets association report failed") from None


def clinvar_clinical_report(
    alleles: tuple[AlleleEvidence, ...], *, source: PinnedPublicSource, reference: PinnedPublicSource,
) -> Iterator[AssociationProof]:
    """Preserve matched pinned VCF aggregates; never invent assertion links.

    The independent VCF/TBI pins and complete dated header are checked even for
    no-match queries. Scanning the compressed VCF avoids trusting tabix coverage
    as proof of condition-specific assertions. No aggregate is called benign by
    absence, and no ClinVar aggregate is sufficient for admission by this rule.
    """
    try:
        alleles = _inputs(alleles)
        with ExitStack() as stack:
            files = _Files(stack)
            maps, fasta, ref_proof = _reference(files, reference)
            data, manifest = files.source(source)
            contract = cv.ClinvarContract.from_mapping(data)
            if contract.release != "2026-08-22" or contract.assembly != "GRCh38":
                raise ValueError
            reference_contract = load_reference_identity(reference.contract_path)
            contract.reference_binding.resolve(reference_contract, required_role="public_annotation_benchmark")
            cv._verify_manifest_contract(manifest, contract, source.contract_sha256)
            for artifact in contract.artifacts:
                output = next(o for o in manifest.outputs if o.path == artifact.path)
                files.open(source.root / artifact.path, output.sha256, output.size_bytes)
                cv._verify_artifact(source.root / artifact.path, artifact)
            artifact = contract.artifact("vcf")
            path = source.root / artifact.path
            cv._read_vcf_header(path, contract)
            provenance = _provenance("clinvar-grch38-2026-08-22-vcf", [{"contract_sha256": source.contract_sha256,
                "manifest_sha256": source.manifest_sha256, "artifacts": [_plain(a) for a in contract.artifacts],
                "scope": "site-aggregates-only; condition-assertion-link-unresolved"}], ref_proof)
            lookups, reasons = {}, {}
            for allele in alleles:
                sequence, reason = _sequence(allele, maps, fasta)
                if sequence is None or not _DNA.fullmatch(sequence[2]) or not _DNA.fullmatch(sequence[3]):
                    reasons[allele.allele_id] = reason or "ambiguous-sequence-exact-link-unresolved"
                    continue
                chrom, pos, ref, alt = sequence
                key = (_public_query_chromosome(chrom), pos, ref, alt)
                lookups.setdefault(key, []).append(allele)
            matched = set()
            for ordinal, text in files.lines(path):
                if text.startswith("#"):
                    continue
                row = parse_clinvar_vcf_row(text)
                fields = row["columns"]
                key = (fields[0], int(fields[1]), fields[3], fields[4])
                if key not in lookups:
                    continue
                record = _source_record("ClinVar VCF", contract.release, artifact.path, ordinal,
                                        artifact.sha256, source.manifest_sha256, row)
                for allele in lookups[key]:
                    matched.add(allele.allele_id)
                    yield _proof(allele, "explicit_clinical_link", "unresolved",
                                 ("exact-allele-aggregate-retained; condition-assertion-link-unresolved",), (record,), provenance)
            for allele in alleles:
                if allele.allele_id not in matched:
                    yield _proof(allele, "explicit_clinical_link", "unresolved",
                                 (reasons.get(allele.allele_id, "no-matched-source-assertion; clinical-significance-unknown"),), (), provenance)
    except Exception:
        raise AdmissionLinkError("offline ClinVar association report failed") from None


def association_report_as_json(proofs: tuple[AssociationProof, ...]) -> str:
    """Explicit controlled serialization; not an authenticity or completion token."""
    return _json({"schema": SCHEMA, "classification": "controlled", "llm_policy": "no-llm", "proofs": proofs}) + "\n"
