"""Independently replay small publisher-bound GRCh38 locus classifications.

This resource classifies reference spans. It supplies no patient sex, dosage,
carriage, disease requirement, partner exclusion, or sequence equivalence.
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import date
import hashlib
from html.parser import HTMLParser
from pathlib import Path, PurePosixPath
import re
from weakref import WeakKeyDictionary

from genosapiens.allele_identity import ReferenceRole
from genosapiens.bounded_json import load_bounded_json, read_bounded_bytes
from genosapiens.evidence_model import SourcedRecord
from genosapiens.evidence_serialization import _canonical, _loads
from genosapiens.genotype_model import LOCUS_SCHEMA, ModelLocusRegion
from genosapiens.manifest import (
    PublicResourceManifest, ResourceOutput, RetentionPolicy, SourceIdentity,
    Transformation, publish_manifest,
)
from genosapiens.reference_concordance import (
    ReferenceIdentityContract, ReferenceConstructionReport, verify_reference_manifest_contract,
    verify_reference_construction,
)
from genosapiens.storage import CapacityPlan, require_capacity


SCHEMA = "genosapiens.grch38-model-locus-resource"
TRANSFORM = {
    "schema": "genosapiens.grch38-model-locus-projection",
    "coordinates": "publisher-1-based-closed; no endpoint adjustment",
    "classification": "NCBI-human-assembled-molecule; NHGRI-numbered-autosome; NCBI-PAR-with-Ensembl-116-agreement",
    "reference_mapping": "accepted-GCA_000001405.15-contract-and-construction; source-UCSC-name; explicit-role-prefix-removal; exact-FAI-length; public-validation-single-leading-chr-source-alias-with-bound-reference-target",
    "span": "partition-X-Y-at-source-PAR-boundaries; unknown-other-sequences",
    "scope": "locus-classification-only; no-sex-dosage-carriage-count-or-partner-exclusion",
}
TRANSFORM_SHA256 = hashlib.sha256((_canonical(TRANSFORM) + "\n").encode()).hexdigest()
_BASE = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/"
_ARCHIVE = "https://jun2026.archive.ensembl.org/"
SOURCE_URLS = {
    "assembly_report": _BASE + "GCA_000001405.15_GRCh38_assembly_report.txt",
    "assembly_regions": _BASE + "GCA_000001405.15_GRCh38_assembly_regions.txt",
    "ensembl_par": _ARCHIVE + "info/docs/genebuild/human_PARS.html",
    "ensembl_coordinates": _ARCHIVE + "info/data/index.html",
    "ensembl_permission": _ARCHIVE + "info/about/legal/disclaimer.html",
    "nhgri_autosome": "https://www.genome.gov/genetics-glossary/Autosome",
    "nhgri_permission": "https://www.genome.gov/about-nhgri/Policies-Guidance/Copyright",
    "ncbi_permission": "https://www.ncbi.nlm.nih.gov/home/about/policies/",
}
MAX_BYTES = 2 * 1024 * 1024
LICENSE = "NCBI molecular-data unrestricted-use policy; NHGRI public domain; Ensembl project data unrestricted; acknowledgements required"
_OWNERS = WeakKeyDictionary()


class ModelLocusResourceError(ValueError):
    """Static resource diagnostics, never untrusted source or case payloads."""


def _bytes(value):
    return (_canonical(value) + "\n").encode("ascii")


def _sha(content):
    return hashlib.sha256(content).hexdigest()


def _relative(value, prefix):
    p = PurePosixPath(value)
    if p.is_absolute() or ".." in p.parts or not p.parts or not p.is_relative_to(prefix):
        raise ModelLocusResourceError("model-locus path is invalid")
    return p.as_posix()


def _bound(path, digest, *, root, size=None):
    if type(digest) is not str or re.fullmatch("[0-9a-f]{64}", digest) is None:
        raise ModelLocusResourceError("model-locus identity is invalid")
    content = read_bounded_bytes(path, reference_root=root, maximum_bytes=MAX_BYTES,
                                 expected_size_bytes=size)
    if _sha(content) != digest:
        raise ModelLocusResourceError("model-locus resource binding is stale")
    return content


def _descriptor(path, repository_root, *, expected_sha256=None, publishing=False):
    data, content = load_bounded_json(path, reference_root=repository_root)
    if expected_sha256 is not None and _sha(content) != expected_sha256:
        raise ModelLocusResourceError("model-locus descriptor is stale")
    if set(data) != {"schema", "classification", "resource_id", "release", "transformation_sha256",
                     "sources", "reference", "output_path", "manifest_path", "manifest_sha256"}:
        raise ModelLocusResourceError("model-locus descriptor fields are invalid")
    if (data["schema"] != SCHEMA or data["classification"] != "public"
            or data["transformation_sha256"] != TRANSFORM_SHA256
            or data["resource_id"] != "grch38-authoritative-model-loci"
            or data["release"] != "grch38-ensembl-116-nhgri-retrieved-2026-09-16"
            or set(data["sources"]) != set(SOURCE_URLS)):
        raise ModelLocusResourceError("model-locus descriptor contract is invalid")
    for role, source in data["sources"].items():
        if (set(source) != {"url", "path", "sha256", "size_bytes", "release"}
                or source["url"] != SOURCE_URLS[role]
                or type(source["size_bytes"]) is not int
                or not 0 < source["size_bytes"] <= MAX_BYTES
                or source["release"] != ("GRCh38-GCA_000001405.15" if role.startswith("assembly_")
                    else "Ensembl-116-June-2026" if role.startswith("ensembl_")
                    else "retrieved-2026-09-16")):
            raise ModelLocusResourceError("model-locus publisher binding is invalid")
        _relative(source["path"], "releases/model-loci")
    if set(data["reference"]) != {"identity_path", "identity_sha256", "manifest_path", "manifest_sha256"}:
        raise ModelLocusResourceError("model-locus reference contract is invalid")
    _relative(data["reference"]["identity_path"], "config/reference-identities")
    _relative(data["reference"]["manifest_path"], "manifests/grch38")
    _relative(data["output_path"], "derived/model-loci")
    _relative(data["manifest_path"], "manifests/model-loci")
    if not publishing and re.fullmatch("[0-9a-f]{64}", data["manifest_sha256"] or "") is None:
        raise ModelLocusResourceError("model-locus publication identity is missing")
    return data, content


class _Text(HTMLParser):
    def __init__(self):
        super().__init__()
        self.parts, self.skip = [], 0

    def handle_starttag(self, tag, attrs):
        if tag in ("script", "style"):
            self.skip += 1

    def handle_endtag(self, tag):
        if tag in ("script", "style"):
            self.skip = max(0, self.skip - 1)

    def handle_data(self, value):
        if not self.skip:
            self.parts.append(value)


def _text(content):
    reader = _Text()
    reader.feed(content.decode("utf-8"))
    return " ".join(" ".join(reader.parts).split())


def _source_facts(sources):
    """Derive publisher facts from retained bytes, never caller projection rows."""
    report, regions = (sources[k].decode("ascii") for k in ("assembly_report", "assembly_regions"))
    for value in (report, regions):
        headers = dict((m[1].strip(), m[2].strip()) for line in value.splitlines()
                       if (m := re.fullmatch(r"# ([^:]+):\s*(.*)", line)))
        if (headers.get("Assembly name") != "GRCh38" or headers.get("Taxid") != "9606"
                or headers.get("GenBank assembly accession") != "GCA_000001405.15"
                or headers.get("Submitter") != "Genome Reference Consortium"):
            raise ModelLocusResourceError("model-locus source assembly is invalid")
    human = _text(sources["nhgri_autosome"])
    if not all(anchor in human for anchor in ("numbered chromosomes", "22 pairs of autosomes", "XX or XY")):
        raise ModelLocusResourceError("human chromosome classification source is invalid")
    permissions = {k: _text(sources[k]) for k in ("ncbi_permission", "nhgri_permission", "ensembl_permission")}
    if ("no restrictions on the use or distribution" not in permissions["ncbi_permission"]
            or "public domain, unless otherwise indicated" not in permissions["nhgri_permission"]
            or "available without restriction" not in permissions["ensembl_permission"]):
        raise ModelLocusResourceError("model-locus source permission is invalid")
    for k in ("ensembl_par", "ensembl_coordinates", "ensembl_permission"):
        if "Ensembl release 116" not in _text(sources[k]):
            raise ModelLocusResourceError("model-locus Ensembl release is invalid")
    if "Ensembl coordinates always have a one-based start" not in _text(sources["ensembl_coordinates"]):
        raise ModelLocusResourceError("model-locus coordinate convention is invalid")
    molecules = {}
    for line in report.splitlines():
        if not line or line.startswith("#"):
            continue
        row = line.split("\t")
        if len(row) != 10:
            raise ModelLocusResourceError("model-locus assembly row is malformed")
        if row[1] != "assembled-molecule":
            continue
        molecule = row[2]
        if molecule in molecules or int(row[8]) <= 0:
            raise ModelLocusResourceError("model-locus molecule is invalid")
        if molecule in {str(i) for i in range(1, 23)} | {"X", "Y"}:
            if row[3] != "Chromosome" or row[7] != "Primary Assembly" or row[9] != "chr" + molecule:
                raise ModelLocusResourceError("model-locus primary chromosome is invalid")
        elif molecule == "MT":
            if row[3] != "Mitochondrion" or row[7] != "non-nuclear" or row[9] != "chrM":
                raise ModelLocusResourceError("model-locus mitochondrial source is invalid")
        else:
            raise ModelLocusResourceError("model-locus assembled molecule is unsupported")
        molecules[molecule] = row
    if set(molecules) != {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}:
        raise ModelLocusResourceError("model-locus source molecule inventory is incomplete")
    par = {"X": [], "Y": []}
    for line in regions.splitlines():
        if not line or line.startswith("#"):
            continue
        row = line.split("\t")
        if len(row) != 8:
            raise ModelLocusResourceError("model-locus region row is malformed")
        if row[4] == "PAR":
            if row[1] not in par or row[0] not in ("PAR#1", "PAR#2"):
                raise ModelLocusResourceError("model-locus PAR source is invalid")
            par[row[1]].append(row)
    documented = re.findall(r"chromosome:GRCh38:Y:(\d+)\s*-\s*(\d+) is shared with X:\s*(\d+)\s*-\s*(\d+) \(PAR([12])\)", _text(sources["ensembl_par"]))
    if len(documented) != 2 or len({r[4] for r in documented}) != 2:
        raise ModelLocusResourceError("model-locus PAR documentation is incomplete")
    for contig, rows in par.items():
        rows.sort(key=lambda r: int(r[2]))
        if len(rows) != 2 or {r[0] for r in rows} != {"PAR#1", "PAR#2"}:
            raise ModelLocusResourceError("model-locus PAR inventory is incomplete")
        previous = 0
        for row in rows:
            start, end = int(row[2]), int(row[3])
            match = next(r for r in documented if "PAR#" + r[4] == row[0])
            expected = tuple(map(int, match[:2] if contig == "Y" else match[2:4]))
            if (not previous < start <= end <= int(molecules[contig][8])
                    or (start, end) != expected):
                raise ModelLocusResourceError("model-locus PAR geometry conflicts")
            previous = end
    return molecules, par


def _reference(data, repository, root):
    ref = data["reference"]
    identity_path = repository / ref["identity_path"]
    identity_content = _bound(identity_path, ref["identity_sha256"], root=repository)
    contract = ReferenceIdentityContract.from_mapping(_loads(identity_content))
    manifest_path = root / ref["manifest_path"]
    raw = _loads(_bound(manifest_path, ref["manifest_sha256"], root=root))
    manifest = PublicResourceManifest.from_mapping(raw)
    verify_reference_manifest_contract(contract, manifest, contract_sha256=ref["identity_sha256"])
    if (contract.primary_assembly_accession != "GCA_000001405.15"
            or contract.transformation.source_prefix != "chr"):
        raise ModelLocusResourceError("model-locus reference assembly is unsupported")
    outputs = {o.path: o for o in manifest.outputs}
    report_path = contract.expected_output("construction_report").path
    original_report = ReferenceConstructionReport.from_mapping(_loads(_bound(
        root / report_path, outputs[report_path].sha256, root=root, size=outputs[report_path].size_bytes)))
    confirmed = verify_reference_construction(contract, manifest, root, contract_sha256=ref["identity_sha256"])
    if confirmed != original_report:
        raise ModelLocusResourceError("model-locus reference verification did not bind the authenticated bytes")
    indexes = {}
    for role in contract.reference_roles:
        path = contract.expected_output(role.fai_output_role).path
        content = _bound(root / path, outputs[path].sha256, root=root, size=outputs[path].size_bytes)
        entries = {}
        for line in content.decode("ascii").splitlines():
            fields = line.split("\t")
            if len(fields) != 5 or not fields[0] or fields[0] in entries:
                raise ModelLocusResourceError("model-locus reference index is invalid")
            length, offset, bases, width = map(int, fields[1:])
            if length <= 0 or offset < 0 or not 0 < bases <= width:
                raise ModelLocusResourceError("model-locus reference index geometry is invalid")
            entries[fields[0]] = length
        if not entries:
            raise ModelLocusResourceError("model-locus reference index is empty")
        indexes[role.role] = entries
    return contract, indexes


def _projection(data, repository, root):
    originals = {role: _bound(root / source["path"], source["sha256"], root=root,
                              size=source["size_bytes"]) for role, source in data["sources"].items()}
    molecules, par = _source_facts(originals)
    contract, indexes = _reference(data, repository, root)
    rows = []
    for role in contract.reference_roles:
        if role.naming not in ("unprefixed", "chr-prefixed", "resource-native"):
            raise ModelLocusResourceError("model-locus reference naming is unsupported")
        for molecule, original in molecules.items():
            source_name, length = original[9], int(original[8])
            name = source_name[3:] if role.naming == "unprefixed" else source_name
            if name not in indexes[role.role]:
                continue  # No resemblance-based classification of other contigs.
            if indexes[role.role][name] != length:
                raise ModelLocusResourceError("model-locus reference length differs from publisher")
            parts = []
            if molecule in ("X", "Y"):
                cursor = 1
                for region in par[molecule]:
                    start, end = int(region[2]), int(region[3])
                    if cursor < start:
                        parts.append((cursor, start - 1, molecule.lower() + "_nonpar"))
                    parts.append((start, end, "par"))
                    cursor = end + 1
                if cursor <= length:
                    parts.append((cursor, length, molecule.lower() + "_nonpar"))
            else:
                parts = [(1, length, "mitochondrial" if molecule == "MT" else "autosome")]
            names = [(name, "reference-dictionary-name")]
            # Public concordance fetches through this declared mapping while
            # native allele identities retain the original source name. Publish
            # exact literal rows for both names; consumers still match all four
            # reference/contig bindings without performing alias inference.
            if role.role == "public_annotation_benchmark" and name.startswith("chr"):
                source_alias = name[3:]
                if not source_alias or source_alias in indexes[role.role]:
                    raise ModelLocusResourceError("model-locus public source alias is ambiguous")
                names.append((source_alias, "public-validation-remove-one-leading-chr"))
            for source_contig, recipe in names:
                for start, end, kind in parts:
                    projection = {"schema": LOCUS_SCHEMA, "coordinate_system": "1-based-closed",
                        "reference_role": role.role, "reference_identity_id": role.identity_id,
                        "reference_manifest_sha256": data["reference"]["manifest_sha256"],
                        "contig": source_contig, "start": start, "end": end, "kind": kind}
                    rows.append({"projection": projection, "assembly_source_row": original,
                                 "par_source_rows": par.get(molecule, []),
                                 "reference_name_mapping": {"source_contig": source_contig,
                                     "reference_contig": name, "recipe": recipe}})
    if not rows:
        raise ModelLocusResourceError("model-locus publication has no supported reference regions")
    rows.sort(key=lambda r: (r["projection"]["reference_role"], r["projection"]["contig"], r["projection"]["start"]))
    return {"schema": SCHEMA, "classification": "public", "transformation_sha256": TRANSFORM_SHA256,
            "sources": data["sources"], "reference": data["reference"], "regions": rows}


def _parameters(data):
    return _sha(_bytes({k: data[k] for k in ("resource_id", "release", "transformation_sha256", "sources", "reference", "output_path")}))


@dataclass(frozen=True, slots=True, weakref_slot=True, repr=False, eq=False)
class ModelLocusResource:
    descriptor_sha256: str
    manifest_sha256: str
    output_sha256: str
    regions: tuple[ModelLocusRegion, ...]

    def __repr__(self):
        return "<authenticated public model-locus resource>"


def authenticated_model_loci(resource):
    if type(resource) is not ModelLocusResource or _OWNERS.get(resource) != _sha(_bytes(resource)):
        raise ModelLocusResourceError("model-locus resource must be independently loaded in this execution")
    return resource.regions


def load_model_locus_resource(descriptor_path, *, reference_root, repository_root, expected_descriptor_sha256):
    """Authenticate all sources, reference mapping and deterministic projection."""
    try:
        if type(expected_descriptor_sha256) is not str or re.fullmatch("[0-9a-f]{64}", expected_descriptor_sha256) is None:
            raise ModelLocusResourceError("model-locus descriptor identity is required")
        root, repo = Path(reference_root).absolute(), Path(repository_root).absolute()
        data, descriptor_bytes = _descriptor(Path(descriptor_path).absolute(), repo,
                                              expected_sha256=expected_descriptor_sha256)
        manifest_path = root / data["manifest_path"]
        raw = _loads(_bound(manifest_path, data["manifest_sha256"], root=root))
        manifest = PublicResourceManifest.from_mapping(raw)
        expected_paths = {s["path"] for s in data["sources"].values()} | {data["output_path"]}
        assembly = data["sources"]["assembly_report"]
        if (manifest.resource_id != data["resource_id"] or manifest.release != data["release"]
                or manifest.assembly != "GRCh38" or manifest.license != LICENSE
                or manifest.source != SourceIdentity(assembly["url"], assembly["sha256"])
                or manifest.transformation != Transformation(TRANSFORM["schema"], _parameters(data))
                or {o.path for o in manifest.outputs} != expected_paths
                or any(o.retention != RetentionPolicy("durable", None) for o in manifest.outputs)):
            raise ModelLocusResourceError("model-locus publication contract differs")
        for output in manifest.outputs:
            _bound(root / output.path, output.sha256, root=root, size=output.size_bytes)
        expected = _projection(data, repo, root)
        output = next(o for o in manifest.outputs if o.path == data["output_path"])
        content = _bound(root / output.path, output.sha256, root=root, size=output.size_bytes)
        if content != _bytes(expected):
            raise ModelLocusResourceError("model-locus projection does not replay publisher sources")
        regions = []
        for i, row in enumerate(expected["regions"]):
            p = row["projection"]
            source = SourcedRecord("GRC/NCBI GRCh38; Ensembl 116; NHGRI chromosome classification",
                data["release"], "model-locus-" + str(i), assembly["sha256"], data["manifest_sha256"],
                {"model_locus_region": p, "publisher_source_rows": row,
                 "publisher_sources": data["sources"], "transformation_sha256": TRANSFORM_SHA256})
            regions.append(ModelLocusRegion(ReferenceRole(p["reference_role"]), p["reference_identity_id"],
                p["reference_manifest_sha256"], p["contig"], p["start"], p["end"], p["kind"], source))
        resource = ModelLocusResource(_sha(descriptor_bytes), data["manifest_sha256"], output.sha256, tuple(regions))
        _OWNERS[resource] = _sha(_bytes(resource))
        return resource
    except Exception:
        raise ModelLocusResourceError("model-locus source, reference or publication qualification failed") from None


def qualify_model_locus_resource(descriptor_path, *, reference_root, repository_root, image_digest, retrieved_on):
    """Publish a new replayed projection; source acquisition is separate."""
    root, repo = Path(reference_root).absolute(), Path(repository_root).absolute()
    data, _ = _descriptor(Path(descriptor_path).absolute(), repo, publishing=True)
    date.fromisoformat(retrieved_on)
    require_capacity(root / "derived/model-loci", CapacityPlan(operation="qualify small model-locus projection",
        source_bytes=0, output_bytes=MAX_BYTES, scratch_bytes=MAX_BYTES))
    projected = _bytes(_projection(data, repo, root))
    from genosapiens.track1_public_validation_coordinator import publish_validation_artifact
    # Existing durable immutable writer rejects symlink ancestry and byte changes.
    publish_validation_artifact(root / data["output_path"], projected)
    outputs = tuple(ResourceOutput(s["path"], s["size_bytes"], s["sha256"], RetentionPolicy("durable", None))
                    for _, s in sorted(data["sources"].items())) + (
        ResourceOutput(data["output_path"], len(projected), _sha(projected), RetentionPolicy("durable", None)),)
    assembly = data["sources"]["assembly_report"]
    manifest = PublicResourceManifest(data["resource_id"], data["release"],
        SourceIdentity(assembly["url"], assembly["sha256"]), LICENSE, "GRCh38", retrieved_on,
        Transformation(TRANSFORM["schema"], _parameters(data)), image_digest, outputs)
    publication = publish_manifest(manifest, root / data["manifest_path"], root)
    return publication
