"""Private schema compatibility only, never case acceptance or an LLM exception.

Public contracts are operator-trusted configuration, not caller-supplied allele
containers. Run offline with public resources read-only. No body, locus, count,
distribution, or genotype override is exposed by the generation API.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass, field
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import stat
from typing import Mapping

from genosapiens.bounded_json import BoundedJSONError, load_bounded_json
from genosapiens.challenge_fixture import (
    ChallengeFixtureContract, ChallengeFixtureError, _bgzf_block, _checksums,
    _load_aliases, _profile_fields, _select_public_records, _source_identity,
    _tabix_bytes,
)
from genosapiens.giab import GiabError, GiabResourceContract
from genosapiens.governance import (
    ArtifactClassification, ArtifactRecord, CaseInventory, GovernanceError,
    validate_relative_classification,
)
from genosapiens.manifest import ManifestError, PublicResourceManifest
from genosapiens.reference_concordance import (
    ReferenceConcordanceError, ReferenceIdentityContract, open_installed_reference,
)
from genosapiens.tabix import TabixError, TabixIndex
from genosapiens.vcf_header import (
    PUBLIC_ASSEMBLY_LABEL, SkeletonCompatibility, VcfHeaderError,
    compare_skeleton_dictionary, read_schema_header,
)


SCHEMA = "genosapiens.controlled-compatibility-fixture"
PROFILE_SCHEMA = "genosapiens.compatibility-genotypes"
SUPPORTED_FORMATS = ("GT", "DP", "AD", "GQ", "PL")
BGZF_PAYLOAD_SIZE = 60 * 1024
POLICY = {"classification": "controlled", "llm_policy": "no-llm", "retention": "delete"}


class CompatibilityFixtureError(RuntimeError):
    """Safe error codes only. Dictionary details, when present, remain controlled."""

    def __init__(self, code: str, compatibility: SkeletonCompatibility | None = None):
        super().__init__(code)
        self.compatibility = compatibility


@dataclass(frozen=True, slots=True, repr=False)
class CompatibilityFixtureResult:
    fixture_vcf: str
    fixture_tbi: str
    report_path: str
    inventory_path: str
    inventory: CaseInventory
    classification: str = field(default="controlled", init=False)
    llm_policy: str = field(default="no-llm", init=False)
    retention: str = field(default="delete", init=False)


def _json_bytes(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def load_compatibility_inventory(
    reference_root: str | Path, inventory_path: str | Path
) -> CaseInventory:
    """Recover the deletion ledger AND its own identity after process interruption."""
    try:
        root = Path(reference_root).resolve(strict=True)
        path = Path(inventory_path).resolve(strict=True)
        relative = path.relative_to(root).as_posix()
        validate_relative_classification(relative, ArtifactClassification.CONTROLLED)
        data, content = load_bounded_json(path)
        if not isinstance(data, Mapping) or any(data.get(key) != value for key, value in POLICY.items()):
            raise CompatibilityFixtureError("compatibility-inventory-invalid")
        inventory = CaseInventory.from_mapping(data)
        expected = {"schema.header.vcf", "compatibility.vcf.gz", "compatibility.vcf.gz.tbi", "report.json"}
        if (
            path.name != "inventory.json" or len(inventory.artifacts) != len(expected)
            or {Path(item.relative_path).name for item in inventory.artifacts} != expected
            or any(item.classification is not ArtifactClassification.CONTROLLED
                   or Path(item.relative_path).parent != Path(relative).parent
                   for item in inventory.artifacts)
        ):
            raise CompatibilityFixtureError("compatibility-inventory-invalid")
        digest = hashlib.sha256(content).hexdigest()
        ledger = ArtifactRecord(
            artifact_id=hashlib.sha256(f"controlled\0{relative}\0{digest}".encode()).hexdigest()[:24],
            relative_path=relative, classification=ArtifactClassification.CONTROLLED,
            kind="compatibility-inventory", sha256=digest, size_bytes=len(content), retention="delete",
        )
        return CaseInventory(inventory.case_id, (*inventory.artifacts, ledger))
    except (BoundedJSONError, GovernanceError, OSError, ValueError, TypeError, RuntimeError):
        raise CompatibilityFixtureError("compatibility-inventory-invalid") from None


def _write_private_file(directory_fd: int, name: str, content: bytes, created: list[str]) -> None:
    descriptor = os.open(
        name, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
        0o600, dir_fd=directory_fd,
    )
    created.append(name)
    with os.fdopen(descriptor, "wb") as output:
        output.write(content)
        output.flush()
        os.fsync(output.fileno())


def _public_path(root: Path, relative: str) -> Path:
    validate_relative_classification(relative, ArtifactClassification.PUBLIC)
    path = root / relative
    # Check before any resource read; public-looking symlinks are not provenance.
    if any(part.is_symlink() for part in (path, *path.parents) if part != root):
        raise CompatibilityFixtureError("public-path-invalid")
    resolved = path.resolve(strict=True)
    validate_relative_classification(
        resolved.relative_to(root).as_posix(), ArtifactClassification.PUBLIC
    )
    if not resolved.is_file():
        raise CompatibilityFixtureError("public-path-invalid")
    return resolved


def _encode_vcf(header: bytes, lines: tuple[bytes, ...]) -> tuple[bytes, bytes]:
    """Chunk both header and records; translate logical offsets to BGZF offsets."""
    payload = header + b"".join(lines)
    blocks = []
    compressed_offsets = []
    compressed_size = 0
    for start in range(0, len(payload), BGZF_PAYLOAD_SIZE):
        block = _bgzf_block(payload[start:start + BGZF_PAYLOAD_SIZE])
        compressed_offsets.append(compressed_size)
        blocks.append(block)
        compressed_size += len(block)
    # Also covers a record ending exactly at a block boundary.
    compressed_offsets.append(compressed_size)
    offsets = {}
    start = len(header)
    for line in lines:
        contig = line.split(b"\t", 1)[0].decode("ascii")
        end = start + len(line)
        start_block, start_within = divmod(start, BGZF_PAYLOAD_SIZE)
        end_block, end_within = divmod(end, BGZF_PAYLOAD_SIZE)
        begin_virtual = (compressed_offsets[start_block] << 16) | start_within
        end_virtual = (compressed_offsets[end_block] << 16) | end_within
        previous = offsets.get(contig)
        offsets[contig] = (previous[0] if previous else begin_virtual, end_virtual)
        start = end
    # Five fixed public profiles need at most five reference chunks. Bin zero
    # deliberately covers whole contigs, with no linear pruning, for this tiny body.
    return b"".join(blocks) + _bgzf_block(b""), _tabix_bytes(offsets)


def _validate_record(line: bytes, contigs: Mapping[str, int], formats: tuple[str, ...],
                     filters: set[str]) -> None:
    """Validate the emitted sequence-only, diploid subset, including vector sizes."""
    fields = line.decode("ascii").rstrip("\n").split("\t")
    if len(fields) != 10:
        raise CompatibilityFixtureError("generated-vcf-syntax-invalid")
    chrom, pos, identifier, ref, alt, qual, filter_value, info, fmt, sample = fields
    alts = alt.split(",")
    if (
        chrom not in contigs or not re.fullmatch(r"[1-9][0-9]*", pos)
        or int(pos) + len(ref) - 1 > min(contigs[chrom], 1 << 29)
        or any(not re.fullmatch(r"[ACGTN]+", allele) for allele in [ref, *alts])
        or ref in alts or len(alts) != len(set(alts)) or identifier != "." or info != "."
        or (qual != "." and not re.fullmatch(r"[0-9]+(?:\.[0-9]+)?", qual))
        or filter_value not in {".", "PASS", *filters}
        or fmt != (":".join(formats) or ".")
        or (not formats and sample != ".")
    ):
        raise CompatibilityFixtureError("generated-vcf-syntax-invalid")
    values = dict(zip(formats, sample.split(":"), strict=True)) if formats else {}
    for name, value in values.items():
        if name == "GT":
            if not re.fullmatch(r"[0-9]+[/|][0-9]+", value) or any(
                int(allele) > len(alts) for allele in re.split(r"[/|]", value)
            ):
                raise CompatibilityFixtureError("generated-genotype-invalid")
            continue
        if value == ".":
            continue
        vector = value.split(",")
        expected = len(alts) + 1 if name == "AD" else (
            (len(alts) + 1) * (len(alts) + 2) // 2 if name == "PL" else 1
        )
        if len(vector) != expected or any(
            item != "." and not re.fullmatch(r"[0-9]+", item) for item in vector
        ):
            raise CompatibilityFixtureError("generated-genotype-invalid")
    if "AD" in values and "DP" in values and "." not in values["AD"] + values["DP"]:
        if sum(int(item) for item in values["AD"].split(",")) != int(values["DP"]):
            raise CompatibilityFixtureError("generated-depth-invalid")
    if "PL" in values and values["PL"] != ".":
        likelihoods = [int(item) for item in values["PL"].split(",")]
        if min(likelihoods) != 0:
            raise CompatibilityFixtureError("generated-likelihood-invalid")
        if "GT" in values:
            a, b = sorted(int(item) for item in re.split(r"[/|]", values["GT"]))
            if likelihoods[b * (b + 1) // 2 + a] != 0:
                raise CompatibilityFixtureError("generated-likelihood-invalid")


def generate_compatibility_fixture(
    *,
    schema_header_path: str | Path,
    public_fixture_contract_path: str | Path,
    giab_contract_path: str | Path,
    reference_identity_path: str | Path,
    reference_root: str | Path,
    output_directory: str | Path,
    image_digest: str,
) -> CompatibilityFixtureResult:
    """Generate in a NEW private controlled directory using trusted public contracts.

    Preserves the whole sanitized header byte-for-byte. Incompatible dictionaries
    (including order-only differences) raise with a controlled comparison attached;
    callers must not log that comparison. No adjudication override is accepted.
    ``report.json`` carries per-allele provenance and declaration exercise counts.
    ``inventory.json`` inventories the products; load it using
    ``load_compatibility_inventory`` to include the ledger itself in deletion,
    including after interruption. The return value uses that same recovery path.
    Nothing returned by this function is an LLM-readable exception.
    """
    created: list[str] = []
    target_created = False
    parent_fd = directory_fd = -1
    target = Path(output_directory)
    try:
        # Validate the ENTIRE bounded input before even opening public resources.
        header = read_schema_header(schema_header_path)
        if not re.fullmatch(r"sha256:[0-9a-f]{64}", image_digest):
            raise CompatibilityFixtureError("image-binding-invalid")
        root = Path(reference_root).resolve(strict=True)
        relative = target.absolute().relative_to(root).as_posix()
        validate_relative_classification(relative, ArtifactClassification.CONTROLLED)
        if target.exists() or any(part.is_symlink() for part in (target, *target.parents)):
            raise CompatibilityFixtureError("output-directory-invalid")
        target = root / relative

        fixture_data, fixture_content = load_bounded_json(public_fixture_contract_path)
        giab_data, giab_content = load_bounded_json(giab_contract_path)
        reference_data, reference_content = load_bounded_json(reference_identity_path)
        fixture = ChallengeFixtureContract.from_mapping(fixture_data)
        giab = GiabResourceContract.from_mapping(giab_data)
        reference = ReferenceIdentityContract.from_mapping(reference_data)
        if (
            header.assembly_label != PUBLIC_ASSEMBLY_LABEL
            or reference.primary_assembly_accession != PUBLIC_ASSEMBLY_LABEL
            or giab.sample != "HG002" or giab.assembly != "GRCh38"
            or fixture.assembly != reference.assembly
            or fixture.source_vcf_role != "truth_vcf"
        ):
            raise CompatibilityFixtureError("grch38-public-source-binding-invalid")
        fixture.reference_binding.resolve(reference, required_role="challenge_alignment")
        giab.reference_binding.resolve(reference, required_role="public_annotation_benchmark")
        public_source = giab.artifact("truth_vcf")
        source_path = _public_path(root, public_source.path)
        source_manifest_path = _public_path(root, fixture.source_manifest_path)
        source_manifest_data, source_manifest_content = load_bounded_json(source_manifest_path)
        source_manifest = PublicResourceManifest.from_mapping(source_manifest_data)
        if source_manifest.release != giab.release or source_manifest.assembly != "GRCh38":
            raise CompatibilityFixtureError("public-source-manifest-invalid")
        _, source_sha256 = _source_identity(fixture, giab, root)

        # Preflight every path the existing public construction validator can open.
        for item in (*reference.source_artifacts, *reference.expected_outputs):
            _public_path(root, item.path)
        identity_sha256 = hashlib.sha256(
            Path(reference_identity_path).read_bytes()
        ).hexdigest()
        manifest_path = _public_path(
            root,
            (
                f"manifests/grch38/{reference.resource_id}."
                f"{identity_sha256[:12]}.json"
            ),
        )
        with open_installed_reference(
            reference, identity_path=reference_identity_path,
            reference_manifest_path=manifest_path, reference_root=root,
        ) as (manifest_sha256, fetch, contigs):
            compatibility = compare_skeleton_dictionary(header, contigs, manifest_sha256)
            if compatibility.status != "compatible":
                raise CompatibilityFixtureError("skeleton-" + compatibility.status, compatibility)
            alias_path, aliases = _load_aliases(reference, root)
            alias_data, alias_content = load_bounded_json(alias_path)
            inverse = alias_data.get("alignment_to_tool")
            if (
                set(aliases.values()) != set(contigs) or len(aliases) != len(contigs)
                or inverse != {value: name for name, value in aliases.items()}
                or any(value != name.removeprefix("chr") for name, value in aliases.items())
            ):
                raise CompatibilityFixtureError("public-aliases-invalid")
            records = _select_public_records(source_path, "HG002")
            order = {name: index for index, (name, _length) in enumerate(header.contigs)}
            if any(aliases.get(record.contig) != record.contig.removeprefix("chr")
                   or aliases.get(record.contig) not in order for record in records):
                raise CompatibilityFixtureError("public-allele-alias-invalid")
            records = tuple(sorted(records, key=lambda item: (
                order[aliases[item.contig]], item.position, item.profile,
            )))
            declared_formats = {item.identifier for item in header.declarations if item.category == "FORMAT"}
            formats = tuple(name for name in SUPPORTED_FORMATS if name in declared_formats)
            filters = {item.identifier for item in header.declarations if item.category == "FILTER"}
            exercised: Counter[tuple[str, str]] = Counter()
            populated: Counter[tuple[str, str]] = Counter()
            lines = []
            provenance = []
            for record in records:
                chrom = aliases[record.contig]
                if fetch(chrom, record.position, len(record.ref)) != record.ref:
                    raise CompatibilityFixtureError("public-allele-ref-mismatch")
                quality, filter_value, profile_format, profile_sample = _profile_fields(record)
                values = dict(zip(profile_format.split(":"), profile_sample.split(":"), strict=True))
                if filter_value != "PASS" and filter_value not in filters:
                    filter_value = "."
                line = (
                    f"{chrom}\t{record.position}\t.\t{record.ref}\t{','.join(record.alts)}\t"
                    f"{quality}\t{filter_value}\t.\t{':'.join(formats) or '.'}\t"
                    f"{':'.join(values[name] for name in formats) or '.'}\n"
                ).encode("ascii")
                _validate_record(line, contigs, formats, filters)
                lines.append(line)
                exercised[("contig", chrom)] += 1
                populated[("contig", chrom)] += 1
                exercised[("FILTER", filter_value)] += 1
                populated[("FILTER", filter_value)] += 1
                for name in formats:
                    exercised[("FORMAT", name)] += 1
                    populated[("FORMAT", name)] += int(any(char != "." for char in values[name].split(",")))
                provenance.append({
                    "record_sha256": hashlib.sha256(line).hexdigest(),
                    "chrom": chrom, "position": record.position, "ref": record.ref,
                    "alts": list(record.alts), "profile": record.profile,
                    "allele_sources": [{
                        "kind": "public_giab_hg002", "resource_id": giab.resource_id,
                        "release": giab.release, "source_chrom": record.contig,
                        "source_position": record.position, "source_ref": record.ref,
                        "source_alt": alt, "source_record_sha256": record.source_record_sha256,
                        "source_vcf_sha256": source_sha256,
                    } for alt in record.alts],
                })
            if _checksums(source_path) != (source_sha256, public_source.md5, public_source.size_bytes):
                raise CompatibilityFixtureError("public-source-changed")
            if source_manifest_path.read_bytes() != source_manifest_content:
                raise CompatibilityFixtureError("public-source-manifest-changed")
            if alias_path.read_bytes() != alias_content:
                raise CompatibilityFixtureError("public-aliases-changed")
            vcf, tbi = _encode_vcf(header.content, tuple(lines))
            report = {
                "schema": SCHEMA, **POLICY, "sample": "SYNTHETIC",
                "header_sha256": header.sha256, "compatibility": asdict(compatibility),
                "reference_accepted": False, "primary_rescue_authorized": False,
                "original_reference_basename_verified": False,
                "challenge_body_used": False, "patient_statistics_used": False,
                "public_caller_distributions_used": False,
                "profile_schema": PROFILE_SCHEMA, "selection_policy": fixture.selection_policy,
                "image_digest": image_digest,
                "reference_binding": asdict(fixture.reference_binding),
                "public_reference_binding": asdict(giab.reference_binding),
                "public_inputs": [
                    {"role": role, "sha256": digest, "classification": "public"}
                    for role, digest in (
                        ("selection_contract", hashlib.sha256(fixture_content).hexdigest()),
                        ("giab_contract", hashlib.sha256(giab_content).hexdigest()),
                        ("reference_contract", hashlib.sha256(reference_content).hexdigest()),
                        ("reference_manifest", manifest_sha256),
                        ("giab_manifest", hashlib.sha256(source_manifest_content).hexdigest()),
                        ("truth_vcf", source_sha256),
                        ("contig_aliases", hashlib.sha256(alias_content).hexdigest()),
                    )
                ],
                "public_source": {
                    "resource_id": giab.resource_id, "release": giab.release,
                    "url": public_source.url, "path": public_source.path,
                    "md5": public_source.md5, "sha256": source_sha256,
                    "size_bytes": public_source.size_bytes,
                    "manifest_path": fixture.source_manifest_path,
                },
                "records": provenance,
                "declarations": [{
                    "category": item.category, "id": item.identifier,
                    "attributes": dict(item.attributes),
                    "exercised_records": exercised[(item.category, item.identifier)],
                    "populated_records": populated[(item.category, item.identifier)],
                    "reason": "exercised" if exercised[(item.category, item.identifier)] else (
                        "no-selected-public-allele" if item.category == "contig" else
                        "unsupported-semantics"
                    ),
                } for item in header.declarations],
                "validation": {
                    "records": len(lines), "public_alleles": sum(len(item.alts) for item in records),
                    "ref_checked": len(lines), "ref_mismatches": 0,
                    "syntax_valid": True, "complete_header_preserved": True,
                    "dictionary_order_sorted": True, "index_queries_verified": len(lines),
                },
            }

        contents = {"schema.header.vcf": header.content, "compatibility.vcf.gz": vcf,
                    "compatibility.vcf.gz.tbi": tbi}
        report["outputs"] = [{"name": name, "sha256": hashlib.sha256(content).hexdigest(),
                              "size_bytes": len(content), **POLICY} for name, content in contents.items()]
        contents["report.json"] = _json_bytes(report)
        run_id = "compatibility-" + hashlib.sha256(contents["report.json"]).hexdigest()[:24]
        artifacts = []
        for name, content in contents.items():
            path = (target / name).relative_to(root).as_posix()
            digest = hashlib.sha256(content).hexdigest()
            artifacts.append(ArtifactRecord(
                artifact_id=hashlib.sha256(f"controlled\0{path}\0{digest}".encode()).hexdigest()[:24],
                relative_path=path, classification=ArtifactClassification.CONTROLLED,
                kind="compatibility-" + name, sha256=digest,
                size_bytes=len(content), retention="delete",
            ))
        inventory = CaseInventory(run_id, tuple(artifacts))
        inventory_content = _json_bytes({**asdict(inventory), **POLICY})

        # Reserve a new private directory; never reuse, replace, or chmod user files.
        parent_fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
        for index, component in enumerate(Path(relative).parts[:-1]):
            try:
                os.mkdir(component, mode=0o700, dir_fd=parent_fd)
            except FileExistsError:
                pass
            next_fd = os.open(component, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW, dir_fd=parent_fd)
            os.close(parent_fd)
            parent_fd = next_fd
            metadata = os.fstat(parent_fd)
            if not (index == 0 and component == "scratch") and (
                metadata.st_mode & 0o077 or metadata.st_uid != os.geteuid()
            ):
                raise CompatibilityFixtureError("output-directory-not-private")
        os.mkdir(target.name, mode=0o700, dir_fd=parent_fd)
        target_created = True
        directory_fd = os.open(target.name, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW, dir_fd=parent_fd)
        opened_target = Path(f"/proc/self/fd/{directory_fd}")
        if opened_target.resolve(strict=True) != target:
            raise CompatibilityFixtureError("output-directory-changed")
        # Persist the deletion ledger first, so an interrupted write is inventoried.
        for name, content in {"inventory.json": inventory_content, **contents}.items():
            if name != "report.json":
                _write_private_file(directory_fd, name, content, created)
        index = TabixIndex.load(opened_target / "compatibility.vcf.gz.tbi")
        index.require_vcf_layout()
        if tuple(item.name for item in index.references) != tuple(dict.fromkeys(item["chrom"] for item in provenance)):
            raise CompatibilityFixtureError("generated-index-invalid")
        with gzip.open(opened_target / "compatibility.vcf.gz", "rb") as source:
            if source.read() != header.content + b"".join(lines):
                raise CompatibilityFixtureError("generated-vcf-roundtrip-invalid")
        for line, record in zip(lines, provenance, strict=True):
            if line.decode("ascii").rstrip("\n") not in index.query(
                opened_target / "compatibility.vcf.gz", record["chrom"], record["position"], record["position"],
            ):
                raise CompatibilityFixtureError("generated-index-query-failed")
        for artifact in inventory.artifacts:
            name = Path(artifact.relative_path).name
            if name == "report.json":
                continue
            path = opened_target / name
            if _checksums(path)[::2] != (artifact.sha256, artifact.size_bytes):
                raise CompatibilityFixtureError("generated-artifact-changed")
            if stat.S_IMODE(path.stat().st_mode) != 0o600:
                raise CompatibilityFixtureError("generated-artifact-not-private")
        if opened_target.resolve(strict=True) != target:
            raise CompatibilityFixtureError("output-directory-changed")
        full_inventory = load_compatibility_inventory(root, opened_target / "inventory.json")
        # The success report is the last product; a partial run cannot claim a pass.
        _write_private_file(directory_fd, "report.json", contents["report.json"], created)
        os.fsync(directory_fd)
        os.fsync(parent_fd)
        return CompatibilityFixtureResult(
            str(target / "compatibility.vcf.gz"), str(target / "compatibility.vcf.gz.tbi"),
            str(target / "report.json"), str(target / "inventory.json"), full_inventory,
        )
    except (CompatibilityFixtureError, VcfHeaderError, BoundedJSONError,
            ChallengeFixtureError, GiabError, ReferenceConcordanceError,
            ManifestError, GovernanceError, TabixError, OSError, ValueError,
            TypeError, KeyError, IndexError, AttributeError, EOFError, RuntimeError) as error:
        try:
            for name in reversed(created):
                try:
                    os.unlink(name, dir_fd=directory_fd)
                except FileNotFoundError:
                    pass
            if target_created:
                os.rmdir(target.name, dir_fd=parent_fd)
        except OSError:
            raise CompatibilityFixtureError("compatibility-cleanup-incomplete") from None
        if isinstance(error, CompatibilityFixtureError):
            raise
        raise CompatibilityFixtureError("compatibility-generation-invalid") from None
    finally:
        if directory_fd >= 0:
            os.close(directory_fd)
        if parent_fd >= 0:
            os.close(parent_fd)
