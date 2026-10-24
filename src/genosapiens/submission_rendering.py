"""Name-only views of normalized alleles, not submission eligibility checks."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict, dataclass, field
import hashlib
import json
from pathlib import Path
from types import MappingProxyType

from genosapiens.bounded_json import (
    BoundedJSONError,
    load_bounded_json,
    read_bounded_bytes,
)
from genosapiens.governance import ArtifactClassification
from genosapiens.grch38_reference import CONTIG_ALIAS_SCHEMA
from genosapiens.llm_guard import LLMPolicy
from genosapiens.manifest import ManifestError, PublicResourceManifest
from genosapiens.reference_concordance import (
    ReferenceConcordanceError,
    ReferenceIdentityContract,
    ReferenceRoleBinding,
    _parse_fai_bytes,
    verify_reference_manifest_contract,
)
from genosapiens.variant_normalization import NormalizedAllele


SUBMISSION_CONTIG_VIEW_SCHEMA = "genosapiens.submission-contig-view"
CANONICAL_CONTIGS = (*(str(number) for number in range(1, 23)), "X", "Y", "M")


class SubmissionRenderingError(RuntimeError):
    """Raised without variant, reference-path, or underlying exception content."""


@dataclass(frozen=True, slots=True, repr=False)
class SubmissionContigView:
    """Retain the normalized identity and its source-link key without copying links.

    POS/REF/ALT and the source reference role/identity remain in normalized_allele.
    Only contig is a rendered name; it does not assign a new reference identity.
    """

    normalized_allele: NormalizedAllele
    contig: str
    reference_manifest_sha256: str
    contig_aliases_sha256: str
    schema: str = field(default=SUBMISSION_CONTIG_VIEW_SCHEMA, init=False)
    classification: ArtifactClassification = field(
        default=ArtifactClassification.CONTROLLED, init=False
    )
    llm_policy: LLMPolicy = field(default=LLMPolicy.NO_LLM, init=False)
    name_only: bool = field(default=True, init=False)
    submission_reference_validated: bool = field(default=False, init=False)
    submission_ready: bool = field(default=False, init=False)


@dataclass(frozen=True, slots=True)
class SubmissionContigRenderer:
    """A loaded name map for one explicitly bound source reference role.

    Create with load_submission_contig_renderer, once per installed bundle/role.
    This snapshot needs no file access when rendering an allele.
    """

    reference_binding: ReferenceRoleBinding
    reference_manifest_sha256: str
    contig_aliases_sha256: str
    _names: Mapping[str, tuple[str, int]] = field(repr=False)
    _canonical_aliases: Mapping[str, tuple[str, int]] = field(repr=False)
    _reference_aliases: Mapping[str, str] = field(
        default_factory=dict, repr=False,
    )
    _source_names: Mapping[str, tuple[str, int]] = field(
        default_factory=dict, repr=False,
    )

    @property
    def source_names(self) -> Mapping[str, tuple[str, int]]:
        """Authenticated source coordinates, including non-submission contigs.

        The loader builds this map once from the bound role's FAI. Keep access
        constant-time: source routing uses it for every occurrence. The fallback
        preserves explicitly constructed renderers' existing name contracts.
        """
        return MappingProxyType(self._source_names) if self._source_names else self.accepted_names

    @property
    def accepted_names(self) -> Mapping[str, tuple[str, int]]:
        """Return authenticated native names plus accepted canonical aliases."""

        return MappingProxyType({**self._names, **self._canonical_aliases})

    @property
    def reference_aliases(self) -> Mapping[str, str]:
        """Return authenticated input names mapped to this role's FAI names.

        Every FAI-native name maps to itself.  Additional aliases come only
        from the validated contig-alias artifact; callers must not infer them
        from a spelling convention.
        """

        return MappingProxyType(dict(self._reference_aliases))

    def render(
        self,
        allele: NormalizedAllele,
        *,
        submission_contig_prefix: str = "chr",
    ) -> SubmissionContigView:
        if submission_contig_prefix not in ("", "chr"):
            raise SubmissionRenderingError("submission contig prefix is invalid")
        if (
            not isinstance(allele, NormalizedAllele)
            or allele.reference_role != self.reference_binding.role
            or allele.reference_identity_id != self.reference_binding.identity_id
        ):
            raise SubmissionRenderingError("normalized allele reference binding is invalid")
        entry = self.accepted_names.get(allele.contig)
        if entry is None:
            raise SubmissionRenderingError("contig has no canonical rendering")
        canonical, length = entry
        if allele.position + len(allele.reference) - 1 > length:
            raise SubmissionRenderingError("normalized allele span is out of bounds")
        return SubmissionContigView(
            normalized_allele=allele,
            contig=submission_contig_prefix + canonical,
            reference_manifest_sha256=self.reference_manifest_sha256,
            contig_aliases_sha256=self.contig_aliases_sha256,
        )


def load_submission_contig_renderer(
    *,
    identity_path: str | Path,
    reference_manifest_path: str | Path,
    reference_root: str | Path,
    reference_binding: ReferenceRoleBinding,
) -> SubmissionContigRenderer:
    """Validate installed v2 naming metadata, not genome sequence or case acceptance.

    The contract and manifest must be the frozen inputs selected by the caller.
    Checksums cover the exact alias and source-role FAI bytes parsed here; full
    genome validation belongs to bundle installation and normalization, not this
    naming boundary. The public-role FAI must match its own declared name inventory;
    only shared canonical aliases require lengths matching the challenge roles.
    """

    try:
        identity_data, identity_content = load_bounded_json(identity_path)
        if not isinstance(identity_data, Mapping):
            raise SubmissionRenderingError("reference identity is invalid")
        contract = ReferenceIdentityContract.from_mapping(identity_data)
        if not isinstance(reference_binding, ReferenceRoleBinding):
            raise SubmissionRenderingError("source reference binding is required")
        role = reference_binding.resolve(contract, required_role=reference_binding.role)

        root = Path(reference_root).resolve(strict=True)
        manifest_path = Path(reference_manifest_path)
        contract_sha256 = hashlib.sha256(identity_content).hexdigest()
        manifest_directory = root / "manifests/grch38"
        allowed_manifests = {
            manifest_directory / f"{contract.resource_id}.json",
            manifest_directory
            / f"{contract.resource_id}.{contract_sha256[:12]}.json",
        }
        if (
            manifest_path.is_symlink()
            or manifest_path.resolve(strict=True)
            not in {path.resolve() for path in allowed_manifests}
            or not manifest_path.resolve().is_relative_to(root)
        ):
            raise SubmissionRenderingError("reference manifest location is invalid")
        manifest_data, manifest_content = load_bounded_json(manifest_path)
        if not isinstance(manifest_data, Mapping):
            raise SubmissionRenderingError("reference manifest is invalid")
        manifest = PublicResourceManifest.from_mapping(manifest_data)
        verify_reference_manifest_contract(
            contract,
            manifest,
            contract_sha256=contract_sha256,
        )

        outputs = {output.path: output for output in manifest.outputs}
        contents = {}
        aliases = None
        for output_role in ("contig_aliases", role.fai_output_role):
            relative = contract.expected_output(output_role).path
            output = outputs[relative]
            path = root / relative
            if path.is_symlink() or not path.resolve(strict=True).is_relative_to(root):
                raise SubmissionRenderingError("reference naming artifact location is invalid")
            if output.retention.policy != "durable" or output.retention.delete_after is not None:
                raise SubmissionRenderingError("reference naming artifact is not durable")
            if output_role == "contig_aliases":
                aliases, content = load_bounded_json(path)
            else:
                content = read_bounded_bytes(path)
            if (
                len(content) != output.size_bytes
                or hashlib.sha256(content).hexdigest() != output.sha256
            ):
                raise SubmissionRenderingError("reference naming artifact checksum mismatch")
            contents[output_role] = content

        if not isinstance(aliases, Mapping):
            raise SubmissionRenderingError("contig alias metadata is invalid")
        expected_roles = {
            item.role: {
                "identity_id": item.identity_id,
                "naming": item.naming,
                "sequence_relationship": item.sequence_relationship,
            }
            for item in contract.reference_roles
        }
        alignment_to_tool = aliases.get("alignment_to_tool")
        tool_to_alignment = aliases.get("tool_to_alignment")
        lengths = aliases.get("sequence_lengths")
        if (
            aliases.get("schema") != CONTIG_ALIAS_SCHEMA
            or aliases.get("contract_id") != contract.contract_id
            or type(aliases.get("sequence_count")) is not int
            or aliases["sequence_count"] != contract.sequence_count
            or aliases.get("reference_roles") != expected_roles
            or aliases.get("canonical_submission_contigs") != list(CANONICAL_CONTIGS)
            or not isinstance(alignment_to_tool, dict)
            or not isinstance(tool_to_alignment, dict)
            or not isinstance(lengths, dict)
        ):
            raise SubmissionRenderingError("contig alias metadata is invalid")
        if (
            len(lengths) != contract.sequence_count
            or any(
                not name
                or any(character.isspace() for character in name)
                or type(length) is not int
                or length <= 0
                for name, length in lengths.items()
            )
            or set(alignment_to_tool) != set(lengths)
            or len(tool_to_alignment) != len(lengths)
            or any(
                not tool
                or any(character.isspace() for character in tool)
                or tool.removeprefix(contract.transformation.source_prefix) != alignment
                for tool, alignment in tool_to_alignment.items()
            )
            or alignment_to_tool != {
                alignment: tool for tool, alignment in tool_to_alignment.items()
            }
            or any(
                lengths.get(name) != length
                for name, length in contract.canonical_contig_lengths.items()
            )
        ):
            raise SubmissionRenderingError("contig aliases or lengths are inconsistent")

        index = _parse_fai_bytes(contents[role.fai_output_role])
        canonical_aliases: dict[str, tuple[str, int]] = {}
        if role.role == "public_annotation_benchmark":
            public_names = aliases.get("public_native_sequence_md5")
            # Names establish coordinate aliases only. Never compare public and
            # challenge M5 values here or treat their equality as allele validation.
            if not isinstance(public_names, dict) or set(public_names) != set(index):
                raise SubmissionRenderingError("public reference has no declared name mapping")
            native_names = {
                name: canonical for name, canonical in tool_to_alignment.items()
                if name in index and canonical in CANONICAL_CONTIGS
            }
            # Public fixtures commonly use the canonical unprefixed aliases from
            # the same authenticated map. Their REF alleles remain subject to the
            # separate public-validation concordance check; this is name-only.
            source_names = native_names
            source_lengths = {name: index[name].length for name in native_names}
            canonical_aliases = {
                canonical: (canonical, index[native].length)
                for native, canonical in native_names.items()
            }
        else:
            source_names = (
                {name: name for name in lengths}
                if role.role == "challenge_alignment"
                else tool_to_alignment
            )
            if set(index) != set(source_names):
                raise SubmissionRenderingError("source reference dictionary has no safe name mapping")
            source_lengths = {name: index[name].length for name in source_names}
        if any(
            source_lengths[name] != lengths[canonical]
            for name, canonical in source_names.items()
        ):
            raise SubmissionRenderingError("source reference dictionary has no safe name mapping")

        reference_aliases = {name: name for name in index}
        for native, canonical in tool_to_alignment.items():
            if native not in index:
                continue
            existing = reference_aliases.get(canonical)
            if existing is not None and existing != native:
                raise SubmissionRenderingError(
                    "source reference dictionary has ambiguous aliases"
                )
            reference_aliases[canonical] = native

        return SubmissionContigRenderer(
            reference_binding=reference_binding,
            reference_manifest_sha256=hashlib.sha256(manifest_content).hexdigest(),
            contig_aliases_sha256=hashlib.sha256(contents["contig_aliases"]).hexdigest(),
            _names=MappingProxyType({
                name: (canonical, lengths[canonical])
                for name, canonical in source_names.items()
                if canonical in CANONICAL_CONTIGS
            }),
            _canonical_aliases=MappingProxyType(canonical_aliases),
            _reference_aliases=MappingProxyType(reference_aliases),
            _source_names=MappingProxyType({
                **{name: (source_names.get(name, name), entry.length)
                   for name, entry in index.items()},
                **canonical_aliases,
            }),
        )
    except (
        BoundedJSONError,
        ManifestError,
        ReferenceConcordanceError,
        OSError,
        RuntimeError,
        TypeError,
        ValueError,
        KeyError,
        AttributeError,
    ):
        # Do not chain exceptions whose messages can contain paths or input data.
        raise SubmissionRenderingError("installed contig rendering binding is invalid") from None


def submission_contig_view_as_json(view: SubmissionContigView) -> str:
    """Serialize a controlled view; source mappings remain in the normalization result."""
    return json.dumps(asdict(view), indent=2, sort_keys=True, ensure_ascii=True) + "\n"
