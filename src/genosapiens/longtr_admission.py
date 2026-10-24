"""Pure interpretation of supplied LongTR records, not a caller or authorization.

The inventory authenticates original input/QC/reference/profile bindings before
calling this adapter. Original calls remain source evidence; parent and member
views reuse AlleleEvidence/IdentityMapping. No source field is rewritten.
"""

from collections.abc import Mapping
from dataclasses import dataclass, replace
import hashlib
import re
from urllib.parse import quote

from genosapiens.allele_identity import GenomicInterval, RepeatAllele, RepeatSize, RepeatSizeUnit, SequenceAllele
from genosapiens.evidence_model import AlleleEvidence, CallEvidence, IdentityMapping, IdentityStage
from genosapiens.evidence_serialization import _canonical
from genosapiens.variant_normalization import VARIANT_NORMALIZATION_SCHEMA
from genosapiens.vcf_qc import technical_flags_for, quality_acceptance_for


COMMIT = "fd73e0d45f239601a12854275ab4efa588010caf"
FORMAT = "longtr-fd73e0d-vcf4.1/v1"
INFO = {"START": ("1", "Integer"), "END": ("1", "Integer"), "MOTIF": (".", "String"),
        "PERIOD": (".", "Integer"), "INEXACT_ALLELE": ("A", "Integer"), "BPDIFFS": ("A", "Integer")}
FORMATS = {"GT": ("1", "String"), "GB": ("1", "String"), "Q": ("1", "Float"), "DP": ("1", "Integer")}
POLICY = {
    "format": FORMAT, "commit": COMMIT, "info": INFO, "formats": FORMATS,
    "source_sha256": {
        "src/genotyper.cpp": "0cd6889e53912f681c25a030a9921c69f10432f817c79579e371f13e4253307f",
        "src/seq_stutter_genotyper.cpp": "d651b4fbe27b1d52df9bfea0b0c61b0c96c4ce90e2cfa3469349e09e9c59ef71",
        "src/region.cpp": "d420f5c8facf1c70fa391a9427057db39d140c7319f671ff9d98bad9fa85c284",
        "src/fasta_reader.cpp": "ce78f0a0f167470e0402b5ac17f6a2fef0ff0946e0c43795cb7f7a0f4d6f37f0",
    },
    "locus": "source-ID-or-explicitly-derived-coordinate-key; no-gene-name-inference",
    "geometry": "START-END-closed-repeat-region; POS-REF-span-is-distinct",
    "sizes": "GB-source-order; base-pairs; only-region-contained-no-flank-change; no-confidence-imputation",
    "identity": "concrete-source-ALT; exact-source-flag; consistent-BPDIFFS-and-known-GB; normalized-view-kept",
    "parent": "one-per-record-sample; shared-lineage; never-an-extra-individual-allele",
    "phase": "source-local-GT-order-only; no-cross-record-block-from-pipe-or-posterior",
    "symbolic": "retain-locus; DEL-index-reordering-not-resolved",
    "limitations": "complex-motif-or-invalid-geometry-retains-unresolved-source-record; no-proxy",
}
POLICY_SHA256 = hashlib.sha256(_canonical(POLICY).encode("ascii")).hexdigest()


def _hash(value):
    return hashlib.sha256(_canonical(value).encode("ascii")).hexdigest()


def validate_header(headers, definitions):
    """Validate this explicitly selected format, not merely a generic VCF flag."""
    if (headers[0] != "##fileformat=VCFv4.1"
            or any(definitions.get(("INFO", k)) != v for k, v in INFO.items())
            or any(definitions.get(("FORMAT", k)) != v for k, v in FORMATS.items())
            or sum(h.startswith("##reference=") for h in headers) != 1
            or not any(h.startswith("##contig=") for h in headers)):
        raise ValueError("LongTR header differs from selected source format")


def _integers(value, separator, *, signed=False):
    if value in (None, "."):
        return None
    pattern = r"-?[0-9]+" if signed else r"[0-9]+"
    values = value.split(separator)
    if any(v != "." and re.fullmatch(pattern, v) is None for v in values):
        raise ValueError
    return tuple(None if v == "." else int(v) for v in values)


@dataclass(frozen=True, repr=False)
class LongtrInterpretation:
    parent: RepeatAllele
    members: tuple[RepeatAllele | None, ...]
    reasons: tuple[tuple[str, ...], ...]
    sample_flags: tuple[str, ...]
    locus_id_source: str
    correspondence: bool


def interpret(source, sample, *, role, reference_id):
    """Interpret source bytes already replayed by inventory, with bounded scope.

    Unknown GT does not remove a concrete exact allele's identity. Conflicting
    GB/GT/sequence correspondence cannot support that sample's member properties.
    GB and all read/posterior fields remain verbatim in the original source.
    """
    fields = source.source_fields
    columns, info = fields["columns"], fields["info"]
    start, end = int(info["START"]), int(info["END"])
    pos, ref, alts = int(columns[1]), columns[3], tuple(columns[4].split(","))
    motif = info["MOTIF"]
    if (re.fullmatch(r"[ACGT]+", motif) is None or _integers(info["PERIOD"], ",") != (len(motif),)
            or not pos <= start <= end <= pos + len(ref) - 1):
        raise ValueError("LongTR locus geometry or motif interpretation unavailable")
    interval = GenomicInterval(columns[0], start, end)
    named = columns[2] != "." and ";" not in columns[2]
    locus = columns[2] if named else "longtr-region:" + quote(columns[0], safe="") + f":{start}-{end}"
    values = dict(zip(fields["format_keys"], fields["sample_values"][sample.sample_name].split(":")))
    flags = () if values.get("FILTER") in (None, ".", "PASS") else ("longtr-sample-filtered",)
    sequence = (ref, *alts)
    gt = tuple(None if t == "." else int(t) for t in re.split(r"[/|]", sample.genotype or "."))
    if values.get("GT", ".") != (sample.genotype or "."):
        raise ValueError("LongTR sample binding differs")
    # Exclude only unsupported sequence claims. Keep every original measurement.
    reasons = [[] for _ in alts]
    inexact, diffs = None, None
    try:
        inexact = _integers(info.get("INEXACT_ALLELE"), ",")
        diffs = _integers(info.get("BPDIFFS"), ",", signed=True)
        if alts != (".",) and (inexact is None or diffs is None or len(inexact) != len(alts)
                               or len(diffs) != len(alts) or any(v not in (0, 1, None) for v in inexact)):
            raise ValueError
    except (ValueError, TypeError):
        inexact = diffs = None
    gb = None
    correspondence = sample.ploidy in (None, 1, 2)
    try:
        gb = _integers(values.get("GB"), "|", signed=True)
        if gb is not None and sample.genotype not in (None, "."):
            if len(gb) != len(gt):
                raise ValueError
            for index, difference in zip(gt, gb, strict=True):
                if index is not None and difference is not None and re.fullmatch(r"[ACGT]+", sequence[index]):
                    if difference != len(sequence[index]) - len(ref):
                        raise ValueError
    except (ValueError, TypeError, IndexError):
        correspondence = False
    left, right = start - pos, pos + len(ref) - 1 - end

    def repeat_length(text):
        # LongTR can retain flank variants. They do not define a repeat-only size.
        if left > 1 or right != 0 or re.fullmatch(r"[ACGT]+", text) is None or len(text) < left + right:
            return None
        if left and text[:left] != ref[:left] or right and text[-right:] != ref[-right:]:
            return None
        return len(text) - left - right

    sizes = None
    if gb is not None:
        sizes = tuple(RepeatSize(
            end - start + 1 + v if v is not None and end - start + 1 + v >= 0
            and left <= 1 and right == 0 and all(repeat_length(s) is not None for s in sequence if s != ".") else None,
            RepeatSizeUnit.BASE_PAIRS) for v in gb)
    parent = RepeatAllele(role, reference_id, locus, interval, motif, sizes, sample.genotype, sample.ploidy)
    members = []
    for n, alt in enumerate(alts):
        if alt == ".":
            if not correspondence:
                reasons[n].append("longtr-genotype-size-correspondence-unresolved")
            members.append(None)
            continue
        if "<DEL>" in alts:
            reasons[n].append("longtr-symbolic-deletion-correspondence-unresolved")
        if re.fullmatch(r"[ACGT]+", alt) is None or alt == ref:
            reasons[n].append("longtr-concrete-member-identity-unavailable")
        if inexact is None or inexact[n] != 0:
            reasons[n].append("longtr-exact-allele-flag-unavailable-or-approximate")
        if diffs is None or diffs[n] != len(alt) - len(ref):
            reasons[n].append("longtr-allele-length-correspondence-unresolved")
        if not correspondence:
            reasons[n].append("longtr-genotype-size-correspondence-unresolved")
        if reasons[n]:
            members.append(None)
        else:
            members.append(replace(parent, allele_sizes=(RepeatSize(repeat_length(alt), RepeatSizeUnit.BASE_PAIRS),),
                exact_sequence=SequenceAllele(role, reference_id, columns[0], pos, ref, alt)))
    return LongtrInterpretation(parent, tuple(members), tuple(tuple(r) for r in reasons), flags,
                               "VCF/ID" if named else "derived-source-region-key", correspondence)


def metadata(allele):
    native = next(i for i in allele.identities if i.identity_id == allele.native_identity_id)
    value = (native.source.source_fields or {}).get("native_interpretation") or {}
    return value if isinstance(value, Mapping) and value.get("repeat_identity") == FORMAT else {}


def is_parent(allele):
    return metadata(allele).get("repeat_representation") == "locus_call"


def is_projection(allele):
    """Check the source/member relation; upstream original-input replay is required."""
    data = metadata(allele)
    if not data or data.get("repeat_representation") != "allele_projection":
        return False
    native = next(i for i in allele.identities if i.identity_id == allele.native_identity_id)
    if data.get("protocol_sha256") != POLICY_SHA256 or len(allele.calls or ()) != 1:
        raise ValueError("invalid LongTR projection binding")
    call, = allele.calls
    if (call.sample is None or call.sample.sample_name != data.get("sample_name") or call.source != native.source
            or call.source.source_release != "VCFv4.1"):
        raise ValueError("invalid LongTR projection source")
    profile = data.get("source_format") or {}
    expected = {"repeat_identity": FORMAT, "input_sha256": native.source.source_sha256,
        "reference_role": native.identity.reference_role.value, "reference_identity_id": native.identity.reference_identity_id,
        "reference_manifest_sha256": native.reference_manifest_sha256}
    if (any(profile.get(k) != v for k, v in expected.items())
            or data.get("parent_call_id") != "longtr-locus-" + _hash((
                native.source.source_sha256, native.source.source_record_id, call.sample.sample_name))):
        raise ValueError("invalid LongTR projection input or reference")
    parent = next(i for i in allele.identities if i.identity_id == data["parent_identity_id"])
    facts = interpret(call.source, call.sample, role=native.identity.reference_role,
                      reference_id=native.identity.reference_identity_id)
    n = call.source.source_fields["allele_number"]
    if (not 1 <= n <= len(facts.members) or facts.members[n - 1] != native.identity
            or parent.source.source_sha256 != native.source.source_sha256
            or parent.source.source_record_id != native.source.source_record_id
            or parent.identity != facts.parent or parent.reference_manifest_sha256 != native.reference_manifest_sha256
            or data.get("parent_identity_sha256") != _hash(parent.identity)
            or parent.identity_id not in native.derived_from_identity_ids):
        raise ValueError("invalid LongTR parent or member correspondence")
    return True


def build_views(row, sources, *, profile, profile_sha256, role, reference_id, reference_sha256, normalized, mappings):
    """Return complete parent/member views and per-ALT accounting for one record."""
    parents, children, by_source, reasons = [], [], {s.source_allele_id: [] for s in sources}, {s.source_allele_id: [] for s in sources}
    for sample in row.qc.sample_calls:
        facts = interpret(row.source, sample, role=role, reference_id=reference_id)
        parent_id = "longtr-locus-" + _hash((row.source.source_sha256, row.source.source_record_id, sample.sample_name))
        parent_native_id = "native-" + parent_id
        child_ids = tuple("longtr-member-" + _hash((parent_id, n + 1)) if m is not None else None
                          for n, m in enumerate(facts.members))
        common = {"adapter": FORMAT, "repeat_identity": FORMAT, "protocol_sha256": POLICY_SHA256,
            "source_format": profile, "source_format_sha256": profile_sha256, "sample_name": sample.sample_name,
            "parent_call_id": parent_id, "parent_identity_id": parent_native_id,
            "parent_identity_sha256": _hash(facts.parent), "locus_id_source": facts.locus_id_source,
            "interval_semantics": "source-START-END-closed-repeat-region", "phase_scope": "source-record-local",
            "projected_member_ids": tuple(i for i in child_ids if i is not None),
            "source_alt_scope_complete": facts.correspondence and sample.genotype is not None and "." not in sample.genotype
                and all(m is not None or s.alternate == "." for m, s in zip(facts.members, sources)),
            "size_semantics": "base-pairs; source-GB-order; confidence-unknown"}
        parent_sources, parent_calls = [], []
        for n, s in enumerate(sources):
            technical_flags = tuple(dict.fromkeys((
                *technical_flags_for(row.qc, sample, s.allele_number),
                *facts.sample_flags,
            )))
            required = bool(technical_flags)
            fields = {**row.source.source_fields, "source_allele_id": s.source_allele_id, "allele_number": s.allele_number,
                "technical_flags": technical_flags, "read_validation_required": required,
                "quality_acceptance": quality_acceptance_for(row.qc, sample, s.allele_number)}
            parent_source = replace(row.source, source_fields={**fields,
                "native_interpretation": {**common, "repeat_representation": "locus_call"}})
            parent_sources.append(parent_source)
            parent_calls.append(CallEvidence(parent_source, sample,
                None if row.qc.filter_status.value == "not-applied" else row.qc.filters, row.source.source_fields["quality"]))
            by_source[s.source_allele_id].append(parent_id)
            reasons[s.source_allele_id].extend(facts.reasons[n])
            member = facts.members[n]
            if member is None:
                continue
            source = replace(row.source, source_fields={**fields,
                "native_interpretation": {**common, "repeat_representation": "allele_projection",
                                          "parent_identity_id": s.source_allele_id}})
            identities = [IdentityMapping(s.source_allele_id, IdentityStage.SOURCE, facts.parent,
                                          reference_sha256, parent_source)]
            predecessors = [s.source_allele_id]
            if s.source_allele_id in mappings:
                m = mappings[s.source_allele_id]
                a = normalized[m.normalized_allele_id]
                identities.append(IdentityMapping(a.normalized_allele_id, IdentityStage.NORMALIZED,
                    SequenceAllele(role, reference_id, a.contig, a.position, a.reference, a.alternate), reference_sha256,
                    source, (s.source_allele_id,), VARIANT_NORMALIZATION_SCHEMA))
                predecessors.append(a.normalized_allele_id)
            child_id = child_ids[n]
            identities.append(IdentityMapping("native-" + child_id, IdentityStage.NATIVE, member, reference_sha256,
                                             source, tuple(predecessors), FORMAT))
            children.append(AlleleEvidence(child_id, "native-" + child_id, tuple(identities),
                (CallEvidence(source, sample, parent_calls[-1].filters, parent_calls[-1].quality),)))
            by_source[s.source_allele_id].append(child_id)
        source_views = tuple(IdentityMapping(s.source_allele_id, IdentityStage.SOURCE, facts.parent, reference_sha256, src)
                             for s, src in zip(sources, parent_sources, strict=True))
        parents.append(AlleleEvidence(parent_id, parent_native_id,
            (*source_views, IdentityMapping(parent_native_id, IdentityStage.NATIVE, facts.parent, reference_sha256, parent_sources[0],
                tuple(s.source_allele_id for s in sources), FORMAT)), tuple(parent_calls)))
    return tuple(parents + children), by_source, reasons
