"""Authenticated offline gnomAD lookup and attachment, never frequency policy.

The normalized input representation is authoritative. Public-native REF agreement
is an additional check, not permission to repair, renormalize, or proxy an allele.
Reports and enriched evidence are controlled/no-LLM even for public test inputs.
"""

from __future__ import annotations

from contextlib import AbstractContextManager
from dataclasses import dataclass, replace
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import selectors
import shutil
import socket
import subprocess
import time
from typing import ClassVar

from genosapiens.allele_identity import MosaicAllele, ReferenceRole, SequenceAllele
from genosapiens.evidence_model import AlleleEvidence, IdentityStage, PopulationEvidence, SourcedRecord, _typed, _validation_session, _construct_owned
from genosapiens.evidence_serialization import _plain, to_json
from genosapiens.governance import ArtifactClassification
from genosapiens.grch38_reference import IndexedFasta
from genosapiens.llm_guard import LLMPolicy
from genosapiens.manifest import PublicResourceManifest, SourceIdentity, Transformation
from genosapiens import gnomad_all_sites as sites
from genosapiens import gnomad_publication as publication
from genosapiens.gnomad_projection import (
    DATASET_VCF_PREFIXES, GnomadProjectionContract, gnomad_vcf_metadata,
    load_gnomad_projection_contract,
)
from genosapiens.gnomad_validation import _parse_info, _validate_header, load_gnomad_validation_contract
from genosapiens.reference_concordance import (
    ReferenceRoleBinding, load_reference_identity, verify_reference_manifest_contract,
)
from genosapiens.submission_rendering import load_submission_contig_renderer
from genosapiens.variant_normalization import (
    SourceAllele, SourceSmallVariantRecord, _normalize_sequence_allele,
    source_allele_id, source_record_id,
)
from genosapiens.vep_annotation import _local, _preflight_resources
from genosapiens.vep_bundle import _file_identity


SCHEMA = "genosapiens.population-annotation"
JOINT_SOURCE_DATASET = "gnomAD joint sites"
ALL_SITES_SOURCE_DATASET = "gnomAD all-sites AN"
SOURCE_RELEASE = "v4.1"
GNOMAD_SUBSETS = tuple(DATASET_VCF_PREFIXES)
GNOMAD_VARIANT_STRATA = ("adjusted", "ancestry", "grpmax", "fafmax95", "fafmax99")
ALL_SITES_STRATUM = "all_sites"
_SHA = re.compile(r"[0-9a-f]{64}")
_ANCESTRIES = ("afr", "amr", "eas", "mid", "nfe", "sas")


from genosapiens.wgs_failures import WgsFailure


class PopulationAnnotationError(WgsFailure):
    """Static diagnostics only; no paths, alleles, or subprocess output."""

    def __init__(self, message, *, category="malformed_evidence"):
        super().__init__(
            "population_lookup", category,
            operation="population_query", exception_kind="project_failure",
        )
        RuntimeError.__init__(self, message)


class _Controlled:
    __slots__ = ()
    classification: ClassVar[ArtifactClassification] = ArtifactClassification.CONTROLLED
    llm_policy: ClassVar[LLMPolicy] = LLMPolicy.NO_LLM

    def __repr__(self):
        return "<controlled population annotation; no-llm>"


def _json(value):
    return json.dumps(_plain(value), sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n"


def _digest(value):
    return hashlib.sha256(value.encode("utf-8") if isinstance(value, str) else value).hexdigest()


@dataclass(frozen=True, slots=True, repr=False)
class PopulationRequest(_Controlled):
    allele: AlleleEvidence
    normalized_identity_id: str | None


@dataclass(frozen=True, slots=True)
class PopulationQueryLimits:
    """Explicit technical bounds, not biological frequency defaults."""

    max_alleles: int
    max_loci_per_query: int
    max_region_argument_bytes: int
    max_rows_per_query: int
    max_output_bytes: int
    timeout_seconds: int

    def __post_init__(self):
        if any(type(v) is not int or v <= 0 for v in (
            self.max_alleles, self.max_loci_per_query, self.max_region_argument_bytes,
            self.max_rows_per_query, self.max_output_bytes, self.timeout_seconds,
        )):
            raise PopulationAnnotationError("invalid population query limits")


@dataclass(frozen=True, slots=True, repr=False)
class GnomadJointReference(_Controlled):
    projection_contract_path: Path
    publication_contract_path: Path
    validation_contract_path: Path
    manifest_path: Path
    manifest_sha256: str
    release_directory: Path


@dataclass(frozen=True, slots=True, repr=False)
class GnomadAllSitesReference(_Controlled):
    contract_path: Path
    reference_root: Path
    bundle_manifest_sha256: str


@dataclass(frozen=True, slots=True, repr=False)
class PopulationReferenceCompatibility(_Controlled):
    identity_path: Path
    identity_sha256: str
    reference_root: Path
    manifest_path: Path
    manifest_sha256: str


@dataclass(frozen=True, slots=True, repr=False)
class PopulationAnnotationResult(_Controlled):
    alleles: tuple[AlleleEvidence, ...]
    report: SourcedRecord


@dataclass(frozen=True, slots=True)
class _VerifiedPopulationFile:
    sha256: str
    size_bytes: int
    stat_identity: tuple[int, int, int, int, int]


def _stat_identity(value: os.stat_result) -> tuple[int, int, int, int, int]:
    return (value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns)


class _PopulationFileCache:
    """Trust qualified published outputs; authenticate other bytes once per stage.

    Published-output reuse checks the pinned manifest and file metadata, not fresh
    content bytes. Same-size modifications that evade metadata checks can therefore
    remain undetected until an actual query exposes them.
    """

    def __init__(self) -> None:
        self.files: dict[Path, _VerifiedPopulationFile] = {}
        self.objects: dict[tuple[object, ...], object] = {}

    def prepared(self, key: tuple[object, ...], build):
        """Keep immutable parsed metadata within one owner session."""
        if key not in self.objects:
            self.objects[key] = build()
        return self.objects[key]

    def verify(self, path: Path, sha256: str | None = None,
               size_bytes: int | None = None,
               *, published_by: str | None = None) -> Path:
        path = _local(path)
        _preflight_resources(files=(path,))
        try:
            before = _stat_identity(path.stat())
        except OSError:
            raise PopulationAnnotationError("population resource byte binding mismatch") from None
        cached = self.files.get(path)
        if cached is None:
            if published_by is not None:
                if (type(published_by) is not str or _SHA.fullmatch(published_by) is None
                        or type(sha256) is not str or _SHA.fullmatch(sha256) is None
                        or type(size_bytes) is not int or size_bytes <= 0
                        or before[2] != size_bytes):
                    raise PopulationAnnotationError("population resource publication binding mismatch")
                cached = _VerifiedPopulationFile(sha256, size_bytes, before)
            else:
                digest, observed_size = _file_identity(path)
                try:
                    after = _stat_identity(path.stat())
                except OSError:
                    raise PopulationAnnotationError("population resource byte binding mismatch") from None
                if before != after or observed_size != before[2]:
                    raise PopulationAnnotationError("population resource byte binding mismatch")
                cached = _VerifiedPopulationFile(digest, observed_size, before)
            self.files[path] = cached
        elif before != cached.stat_identity:
            raise PopulationAnnotationError("population resource byte binding mismatch")
        if (sha256 is not None and (
                type(sha256) is not str or _SHA.fullmatch(sha256) is None
                or sha256 != cached.sha256)
                or size_bytes is not None and size_bytes != cached.size_bytes):
            raise PopulationAnnotationError("population resource byte binding mismatch")
        return path

    def assert_stable(self) -> None:
        for path, cached in tuple(self.files.items()):
            try:
                _preflight_resources(files=(path,))
                current = _stat_identity(path.stat())
            except (OSError, ValueError, RuntimeError):
                raise PopulationAnnotationError("population resource byte binding mismatch") from None
            if current != cached.stat_identity:
                raise PopulationAnnotationError("population resource byte binding mismatch")

    def report_files(self, selected: dict[Path, tuple[str, int]]) -> None:
        """Ensure report bindings agree with the stage cache without rereading bytes."""
        for path, (sha256, size_bytes) in selected.items():
            cached = self.files.get(path)
            if cached is None or (cached.sha256, cached.size_bytes) != (sha256, size_bytes):
                raise PopulationAnnotationError("population resource byte binding mismatch")


class PopulationResourceSession(AbstractContextManager["PopulationResourceSession"]):
    """One offline population-resource lifetime shared by bounded input batches."""

    def __init__(
        self, *, joint: GnomadJointReference, all_sites: GnomadAllSitesReference,
        reference: PopulationReferenceCompatibility, runtime_image: str,
        limits: PopulationQueryLimits,
    ) -> None:
        self.joint = joint
        self.all_sites = all_sites
        self.reference = reference
        self.runtime_image = runtime_image
        self.limits = limits
        self._files = _PopulationFileCache()
        self._closed = False
        # Physical work counters are telemetry, never scientific report identity.
        self.command_invocations = {'bcftools': 0, 'tabix': 0}
        self.work_counts = {'input_encodes': 0, 'decorated_alleles': 0, 'report_encodes': 0,
                            'screening_rows': 0, 'physical_queries': 0, 'proofs': 0}

    def query(self, requests: tuple[PopulationRequest, ...], *, _case_cache=None) -> PopulationAnnotationResult:
        if self._closed:
            raise PopulationAnnotationError("population resource session is closed")
        self._files.assert_stable()
        result = _run_population_annotation(
            requests, joint=self.joint, all_sites=self.all_sites, reference=self.reference,
            runtime_image=self.runtime_image, limits=self.limits, file_cache=self._files,
            case_cache=_case_cache, command_metrics=self.command_invocations,
            work_counts=self.work_counts,
        )
        self._files.assert_stable()
        return result

    def __exit__(self, *_args: object) -> None:
        try:
            self._files.assert_stable()
        finally:
            self._closed = True

    def query_with_noncontribution(self, requests, *, _case_cache=None):
        """Explicit supplemental proof output; ordinary query semantics stay exact."""
        if self._closed:
            raise PopulationAnnotationError('population resource session is closed')
        self._files.assert_stable()
        proofs = []
        result = _run_population_annotation(requests, joint=self.joint, all_sites=self.all_sites,
            reference=self.reference, runtime_image=self.runtime_image, limits=self.limits,
            file_cache=self._files, case_cache=_case_cache, noncontribution_output=proofs,
            command_metrics=self.command_invocations, work_counts=self.work_counts)
        self._files.assert_stable()
        return result, tuple(proofs)

    def _screening_facts(self, requests, *, _case_cache=None):
        """Private adaptation output; no decorated alleles or full report."""
        if self._closed:
            raise PopulationAnnotationError('population resource session is closed')
        self._files.assert_stable()
        rows, proofs = [], []
        _run_population_annotation(requests, joint=self.joint, all_sites=self.all_sites,
            reference=self.reference, runtime_image=self.runtime_image, limits=self.limits,
            file_cache=self._files, case_cache=_case_cache, noncontribution_output=proofs,
            screening_output=rows, command_metrics=self.command_invocations, work_counts=self.work_counts)
        self._files.assert_stable()
        return _PopulationScreeningResult(tuple(rows), tuple(proofs))


@dataclass(frozen=True, slots=True, repr=False)
class _PopulationScreeningRow:
    allele_id: str
    normalized_identity_id: str | None
    input_sha256: str
    population: tuple[PopulationEvidence, ...]
    query_receipts: tuple[bytes, ...]


@dataclass(frozen=True, slots=True, repr=False)
class _PopulationScreeningResult:
    rows: tuple[_PopulationScreeningRow, ...]
    proofs: tuple


@dataclass(frozen=True, slots=True, repr=False)
class _PopulationQueryText:
    text: str
    receipt: bytes


def _raw_population_fact(fact, *, scope, region):
    """Decode a bound raw cache fact and its original physical query receipt."""
    from collections.abc import Mapping
    if (fact.status not in ('positive', 'no_match') or not isinstance(fact.value, Mapping)
            or set(fact.value) != {'lines', 'receipt'} or type(fact.value['lines']) is not tuple
            or any(type(line) is not str for line in fact.value['lines'])
            or bool(fact.value['lines']) != (fact.status == 'positive')
            or type(fact.value['receipt']) is not str):
        raise ValueError
    receipt = fact.value['receipt'].encode()
    data = json.loads(receipt)
    if (_json(data).encode() != receipt or set(data) != {'schema', 'version', 'tool', 'tool_sha256',
            'runtime_image', 'argv', 'regions', 'resource_scope', 'stdout_sha256', 'stdout_bytes', 'exit_status'}
            or data['schema'] != 'genosapiens.population-physical-query' or type(data['version']) is not int
            or data['version'] != 1 or data['resource_scope'] != scope
            or data['runtime_image'] != scope['runtime_image']
            or data['tool'] not in scope['tools'] or data['tool_sha256'] != scope['tools'][data['tool']]
            or type(data['argv']) is not list or any(type(v) is not str for v in data['argv'])
            or type(data['regions']) is not list or not data['regions']
            or type(data['stdout_sha256']) is not str or _SHA.fullmatch(data['stdout_sha256']) is None
            or type(data['stdout_bytes']) is not int or data['stdout_bytes'] < 0
            or type(data['exit_status']) is not int or data['exit_status'] != 0):
        raise ValueError
    expected_args = ([scope['artifact'], *data['regions']] if data['tool'] == 'tabix' else
        ['query', '--regions', ','.join(data['regions']), '--format',
         '%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER\\t%INFO\\n', scope['artifact']])
    if (not data['argv'] or not Path(data['argv'][0]).is_absolute()
            or data['argv'][0] != scope['tool_paths'][data['tool']] or data['argv'][1:] != expected_args):
        raise ValueError
    def span(value):
        match = re.fullmatch(r'([^:]+):([0-9]+)-([0-9]+)', value)
        if match is None: raise ValueError
        contig, start, end = match.groups()
        if not 0 < int(start) <= int(end): raise ValueError
        return contig, int(start), int(end)
    contig, start, end = span(region)
    if not any(c == contig and lo <= start <= end <= hi for c, lo, hi in map(span, data['regions'])):
        raise ValueError
    return _PopulationQueryText(''.join(line+'\n' for line in fact.value['lines']), receipt)


def _command(executable, arguments, limits):
    """Drain both pipes with a wall-clock and byte bound; never return a prefix."""
    import signal
    from genosapiens.wgs_query_limits import _check_query_cancellation
    process = None
    try:
        _check_query_cancellation()
        process = subprocess.Popen((str(executable), *arguments), stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, stdin=subprocess.DEVNULL,start_new_session=True)
        output, errors = bytearray(), bytearray()
        deadline = time.monotonic() + limits.timeout_seconds
        with selectors.DefaultSelector() as selector:
            selector.register(process.stdout, selectors.EVENT_READ, output)
            selector.register(process.stderr, selectors.EVENT_READ, errors)
            while selector.get_map() or process.poll() is None:
                _check_query_cancellation()
                remaining = deadline - time.monotonic()
                if remaining <= 0:
                    raise PopulationAnnotationError("population query timed out", category="timeout")
                for key, _ in selector.select(min(remaining,.1)):
                    chunk = os.read(key.fd, 65536)
                    if not chunk:
                        selector.unregister(key.fileobj)
                    else:
                        key.data.extend(chunk)
                        if len(output) + len(errors) > limits.max_output_bytes:
                            raise PopulationAnnotationError("population query exceeded byte bound", category="response_capacity")
            if process.wait(timeout=max(0.001, deadline - time.monotonic())) != 0 or errors:
                raise PopulationAnnotationError("population query failed or emitted diagnostics")
        return output.decode("ascii", errors="strict")
    except (OSError, UnicodeError, subprocess.TimeoutExpired):
        raise PopulationAnnotationError("population query execution failed") from None
    finally:
        if process is not None:
            import sys
            primary=sys.exc_info()[1]
            try:
                os.killpg(process.pid,signal.SIGKILL)
            except ProcessLookupError:
                pass
            process.wait()
            process.stdout.close()
            process.stderr.close()
            # The admission process may be an active subreaper. Reap only this
            # query's adopted descendants, never children of another owner.
            deadline=time.monotonic()+5
            while True:
                try:
                    child,_=os.waitpid(-process.pid,os.WNOHANG)
                except ChildProcessError:
                    break
                if child==0:
                    if time.monotonic()>=deadline:
                        cleanup=PopulationAnnotationError('population query descendants did not stop')
                        if primary is None:raise cleanup
                        if isinstance(primary,WgsFailure):primary.add_secondary(cleanup)
                        break
                    time.sleep(.01)


def _lines(content, limits):
    if len(content.encode("ascii")) > limits.max_output_bytes:
        raise PopulationAnnotationError("population query exceeded byte bound")
    lines = content.splitlines()
    if len(lines) > limits.max_rows_per_query:
        raise PopulationAnnotationError("population query exceeded row bound")
    if content and (not content.endswith("\n") or any(not line for line in lines)):
        raise PopulationAnnotationError("incomplete population query output")
    return lines


def _locus_batches(alleles, limits):
    """Return deterministic same-contig region chunks within both argument bounds."""
    by_contig = {}
    for allele in alleles:
        if allele is not None:
            by_contig.setdefault(allele.contig, set()).add(allele.position)
    batches = []
    for contig in sorted(by_contig):
        regions = []
        encoded_bytes = 0
        for position in sorted(by_contig[contig]):
            region = f"{contig}:{position}-{position}"
            added_bytes = len(region.encode("ascii")) + (1 if regions else 0)
            if added_bytes > limits.max_region_argument_bytes:
                raise PopulationAnnotationError("population query region exceeded argument bound")
            if (len(regions) == limits.max_loci_per_query
                    or encoded_bytes + added_bytes > limits.max_region_argument_bytes):
                batches.append((contig, tuple(regions)))
                regions = []
                encoded_bytes = 0
                added_bytes = len(region.encode("ascii"))
            regions.append(region)
            encoded_bytes += added_bytes
        if regions:
            batches.append((contig, tuple(regions)))
    return tuple(batches)


def _batch_record(*, tool, source, contig, regions, content):
    return {
        "batch_id": f"{tool}:{source}:{contig}:{regions[0]}:{regions[-1]}",
        "tool": tool,
        "source": source,
        "contig": contig,
        "regions": regions,
        "region_argument_bytes": len(",".join(regions).encode("ascii")),
        "returned_row_count": len(content.splitlines()),
        "output_sha256": _digest(content),
        "returned_records": tuple(content.splitlines()),
    }


def _raw_info(info):
    return {} if info == "." else dict(item.split("=", 1) for item in info.split(";"))


def _population_info(info, projection, metadata=None):
    # VCF permits a single dot for an entirely missing fixed-length vector.
    # Treat it as unavailable, not a one-element ancestry array or six zeros.
    # Check duplicate keys before omitting these tokens from the existing parser.
    metadata = gnomad_vcf_metadata(projection)["info"] if metadata is None else metadata
    tokens, seen = [], set()
    for token in (() if info == "." else info.split(";")):
        key, value = token.split("=", 1)
        if key in seen or key not in metadata or not value:
            raise ValueError
        seen.add(key)
        if value == "." and metadata[key]["Number"] not in ("1", "."):
            continue
        tokens.append(token)
    return _parse_info(";".join(tokens) or ".", projection)


def _frequency(ac, an):
    return ac / an if ac is not None and an is not None and an > 0 and 0 <= ac <= an else None


def _faf(value):
    if value is None or not math.isfinite(value):
        return None
    if not 0 <= value <= 1:
        raise ValueError
    return value


def parse_gnomad_record(
    line: str, *, expected: SequenceAllele, projection: GnomadProjectionContract,
    source: SourcedRecord, info_metadata=None,
) -> tuple[PopulationEvidence, ...]:
    """Pure parser for one exact six-column bcftools query record.

    This does NOT authenticate a resource. Only run_population_annotation and its
    replay validator authenticate files. Raw lexical INFO values preserve missing
    entries and nonfinite statistics without putting NaN/Inf into evidence codecs.
    """
    try:
        if type(source) is not SourcedRecord or type(line) is not str or "\n" in line.removesuffix("\n"):
            raise ValueError
        row = line.rstrip("\n").split("\t")
        if (len(row) != 6 or row[:4] != [expected.contig, str(expected.position),
                                        expected.reference, expected.alternate]
                or any(ord(c) < 32 or ord(c) > 126 for field in row for c in field)):
            raise ValueError
        metadata = (
            gnomad_vcf_metadata(projection)["info"]
            if info_metadata is None else info_metadata
        )
        values = _population_info(row[5], projection, metadata)
        raw = _raw_info(row[5])
        for key, value in values.items():
            items = value if isinstance(value, list) else [value]
            if metadata[key]["Type"] == "Integer" and any(v is not None and v < 0 for v in items):
                raise ValueError
            if (key == "DISC" or key.startswith("REGION_")) and value not in (None, 0, 1):
                raise ValueError
        disc = values["DISC"]
        common = {key: raw.get(key) for key in metadata if key.startswith(("DISC", "REGION_"))}
        records = []
        for dataset in projection.datasets:
            prefix = DATASET_VCF_PREFIXES[dataset.name]
            filters = None
            if dataset.has_filters and values[prefix + "_FILTERS_DEFINED"] == 1:
                filters = tuple(values[prefix + "_FILTERS"])
            strata = [("adjusted", None, None), *(
                ("ancestry", ancestry, i) for i, ancestry in enumerate(projection.major_genetic_ancestries)
            ), ("grpmax", values[prefix + "_GRPMAX_ANC"], None),
                ("fafmax95", values[prefix + "_FAFMAX95_ANC"], None),
                ("fafmax99", values[prefix + "_FAFMAX99_ANC"], None)]
            for stratum, ancestry, index in strata:
                keys = []

                def field(suffix):
                    key = prefix + "_" + suffix
                    keys.append(key)
                    value = values[key]
                    return value[index] if index is not None and value is not None else value

                ac = an = hom = hemi = faf95 = faf99 = None
                if stratum in ("adjusted", "ancestry", "grpmax"):
                    def count(suffix):
                        return field("GRPMAX_" + suffix if stratum == "grpmax" else
                                     suffix + ("_ANC" if stratum == "ancestry" else ""))
                    ac, an, hom = count("AC"), count("AN"), count("NHOMALT")
                    if stratum == "grpmax":
                        field("GRPMAX_AF")  # Retain publisher AF even when it disagrees with AC/AN.
                        field("GRPMAX_ANC")
                    else:
                        hemi = count("NHEMI")
                        faf95, faf99 = _faf(count("FAF95")), _faf(count("FAF99"))
                else:
                    suffix = "FAFMAX95" if stratum == "fafmax95" else "FAFMAX99"
                    value = _faf(field(suffix))
                    field(suffix + "_ANC")
                    if stratum == "fafmax95":
                        faf95 = value
                    else:
                        faf99 = value
                if stratum == "ancestry" and ancestry not in dataset.faf_ancestries and (faf95 is not None or faf99 is not None):
                    raise ValueError
                af = _frequency(ac, an)
                count_status = ("not_applicable" if stratum.startswith("fafmax") else
                                "missing_ac_or_an" if ac is None or an is None else
                                "inconsistent_ac_an" if ac > an else
                                "zero_an" if an == 0 else "valid_ac_an")
                filtered = any(f not in (None, "PASS") for f in filters or ())
                state = ("unknown_counts" if af is None else "observed_alternate" if ac > 0 else
                         "ac0_filtered_observation" if filtered else "observed_zero_alternate")
                if count_status == "not_applicable":
                    state = "published_faf_maximum" if faf95 is not None or faf99 is not None else "unknown_faf_maximum"
                identifier = f"{row[0]}:{row[1]}:{row[2]}:{row[3]}:{dataset.name}:{stratum}:{ancestry or 'overall'}"
                fields = {
                    "schema": SCHEMA, "stratum": stratum, "ancestry_index": index,
                    "ancestry_order": list(projection.major_genetic_ancestries),
                    "state": state, "exact_variant_record": True,
                    "count_status": count_status,
                    "raw_fields": {key: raw.get(key) for key in keys},
                    "raw_context": common, "projection_filter": row[4],
                    "source_filters_defined": values.get(prefix + "_FILTERS_DEFINED"),
                    "source_filters_raw": raw.get(prefix + "_FILTERS"),
                    "query_record_sha256": _digest(line.rstrip("\n") + "\n"),
                    "annotation_binding": source.source_fields,
                }
                records.append(_construct_owned(PopulationEvidence,
                    _construct_owned(SourcedRecord, source.source_dataset, source.source_release,
                        identifier, source.source_sha256, source.public_manifest_sha256, fields),
                    subset=dataset.name, ancestry=ancestry, stratum=stratum,
                    allele_count=ac, allele_number=an, allele_frequency=af,
                    filtering_allele_frequency=faf95, filtering_allele_frequency_99=faf99,
                    homozygote_count=hom, hemizygote_count=hemi,
                    callable=None if an is None or count_status == "inconsistent_ac_an" else an > 0,
                    filters=filters,
                    discordant=None if disc is None else bool(disc),
                ))
        return tuple(records)
    except MemoryError:
        raise
    except Exception:
        raise PopulationAnnotationError("invalid gnomAD population record") from None


def _selected(request):
    if type(request) is not PopulationRequest or type(request.allele) is not AlleleEvidence:
        raise ValueError
    request.allele.__post_init__()
    mappings = {m.identity_id: m for m in request.allele.identities}
    native = mappings[request.allele.native_identity_id]
    identity = native.identity
    if type(identity) not in (SequenceAllele, MosaicAllele):
        return None, "unsupported_native_class"
    sequence = identity.sequence if type(identity) is MosaicAllele else identity
    if request.normalized_identity_id is None:
        # An omitted selector cannot hide existing normalization, even when its
        # mapping is unlinked or contradictory and needs explicit rejection.
        if any(mapping.stage is IdentityStage.NORMALIZED for mapping in mappings.values()):
            raise PopulationAnnotationError("population normalized identity selector is required")
        return None, "no_authoritative_normalized_identity"
    selected = mappings.get(request.normalized_identity_id)
    if (selected is None or selected.stage is not IdentityStage.NORMALIZED
            or type(selected.identity) is not SequenceAllele or selected.identity != sequence
            or selected.reference_manifest_sha256 != native.reference_manifest_sha256):
        raise ValueError
    # Follow the supplied DAG; an unrelated normalized mapping cannot authorize a join.
    pending, ancestors = list(native.derived_from_identity_ids), set()
    while pending:
        key = pending.pop()
        if key not in ancestors:
            ancestors.add(key)
            pending.extend(mappings[key].derived_from_identity_ids)
    if selected.identity_id not in ancestors:
        raise ValueError
    if any(re.fullmatch(r"[ACGT]+", v) is None for v in (sequence.reference, sequence.alternate)):
        return None, "ambiguous_sequence"
    return selected, None


def _coverage_records(coverage, raw_rows, *, matched, sources, bundle_sha):
    records = []
    for subset in GNOMAD_SUBSETS:
        value = getattr(coverage, subset) if subset != "joint" else None
        an = coverage.joint_an if subset == "joint" else None if value is None else value.an
        state = ("coverage_only_observed_record" if matched else "callable_zero_alternate"
                 if an is not None and an > 0 else "unknown_no_callable_data")
        source = sources[subset] if subset != "joint" else SourcedRecord(
            ALL_SITES_SOURCE_DATASET, SOURCE_RELEASE, coverage.locus, bundle_sha, bundle_sha, None,
        )
        if type(source) is not SourcedRecord:
            raise ValueError('invalid population source record')
        fields = {
            "schema": SCHEMA, "state": state, "locus": coverage.locus,
            "exact_variant_record": matched, "coverage": _plain(value) if subset != "joint" else _plain(coverage),
            "raw_rows": raw_rows[subset] if subset != "joint" else raw_rows,
            "an_derivation": "exomes_plus_genomes_both_known" if subset == "joint" else "published_locus_AN",
            "source_bindings": {name: _plain(s) for name, s in sources.items()} if subset == "joint" else source.source_fields,
            "coverage_scope": "locus_only_not_ancestry_or_full_variant_span",
            "zero_observation_scope": (
                "exact_allele_absent_from_joint_projection_at_callable_locus"
                if not matched and an is not None and an > 0 else None
            ),
        }
        records.append(_construct_owned(PopulationEvidence,
            _construct_owned(SourcedRecord, source.source_dataset, source.source_release,
                coverage.locus + ":" + subset, source.source_sha256,
                source.public_manifest_sha256, fields),
            subset=subset, stratum=ALL_SITES_STRATUM, allele_number=an,
            callable=None if an is None else an > 0,
        ))
    return tuple(records)


def validate_joint_contract_bindings(projection, contract, validation, projection_sha256):
    """Use the consumer's current source/projection/validation relationships."""
    if (projection.source.assembly != "GRCh38" or projection.source.release != "v4.1"
            or contract.assembly != "GRCh38" or contract.release != "v4.1"
            or projection.major_genetic_ancestries != _ANCESTRIES
            or validation.ancestry_order != _ANCESTRIES
            or validation.projection_contract_sha256 != projection_sha256
            or contract.projection_id != projection.projection_id
            or validation.projection_id != projection.projection_id
            or validation.validation_id != contract.pilot_validation_id):
        raise PopulationAnnotationError("population current contract bindings mismatch")


def validate_joint_shard_header(header, projection, contig, reference_length):
    """Share the actual consumer's INFO and public-reference dictionary checks."""
    _validate_header(header, projection)
    dictionary = re.findall(r"^##contig=<ID=([^,>]+),length=([0-9]+)(?:,[^>]*)?>$", header, re.M)
    declared = dict(dictionary)
    if len(dictionary) != len(declared) or declared.get(contig) != str(reference_length):
        raise PopulationAnnotationError("population shard reference dictionary mismatch")


def _cached_region_text(cache, *, role, scope, contig, regions, fetch, limits):
    """Reuse per-point raw facts and rebuild the original full-batch response.

    Joint queries include upstream deletions overlapping a requested point.
    Each point retains those rows; their union is reconciled by the existing
    complete-response validator before any new fact is committed.
    """
    if cache is None:
        fetched = fetch(regions)
        return fetched.text, lambda _text:None, {region:fetched.receipt for region in regions}
    known,missing = {},[]
    def key(region):
        return dict(scope,contig=contig,region=region)
    def belongs(line,region):
        fields=line.split('\t')
        point=int(region.rsplit(':',1)[1].split('-',1)[0])
        if len(fields)<3 or fields[0]!=contig or not fields[1].isdecimal():
            raise ValueError
        position=int(fields[1])
        return (position<=point<position+len(fields[2]) if role=='population_joint' else position==point)
    for region in regions:
        fact=cache.get(role,key(region))
        if fact is None:
            missing.append(region)
        else:
            raw = _raw_population_fact(fact, scope=scope, region=region)
            if any(not belongs(line,region) for line in raw.text.splitlines()):
                raise ValueError
            known[region]=raw
    fetched = None if not missing else fetch(tuple(missing))
    fresh=() if fetched is None else _lines(fetched.text,limits)
    if len(fresh)!=len(set(fresh)) or any(not any(belongs(line,r) for r in missing) for line in fresh):
        raise ValueError
    lines=set(fresh)
    for raw in known.values():lines.update(raw.text.splitlines())
    text=''.join(line+'\n' for line in sorted(lines))
    _lines(text,limits)
    def save(validated_text):
        validated=validated_text.splitlines()
        for region in missing:
            values=tuple(line for line in validated if belongs(line,region))
            cache.put(role,key(region),status='positive' if values else 'no_match',
                value=dict(lines=values, receipt=fetched.receipt.decode()))
    receipts = {region:raw.receipt for region,raw in known.items()}
    receipts.update((region,fetched.receipt) for region in missing)
    return text,save,receipts


def _span_query_results(spans, *, fetch, limits, cache=None, scope=None):
    """Bounded physical queries; retain exact logical per-span evidence.

    Overlapping regions are merged before tabix so a source row cannot appear
    twice merely because two requested spans overlap. Source duplicates and
    contradictory rows still fail. Cached facts keep the original per-span key.
    """
    pending = []
    cost = argument_bytes = 0

    def region(span):
        contig, start, end = span
        if (type(contig) is not str or not re.fullmatch(r'(?:chr)?(?:[1-9]|1[0-9]|2[0-2])', contig)
                or type(start) is not int or type(end) is not int or not 0 < start <= end):
            raise ValueError
        return f'{contig}:{start}-{end}'

    def key(span):
        return dict(scope or {}, noncontribution_span=region(span))

    def validate(text, span):
        contig, start, end = span
        lines = _lines(text, limits)
        # This strict owner checks duplicates, geometry and every source field,
        # including callable rows that do not yield a noncontribution proof.
        sites.noncontribution_query_basis(text, contig=contig, start=start, end=end,
            status='matched' if text else 'no_match')
        return tuple(sorted(lines, key=lambda line: int(line.split('\t')[1])))

    def execute(batch):
        merged = []
        for contig, start, end in batch:
            if merged and merged[-1][0] == contig and start <= merged[-1][2] + 1:
                old = merged[-1]; merged[-1] = (contig, old[1], max(old[2], end))
            else:
                merged.append((contig, start, end))
        fetched = fetch(tuple(region(span) for span in merged))
        text = fetched.text
        lines = _lines(text, limits)
        indexed = {}
        for line in lines:
            contig, position, _ = sites._parse_derived_record(line)
            if ((contig, position) in indexed or
                    not any(contig == c and start <= position <= end for c, start, end in merged)):
                raise ValueError
            indexed[contig, position] = line
        for span in batch:
            contig, start, end = span
            content = ''.join(line+'\n' for (c, pos), line in sorted(indexed.items())
                if c == contig and start <= pos <= end)
            values = validate(content, span)
            if cache is not None:
                cache.put('population_all_sites', key(span),
                    status='positive' if values else 'no_match',
                    value=dict(lines=values, receipt=fetched.receipt.decode()))
            yield span, _PopulationQueryText(content, fetched.receipt)

    for span in sorted(set(spans)):
        encoded = region(span)
        if len(encoded.encode('ascii')) > limits.max_region_argument_bytes:
            raise PopulationAnnotationError('population query region exceeded argument bound')
        fact = None if cache is None else cache.get('population_all_sites', key(span))
        if fact is not None:
            raw = _raw_population_fact(fact, scope=scope, region=encoded)
            values = validate(raw.text, span)
            yield span, _PopulationQueryText(''.join(line+'\n' for line in values), raw.receipt)
            continue
        added = len(encoded.encode('ascii')) + bool(pending)
        width = span[2] - span[1] + 1
        if pending and (span[0] != pending[0][0] or len(pending) >= limits.max_loci_per_query
                or argument_bytes + added > limits.max_region_argument_bytes
                or cost + width > limits.max_rows_per_query):
            yield from execute(pending)
            pending = []; cost = argument_bytes = 0
            added = len(encoded.encode('ascii'))
        pending.append(span); cost += width; argument_bytes += added
    if pending:
        yield from execute(pending)


@_validation_session(isolated=True)
def _cached_population_adaptation(cache, *, key, produce, source=None):
    """Adapt invariant source facts once; bind the current batch on each use."""
    if cache is None:
        return produce()
    found=cache.population_adaptation(key)
    cold=found is None
    if found is None:
        cache.metrics['population_adaptations']+=1
        rows=produce()
        status=('positive' if key['kind']=='joint' or key.get('matched') else
            'no_match' if any(row.callable for row in rows) else 'uncertain')
        found=cache._put_population_adaptation(key,status=status,rows=rows)
    fact,rows=found
    expected=('positive' if key['kind']=='joint' or key.get('matched') else
        'no_match' if any(row.callable for row in rows) else 'uncertain')
    if fact.status!=expected:
        raise ValueError
    if source is not None and not cold:
        # Batch IDs/query digests belong to the current report, not the reusable
        # population counts. A cold producer already used this exact source;
        # only warm facts need rebinding. Preserve every other source field.
        rows=tuple(replace(row,source=replace(row.source,source_fields={
            **row.source.source_fields,'annotation_binding':source.source_fields})) for row in rows)
    return rows


def _run_population_annotation(
    requests: tuple[PopulationRequest, ...], *, joint: GnomadJointReference,
    all_sites: GnomadAllSitesReference, reference: PopulationReferenceCompatibility,
    runtime_image: str, limits: PopulationQueryLimits, file_cache: _PopulationFileCache,
    case_cache=None, noncontribution_output=None, command_metrics=None,
    screening_output=None, work_counts=None,
) -> PopulationAnnotationResult:
    """Attach a complete bounded batch using read-only files in an offline container.

    No public callback, download/bootstrap path, write, or network fallback exists.
    Only selected BCF/CSI shards are opened. Completed publication manifests
    supply trusted output digests; selected data files are checked for safe paths,
    exact size, stable metadata and query semantics without a fresh full hash.
    runtime_image is an externally verified OCI pin, not self-certified here.
    """
    try:
        if (type(requests) is not tuple or type(joint) is not GnomadJointReference
                or type(all_sites) is not GnomadAllSitesReference
                or type(reference) is not PopulationReferenceCompatibility
                or type(limits) is not PopulationQueryLimits
                or type(runtime_image) is not str
                or re.fullmatch(r"[A-Za-z0-9_./:-]+@sha256:[0-9a-f]{64}", runtime_image) is None):
            raise ValueError
        limits.__post_init__()
        if (len(requests) > limits.max_alleles
                or len({r.allele.allele_id for r in requests}) != len(requests)):
            raise ValueError
        if any(name != "lo" for _, name in socket.if_nameindex()):
            raise PopulationAnnotationError("population annotation requires an isolated network namespace")
        selected = [_selected(r) for r in requests]
        # Batch-local encodings serve both normalization and report identity.
        # The resource session must never retain request/evidence objects.
        encoded_inputs = tuple(to_json(request.allele) for request in requests)
        if work_counts is not None: work_counts['input_encodes'] += len(requests)
        files: dict[Path, tuple[str, int]] = {}

        def verify(path, sha=None, size=None, *, published_by=None):
            path = file_cache.verify(path, sha, size, published_by=published_by)
            cached = file_cache.files[path]
            binding = (cached.sha256, cached.size_bytes)
            if path in files and files[path] != binding:
                raise PopulationAnnotationError("population resources changed during lookup")
            files[path] = binding
            if path.suffix == ".json" and cached.size_bytes > limits.max_output_bytes:
                raise PopulationAnnotationError("population metadata exceeded byte bound")
            return path

        def document(path, sha=None):
            path = verify(path, sha)
            if files[path][1] > limits.max_output_bytes:
                raise PopulationAnnotationError("population metadata exceeded byte bound")
            return file_cache.prepared(
                ("document", path),
                lambda: json.loads(path.read_text(encoding="utf-8")),
            )

        def manifest(path, sha=None):
            path = verify(path, sha)
            return file_cache.prepared(
                ("manifest", path),
                lambda: PublicResourceManifest.from_mapping(document(path, sha)),
            )

        adapter_path = verify(Path(__file__))
        projection_path = verify(joint.projection_contract_path)
        publication_path = verify(joint.publication_contract_path)
        validation_path = verify(joint.validation_contract_path)
        projection = file_cache.prepared(
            ("projection", projection_path),
            lambda: load_gnomad_projection_contract(projection_path),
        )
        contract = file_cache.prepared(
            ("publication_contract", publication_path),
            lambda: publication.load_gnomad_publication_contract(publication_path),
        )
        validation = file_cache.prepared(
            ("validation_contract", validation_path),
            lambda: load_gnomad_validation_contract(validation_path),
        )
        projection_info = file_cache.prepared(
            ("projection_info", projection_path),
            lambda: gnomad_vcf_metadata(projection)["info"],
        )
        file_cache.prepared(
            ("joint_contract_bindings", files[projection_path][0],
             files[publication_path][0], files[validation_path][0]),
            lambda: validate_joint_contract_bindings(
                projection, contract, validation, files[projection_path][0],
            ),
        )
        joint_manifest = manifest(joint.manifest_path, joint.manifest_sha256)
        transformation = file_cache.prepared(
            ("transformation", files[publication_path][0],
             files[projection_path][0], files[validation_path][0]),
            lambda: publication._transformation_parameters_sha256(
                publication_path, projection_path, validation_path, validation,
            ),
        )
        file_cache.prepared(
            ("joint_manifest_contract", joint.manifest_sha256,
             files[publication_path][0], files[projection_path][0], transformation),
            lambda: publication._verify_manifest_contract(
                joint_manifest, contract, projection, transformation,
            ),
        )
        joint_outputs = {o.path: o for o in joint_manifest.outputs}
        if sum(o.size_bytes for o in joint_manifest.outputs) > contract.aggregate_output_cap_bytes:
            raise ValueError

        # Authenticate reference naming metadata and actual sequence bytes. Only
        # explicitly selected roles are opened, not the full reconstructed bundle.
        identity_path = verify(reference.identity_path, reference.identity_sha256)
        identity = file_cache.prepared(
            ("reference_identity", identity_path),
            lambda: load_reference_identity(identity_path),
        )
        if identity.primary_assembly_accession != "GCA_000001405.15":
            raise ValueError
        reference_manifest = manifest(reference.manifest_path, reference.manifest_sha256)
        file_cache.prepared(
            ("reference_manifest_contract", files[identity_path][0],
             reference.manifest_sha256, reference.identity_sha256),
            lambda: verify_reference_manifest_contract(
                identity, reference_manifest,
                contract_sha256=reference.identity_sha256,
            ),
        )
        reference_outputs = {o.path: o for o in reference_manifest.outputs}
        reference_root = _local(reference.reference_root)
        alias_path = identity.expected_output("contig_aliases").path
        alias_output = reference_outputs[alias_path]
        verify(reference_root / alias_path, alias_output.sha256, alias_output.size_bytes,
               published_by=reference.manifest_sha256)
        roles = {ReferenceRole.PUBLIC_ANNOTATION_BENCHMARK.value}
        roles.update(m.identity.reference_role.value for m, reason in selected if m is not None)
        names, fastas = {}, {}
        for role_name in sorted(roles):
            role = identity.reference_role(role_name)
            paths = []
            for output_role in (role.fasta_output_role, role.fai_output_role):
                output = reference_outputs[identity.expected_output(output_role).path]
                paths.append(verify(
                    reference_root / output.path, output.sha256, output.size_bytes,
                    published_by=reference.manifest_sha256,
                ))
            fastas[role_name] = file_cache.prepared(
                ("fasta", *paths), lambda: IndexedFasta(*paths),
            )
            names[role_name] = file_cache.prepared(
                ("renderer", role_name, identity_path, reference.manifest_path),
                lambda: load_submission_contig_renderer(
                    identity_path=identity_path, reference_manifest_path=reference.manifest_path,
                    reference_root=reference_root,
                    reference_binding=ReferenceRoleBinding(role_name, role.identity_id),
                ),
            )
        public_role = identity.reference_role(ReferenceRole.PUBLIC_ANNOTATION_BENCHMARK.value)
        public_fasta = fastas[public_role.role]
        canonical_public = {canonical: name for name, (canonical, _) in names[public_role.role]._names.items()}
        resource_alleles = []
        reasons = []
        for request, (mapping, reason), encoded_input in zip(
                requests, selected, encoded_inputs, strict=True):
            resource_allele = None
            if mapping is not None:
                sequence = mapping.identity
                role = identity.reference_role(sequence.reference_role.value)
                if (mapping.reference_manifest_sha256 != reference.manifest_sha256
                        or sequence.reference_identity_id != role.identity_id):
                    raise PopulationAnnotationError("population source reference binding mismatch")
                name = names[role.role].accepted_names.get(sequence.contig)
                canonical = None if name is None else name[0]
                contig = canonical_public.get(canonical)
                if canonical == "M":
                    reason = "mitochondrial_outside_joint_product"
                elif contig not in contract.canonical_contigs:
                    reason = "contig_outside_joint_product"
                else:
                    source_fasta = fastas[role.role]
                    source_names = tuple(native for native, (key, _) in names[role.role]._names.items()
                                         if key == canonical)
                    if len(source_names) != 1:
                        raise PopulationAnnotationError("population source contig mapping is ambiguous")
                    source_contig, = source_names
                    if (source_fasta.fetch(source_contig, sequence.position, len(sequence.reference)) != sequence.reference
                            or public_fasta.fetch(contig, sequence.position, len(sequence.reference)) != sequence.reference):
                        raise PopulationAnnotationError("population reference allele mismatch")
                    # Verify, never replace, the selected normalization. Bound any
                    # upstream reference walk on an incorrectly labelled input.
                    fetched = 0

                    def fetch(c, p, n):
                        nonlocal fetched
                        if c != sequence.contig:
                            raise PopulationAnnotationError("population normalization changed source contig")
                        fetched += n
                        if fetched > limits.max_output_bytes:
                            raise PopulationAnnotationError("normalization verification exceeded bound")
                        return source_fasta.fetch(source_contig, p, n)

                    source_record = SourceSmallVariantRecord(1, sequence.contig, sequence.position,
                                                            sequence.reference, (sequence.alternate,))
                    input_sha = _digest(encoded_input)
                    record_id = source_record_id(input_sha, source_record)
                    source_allele = SourceAllele(source_allele_id(record_id, 1), record_id, input_sha,
                                                source_record, 1, sequence.alternate)
                    normalized = _normalize_sequence_allele(
                        source_allele, reference_identity_id=sequence.reference_identity_id,
                        reference_role=sequence.reference_role.value, fetch_reference=fetch,
                    )
                    if (normalized.position, normalized.reference, normalized.alternate) != (
                        sequence.position, sequence.reference, sequence.alternate,
                    ):
                        raise PopulationAnnotationError("population input is not authoritative normalized sequence")
                    resource_allele = SequenceAllele(ReferenceRole.PUBLIC_ANNOTATION_BENCHMARK,
                        public_role.identity_id, contig, sequence.position, sequence.reference, sequence.alternate)
            resource_alleles.append(resource_allele)
            reasons.append(reason)

        tools = {}
        for tool in ("bcftools", "tabix"):
            executable = shutil.which(tool)
            if executable is None:
                raise PopulationAnnotationError("population runtime is missing bcftools or tabix")
            tools[tool] = verify(executable)

        def command(tool, *args):
            if command_metrics is not None:
                command_metrics[tool] += 1
            return _command(tools[tool], args, limits)

        shards, headers = {}, {}
        for contig in sorted({a.contig for a in resource_alleles if a is not None}):
            name = publication._expected_shard_name(contig)
            path = _local(joint.release_directory) / name
            for filename in (name, name + ".csi"):
                output = joint_outputs[contract.output_directory + "/" + filename]
                if output.retention.policy != "durable" or output.retention.delete_after is not None:
                    raise ValueError
                verified = verify(
                    path.parent / filename, output.sha256, output.size_bytes,
                    published_by=joint.manifest_sha256,
                )
                with gzip.open(verified, "rb") as handle:
                    if handle.read(4 if filename.endswith(".csi") else 5) != (
                        b"CSI\x01" if filename.endswith(".csi") else b"BCF\x02\x02"
                    ):
                        raise ValueError
            def prepare_joint_header():
                header = command("bcftools", "view", "--header-only", "--no-version", str(path))
                validate_joint_shard_header(
                    header, projection, contig, public_fasta.index[contig].length,
                )
                index_stats = command("bcftools", "index", "--stats", str(path))
                stats = index_stats.rstrip("\n").split("\t")
                if (len(stats) != 3 or stats[0] != contig or not stats[2].isdecimal()
                        or int(stats[2]) <= 0):
                    raise ValueError
                return {"header_sha256": _digest(header), "index_stats": index_stats}

            shards[contig] = path
            headers[contig] = file_cache.prepared(
                ("joint_header", path, projection_path, contig,
                 public_fasta.index[contig].length),
                prepare_joint_header,
            )

        all_contract_path = verify(all_sites.contract_path)
        all_contract = file_cache.prepared(
            ("all_sites_contract", all_contract_path),
            lambda: sites.load_gnomad_all_sites_contract(all_contract_path),
        )
        all_root = _local(all_sites.reference_root)
        bundle = document(all_root / all_contract.bundle_manifest_path, all_sites.bundle_manifest_sha256)
        source_hashes, transforms, coverage_paths, coverage_sources, indexed_contigs = {}, [], {}, {}, {}
        entries = file_cache.prepared(
            ("all_sites_bundle", all_contract_path, all_sites.bundle_manifest_sha256),
            lambda: sites.validate_bundle_metadata(
                bundle, all_contract, files[all_contract_path][0],
            ),
        )
        for source, (entry, transform) in zip(all_contract.sources, entries, strict=True):
            source_hashes[source.name] = entry["source_manifest_sha256"]
            transforms.append(transform)
            current = manifest(all_root / entry["source_manifest_path"], entry["source_manifest_sha256"])
            file_cache.prepared(
                ("all_sites_source_contract", entry["source_manifest_sha256"],
                 files[all_contract_path][0], source.name),
                lambda: sites._verify_source_manifest_contract(
                    current, all_contract, source, files[all_contract_path][0],
                ),
            )
            for output in current.outputs:
                if output.retention.policy != "durable" or output.retention.delete_after is not None:
                    raise ValueError
                verify(
                    all_root / output.path, output.sha256, output.size_bytes,
                    published_by=entry["source_manifest_sha256"],
                )
            path = all_root / all_contract.output_directory / source.output_name
            def prepare_all_sites_contigs():
                with gzip.open(path, "rb") as handle:
                    if (handle.readline(limits.max_output_bytes + 1).decode("ascii")
                            .rstrip("\n") != sites.DERIVED_HEADER):
                        raise ValueError
                with gzip.open(str(path) + ".tbi", "rb") as handle:
                    if handle.read(4) != b"TBI\x01":
                        raise ValueError
                contigs = tuple(command("tabix", "--list-chroms", str(path)).splitlines())
                if contigs != transform.contigs:
                    raise ValueError
                return contigs

            contigs = file_cache.prepared(
                ("all_sites_contigs", path, entry["source_manifest_sha256"]),
                prepare_all_sites_contigs,
            )
            indexed_contigs[source.name] = contigs
            coverage_paths[source.name] = path
            coverage_sources[source.name] = SourcedRecord(
                ALL_SITES_SOURCE_DATASET, SOURCE_RELEASE, source.name, files[path][0], entry["source_manifest_sha256"],
                {"artifact": str(path), "index_sha256": files[Path(str(path) + ".tbi")][0],
                 "bundle_manifest_sha256": all_sites.bundle_manifest_sha256,
                 "contract_sha256": files[all_contract_path][0]},
            )
        if bundle != _plain(sites._bundle_document(
            contract=all_contract, contract_sha256=files[all_contract_path][0], image_digest=bundle["image_digest"],
            retrieved_on=bundle["retrieved_on"], transforms=transforms, source_manifest_sha256s=source_hashes,
        )) or re.fullmatch(r"sha256:[0-9a-f]{64}", bundle["image_digest"]) is None:
            raise ValueError

        fact_scope = {
            'adapter_sha256':files[adapter_path][0], 'limits':_plain(limits),
            'runtime_image':runtime_image, 'tools':{name:files[path][0] for name,path in tools.items()},
            'tool_paths':{name:str(path) for name,path in tools.items()},
            'joint_manifest_sha256':joint.manifest_sha256,
            'all_sites_manifest_sha256':all_sites.bundle_manifest_sha256,
            'reference_manifest_sha256':reference.manifest_sha256,
            'reference_identity_sha256':reference.identity_sha256,
            'projection_contract_sha256':files[projection_path][0],
            'all_sites_contract_sha256':files[all_contract_path][0],
        }
        def query(tool, scope, regions, *args):
            # Hash the actual complete stdout before logical sorting/splitting.
            content = command(tool, *args)
            receipt = _json(dict(schema='genosapiens.population-physical-query', version=1,
                tool=tool, tool_sha256=files[tools[tool]][0], runtime_image=runtime_image,
                argv=(str(tools[tool]), *args), regions=regions, resource_scope=scope,
                stdout_sha256=_digest(content), stdout_bytes=len(content.encode()), exit_status=0)).encode()
            if work_counts is not None: work_counts['physical_queries'] += 1
            return _PopulationQueryText(content, receipt)

        query_batches = []
        joint_rows = {}
        joint_batches = {}
        receipts_by_locus = {}
        coverage_rows = {name: {} for name in coverage_paths}
        coverage_batches = {name: {} for name in coverage_paths}
        for contig, regions in _locus_batches(resource_alleles, limits):
            requested_positions = {
                int(region.rsplit(":", 1)[1].split("-", 1)[0]) for region in regions
            }
            scope = dict(fact_scope,artifact=str(shards[contig]),artifact_sha256=files[shards[contig]][0],
                index_sha256=files[Path(str(shards[contig])+'.csi')][0])
            def fetch_joint(selected_regions):
                return query(
                "bcftools", scope, selected_regions, "query", "--regions", ",".join(selected_regions),
                "--format", "%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER\\t%INFO\\n",
                str(shards[contig]),
            )
            text,save_joint,receipts = _cached_region_text(case_cache,role='population_joint',scope=scope,
                contig=contig,regions=regions,fetch=fetch_joint,limits=limits)
            for region,receipt in receipts.items():
                locus = region.split('-')[0]
                receipts_by_locus.setdefault(locus, set()).add(receipt)
            lines = _lines(text, limits)
            seen = set()
            for line in lines:
                row = line.split("\t")
                if (
                    len(row) != 6
                    or row[0] != contig
                    or not row[1].isdecimal()
                    or any(re.fullmatch(r"[ACGT]+", value) is None for value in row[2:4])
                    or row[2] == row[3]
                ):
                    raise ValueError
                position = int(row[1])
                # Region jumps may return an upstream deletion overlapping a
                # requested point. Account for it, but attach only exact POS.
                if not any(position <= requested < position + len(row[2])
                           for requested in requested_positions):
                    raise ValueError
                key = tuple(row[:4])
                if key in seen:
                    raise ValueError
                seen.add(key)
                _population_info(row[5], projection, projection_info)
            text = "".join(line + "\n" for line in sorted(lines))
            save_joint(text)
            batch = _batch_record(
                tool="bcftools", source=JOINT_SOURCE_DATASET, contig=contig,
                regions=regions, content=text,
            )
            query_batches.append(batch)
            for line in text.splitlines():
                row = line.split("\t")
                position = int(row[1])
                if position in requested_positions:
                    joint_rows.setdefault(f"{contig}:{position}", []).append(line)
            for position in requested_positions:
                joint_batches[f"{contig}:{position}"] = batch

            for name, path in coverage_paths.items():
                if contig not in indexed_contigs[name]:
                    batch = _batch_record(
                        tool="tabix", source=name, contig=contig,
                        regions=regions, content="",
                    )
                    batch["status"] = "contig-not-indexed"
                    query_batches.append(batch)
                    for position in requested_positions:
                        coverage_rows[name][f"{contig}:{position}"] = ""
                        coverage_batches[name][f"{contig}:{position}"] = batch
                    continue
                scope = dict(fact_scope,source=name,artifact=str(path),artifact_sha256=files[path][0],
                    index_sha256=files[Path(str(path)+'.tbi')][0])
                text,save_coverage,receipts = _cached_region_text(case_cache,role='population_all_sites',scope=scope,
                    contig=contig,regions=regions,
                    fetch=lambda selected_regions:query('tabix',scope,selected_regions,str(path),*selected_regions),limits=limits)
                for region,receipt in receipts.items():
                    receipts_by_locus[region.split('-')[0]].add(receipt)
                lines = _lines(text, limits)
                seen_positions = set()
                for line in lines:
                    parsed_contig, position, _ = sites._parse_derived_record(line)
                    if (parsed_contig != contig or position not in requested_positions
                            or position in seen_positions):
                        raise ValueError
                    seen_positions.add(position)
                text = "".join(line + "\n" for line in sorted(
                    lines, key=lambda item: (int(item.split("\t", 2)[1]), item)
                ))
                save_coverage(text)
                batch = _batch_record(
                    tool="tabix", source=name, contig=contig,
                    regions=regions, content=text,
                )
                query_batches.append(batch)
                for line in text.splitlines():
                    parsed_contig, position, _ = sites._parse_derived_record(line)
                    coverage_rows[name][f"{parsed_contig}:{position}"] = line + "\n"
                for position in requested_positions:
                    locus = f"{contig}:{position}"
                    coverage_rows[name].setdefault(locus, "")
                    coverage_batches[name][locus] = batch

        coverage_cache = {}
        for locus in sorted(joint_batches):
            coverage, raw_rows = {}, {}
            for name in coverage_paths:
                raw = coverage_rows[name][locus]
                lines = raw.splitlines()
                value = None
                if lines:
                    contig, position, value = sites._parse_derived_record(lines[0])
                    if f"{contig}:{position}" != locus:
                        raise ValueError
                coverage[name], raw_rows[name] = value, raw
            coverage_cache[locus] = (sites.AllSitesLocusCoverage(locus, **coverage), raw_rows)

        alleles, assessments, proof_candidates = [], [], []
        screening_records, receipts_by_allele = [], {}
        for request, allele, reason in zip(requests, resource_alleles, reasons, strict=True):
            evidence = ()
            assessment = {"allele_id": request.allele.allele_id, "native_identity_id": request.allele.native_identity_id,
                          "normalized_identity_id": request.normalized_identity_id,
                          "resource_allele": _plain(allele), "reason": reason, "status": "unsupported"}
            if allele is None and case_cache is not None:
                # Unsupported is an explicit source fact, never a no-match.
                key = dict(fact_scope, kind='unsupported', request=_plain(request), reason=reason)
                fact = case_cache.get('population_adapted', key)
                if fact is None:
                    case_cache.put('population_adapted', key, status='unsupported', value=[])
                elif fact.status != 'unsupported' or fact.value != ():
                    raise ValueError
            if allele is not None:
                locus = f"{allele.contig}:{allele.position}"
                matches = [line for line in joint_rows.get(locus, ()) if line.split("\t")[:4] ==
                           [allele.contig, str(allele.position), allele.reference, allele.alternate]]
                if len(matches) > 1:
                    raise ValueError
                path = shards[allele.contig]
                joint_batch = joint_batches[locus]
                source = SourcedRecord(JOINT_SOURCE_DATASET, SOURCE_RELEASE, locus, files[path][0], joint.manifest_sha256,
                    {"artifact": str(path), "index_sha256": files[Path(str(path) + ".csi")][0],
                     "projection_contract_sha256": files[projection_path][0],
                     "query_batch_id": joint_batch["batch_id"],
                     "query_sha256": joint_batch["output_sha256"]})
                if matches:
                    key=dict(fact_scope,kind='joint',resource_allele=_plain(allele),raw_line=matches[0],
                        source={**_plain(source),'source_fields':{k:v for k,v in source.source_fields.items()
                            if k not in ('query_batch_id','query_sha256')}})
                    evidence = _cached_population_adaptation(case_cache,key=key,source=source,
                        produce=lambda:parse_gnomad_record(matches[0],expected=allele,projection=projection,
                            source=source,info_metadata=projection_info))
                coverage, raw_rows = coverage_cache[locus]
                evidence += _cached_population_adaptation(case_cache,
                    key=dict(fact_scope,kind='coverage',coverage=_plain(coverage),raw_rows=raw_rows,
                        matched=bool(matches),sources=_plain(coverage_sources)),
                    produce=lambda:_coverage_records(coverage,raw_rows,matched=bool(matches),sources=coverage_sources,
                        bundle_sha=all_sites.bundle_manifest_sha256))
                assessment.update(status="matched" if matches else "no_match", query_locus=locus,
                                  query_batch_id=joint_batch["batch_id"],
                                  query_sha256=joint_batch["output_sha256"],
                                  all_sites_batch_ids={name: coverage_batches[name][locus]["batch_id"]
                                                       for name in coverage_paths},
                                  all_sites_rows=raw_rows)
            population = (request.allele.population or ()) + evidence
            receipts_by_allele[request.allele.allele_id] = (set() if allele is None
                else set(receipts_by_locus[f'{allele.contig}:{allele.position}']))
            if noncontribution_output is not None and allele is not None:
                from genosapiens.common_allele_filter import PopulationNonContributionProof
                adjusted_exomes = tuple(row for row in evidence if row.subset == 'exomes' and row.stratum == 'adjusted')
                if (len(adjusted_exomes) == 1 and allele.contig in indexed_contigs['exomes']
                        and re.fullmatch(r'(?:chr)?(?:[1-9]|1[0-9]|2[0-2])', allele.contig)):
                    proof_candidates.append((request, allele, adjusted_exomes[0],
                        raw_rows['exomes'], coverage_batches['exomes'][locus]))
            if screening_output is not None:
                screening_records.append(population)
            else:
                attached = replace(request.allele, population=population) if evidence else request.allele
                alleles.append(attached)
                assessment["appended_evidence_count"] = len(evidence)
                assessment["output_evidence_sha256"] = _digest(to_json(attached))
                assessments.append(assessment)
                if work_counts is not None:
                    work_counts['decorated_alleles'] += bool(evidence)
                    work_counts['report_encodes'] += 1

        # Keep only compact proofs/results across query batches, never all span
        # output text. The logical evidence order remains the request order.
        span_candidates = {}
        proof_results = [None] * len(proof_candidates)
        span_batches = [None] * len(proof_candidates)
        def make_proof(request, allele, subset, text, batch):
            mapping = next(m for m in request.allele.identities if m.identity_id == request.normalized_identity_id)
            sequence = mapping.identity
            start, end = allele.position, allele.position+len(allele.reference)-1
            status = ('matched' if text else 'no_match') if batch.get('status') is None else batch['status']
            basis = sites.noncontribution_query_basis(text, contig=allele.contig, start=start, end=end, status=status)
            if basis is not None:
                raw_context = subset.source.source_fields.get('raw_context', {})
                fields = coverage_sources['exomes'].source_fields
                return PopulationNonContributionProof(
                    subset='exomes', assembly='GRCh38', reference_identity_id=sequence.reference_identity_id,
                    reference_manifest_sha256=reference.manifest_sha256, normalized_allele_id=request.normalized_identity_id,
                    contig=sequence.contig, position=sequence.position, reference=sequence.reference,
                    alternate=sequence.alternate, checked_start=start, checked_end=end,
                    coordinate_convention='one_based_inclusive', all_sites_resource_sha256=coverage_sources['exomes'].source_sha256,
                    all_sites_contract_sha256=fields['contract_sha256'], all_sites_index_sha256=fields['index_sha256'],
                    query_identity_sha256=_digest(_json(dict(batch=batch, allele=_plain(allele)))),
                    query_status=status, query_output_sha256=batch['output_sha256'], proof_basis=basis,
                    raw_variant_ac=subset.allele_count, raw_variant_an=subset.allele_number,
                    source_fields_sha256=_digest(_json(dict(raw_context))))


        for index, (request, allele, subset, anchor_text, anchor_batch) in enumerate(proof_candidates):
            span = (allele.contig, allele.position, allele.position+len(allele.reference)-1)
            if span[1] == span[2]:
                proof_results[index] = make_proof(request, allele, subset, anchor_text, anchor_batch)
            else:
                span_candidates.setdefault(span, []).append((index, request, allele, subset))
        if span_candidates:
            path = coverage_paths['exomes']
            scope = dict(fact_scope, source='exomes', artifact=str(path), artifact_sha256=files[path][0],
                index_sha256=files[Path(str(path)+'.tbi')][0])
            for (contig, start, end), raw in _span_query_results(tuple(span_candidates),
                    fetch=lambda regions: query('tabix', scope, regions, str(path), *regions),
                    limits=limits, cache=case_cache, scope=scope):
                batch = _batch_record(tool='tabix', source='exomes', contig=contig,
                    regions=(f'{contig}:{start}-{end}',), content=raw.text)
                for index, request, allele, subset in span_candidates[contig, start, end]:
                    receipts_by_allele[request.allele.allele_id].add(raw.receipt)
                    span_batches[index] = batch
                    proof_results[index] = make_proof(request, allele, subset, raw.text, batch)
        query_batches.extend(batch for batch in span_batches if batch is not None)
        if noncontribution_output is not None:
            noncontribution_output.extend(proof for proof in proof_results if proof is not None)
            if work_counts is not None: work_counts['proofs'] += sum(proof is not None for proof in proof_results)

        # Later batches reuse authenticated bytes only while the same regular
        # files remain in place. This is a metadata check, not another full hash.
        file_cache.assert_stable()
        file_cache.report_files(files)
        if screening_output is not None:
            screening_output.extend(_PopulationScreeningRow(request.allele.allele_id,
                request.normalized_identity_id, _digest(encoded), population,
                tuple(sorted(receipts_by_allele[request.allele.allele_id])))
                for request,encoded,population in zip(requests,encoded_inputs,screening_records,strict=True))
            if work_counts is not None: work_counts['screening_rows'] += len(requests)
            return None
        binding = {
            "schema": SCHEMA, "classification": "controlled", "llm_policy": "no-llm",
            "runtime_image": runtime_image, "limits": _plain(limits),
            "adapter_sha256": files[adapter_path][0],
            "input_sha256": _digest(_json([{"allele": json.loads(encoded),
                                           "normalized_identity_id": r.normalized_identity_id}
                                          for r, encoded in zip(requests, encoded_inputs, strict=True)])),
            "resources": {"joint": _plain(joint), "all_sites": _plain(all_sites), "reference": _plain(reference)},
            "files": [{"path": str(p), "sha256": sha, "size_bytes": size} for p, (sha, size) in sorted(files.items())],
            "scope": "selected_joint_shards;both_all_sites_products;selected_reference_roles",
            "headers": headers, "query_batches": query_batches, "assessments": assessments,
            "query_summary": {
                "requested_alleles": len(requests),
                "supported_alleles": sum(allele is not None for allele in resource_alleles),
                "unique_loci": len(joint_batches),
                "bcftools_invocations": sum(item["tool"] == "bcftools" for item in query_batches),
                "tabix_invocations": sum(item["tool"] == "tabix" and item.get("status") is None
                                          for item in query_batches),
            },
            "tools": {name: {"path": str(p), "sha256": files[p][0]} for name, p in tools.items()},
        }
        # Paths are local execution bindings, not JSON-native dataclass leaves.
        binding["resources"] = {name: {k: str(v) if isinstance(v, Path) else v for k, v in value.items()}
                                for name, value in binding["resources"].items()}
        report = SourcedRecord("genosapiens.population-annotation", "current", "batch", _digest(_json(binding)),
                               joint.manifest_sha256, binding)
        return PopulationAnnotationResult(tuple(alleles), report)
    except (PopulationAnnotationError, MemoryError):
        raise
    except MemoryError:
        raise
    except Exception:
        raise PopulationAnnotationError("offline population annotation validation failed", category="unexpected") from None


def run_population_annotation(
    requests: tuple[PopulationRequest, ...], *, joint: GnomadJointReference,
    all_sites: GnomadAllSitesReference, reference: PopulationReferenceCompatibility,
    runtime_image: str, limits: PopulationQueryLimits,
) -> PopulationAnnotationResult:
    """Run one bounded batch in a new authenticated resource session.

    Multi-batch callers should hold :class:`PopulationResourceSession` for the
    stage so selected large files are hashed once and checked for stable identity
    around every query.
    """
    with PopulationResourceSession(
        joint=joint, all_sites=all_sites, reference=reference,
        runtime_image=runtime_image, limits=limits,
    ) as session:
        return session.query(requests)


def validate_population_annotation(
    result: PopulationAnnotationResult, requests: tuple[PopulationRequest, ...], *,
    joint: GnomadJointReference, all_sites: GnomadAllSitesReference,
    reference: PopulationReferenceCompatibility, runtime_image: str, limits: PopulationQueryLimits,
) -> None:
    """Consumer authentication: rerun queries and verify bytes, not an asserted pass.

    The caller supplies trusted immutable pins independently of the report. Replay
    compares the full report and all enriched evidence, including unrelated fields.
    """
    try:
        if type(result) is not PopulationAnnotationResult:
            raise ValueError
        replay = run_population_annotation(requests, joint=joint, all_sites=all_sites, reference=reference,
                                           runtime_image=runtime_image, limits=limits)
        if result != replay:
            raise ValueError
    except Exception:
        raise PopulationAnnotationError("population annotation replay validation failed") from None


def population_annotation_report_as_json(result: PopulationAnnotationResult) -> str:
    """Explicit controlled serialization; not a validation or declassification."""
    return _json(result.report)
