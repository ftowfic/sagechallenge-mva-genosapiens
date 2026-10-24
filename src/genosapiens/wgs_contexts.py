"""Disk-backed global context closure for bounded WGS admission batches."""

from __future__ import annotations

from contextlib import closing
from dataclasses import InitVar, dataclass, fields, replace
import hashlib
import json
import os
from pathlib import Path
import sqlite3
import sys

from genosapiens.wgs_failures import failure_boundary, static_failure, run_cleanup, WgsFailure
import tempfile
import time
from types import MappingProxyType
from typing import Callable, Iterable, Iterator, Mapping

from genosapiens.wgs_context_payloads import Encoder, Payloads, expand, ContextPayloadError
from genosapiens.wgs_admission_projection import AdmissionProjectionReader, AdmissionProjectionError
from genosapiens.bounded_wgs_config import BoundedWgsExecutionConfig
from genosapiens.wgs_context_runtime import ContextWriterProfile
from genosapiens.candidate_admission import (
    AdmissionDecision,
    CandidateContextSeed,
    ShardedCandidateAdmission,
    _bytes,
    _hash,
    _requirement_key,
    _record_bytes,
    _unique,
)
from genosapiens.evidence_model import AlleleReference, _typed
from genosapiens.evidence_model import (
    AlleleEvidence,
    DiseaseGeneEvidence,
    IdentityStage,
    PublicSourcedEvidence,
    SourcedRecord,
)
from genosapiens.clinical_annotation import ClinicalAnnotationResult
from genosapiens.wgs_failures import WgsFailure
from genosapiens.evidence_serialization import _plain
from genosapiens.genotype_hypotheses import (
    ConstructedHypothesis,
    GenotypeHypothesesError,
    HypothesisPolicy,
    PartnerSearchAssessment,
    evaluate_complete_context,
    hypothesis_policy_bytes,
    _pairs,
    _requires_pairs,
)
from genosapiens.genotype_properties import (
    SuppliedGenotypeProperties,
    assess_ledger_genotypes,
)
from genosapiens.population_annotation import (
    PopulationAnnotationResult,
    PopulationRequest,
    PopulationResourceSession,
)
from genosapiens.population_policy import (
    PopulationPolicyResult,
    PreparedPopulationPolicy,
    evaluate_population_policy,
)
from genosapiens.partner_relevance import (
    ContextRelevanceAssessment,
    SourceContextRelevance,
)
from genosapiens.wgs_admission import (
    LedgerAdmissionInput,
)


CONTEXT_INDEX_SCHEMA = "genosapiens.wgs-context-index"
_CLOSED_CONTEXT_SEAL = object()
CLOSED_CONTEXT_VERSION = 1
CLOSED_MEMBER_PAGE_BYTES = 16 * 1024**2
CLOSED_RECORD_BYTES = 64 * 1024**2


class WgsContextError(RuntimeError):
    """Static global-context closure failure."""

    operation = "context_annotate"


_CONTEXT_BYTES_SEAL = object()


_DECISION_FIXED_FIELDS = tuple(
    field.name for field in fields(AdmissionDecision)
    if field.name not in ("evidence_refs", "source_refs")
)


def _decision_reference_bytes(values):
    return [(_record_bytes(value), value) for value in values]


def _reference_list_extra_bytes(items):
    # Each key contains one trailing newline, absent inside a JSON list.
    return sum(len(key)-1 for key, _ in items) + max(0, len(items)-1)


def _merge_reference_bytes(target, items, total):
    for key, value in items:
        if key not in target:
            total += len(key)-1 + bool(target)
        target[key] = value  # _unique retains the last equal canonical key.
    return total


@dataclass(frozen=True, slots=True, repr=False)
class _ContextByteOwner:
    seed: object
    members: tuple
    decisions: tuple
    sources: tuple
    total: int
    seal: object


@dataclass(frozen=True, slots=True, repr=False)
class CompleteContext:
    ordinal: int
    seed: CandidateContextSeed
    members: tuple[LedgerAdmissionInput, ...]
    canonical_member_bytes: int
    admission_decisions: tuple[AdmissionDecision, ...] = ()
    public_sources: tuple[PublicSourcedEvidence, ...] = ()
    _byte_owner: InitVar[object] = None

    def __post_init__(self, _byte_owner) -> None:
        if _byte_owner is None:
            total = len(_bytes(self.seed)) + sum(len(v.as_bytes()) for v in self.members) \
                + sum(len(_bytes(v)) for v in self.admission_decisions) \
                + sum(len(_bytes(v)) for v in self.public_sources)
        elif (type(_byte_owner) is _ContextByteOwner and _byte_owner.seal is _CONTEXT_BYTES_SEAL
                and _byte_owner.seed is self.seed and _byte_owner.members is self.members
                and _byte_owner.decisions is self.admission_decisions and _byte_owner.sources is self.public_sources):
            total = _byte_owner.total
        else:
            raise WgsContextError('complete context byte owner is invalid')
        references = tuple(
            AlleleReference(value.allele.allele_id, value.allele.native_identity_id)
            for value in self.members
        )
        if (
            not isinstance(self.ordinal, int)
            or isinstance(self.ordinal, bool)
            or self.ordinal <= 0
            or type(self.seed) is not CandidateContextSeed
            or type(self.members) is not tuple
            or any(type(value) is not LedgerAdmissionInput for value in self.members)
            or type(self.admission_decisions) is not tuple
            or any(
                type(value) is not AdmissionDecision
                for value in self.admission_decisions
            )
            or type(self.public_sources) is not tuple
            or any(
                type(value) is not PublicSourcedEvidence
                for value in self.public_sources
            )
            or tuple(value.decision_id for value in self.admission_decisions)
            != tuple(sorted({
                value.decision_id for value in self.admission_decisions
            }))
            or tuple(value.evidence_id for value in self.public_sources)
            != tuple(sorted({value.evidence_id for value in self.public_sources}))
            or self.seed.allele_refs != references
            or self.canonical_member_bytes != total
        ):
            raise WgsContextError("complete context is invalid")

    @property
    def context_id(self) -> str:
        return self.seed.context.context_id


@dataclass(frozen=True, slots=True, repr=False)
class ClosedPublicRecord:
    evidence_id: str
    public_manifest_sha256: str
    record: SourcedRecord | None


class ClosedContext:
    """Full global membership, with bounded access to its retained evidence.

    Valid only while iter_closed_contexts owns the sealed admission provider and
    global index. This operation-local handle is never a serialized trust flag.
    Total logical input bytes may exceed the old rich-context resident bound.
    """
    def __init__(self, ordinal, seed, connection, reader, config, *, scratch_directory, payloads=None, _seal=None):
        if _seal is not _CLOSED_CONTEXT_SEAL:
            raise WgsContextError("closed context requires global source closure")
        self.ordinal, self.seed = ordinal, seed
        self._connection, self._reader, self._config = connection, reader, config
        self._scratch_directory = Path(scratch_directory)
        self._payloads = payloads
        self.context_id = seed.context.context_id
        self._baseline = connection.execute(
            'SELECT 1 FROM baseline_seeds WHERE context_id=?', (self.context_id,)).fetchone() is not None
        count, total = connection.execute(
            "SELECT count(*),coalesce(sum(i.input_size),0) FROM members m JOIN inputs i USING(allele_id) WHERE m.context_id=?",
            (self.context_id,)).fetchone()
        if count > config.contexts.max_member_identities:
            raise WgsContextError("complete context exceeds member capacity")
        self.member_count, self.input_bytes = count, total
        # References are bounded metadata (at most 100,000), not decoded rich
        # members. Preserve the original retained-ordinal membership order.
        references = tuple(AlleleReference(key, native) for key, native in connection.execute(
            "SELECT i.allele_id,i.native_identity_id FROM members m JOIN inputs i USING(allele_id) "
            "WHERE m.context_id=? ORDER BY i.retained_ordinal", (self.context_id,)))
        if set(references) != set(seed.allele_refs):
            raise WgsContextError("closed context membership disagrees with merged seeds")
        self.seed = replace(seed, allele_refs=references)
        if len(_bytes(self.seed)) > config.contexts.max_context_bytes:
            raise WgsContextError("closed context metadata exceeds capacity")

    def iter_members(self, *, allele_order=False):
        """Read complete membership in source order or lexical pair order."""
        for member, _ in self.iter_sized_members(allele_order=allele_order):
            yield member

    def iter_sized_members(self, *, allele_order=False):
        """Return members with the canonical input sizes verified by their owner.

        get_inputs authenticates each requested reference, byte length, digest
        and strict admission codec (or reuses its exact private owned result).
        Consumers can size a batch without serializing that same input again.
        """
        def fetch(page):
            values = self._reader.get_inputs(tuple(page))
            try:
                for ref, member in zip(page, values, strict=True):
                    if (member.retained_ordinal, member.identity_id,
                            None if member.allele is None else member.allele.allele_id) != ref[1:4]:
                        raise WgsContextError("closed context member identity changed")
                    yield member, ref[5]
            except ValueError:
                raise WgsContextError("closed context member coverage changed") from None

        ordering = "i.allele_id" if allele_order else "i.retained_ordinal"
        try:
            cursor = self._connection.execute(
                "SELECT i.projection_binding_sha256,i.retained_ordinal,i.identity_id,i.allele_id,i.input_sha256,i.input_size "
                "FROM members m JOIN inputs i USING(allele_id) WHERE m.context_id=? ORDER BY " + ordering,
                (self.context_id,))
            page, size, count = [], 0, 0
            for row in cursor:
                if row[5] > CLOSED_MEMBER_PAGE_BYTES:
                    raise WgsContextError("closed context member exceeds page capacity")
                if page and (size + row[5] > CLOSED_MEMBER_PAGE_BYTES
                             or len(page) >= self._config.sqlite.cursor_batch_rows):
                    yield from fetch(page)
                    page, size = [], 0
                page.append(tuple(row))
                size += row[5]
                count += 1
            if page:
                yield from fetch(page)
            if count != self.member_count:
                raise WgsContextError("closed context member coverage changed")
        except (sqlite3.Error, AdmissionProjectionError):
            raise WgsContextError("closed context owner is unavailable or changed") from None

    def iter_decisions(self):
        """Validate every row, then construct each merged decision once."""
        try:
            for key in self.seed.decision_ids:
                first = None
                refs = sources = None
                cursor = self._connection.execute(
                    "SELECT d.payload_sha256 FROM decisions d WHERE d.decision_id=? AND "
                    "(? OR EXISTS(SELECT 1 FROM seeds s WHERE s.context_id=? AND s.batch_ordinal=d.batch_ordinal)) "
                    "ORDER BY d.batch_ordinal,d.payload_sha256",
                    (key, self._baseline, self.context_id))
                for (digest,) in cursor:
                    content = self._payloads.expanded(digest)
                    if len(content) > CLOSED_RECORD_BYTES:
                        raise WgsContextError('closed context decision exceeds record capacity')
                    value = _typed(json.loads(content), AdmissionDecision, decode=True)
                    if not self._connection.execute('SELECT 1 FROM members WHERE context_id=? AND allele_id=?',
                            (self.context_id, value.subject_id)).fetchone():
                        continue
                    if first is None:
                        first = value
                        first_size = len(_bytes(first))
                        size = first_size
                    else:
                        if any(getattr(first, name) != getattr(value, name) for name in _DECISION_FIXED_FIELDS):
                            raise WgsContextError('cross-shard admission decision conflicts')
                        if refs is None:
                            old_refs = _decision_reference_bytes(first.evidence_refs)
                            old_sources = _decision_reference_bytes(first.source_refs)
                            fixed_size = (first_size - _reference_list_extra_bytes(old_refs)
                                          - _reference_list_extra_bytes(old_sources))
                            refs, sources = {}, {}
                            refs_size = _merge_reference_bytes(refs, old_refs, 0)
                            sources_size = _merge_reference_bytes(sources, old_sources, 0)
                        refs_size = _merge_reference_bytes(refs, _decision_reference_bytes(value.evidence_refs), refs_size)
                        sources_size = _merge_reference_bytes(sources, _decision_reference_bytes(value.source_refs), sources_size)
                        size = fixed_size+refs_size+sources_size
                    # Preserve the capacity gate after EVERY accepted intermediate row.
                    if size > CLOSED_RECORD_BYTES:
                        raise WgsContextError('closed context decision exceeds record capacity')
                if first is not None:
                    if refs is not None:
                        first = replace(first, evidence_refs=tuple(refs[k] for k in sorted(refs)),
                                        source_refs=tuple(sources[k] for k in sorted(sources)))
                        # Check the incremental accounting against the final canonical record.
                        actual_size = len(_bytes(first))
                        if actual_size != size:
                            raise WgsContextError('merged decision canonical size changed')
                        if actual_size > CLOSED_RECORD_BYTES:
                            raise WgsContextError('closed context decision exceeds record capacity')
                    yield first
        except sqlite3.Error:
            raise WgsContextError('closed context owner is unavailable or changed') from None

    def iter_public_records(self):
        """Close public-owner conflicts on disk before exposing individual rows."""
        connection = self._connection
        table = ("context_public_baseline" if self._baseline else
                 "context_public_" + hashlib.sha256(self.context_id.encode()).hexdigest())
        temporary = store = cursor = None
        if self._baseline and getattr(self._payloads, 'public_closed', False):
            yield from self._public_table_rows(connection, table)
            return
        try:
            # A SQLite TEMP table can spill onto the runtime's small tmpfs.
            # Keep this context's deduplication on explicit disk scratch, in a
            # separate connection so suspended source cursors remain valid.
            if self._baseline:
                store = connection
            else:
                temporary = tempfile.TemporaryDirectory(prefix=".context-public-",dir=self._scratch_directory)
                store = _open_index(Path(temporary.name)/"public.sqlite",self._config)
            store.execute(f"CREATE TABLE {table} (evidence_id TEXT, manifest TEXT, record_id TEXT, content BLOB, "
                "PRIMARY KEY(evidence_id,record_id)) WITHOUT ROWID")
            store.execute("CREATE TABLE fragments(evidence_id TEXT,digest TEXT,PRIMARY KEY(evidence_id,digest)) WITHOUT ROWID")
            def fragments():
                # Walk both existing primary-key indexes; a GROUP BY/join
                # sort could spill all source fragments onto SQLite temp.
                batch = 0
                while True:
                    row = (connection.execute('SELECT batch_ordinal FROM batches WHERE batch_ordinal>? '
                        'ORDER BY batch_ordinal LIMIT 1', (batch,)).fetchone() if self._baseline else
                        connection.execute("SELECT batch_ordinal FROM seeds WHERE context_id=? AND batch_ordinal>? "
                        "ORDER BY batch_ordinal LIMIT 1",(self.context_id,batch)).fetchone())
                    if row is None:
                        return
                    batch, = row
                    previous = ''
                    while True:
                        size = self._config.sqlite.cursor_batch_rows
                        with closing(connection.execute("SELECT evidence_id,source_sha256 FROM public_sources "
                                "WHERE batch_ordinal=? AND evidence_id>? ORDER BY evidence_id LIMIT ?",
                                (batch,previous,size))) as page:
                            rows = page.fetchmany(size)
                        if not rows:
                            break
                        for evidence_id,digest in rows:
                            previous = evidence_id
                            if store.execute("INSERT OR IGNORE INTO fragments VALUES (?,?)",(evidence_id,digest)).rowcount:
                                value = connection.execute("SELECT canonical_source FROM public_payloads WHERE sha256=?",(digest,)).fetchone()
                                if value is None:
                                    raise WgsContextError("closed public source fragment is missing")
                                yield value
            cursor = fragments()
            for (content,) in cursor:
                content = expand(content, self._payloads.get)
                if len(content) > CLOSED_RECORD_BYTES:
                    raise WgsContextError("closed public source fragment exceeds record capacity")
                source = _typed(json.loads(content), PublicSourcedEvidence, decode=True)
                prior = store.execute(f"SELECT manifest FROM {table} WHERE evidence_id=? LIMIT 1", (source.evidence_id,)).fetchone()
                if prior is not None and prior[0] != source.public_manifest_sha256:
                    raise WgsContextError("cross-shard public source evidence conflicts")
                # Preserve empty owners as an explicit row as well.
                if source.records:
                    store.execute(f"DELETE FROM {table} WHERE evidence_id=? AND record_id=''", (source.evidence_id,))
                records = source.records.items() if source.records else (() if prior is not None else (("", None),))
                for key, record in records:
                    content = _bytes(record)
                    if len(content) > CLOSED_RECORD_BYTES:
                        raise WgsContextError("closed public evidence exceeds record capacity")
                    prior = store.execute(f"SELECT content FROM {table} WHERE evidence_id=? AND record_id=?",
                        (source.evidence_id, key)).fetchone()
                    if prior is None:
                        store.execute(f"INSERT INTO {table} VALUES (?,?,?,?)",
                            (source.evidence_id, source.public_manifest_sha256, key, content))
                    elif prior[0] != content:
                        raise WgsContextError("cross-shard public source evidence conflicts")
                    if self._baseline:
                        self._payloads.transactions.add(len(content))
            cursor.close()
            cursor = None
            store.execute('DROP TABLE fragments')
            if self._baseline:
                self._payloads.transactions.commit(final=True)
                self._payloads.public_closed = True
            yield from self._public_table_rows(store, table)
        except sqlite3.Error:
            raise WgsContextError("closed context owner is unavailable or changed") from None
        finally:
            if cursor is not None:
                cursor.close()
            if store is not None and store is not connection:
                store.close()
            if temporary is not None:
                temporary.cleanup()

    @staticmethod
    def _public_table_rows(store, table):
        for evidence_id, manifest, key, content in store.execute(
                f"SELECT evidence_id,manifest,record_id,content FROM {table} ORDER BY evidence_id,record_id"):
            yield ClosedPublicRecord(evidence_id, manifest,
                None if key == "" else _typed(json.loads(content), SourcedRecord, decode=True))

    @property
    def public_owner(self):
        return self._payloads if self._baseline else None


@dataclass(frozen=True, slots=True, repr=False)
class ContextHypothesisOutput:
    """One complete context's bounded genotype and hypothesis output."""

    ordinal: int
    context_id: str
    alleles: tuple[AlleleEvidence, ...]
    genotypes: SuppliedGenotypeProperties
    search: PartnerSearchAssessment
    hypotheses: tuple[ConstructedHypothesis, ...]

    def __post_init__(self) -> None:
        if (
            not isinstance(self.ordinal, int)
            or isinstance(self.ordinal, bool)
            or self.ordinal <= 0
            or not isinstance(self.context_id, str)
            or not self.context_id
            or type(self.alleles) is not tuple
            or any(type(value) is not AlleleEvidence for value in self.alleles)
            or len({value.allele_id for value in self.alleles}) != len(self.alleles)
            or type(self.genotypes) is not SuppliedGenotypeProperties
            or type(self.search) is not PartnerSearchAssessment
            or self.search.context_id != self.context_id
            or type(self.hypotheses) is not tuple
            or any(type(value) is not ConstructedHypothesis for value in self.hypotheses)
            or any(value.admission_context_id != self.context_id for value in self.hypotheses)
            or {value.allele_id for value in self.alleles}
            != {value.allele_id for value in self.genotypes.allele_properties}
        ):
            raise WgsContextError("context hypothesis output is invalid")

    def as_bytes(self) -> bytes:
        return _bytes({
            "schema": "genosapiens.wgs-context-hypotheses",
            "version": 1,
            "classification": "controlled",
            "llm_policy": "no-llm",
            "result": self,
        })


@dataclass(frozen=True, slots=True, repr=False)
class PopulationAssessedContext:
    """One context after unchanged downstream annotation and population policy."""

    context: ContextHypothesisOutput
    population_annotation: PopulationAnnotationResult
    clinical_annotation: ClinicalAnnotationResult | None
    assessments: tuple[PopulationPolicyResult, ...]
    early_protected_allele_ids: tuple[str, ...]
    technical_obligation_allele_ids: tuple[str, ...]

    def __post_init__(self) -> None:
        allele_ids = {value.allele_id for value in self.context.alleles}
        hypotheses = {
            value.hypothesis.hypothesis_id: value.hypothesis
            for value in self.context.hypotheses
        }
        if (
            type(self.context) is not ContextHypothesisOutput
            or type(self.population_annotation) is not PopulationAnnotationResult
            or self.clinical_annotation is not None
            and type(self.clinical_annotation) is not ClinicalAnnotationResult
            or type(self.assessments) is not tuple
            or any(type(value) is not PopulationPolicyResult for value in self.assessments)
            or len({value.hypothesis.hypothesis_id for value in self.assessments})
            != len(self.assessments)
            or set(value.hypothesis.hypothesis_id for value in self.assessments)
            != set(hypotheses)
            or any(
                value.hypothesis != hypotheses[value.hypothesis.hypothesis_id]
                for value in self.assessments
            )
            or {value.allele_id for value in self.population_annotation.alleles}
            != allele_ids
            or self.clinical_annotation is not None
            and {value.allele_id for value in self.clinical_annotation.alleles}
            != allele_ids
            or self.early_protected_allele_ids
            != tuple(sorted(set(self.early_protected_allele_ids)))
            or self.technical_obligation_allele_ids
            != tuple(sorted(set(self.technical_obligation_allele_ids)))
            or not set(self.early_protected_allele_ids) <= allele_ids
            or not set(self.technical_obligation_allele_ids) <= allele_ids
        ):
            raise WgsContextError("population-assessed context is invalid")

    def as_bytes(self) -> bytes:
        return _bytes({
            "schema": "genosapiens.wgs-population-assessed-context",
            "version": 1,
            "classification": "controlled",
            "llm_policy": "no-llm",
            "result": self,
        })


def _open_index(path: Path, config: BoundedWgsExecutionConfig, *, context_wal=False,
                writer_profile=None) -> sqlite3.Connection:
    if writer_profile is not None and (not context_wal or type(writer_profile) is not ContextWriterProfile):
        raise WgsContextError('profile requires the disposable context writer')
    profile = ContextWriterProfile() if writer_profile is None else writer_profile
    descriptor = -1
    connection = None
    created = False
    try:
        descriptor = os.open(
            path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600,
        )
        created = True
        os.fchmod(descriptor, 0o600)
        os.close(descriptor)
        descriptor = -1
        connection = sqlite3.connect(path)
        connection.enable_load_extension(False)
        connection.execute("PRAGMA foreign_keys=ON")
        connection.execute(
            f"PRAGMA cache_size=-{max(1, (profile.cache_bytes if context_wal else config.sqlite.cache_bytes) // 1024)}"
        )
        connection.execute("PRAGMA temp_store=FILE")
        connection.execute("PRAGMA mmap_size=0")
        connection.execute("PRAGMA trusted_schema=OFF")
        mode = 'WAL' if context_wal else 'DELETE'
        if connection.execute('PRAGMA journal_mode='+mode).fetchone()[0].upper() != mode:
            raise WgsContextError('context journal mode is unavailable')
        connection.execute('PRAGMA synchronous='+ (profile.synchronous if context_wal else 'FULL'))
        if context_wal:
            page_size = connection.execute('PRAGMA page_size').fetchone()[0]
            pages = max(1, profile.autocheckpoint_bytes//page_size)
            connection.execute(f'PRAGMA wal_autocheckpoint={pages}')
            expected = {'cache_size': -profile.cache_bytes//1024,
                        'synchronous': 2 if profile.synchronous == 'FULL' else 1,
                        'wal_autocheckpoint': pages, 'foreign_keys': 1,
                        'temp_store': 1, 'mmap_size': 0, 'trusted_schema': 0}
            if any(connection.execute('PRAGMA '+key).fetchone()[0] != value
                   for key, value in expected.items()):
                raise WgsContextError('effective context writer profile differs')
        return connection
    except (OSError, sqlite3.Error, WgsContextError):
        if descriptor >= 0:
            os.close(descriptor)
        if connection is not None:
            connection.close()
        if created:
            for owned in (path, path.with_name(path.name+'-wal'), path.with_name(path.name+'-shm')):
                try:
                    owned.unlink()
                except FileNotFoundError:
                    pass
        raise WgsContextError("cannot create global context index") from None


class _ContextTransactions:
    """Single-connection disposable WAL; never a durable restart boundary.

    Thresholds trigger commits, with at most one bounded 64-MiB part overshoot.
    WAL checkpoints never run on another connection.
    """
    def __init__(self, connection, path, profile=None, *, observation=None):
        self.connection, self.path = connection, Path(path)
        self.profile = ContextWriterProfile() if profile is None else profile
        if type(self.profile) is not ContextWriterProfile:
            raise WgsContextError('invalid context transaction profile')
        self.bytes = 0
        self.changes = connection.total_changes
        self.commits = self.checkpoints = 0
        self.commit_seconds = self.checkpoint_seconds = 0.0
        self.max_supplied_bytes = self.max_changes = self.max_trigger_overshoot_bytes = 0
        self.max_wal_overshoot_bytes = 0
        self.allocation = {'db': 0, 'wal': 0, 'shm': 0}
        self.peak_allocation = dict(self.allocation)
        self.observation=observation
        self.supplied_bytes=0
        self.part_write_seconds=0.0

    def observe(self):
        if self.observation is not None:
            self.observation(dict(profile_sha256=self.profile.sha256,
                supplied_bytes=self.supplied_bytes,total_changes=self.connection.total_changes,
                pending_supplied_bytes=self.bytes,pending_changes=self.connection.total_changes-self.changes,
                commits=self.commits,commit_seconds=self.commit_seconds,
                checkpoints=self.checkpoints,checkpoint_seconds=self.checkpoint_seconds,
                part_write_seconds_inclusive=self.part_write_seconds,
                max_trigger_overshoot_bytes=self.max_trigger_overshoot_bytes,
                sidecar_bytes=dict(self.allocation),peak_sidecar_bytes=dict(self.peak_allocation),
                max_wal_bytes=self.profile.max_wal_bytes))

    def check_allocation(self):
        for key, suffix in (('db', ''), ('wal', '-wal'), ('shm', '-shm')):
            path = self.path.with_name(self.path.name+suffix)
            try:
                size = path.stat().st_size
            except FileNotFoundError:
                size = 0
            self.allocation[key] = size
            self.peak_allocation[key] = max(self.peak_allocation[key], size)
        overshoot = max(0, self.allocation['wal']-self.profile.max_wal_bytes)
        self.max_wal_overshoot_bytes = max(self.max_wal_overshoot_bytes, overshoot)
        if overshoot:
            raise WgsContextError('context WAL allocation exceeds capacity')

    def add(self, size):
        if type(size) is not int or not 0 <= size <= 64*1024**2:
            raise WgsContextError('context transaction unit exceeds capacity')
        self.bytes += size
        self.supplied_bytes += size
        changes = self.connection.total_changes-self.changes
        self.max_supplied_bytes = max(self.max_supplied_bytes, self.bytes)
        self.max_changes = max(self.max_changes, changes)
        self.max_trigger_overshoot_bytes = max(self.max_trigger_overshoot_bytes,
                                               max(0, self.bytes-self.profile.transaction_bytes))
        self.check_allocation()
        if self.bytes >= self.profile.transaction_bytes or changes >= self.profile.transaction_changes:
            self.commit()
        else:
            self.observe()

    def commit(self, *, final=False):
        self.check_allocation()
        started = time.monotonic()
        self.connection.commit()
        self.commit_seconds += time.monotonic()-started
        self.commits += 1
        self.bytes = 0
        self.changes = self.connection.total_changes
        # FULL is restored outside any transaction before exposing contexts.
        if final:
            self.connection.execute('PRAGMA synchronous=FULL')
            if self.connection.execute('PRAGMA synchronous').fetchone() != (2,):
                raise WgsContextError('context closure could not restore FULL')
        self.check_allocation()
        if final or self.allocation['wal'] > self.profile.truncate_bytes:
            started = time.monotonic()
            if self.connection.execute('PRAGMA wal_checkpoint(TRUNCATE)').fetchone() != (0, 0, 0):
                raise WgsContextError('context checkpoint did not complete')
            self.checkpoint_seconds += time.monotonic()-started
            self.checkpoints += 1
            self.check_allocation()
        self.observe()


def _context_keys(connection, size):
    # Close each bounded read cursor before exposing a context that can build
    # its shared public table and checkpoint on this same connection.
    previous = ''
    while True:
        with closing(connection.execute('SELECT context_id FROM contexts WHERE context_id>? '
                'AND EXISTS(SELECT 1 FROM members WHERE members.context_id=contexts.context_id) '
                'ORDER BY context_id LIMIT ?', (previous, size))) as cursor:
            rows = cursor.fetchmany(size)
        if not rows:
            return
        for (previous,) in rows:
            yield previous


def _rows(cursor, size: int):
    while True:
        values = cursor.fetchmany(size)
        if not values:
            return
        yield from values


def _grouping_bytes(seed: CandidateContextSeed) -> bytes:
    context = seed.context
    requirements = {_requirement_key(value) for value in context.requirements}
    if len(requirements) != 1:
        raise WgsContextError("context requirement identity is invalid")
    requirement = next(iter(requirements))
    target = context.locus_id or context.gene_id
    expected = "context-" + _hash([
        context.native_disease_id,
        seed.target_namespace,
        target,
        requirement,
    ])
    if context.context_id != expected:
        raise WgsContextError("context grouping identity is invalid")
    return _bytes({
        "schema": CONTEXT_INDEX_SCHEMA,
        "context_id": context.context_id,
        "disease_id": context.disease_id,
        "native_disease_id": context.native_disease_id,
        "target_namespace": seed.target_namespace,
        "gene_id": context.gene_id,
        "locus_id": context.locus_id,
        "inheritance": context.inheritance,
        "model": context.model,
        "requirement": requirement,
        "locus_bindings": context.locus_bindings,
    })


def _load_seed(content: bytes) -> CandidateContextSeed:
    try:
        data = json.loads(content)
        value = _typed(data, CandidateContextSeed, decode=True)
        if _bytes(value) != content:
            raise ValueError
        return value
    except (ValueError, TypeError, RuntimeError, UnicodeError):
        raise WgsContextError("candidate context seed is invalid") from None


def _merge_seeds(values: tuple[CandidateContextSeed, ...]) -> CandidateContextSeed:
    if not values:
        raise WgsContextError("context seed coverage is empty")
    first = values[0]
    if any(
        value.target_namespace != first.target_namespace
        or value.phenotype_support != first.phenotype_support
        or value.support_ledger != first.support_ledger
        or _grouping_bytes(value) != _grouping_bytes(first)
        for value in values[1:]
    ):
        raise WgsContextError("cross-shard context claims conflict")
    context = replace(
        first.context,
        requirements=_unique(
            requirement for value in values
            for requirement in value.context.requirements
        ),
        evidence_refs=_unique(
            reference for value in values
            for reference in value.context.evidence_refs
        ),
    )
    return replace(
        first,
        context=context,
        routes=tuple(sorted({route for value in values for route in value.routes})),
        allele_refs=_unique(
            reference for value in values for reference in value.allele_refs
        ),
        phenotype_proof_ids=tuple(sorted({
            proof for value in values for proof in value.phenotype_proof_ids
        })),
        decision_ids=tuple(sorted({
            decision for value in values for decision in value.decision_ids
        })),
    )


def _merge_seed_stream(values, limit):
    """Accumulate unique additions once, retaining a single shared base."""
    first = next(values, None)
    if first is None:
        raise WgsContextError('context seed coverage is empty')
    grouping = _grouping_bytes(first)
    maps = {key: {} for key in ('requirements', 'evidence_refs', 'allele_refs')}
    sets = {key: set() for key in ('routes', 'phenotype_proof_ids', 'decision_ids')}
    resident = len(_bytes(first))
    def add(value):
        nonlocal resident
        if (value.target_namespace != first.target_namespace
                or value.phenotype_support != first.phenotype_support
                or value.support_ledger != first.support_ledger or _grouping_bytes(value) != grouping):
            raise WgsContextError('cross-shard context claims conflict')
        for key, target in maps.items():
            owner = value.context if key in ('requirements', 'evidence_refs') else value
            for item in getattr(owner, key):
                content = _bytes(item)
                if content not in target:
                    target[content] = item
                    resident += len(content)
        for key, target in sets.items():
            for item in getattr(value, key):
                if item not in target:
                    target.add(item)
                    resident += len(_bytes(item))
        if resident > limit * 2:
            # Includes the original base plus its initial entries; the exact
            # complete logical size is checked after the single final merge.
            raise WgsContextError('complete context exceeds capacity')
    add(first)
    for value in values:
        add(value)
    merged = {key: tuple(value for _, value in sorted(target.items())) for key, target in maps.items()}
    context = replace(first.context, requirements=merged.pop('requirements'), evidence_refs=merged.pop('evidence_refs'))
    result = replace(first, context=context, **merged, **{key: tuple(sorted(value)) for key,value in sets.items()})
    if len(_bytes(result)) > limit:
        raise WgsContextError('complete context exceeds capacity')
    return result


def _bounded_relevance(value: CompleteContext) -> ContextRelevanceAssessment:
    binding = {
        "schema": "genosapiens.wgs-context-relevance",
        "version": 1,
        "context_sha256": _hash(value.seed.context),
        "seed_sha256": _hash(value.seed),
        "scope": "globally-closed-retained-context-membership",
    }
    binding_sha256 = _hash(binding)
    entries = {}
    for member in value.members:
        allele = member.allele
        if allele is None:
            raise WgsContextError("complete context member has no allele")
        for call in allele.calls or ():
            fields = call.source.source_fields or {}
            source_id = fields.get("source_allele_id")
            record_id = fields.get(
                "source_record_id", call.source.source_record_id
            )
            if not isinstance(source_id, str) or not isinstance(record_id, str):
                raise WgsContextError(
                    "complete context source membership is invalid"
                )
            previous = entries.get(source_id)
            derived = (allele.allele_id,)
            if previous is not None:
                if previous.source_record_id != record_id:
                    raise WgsContextError(
                        "complete context source membership conflicts"
                    )
                derived = tuple(sorted({
                    *previous.derived_allele_ids, allele.allele_id,
                }))
            entries[source_id] = SourceContextRelevance(
                record_id, source_id, derived, "derived", (), "relevant", (),
                (call.source,), value.seed.context.evidence_refs, (), (), (),
                (
                    "globally-closed-admitted-context-member",
                    "source-target-association-replayed-upstream",
                ),
                binding_sha256,
            )
    return ContextRelevanceAssessment(
        value.seed.context,
        tuple(entries[key] for key in sorted(entries)),
        binding,
    )


def annotate_complete_context(
    complete_context: CompleteContext,
    *,
    population_session: PopulationResourceSession,
    execution_config: BoundedWgsExecutionConfig,
    clinical_annotator: Callable[
        [tuple[AlleleEvidence, ...]], ClinicalAnnotationResult
    ] | None = None,
    query_cache=None,
) -> tuple[
    CompleteContext,
    PopulationAnnotationResult,
    ClinicalAnnotationResult | None,
]:
    """Attach downstream population and clinical evidence after global closure."""

    operation = "context_annotation_inputs"
    try:
        if (
            type(complete_context) is not CompleteContext
            or not complete_context.members
            or type(population_session) is not PopulationResourceSession
            or type(execution_config) is not BoundedWgsExecutionConfig
            or clinical_annotator is not None and not callable(clinical_annotator)
        ):
            raise ValueError
        requests = []
        for member in complete_context.members:
            allele = member.allele
            if allele is None or allele.population or allele.clinical:
                raise ValueError
            normalized = tuple(
                value.identity_id for value in allele.identities
                if value.stage is IdentityStage.NORMALIZED
            )
            if len(normalized) > 1:
                raise ValueError
            requests.append(PopulationRequest(
                allele, normalized[0] if normalized else None,
            ))
        batch_size = min(execution_config.public_queries.batch_rows, population_session.limits.max_alleles)
        population_rows, population_reports = [], []
        clinical_rows, clinical_coverage, clinical_bindings = [], [], []
        clinical_stores = {}
        resident_bytes = complete_context.canonical_member_bytes
        def retain(value):
            nonlocal resident_bytes
            resident_bytes += len(_bytes(value))
            if resident_bytes > execution_config.contexts.max_context_bytes:
                raise WgsContextError("annotated complete context exceeds capacity")
        for offset in range(0, len(requests), batch_size):
            batch_requests = tuple(requests[offset:offset + batch_size])
            operation = "context_population_query"
            population_batch = (population_session.query(batch_requests) if query_cache is None
                else query_cache.population_results(batch_requests,session=population_session))
            operation = "context_population_bind"
            expected = {value.allele.allele_id: value.allele for value in batch_requests}
            if len(population_batch.alleles) != len(expected) or {value.allele_id for value in population_batch.alleles} != set(expected):
                raise ValueError
            for allele in population_batch.alleles:
                if replace(allele, population=expected[allele.allele_id].population) != expected[allele.allele_id]:
                    raise ValueError
            retain(population_batch)
            population_rows.extend(population_batch.alleles)
            population_reports.append(population_batch.report)
            if clinical_annotator is not None:
                operation = "context_clinical_query"
                clinical_batch = clinical_annotator(population_batch.alleles)
                operation = "context_clinical_bind"
                expected_clinical = {value.allele_id: value for value in population_batch.alleles}
                if len(clinical_batch.alleles) != len(expected_clinical) or {
                    value.allele_id for value in clinical_batch.alleles
                } != set(expected_clinical) or any(
                    replace(value, clinical=expected_clinical[value.allele_id].clinical) != expected_clinical[value.allele_id]
                    for value in clinical_batch.alleles
                ):
                    raise ValueError
                retain(clinical_batch)
                clinical_rows.extend(clinical_batch.alleles)
                clinical_coverage.extend(clinical_batch.coverage)
                clinical_bindings.append(clinical_batch.bindings)
                for store in clinical_batch.public_sourced_evidence:
                    key = (store.evidence_id, store.public_manifest_sha256)
                    records = clinical_stores.setdefault(key, {})
                    for record_id, record in store.records.items():
                        if records.setdefault(record_id, record) != record:
                            raise ValueError
        operation = "context_annotation_finalize"
        def batch_binding(role, reports):
            fields = {"scope": "complete-context", "batches": tuple(_plain(value) for value in reports)}
            digest = hashlib.sha256(_bytes(fields)).hexdigest()
            return SourcedRecord("genosapiens.wgs-context-" + role, "1", digest, digest, None, fields)
        population_result = PopulationAnnotationResult(tuple(population_rows), batch_binding("population", population_reports))
        original = {
            value.allele.allele_id: value.allele
            for value in complete_context.members
        }
        populated = {
            value.allele_id: value for value in population_result.alleles
        }
        if set(populated) != set(original) or any(
            replace(populated[key], population=original[key].population)
            != original[key]
            for key in original
        ):
            raise ValueError
        clinical_result = None
        enriched = population_result.alleles
        if clinical_annotator is not None:
            clinical_result = ClinicalAnnotationResult(
                tuple(clinical_rows), tuple(clinical_coverage),
                batch_binding("clinical", clinical_bindings),
                tuple(PublicSourcedEvidence(*key, records) for key, records in sorted(clinical_stores.items())),
            )
            enriched = clinical_result.alleles
        enriched_by_id = {value.allele_id: value for value in enriched}
        if set(enriched_by_id) != set(original) or any(
            replace(
                enriched_by_id[key],
                population=original[key].population,
                clinical=original[key].clinical,
            ) != original[key]
            for key in original
        ):
            raise ValueError
        members = tuple(replace(
            value, allele=enriched_by_id[value.allele.allele_id],
        ) for value in complete_context.members)
        total = len(_bytes(complete_context.seed)) + sum(
            len(value.as_bytes()) for value in members
        ) + sum(len(_bytes(value)) for value in complete_context.admission_decisions) \
            + sum(len(_bytes(value)) for value in complete_context.public_sources)
        if total > execution_config.contexts.max_context_bytes:
            raise WgsContextError(
                "annotated complete context exceeds capacity"
            )
        return (
            replace(complete_context,members=members,canonical_member_bytes=total,
                _byte_owner=_ContextByteOwner(complete_context.seed,members,
                    complete_context.admission_decisions,complete_context.public_sources,total,_CONTEXT_BYTES_SEAL)),
            population_result,
            clinical_result,
        )
    except WgsContextError as error:
        error.operation = operation
        raise
    except WgsFailure:
        raise
    except Exception:
        error = WgsContextError("complete context annotation failed")
        error.operation = operation
        raise error from None


def assess_context_population(
    context: ContextHypothesisOutput,
    *,
    population_annotation: PopulationAnnotationResult,
    clinical_annotation: ClinicalAnnotationResult | None,
    prepared_policy: PreparedPopulationPolicy,
    early_protected_allele_ids: tuple[str, ...] = (),
) -> PopulationAssessedContext:
    """Assess hypotheses with the owning orchestration's curation-validated policy."""

    try:
        if type(context) is not ContextHypothesisOutput:
            raise ValueError
        early = tuple(sorted(set(early_protected_allele_ids)))
        if early != early_protected_allele_ids:
            raise ValueError
        by_id = {value.allele_id: value for value in context.alleles}
        assessments = []
        for value in context.hypotheses:
            selected = tuple(
                by_id[key]
                for key in value.hypothesis.evidence_subject_ids
            )
            assessments.append(evaluate_population_policy(
                value.hypothesis,
                selected,
                prepared_policy=prepared_policy,
            ))
        technical = tuple(sorted(
            value.allele_id for value in context.genotypes.allele_properties
            if value.read_validation_required
        ))
        return PopulationAssessedContext(
            context,
            population_annotation,
            clinical_annotation,
            tuple(assessments),
            early,
            technical,
        )
    except Exception:
        raise WgsContextError(
            "complete context population assessment failed"
        ) from None


def generate_context_hypotheses(
    complete_context: CompleteContext,
    *,
    sample_name: str,
    policy: HypothesisPolicy,
    expected_policy_sha256: str,
    admission_policy_sha256: str,
    disease_gene_evidence: tuple[DiseaseGeneEvidence, ...],
    observation_ids: tuple[str, ...],
    outside_context_partner_uncertainty: bool = True,
) -> ContextHypothesisOutput:
    """Evaluate one globally closed context with the shared scientific kernel."""

    try:
        if (
            type(complete_context) is not CompleteContext
            or not complete_context.members
            or type(policy) is not HypothesisPolicy
            or type(disease_gene_evidence) is not tuple
            or any(
                type(value) is not DiseaseGeneEvidence
                for value in disease_gene_evidence
            )
            or type(observation_ids) is not tuple
            or any(not isinstance(value, str) or not value for value in observation_ids)
            or len(set(observation_ids)) != len(observation_ids)
            or type(outside_context_partner_uncertainty) is not bool
        ):
            raise ValueError
        pin = hypothesis_policy_bytes(
            policy, admission_policy_sha256=admission_policy_sha256,
        )
        if hashlib.sha256(pin).hexdigest() != expected_policy_sha256:
            raise ValueError
        parents = {
            value.parent_source_manifest_sha256
            for value in complete_context.members
        }
        references = {
            mapping.reference_manifest_sha256
            for value in complete_context.members
            for mapping in value.allele.identities
            if mapping.identity_id == value.allele.native_identity_id
        }
        if len(parents) != 1 or len(references) != 1:
            raise ValueError
        source_manifest_sha256, = parents
        reference_manifest_sha256, = references
        alleles = {
            value.allele.allele_id: value.allele
            for value in complete_context.members
        }
        keys = tuple(sorted(alleles))
        if len(keys) != len(complete_context.members):
            raise ValueError
        pair_universe, pair_keys = None, ()
        if _requires_pairs(complete_context.seed.context):
            # Keep the frozen policy's guard, but never publish a partial search.
            pair_universe = _pairs(keys, alleles, policy.maximum_pairs_per_context)
            total, pair_keys = pair_universe
            if total != len(pair_keys):
                raise WgsContextError("complete context pair capacity exceeded")
        genotypes = assess_ledger_genotypes(
            complete_context.members,
            sample_name=sample_name,
            source_manifest_sha256=source_manifest_sha256,
            reference_manifest_sha256=reference_manifest_sha256,
            pairs=pair_keys,
        )
        search, hypotheses = evaluate_complete_context(
            context=complete_context.seed.context,
            relevance=_bounded_relevance(complete_context),
            keys=keys,
            alleles=alleles,
            genotypes=genotypes,
            policy=policy,
            disease_gene_evidence=disease_gene_evidence,
            sample_name=sample_name,
            observation_ids=observation_ids,
            input_sha256=source_manifest_sha256,
            inventory_sha256=genotypes.inventory_sha256,
            policy_sha256=expected_policy_sha256,
            pair_universe=pair_universe,
            source_inventory_status="assessed-domain-complete",
            outside_context_partner_uncertainty=
                outside_context_partner_uncertainty,
        )
        return ContextHypothesisOutput(
            complete_context.ordinal,
            complete_context.context_id,
            tuple(alleles[key] for key in keys),
            genotypes,
            search,
            tuple(sorted(
                hypotheses,
                key=lambda value: value.hypothesis.hypothesis_id,
            )),
        )
    except WgsContextError:
        raise
    except (ValueError, TypeError, RuntimeError, KeyError, AttributeError,
            StopIteration, GenotypeHypothesesError):
        raise WgsContextError("complete context hypothesis replay failed") from None


def iter_context_hypotheses(
    contexts: Iterable[CompleteContext],
    **kwargs,
) -> Iterator[ContextHypothesisOutput]:
    """Yield one bounded complete-context result at a time in cursor order."""

    for value in iter_member_contexts(contexts):
        yield generate_context_hypotheses(value, **kwargs)


def iter_member_contexts(
    contexts: Iterable[CompleteContext],
) -> Iterator[CompleteContext]:
    """Keep globally checked contexts that have an allele to evaluate."""

    previous = 0
    for value in contexts:
        if type(value) is not CompleteContext or value.ordinal <= previous:
            raise WgsContextError("complete context stream order is invalid")
        previous = value.ordinal
        if value.members:
            yield value


def _iter_stored_contexts(connection, payloads, admission_reader, execution_config, *, scratch_path, progress, total_contexts,
                          science_selector=None):
    """One shared deterministic merge owner for fresh and checkpointed inputs."""
    from genosapiens.wgs_runtime import analysis_progress
    roster = None
    if science_selector is not None:
        from genosapiens.wgs_science_selection import ScientificContextSelector
        if type(science_selector) is not ScientificContextSelector:
            raise WgsContextError('scientific selector has a foreign owner')
        roster = science_selector.build(connection, payloads, page_rows=execution_config.sqlite.cursor_batch_rows)
    for ordinal, context_id in enumerate(_context_keys(
        connection, execution_config.sqlite.cursor_batch_rows), start=1):
        if roster is not None and not roster.selected(context_id):
            continue
        seed_cursor = connection.execute(
            'SELECT p.canonical_seed FROM ('
            'SELECT seed_sha256,0 AS ordinal FROM baseline_seeds WHERE context_id=? UNION ALL '
            'SELECT seed_sha256,min(batch_ordinal) AS ordinal FROM seeds WHERE context_id=? GROUP BY seed_sha256'
            ') s JOIN seed_payloads p ON p.sha256=s.seed_sha256 ORDER BY s.ordinal,s.seed_sha256',
            (context_id,context_id))
        with failure_boundary('bounded_products', 'context_merge', known_errors=(WgsContextError, ContextPayloadError)):
            seed = _merge_seed_stream(
                (_load_seed(expand(content, payloads.get)) for (content,) in seed_cursor),
                execution_config.contexts.max_context_bytes)
        if not seed.allele_refs:
            raise WgsContextError('nonempty context lost its assessed members')
        merged = ordinal
        analysis_progress(progress, 'context_assembly', 'context_merge',
            counts={'contexts_merged':merged, 'contexts_emitted':merged},
            unit='contexts', completed=merged, total=total_contexts)
        yield ClosedContext(ordinal, seed, connection, admission_reader, execution_config,
                            scratch_directory=scratch_path.parent, payloads=payloads,
                            _seal=_CLOSED_CONTEXT_SEAL)


def iter_closed_contexts(
    admissions: Iterable[ShardedCandidateAdmission],
    *,
    admission_reader: AdmissionProjectionReader,
    scratch_path: str | Path,
    execution_config: BoundedWgsExecutionConfig,
    context_executor=None,
    validation_audit=None,
    writer_profile=None,
    assessment_owner=None,
    profile_reader=None,
    science_selector=None,
) -> Iterator[ClosedContext]:
    """Globally group all admission shards before yielding any complete context."""

    scratch_path = Path(scratch_path)
    if (
        not isinstance(admission_reader,AdmissionProjectionReader)
        or not scratch_path.is_absolute()
        or not isinstance(execution_config, BoundedWgsExecutionConfig)
        or execution_config.contexts.pair_truncation_allowed
    ):
        raise WgsContextError("global context cursor binding is invalid")
    profile = ContextWriterProfile() if writer_profile is None else writer_profile
    connection = _open_index(scratch_path, execution_config, context_wal=True, writer_profile=profile)
    payloads = Payloads(connection)
    transactions = _ContextTransactions(connection, scratch_path, profile)
    payloads.transactions = transactions
    profile_rows=None
    from genosapiens.wgs_runtime import analysis_progress
    progress = getattr(context_executor, 'progress', None)
    merged = 0
    try:
        payloads.create()
        from genosapiens.wgs_context_admission import (
            ContextAdmissionExecutor, SparseContextBatch, EncodedContextPart, _SEAL,
        )
        if context_executor is not None and type(context_executor) is not ContextAdmissionExecutor:
            raise WgsContextError('global context executor is invalid')
        from genosapiens.wgs_association_projection import _discovery_domain
        if profile_reader is None:
            expected_inputs=dict(admission_reader.eligibility_manifest.input_dispositions)['assess']
            profile_sha256=None
        else:
            expected_inputs,_,_,profile_sha256=_discovery_domain(admission_reader._owner,profile_reader)
        if context_executor is not None and context_executor.profile_reader is not profile_reader:
            raise WgsContextError('global context profile owner differs')
        if context_executor is not None:
            transactions.observation=context_executor.observe_writer
        if assessment_owner is not None:
            from genosapiens.wgs_association_projection import AssessmentProjectionOwner
            if (type(assessment_owner) is not AssessmentProjectionOwner
                    or assessment_owner.profile_reader is not profile_reader):
                raise WgsContextError('foreign assessment evidence owner')
        connection.executescript("""
            CREATE TABLE batches(
                batch_ordinal INTEGER PRIMARY KEY,
                parent_source_sha256 TEXT NOT NULL,
                parent_working_sha256 TEXT NOT NULL,
                native_contract_sha256 TEXT NOT NULL,
                policy_sha256 TEXT NOT NULL
            );
            CREATE TABLE inputs(
                retained_ordinal INTEGER PRIMARY KEY,
                allele_id TEXT UNIQUE,
                native_identity_id TEXT,
                projection_binding_sha256 TEXT NOT NULL,
                identity_id TEXT NOT NULL UNIQUE,
                input_sha256 TEXT NOT NULL,
                input_size INTEGER NOT NULL
            );
            CREATE TABLE assessment_subjects(subject_id TEXT PRIMARY KEY,retained_ordinal INTEGER UNIQUE NOT NULL
                REFERENCES inputs(retained_ordinal));
            CREATE TABLE contexts(
                context_id TEXT PRIMARY KEY,
                grouping_bytes BLOB NOT NULL
            );
            CREATE TABLE seed_payloads(sha256 TEXT PRIMARY KEY,canonical_seed BLOB NOT NULL);
            CREATE TABLE baseline_seeds(
                context_id TEXT PRIMARY KEY REFERENCES contexts(context_id),
                seed_sha256 TEXT NOT NULL REFERENCES seed_payloads(sha256)
            );
            CREATE TABLE seeds(
                context_id TEXT NOT NULL REFERENCES contexts(context_id),
                batch_ordinal INTEGER NOT NULL REFERENCES batches(batch_ordinal),
                seed_sha256 TEXT NOT NULL REFERENCES seed_payloads(sha256),
                PRIMARY KEY(context_id,batch_ordinal,seed_sha256)
            );
            CREATE TABLE members(
                context_id TEXT NOT NULL REFERENCES contexts(context_id),
                allele_id TEXT NOT NULL REFERENCES inputs(allele_id),
                PRIMARY KEY(context_id,allele_id)
            );
            CREATE TABLE decisions(
                batch_ordinal INTEGER NOT NULL REFERENCES batches(batch_ordinal),
                decision_id TEXT NOT NULL,
                payload_sha256 TEXT NOT NULL REFERENCES context_payloads(sha256),
                PRIMARY KEY(decision_id,batch_ordinal,payload_sha256)
            );
            CREATE TABLE public_payloads(sha256 TEXT PRIMARY KEY,canonical_source BLOB NOT NULL);
            CREATE TABLE public_sources(
                batch_ordinal INTEGER NOT NULL REFERENCES batches(batch_ordinal),
                evidence_id TEXT NOT NULL,
                source_sha256 TEXT NOT NULL REFERENCES public_payloads(sha256),
                PRIMARY KEY(batch_ordinal,evidence_id)
            );
        """)
        if validation_audit is not None or profile.audit_index == 'always':
            connection.execute('CREATE INDEX decisions_by_batch ON decisions(batch_ordinal,decision_id,payload_sha256)')
        def _store_part(part, batch_ordinal=None, batch_inputs=None):
            if type(part) is not EncodedContextPart or part._seal is not _SEAL:
                raise WgsContextError('context rows have no validated worker owner')
            payloads.put_fragments(part.fragments)
            for key, grouping, content, refs in part.seeds:
                found = connection.execute('SELECT grouping_bytes FROM contexts WHERE context_id=?', (key,)).fetchone()
                if found is None:
                    connection.execute('INSERT INTO contexts VALUES (?,?)', (key, grouping))
                elif found[0] != grouping:
                    raise WgsContextError('cross-shard context grouping conflicts')
                references = dict(refs)
                if (len(references) != len(refs) or batch_ordinal is None and references
                        or batch_ordinal is not None and any(batch_inputs.get(k) != v for k,v in refs)):
                    raise WgsContextError('candidate context membership is invalid')
                digest = hashlib.sha256(content).hexdigest()
                connection.execute('INSERT OR IGNORE INTO seed_payloads VALUES (?,?)', (digest, content))
                if batch_ordinal is None:
                    connection.execute('INSERT INTO baseline_seeds VALUES (?,?)', (key, digest))
                else:
                    connection.execute('INSERT OR IGNORE INTO seeds VALUES (?,?,?)', (key,batch_ordinal,digest))
                    connection.executemany('INSERT OR IGNORE INTO members VALUES (?,?)', ((key,k) for k in references))
            if batch_ordinal is None:
                if part.decisions or part.public_sources:
                    raise WgsContextError('baseline context has batch evidence')
                transactions.add(sum(len(row[1])+len(row[2]) for row in part.seeds)
                                 + sum(len(row[2]) for row in part.fragments))
                return
            connection.executemany('INSERT OR IGNORE INTO decisions VALUES (?,?,?)',
                ((batch_ordinal,key,payloads.put(content)) for key,content in part.decisions))
            for key, content in part.public_sources:
                digest = hashlib.sha256(content).hexdigest()
                prior = connection.execute('SELECT source_sha256 FROM public_sources WHERE batch_ordinal=? AND evidence_id=?',
                    (batch_ordinal,key)).fetchone()
                if prior is not None and prior[0] != digest:
                    raise WgsContextError('public source identity conflicts within batch')
                connection.execute('INSERT OR IGNORE INTO public_payloads VALUES (?,?)',(digest,content))
                connection.execute('INSERT OR IGNORE INTO public_sources VALUES (?,?,?)',(batch_ordinal,key,digest))
            transactions.add(sum(len(row[1])+len(row[2]) for row in part.seeds)
                + sum(len(row[1]) for rows in (part.decisions, part.public_sources) for row in rows)
                + sum(len(row[2]) for row in part.fragments))

        def store_part(part,batch_ordinal=None,batch_inputs=None):
            started=time.monotonic()
            try:
                with failure_boundary('bounded_products', 'context_store', known_errors=(WgsContextError, ContextPayloadError)):
                    return _store_part(part,batch_ordinal,batch_inputs)
            finally:
                transactions.part_write_seconds+=time.monotonic()-started
                run_cleanup([('diagnostic_publish', transactions.observe)], primary=sys.exc_info()[1])

        # The producer completion retains the full original source accounting.
        # Context input rows are written as subjects are actually consumed.
        if context_executor is not None:
            context_executor.reuse_completed_accounting(admission_reader)
        if profile_reader is not None and context_executor is None:
            profile_rows = (row for row in profile_reader.iter_rows()
                if row.decision['disposition'] in {'selected','unresolved'})
        parents = policy = None
        batch_count = input_count = context_rows = previous_ordinal = 0
        for batch_ordinal, admission in enumerate(admissions, start=1):
            sparse = type(admission) is SparseContextBatch
            if type(admission) is not ShardedCandidateAdmission and not sparse:
                raise WgsContextError("candidate admission shard is invalid")
            if sparse:
                admission.require_source()
            values=admission.subjects if sparse else admission.inputs
            if (not values
                    or sparse and (context_executor is None or admission._seal is not _SEAL)):
                raise WgsContextError("candidate admission shard is invalid")
            current_parents = (
                values[0].parent_source_manifest_sha256,
                values[0].parent_working_set_manifest_sha256,
                values[0].native_source_contract_sha256,
            )
            current_policy = admission.binding.get("policy_sha256")
            if parents is None:
                parents, policy = current_parents, current_policy
            elif current_parents != parents or current_policy != policy:
                raise WgsContextError("candidate admission shard parents conflict")
            connection.execute(
                "INSERT INTO batches VALUES (?,?,?,?,?)",
                (batch_ordinal, *current_parents, current_policy),
            )
            batch_inputs = {}
            references=admission_reader.refs_for_subjects if sparse else admission_reader.refs_for_inputs
            refs=references(values,expected_rows_sha256=admission.binding.get('input_rows_sha256'))
            for value,ref in zip(values,refs,strict=True):
                if profile_rows is not None:
                    selected = next(profile_rows, None)
                    if selected is None or selected.original_ref != ref:
                        raise WgsContextError('candidate admission is outside selected profile')
                if value.retained_ordinal <= previous_ordinal:
                    raise WgsContextError('candidate admission assessed row order is invalid')
                previous_ordinal=value.retained_ordinal
                allele_id=value.allele_id if sparse else None if value.allele is None else value.allele.allele_id
                native_id=value.native_identity_id if sparse else None if value.allele is None else value.allele.native_identity_id
                connection.execute('INSERT INTO inputs VALUES (?,?,?,?,?,?,?)',
                    (ref[1], ref[3], native_id, ref[0], ref[2], ref[4], ref[5]))
                if allele_id is not None: batch_inputs[allele_id] = native_id
                connection.execute('INSERT INTO assessment_subjects VALUES (?,?)',(value.subject_id,value.retained_ordinal))
                input_count += 1
                transactions.add(sum(len(v) for v in ref if type(v) is str)+64)
            if sparse:
                with closing(admission.parts) as parts:
                    for part in parts:
                        store_part(part, batch_ordinal, batch_inputs)
                        context_rows += len(part.seeds)
                admission.require_source()
                if validation_audit is not None:
                    validation_audit.admission_projected(admission.subjects, sorted(payloads.expanded(row[0]) for row in connection.execute(
                        'SELECT payload_sha256 FROM decisions WHERE batch_ordinal=?',
                        (batch_ordinal,))),admission_reader=admission_reader)
                batch_count += 1
                if profile.sparse_batch_commit:
                    transactions.commit()
                continue
            seen_contexts = set()
            for seed in admission.contexts:
                if seed.context.context_id in seen_contexts:
                    raise WgsContextError('duplicate context within batch')
                seen_contexts.add(seed.context.context_id)
                grouping = _grouping_bytes(seed)
                found = connection.execute(
                    "SELECT grouping_bytes FROM contexts WHERE context_id=?",
                    (seed.context.context_id,),
                ).fetchone()
                if found is None:
                    connection.execute(
                        "INSERT INTO contexts VALUES (?,?)",
                        (seed.context.context_id, grouping),
                    )
                elif found[0] != grouping:
                    raise WgsContextError("cross-shard context grouping conflicts")
                references = {value.allele_id: value for value in seed.allele_refs}
                if (
                    len(references) != len(seed.allele_refs)
                    or set(references) - set(batch_inputs)
                    or any(
                        reference.identity_id != batch_inputs[allele_id]
                        for allele_id, reference in references.items()
                    )
                ):
                    raise WgsContextError("candidate context membership is invalid")
                encoder = Encoder()
                content = encoder.record(seed)
                payloads.put_fragments(encoder.rows())
                seed_sha = hashlib.sha256(content).hexdigest()
                connection.execute('INSERT OR IGNORE INTO seed_payloads VALUES (?,?)',(seed_sha,content))
                connection.execute('INSERT INTO seeds VALUES (?,?,?)',
                    (seed.context.context_id,batch_ordinal,seed_sha))
                connection.executemany(
                    "INSERT OR IGNORE INTO members VALUES (?,?)",
                    (
                        (seed.context.context_id, reference.allele_id)
                        for reference in seed.allele_refs
                    ),
                )
                context_rows += 1
                transactions.add(len(content)+sum(len(row[2]) for row in encoder.rows()))
            for value in admission.decisions:
                encoder = Encoder()
                content = encoder.record(value)
                payloads.put_fragments(encoder.rows())
                connection.execute('INSERT INTO decisions VALUES (?,?,?)',
                    (batch_ordinal,value.decision_id,payloads.put(content)))
                transactions.add(len(content)+sum(len(row[2]) for row in encoder.rows()))
            for value in admission.public_sources:
                encoder = Encoder()
                content = encoder.record(value)
                payloads.put_fragments(encoder.rows())
                digest = hashlib.sha256(content).hexdigest()
                connection.execute('INSERT OR IGNORE INTO public_payloads VALUES (?,?)',(digest,content))
                connection.execute('INSERT INTO public_sources VALUES (?,?,?)',
                    (batch_ordinal,value.evidence_id,digest))
                transactions.add(len(content)+sum(len(row[2]) for row in encoder.rows()))
            batch_count += 1
            if batch_count % execution_config.sqlite.cursor_batch_rows == 0:
                transactions.commit()
        transactions.commit(final=True)
        if input_count != expected_inputs:
            raise WgsContextError('candidate admission global assessed coverage is incomplete')
        if profile_reader is not None:
            _discovery_domain(admission_reader._owner,profile_reader)
        if context_executor is not None:
            with closing(context_executor.baseline_parts()) as parts:
                for part in parts: store_part(part)
            transactions.commit()
            context_executor.finish()
        if assessment_owner is not None:
            if policy is not None and policy!=assessment_owner.policy_sha256:
                raise WgsContextError('assessment policy closure differs')
            analysis_progress(progress, 'context_assembly', 'assessment_publication')
            with failure_boundary('bounded_products', 'context_assessment_publish'):
                assessment_owner.publish(connection,payloads,admission_reader)
        if profile_rows is not None and next(profile_rows, None) is not None:
            raise WgsContextError('global context profile accounting is incomplete')
        if batch_count == 0:
            if science_selector is not None:
                yield from _iter_stored_contexts(connection, payloads, admission_reader, execution_config,
                    scratch_path=scratch_path, progress=progress, total_contexts=0,
                    science_selector=science_selector)
            if context_executor is not None:
                context_executor.write('complete', force=True)
            return
        if connection.execute("PRAGMA foreign_key_check").fetchone() is not None:
            raise WgsContextError("global context index foreign keys are invalid")
        if connection.execute("PRAGMA integrity_check").fetchone()[0] != "ok":
            raise WgsContextError("global context index integrity is invalid")

        total_contexts = len(context_executor.context_keys) if context_executor is not None else None
        analysis_progress(progress, 'context_assembly', 'context_merge',
            counts={'contexts_obtained':total_contexts,'contexts_merged':0,'contexts_emitted':0})
        yield from _iter_stored_contexts(connection, payloads, admission_reader, execution_config,
            scratch_path=scratch_path, progress=progress, total_contexts=total_contexts,
            science_selector=science_selector)
        if context_executor is not None:
            context_executor.write('complete', force=True)
    except BaseException as error:
        if isinstance(error, (WgsFailure, WgsContextError, KeyboardInterrupt, GeneratorExit)):
            primary = error
        elif isinstance(error, (AdmissionProjectionError, ContextPayloadError)):
            primary = WgsFailure('bounded_products', 'malformed_evidence',
                operation='context_store', exception_kind='project_failure')
        else:
            primary = static_failure(error, 'bounded_products', operation='context_store')
        if context_executor is not None:
            run_cleanup([('diagnostic_publish', lambda: context_executor.write(
                'interrupted' if isinstance(error, GeneratorExit) else 'failed', force=True))], primary=primary)
        raise primary from None
    finally:
        actions = []
        if profile_rows is not None:
            actions.append(('worker_shutdown', profile_rows.close))
        close = getattr(admissions, 'close', None)
        if close is not None:
            actions.append(('worker_shutdown', close))
        actions.append(('context_store', connection.close))
        for path in (scratch_path, scratch_path.with_name(scratch_path.name+'-wal'),
                     scratch_path.with_name(scratch_path.name+'-shm')):
            actions.append(('context_store', lambda path=path: path.unlink(missing_ok=True)))
        run_cleanup(actions, primary=sys.exc_info()[1])


def _materialize_complete_context(closed):
    """Small caller adapter over the same bounded, globally closed owner."""
    budget = closed._config.contexts.max_context_bytes
    total = len(_bytes(closed.seed)) + closed.input_bytes
    if total > budget:
        raise WgsContextError('complete context exceeds capacity')
    members = tuple(closed.iter_members())
    decisions = []
    for value in closed.iter_decisions():
        total += len(_bytes(value))
        if total > budget:
            raise WgsContextError('complete context exceeds capacity')
        decisions.append(value)
    public = []
    evidence_id = manifest = None
    records = {}
    for row in closed.iter_public_records():
        if evidence_id is not None and row.evidence_id != evidence_id:
            public.append(PublicSourcedEvidence(evidence_id, manifest, records))
            records = {}
        evidence_id, manifest = row.evidence_id, row.public_manifest_sha256
        if row.record is not None:
            records[row.record.source_record_id] = row.record
            total += len(_bytes(row.record))
            if total > budget:
                raise WgsContextError('complete context exceeds capacity')
    if evidence_id is not None:
        public.append(PublicSourcedEvidence(evidence_id, manifest, records))
    decisions, public = tuple(decisions), tuple(public)
    total = len(_bytes(closed.seed)) + closed.input_bytes
    total += sum(len(_bytes(value)) for value in (*decisions, *public))
    if total > budget:
        raise WgsContextError('complete context exceeds capacity')
    return CompleteContext(closed.ordinal, closed.seed, members, total, decisions, public,
        _ContextByteOwner(closed.seed,members,decisions,public,total,_CONTEXT_BYTES_SEAL))


def iter_complete_contexts(admissions, *, admission_reader, scratch_path, execution_config):
    """Materialize only small contexts; the bounded engine consumes closed cursors."""
    for closed in iter_closed_contexts(admissions, admission_reader=admission_reader,
            scratch_path=scratch_path, execution_config=execution_config):
        yield _materialize_complete_context(closed)
