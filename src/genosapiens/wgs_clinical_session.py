"""One authenticated bounded clinical source session for WGS analysis.

The canonical public protection index supplies complete negative lookup coverage
and exact raw-source locators.  Matching source rows are replayed from ClinVar
with its tabix index and from the named Open Targets parquet row.  Resource
authentication and DuckDB setup happen once for the stage, not once per
admission shard or complete context.
"""

from __future__ import annotations

from collections import OrderedDict
from contextlib import AbstractContextManager, ExitStack, contextmanager
from dataclasses import dataclass, replace
import hashlib
import json
import os
from pathlib import Path
import shutil
import socket
import sqlite3
from types import MappingProxyType

from genosapiens import admission_links as links, clinical_annotation as clinical
from genosapiens import clinvar as cv
from genosapiens.bounded_wgs_config import BoundedWgsExecutionConfig
from genosapiens.candidate_admission import (
    ShardedClinicalAssociations,
    ShardedClinicalProof,
)
from genosapiens.common_allele_public_indexes import (
    MAX_ROW_BYTES,
    CommonAllelePublicIndexSession,
    PublicIndexQuery,
    VerifiedPublicIndexes,
)
from genosapiens.evidence_model import (
    AlleleEvidence,
    ClinicalEvidence,
    EvidenceReference,
    PublicSourcedEvidence,
    SourcedRecord,
    _typed,
)
from genosapiens.evidence_serialization import _plain
from genosapiens.wgs_failures import WgsFailure, failure_boundary, static_failure
from genosapiens.wgs_query_limits import QueryBudget


class WgsClinicalSessionError(RuntimeError):
    """Static failure while replaying authenticated public clinical sources."""


def _canonical(value: object) -> bytes:
    return (
        json.dumps(
            _plain(value), sort_keys=True, separators=(",", ":"),
            ensure_ascii=True, allow_nan=False,
        )
        + "\n"
    ).encode("ascii")


def _reference(record: SourcedRecord) -> EvidenceReference:
    manifest = record.public_manifest_sha256
    return EvidenceReference(
        "clinical-source-" + manifest,
        record.source_record_id,
        manifest,
    )


@dataclass(frozen=True, slots=True, repr=False)
class _Query:
    allele: AlleleEvidence
    chromosome: str | None
    position: int | None
    reference: str | None
    alternate: str | None
    reason: str | None

def _adapt_clinical_sources(sources):
    clinical_rows = [clinical._projection(record) for record in sources]
    conflicts = {
        _reference(value.source): set() for value in clinical_rows
    }
    for item in clinical_rows:
        row = item.source.source_fields["record"]
        dataset = item.source.source_dataset
        fields = row["info"] if dataset == "ClinVar VCF" else row
        affirmative = clinical.affirmative_source_conflict(dataset, fields)
        if affirmative or len(clinical._polarities(item.significance)) == 2:
            conflicts[_reference(item.source)].add(
                _reference(item.source)
            )
    linked = [
        value for value in clinical_rows
        if value.explicit_condition_link
    ]
    for index, left in enumerate(linked):
        for right in linked[index + 1:]:
            if (left.condition_ids, left.target_id) != (
                right.condition_ids, right.target_id,
            ):
                continue
            left_signs = clinical._polarities(left.significance)
            right_signs = clinical._polarities(right.significance)
            if (
                "pathogenic" in left_signs
                and "benign" in right_signs
                or "benign" in left_signs
                and "pathogenic" in right_signs
            ):
                conflicts[_reference(left.source)].add(
                    _reference(right.source)
                )
                conflicts[_reference(right.source)].add(
                    _reference(left.source)
                )
    clinical_rows = tuple(replace(
        value,
        conflict_refs=tuple(sorted(
            conflicts[_reference(value.source)],
            key=lambda child: (child.evidence_id, child.record_id),
        )),
    ) for value in clinical_rows)
    return clinical_rows



_RETENTION_CACHE_BYTES = 16 * 1024 * 1024


def _retention_binding(session):
    from genosapiens.variant_retention import POLICY_SHA256
    return (session.indexes.manifest_sha256, session.indexes.closure_sha256,
            *((source.root, source.manifest_sha256, source.contract_sha256)
              for source in (session.open_targets, session.clinvar, session.reference)),
            session._index, POLICY_SHA256, _adapt_retention_records)


def _adapt_retention_records(session, records, exact, bindings, reference_identity):
    from genosapiens.retention_sources import RetentionEvidence, clinical_from_index
    from genosapiens.curated_alleles import functional_from_index
    from genosapiens.variant_retention import RetentionError
    from genosapiens.wgs_admission_diagnostics import boundary, tag
    self = session
    with boundary('profile_clinical_adapt', reason='malformed_evidence', known_errors=(RetentionError,)):
        clinical_facts, functional_facts = [], []
        for row in records:
            if row['record_type'] != 'exact':
                continue
            if row['source'] in {'clinvar', 'open_targets'}:
                fact = clinical_from_index(row, reference_identity=reference_identity,
                    reference_sha256=self.reference.manifest_sha256,
                    source_release_sha256=bindings[row['source']])
                if fact.allele != exact:
                    raise tag(WgsClinicalSessionError('selection index exact identity differs'),
                        'profile_clinical_bind', 'exact_identity_mismatch')
                clinical_facts.append(fact)
            elif row['source'] == 'curated_functional':
                fact = functional_from_index(row)
                if fact.allele != exact or fact.source.release_sha256 != bindings['curated']:
                    raise tag(WgsClinicalSessionError('selection curated identity or source differs'),
                        'profile_clinical_bind',
                        'source_binding_mismatch' if fact.source.release_sha256 != bindings['curated']
                        else 'exact_identity_mismatch')
                functional_facts.append(fact)
        return tuple(clinical_facts), tuple(functional_facts)


class _RetentionWindow:
    """Private, ordered access to one session's normalized batch results.

    Raw results retain the existing query byte reservation. Adaptation consumes
    one result at a time, so dense windows never retain a second typed batch.
    Typed entries and the session's serialized cache share its 16-MiB allowance.
    Eviction only repeats adaptation of owned raw rows, never a physical query.
    No entry survives this operation or crosses an external codec boundary.
    """
    def __init__(self, session, queries, requested, records, missing, complete,
                 bindings, reference_identity, budget):
        self.session, self.binding = session, _retention_binding(session)
        self.queries, self.requested, self.records = queries, requested, records
        self.allele_keys = tuple((q.allele.allele_id, q.allele.native_identity_id,
            _canonical(links._native(q.allele))) for q in queries)
        self.missing, self.complete, self.bindings = missing, complete, bindings
        self.reference_identity, self.budget = reference_identity, budget
        self.position = self.bytes = 0
        self.cache = OrderedDict()
        self.closed = False

    def _check(self, session):
        if (self.closed or session is not self.session
                or getattr(session, '_retention_owner', None) is not self
                or _retention_binding(session) != self.binding):
            raise WgsClinicalSessionError('clinical result owner or source binding differs')
        session._active()
        session._index.assert_stable()
        if session._reference_sites is not None:
            session._reference_sites.assert_stable()

    def take(self, session, allele):
        from genosapiens.retention_sources import RetentionEvidence
        self._check(session)
        if (self.position >= len(self.queries) or self.queries[self.position].allele is not allele
                or (allele.allele_id, allele.native_identity_id) != self.allele_keys[self.position][:2]
                or _canonical(links._native(allele)) != self.allele_keys[self.position][2]):
            raise WgsClinicalSessionError('clinical result allele owner differs')
        query, request = self.queries[self.position], self.requested[self.position]
        self.position += 1
        if request is None:
            return RetentionEvidence(None, (), (), self.complete, query.reason)
        key, exact = request
        if key in self.cache:
            self.cache.move_to_end(key)
            return self.cache[key][0]
        # Reserve before constructing typed children. Canonical source fields
        # bound copied text. The fixed allowance covers each small typed graph.
        charge = sum(4 * len(_canonical(row)) + 4096 for row in self.records[key]) + len(key) + 4096
        self.budget.check()
        if charge > self.budget.limits.max_response_bytes:
            from genosapiens.wgs_failures import WgsFailure
            raise WgsFailure('clinical_adaptation', 'response_capacity',
                operation='clinical_response_adapt', exception_kind='project_failure')
        while self.cache and self.bytes + session._retention_cache_bytes + charge > _RETENTION_CACHE_BYTES:
            _, (_, size) = self.cache.popitem(last=False)
            self.bytes -= size
        while session._retention_cache and self.bytes + session._retention_cache_bytes + charge > _RETENTION_CACHE_BYTES:
            old_key, old_content = session._retention_cache.popitem(last=False)
            session._retention_cache_bytes -= len(old_key) + len(old_content)
        clinical, functional = _adapt_retention_records(session, self.records[key], exact,
            self.bindings, self.reference_identity)
        result = RetentionEvidence(exact, clinical, functional, self.complete)
        session.retention_metrics['adaptations'] = session.retention_metrics.get('adaptations', 0) + 1
        # An oversized entry uses only the bounded transient query allowance.
        # It is consumed before the next entry and is never retained in cache.
        if charge <= _RETENTION_CACHE_BYTES:
            self.cache[key] = result, charge
            self.bytes += charge
        session.retention_metrics['typed_peak_bytes'] = max(
            session.retention_metrics.get('typed_peak_bytes', 0), self.bytes, charge)
        return result

    def finish(self):
        self._check(self.session)
        if self.position != len(self.queries):
            raise WgsClinicalSessionError('clinical result window is incomplete')
        # Publish only after all required evidence passed adaptation. Release
        # typed entries first so both representations share the same allowance.
        self.cache.clear()
        self.bytes = 0
        for key in self.missing:
            content = _canonical(self.records[key])
            if len(content) + len(key) > _RETENTION_CACHE_BYTES:
                continue
            self.session._retention_cache[key] = content
            self.session._retention_cache_bytes += len(content) + len(key)
            while len(self.session._retention_cache) > 256 or self.session._retention_cache_bytes > _RETENTION_CACHE_BYTES:
                old_key, old_content = self.session._retention_cache.popitem(last=False)
                self.session._retention_cache_bytes -= len(old_key) + len(old_content)

    def close(self):
        self.closed = True
        self.cache.clear()
        self.records.clear()
        self.queries = self.requested = ()
        self.allele_keys = ()
        self.bytes = 0


class WgsClinicalResourceSession(
    AbstractContextManager["WgsClinicalResourceSession"]
):
    """Reuse qualified index/source bindings with stable readers for one stage."""

    def __init__(
        self,
        *,
        indexes: VerifiedPublicIndexes,
        open_targets: links.PinnedPublicSource,
        clinvar: links.PinnedPublicSource,
        reference: links.PinnedPublicSource,
        execution_config: BoundedWgsExecutionConfig,
        scratch_path: str | Path,
    ) -> None:
        scratch = Path(scratch_path)
        if (
            type(indexes) is not VerifiedPublicIndexes
            or type(open_targets) is not links.PinnedPublicSource
            or type(clinvar) is not links.PinnedPublicSource
            or type(reference) is not links.PinnedPublicSource
            or not isinstance(execution_config, BoundedWgsExecutionConfig)
            or not scratch.is_absolute()
        ):
            raise WgsClinicalSessionError("bounded clinical session is invalid")
        self.indexes = indexes
        self.open_targets = open_targets
        self.clinvar = clinvar
        self.reference = reference
        self.execution_config = execution_config
        self.scratch_path = scratch
        self._stack: ExitStack | None = None
        self._files = None
        self._index = None
        self._maps = None
        self._fasta = None
        self._reference_proof = None
        self._ot_release = None
        self._ot_contract = None
        self._ot_artifacts = None
        self._clinvar_contract = None
        self._clinvar_artifact = None
        self._cache = None
        self._query_ordinal = 0
        self.query_cache = None
        self._reference_sites = None
        self._retention_cache = OrderedDict()
        self._retention_cache_bytes = 0
        self.retention_metrics = dict(query_batches=0, unique_queries=0, cache_hits=0)

    def __enter__(self) -> "WgsClinicalResourceSession":
        if self._stack is not None or self.scratch_path.exists():
            raise WgsClinicalSessionError("bounded clinical session is already active")
        if any(name != "lo" for _, name in socket.if_nameindex()):
            raise WgsClinicalSessionError("bounded clinical session requires offline execution")
        self.scratch_path.mkdir(mode=0o700, parents=True)
        stack = ExitStack()
        try:
            files = links._Files(stack, reuse_published_outputs=True)
            maps, fasta, reference_proof = links._reference(files, self.reference)

            clinvar_data, clinvar_manifest = files.source(self.clinvar)
            clinvar_contract = cv.ClinvarContract.from_mapping(clinvar_data)
            cv._verify_manifest_contract(
                clinvar_manifest, clinvar_contract, self.clinvar.contract_sha256,
            )
            for artifact in clinvar_contract.artifacts:
                output = next(
                    value for value in clinvar_manifest.outputs
                    if value.path == artifact.path
                )
                files.open(
                    self.clinvar.root / artifact.path,
                    output.sha256,
                    output.size_bytes,
                )
            clinvar_artifact = clinvar_contract.artifact("vcf")
            cv._read_vcf_header(
                self.clinvar.root / clinvar_artifact.path,
                clinvar_contract,
            )

            ot_release, ot_contract, ot_artifacts = links._open_targets(
                files, self.open_targets,
                reuse_published_reconciliation=True,
            )
            cache = sqlite3.connect(self.scratch_path / "raw-record-cache.sqlite")
            stack.callback(cache.close)
            cache.execute("PRAGMA journal_mode=DELETE")
            cache.execute("PRAGMA synchronous=FULL")
            cache.execute("PRAGMA temp_store=FILE")
            cache.execute("PRAGMA trusted_schema=OFF")
            cache.execute("PRAGMA cache_size=-65536")
            cache.execute("""
                CREATE TABLE raw_records(
                    projection_sha256 TEXT NOT NULL,
                    record_ordinal INTEGER NOT NULL,
                    canonical_record BLOB NOT NULL,
                    record_sha256 TEXT NOT NULL,
                    PRIMARY KEY(projection_sha256, record_ordinal)
                )
            """)
            cache.execute('CREATE INDEX raw_record_hash ON raw_records(record_sha256)')
            cache.commit()

            self._files = files
            self._index = stack.enter_context(
                CommonAllelePublicIndexSession(self.indexes, limits=self.execution_config.public_queries)
            )
            self._maps = maps
            self._fasta = fasta
            self._reference_proof = reference_proof
            self._ot_release = ot_release
            self._ot_contract = ot_contract
            self._ot_artifacts = ot_artifacts
            self._clinvar_contract = clinvar_contract
            self._clinvar_artifact = clinvar_artifact
            self._cache = cache
            self._stack = stack
            return self
        except (Exception, KeyboardInterrupt) as error:
            primary = static_failure(
                error,
                "resource_session",
                operation="resource_session_open",
            )
            try:
                stack.close()
            except (Exception, KeyboardInterrupt) as cleanup_error:
                primary.add_secondary(static_failure(
                    cleanup_error,
                    "resource_session",
                    operation="resource_session_close",
                ))
            shutil.rmtree(self.scratch_path, ignore_errors=True)
            raise primary from None

    def _active(self) -> None:
        if self._stack is None:
            raise WgsClinicalSessionError("bounded clinical session is inactive")

    def _queries(self, alleles: tuple[AlleleEvidence, ...]) -> tuple[_Query, ...]:
        self._active()
        values = []
        for allele in links._inputs(alleles):
            sequence, reason = links._sequence(allele, self._maps, self._fasta)
            if sequence is None or not links._DNA.fullmatch(sequence[2]) \
                    or not links._DNA.fullmatch(sequence[3]):
                values.append(_Query(
                    allele, None, None, None, None,
                    reason or "ambiguous-sequence-exact-link-unresolved",
                ))
                continue
            chromosome, position, reference, alternate = sequence
            values.append(_Query(
                allele,
                links._public_query_chromosome(chromosome),
                position,
                reference,
                alternate,
                None,
            ))
        return tuple(values)

    def _budget(self, size: int) -> QueryBudget:
        budget = QueryBudget(self.execution_config.public_queries, max_row_bytes=MAX_ROW_BYTES, origin="clinical_replay")
        budget.batch(size)
        return budget

    def _cached(self, projection: object, budget: QueryBudget) -> tuple[SourcedRecord, ...] | None:
        key = hashlib.sha256(_canonical(projection)).hexdigest()
        records = []
        with budget.sqlite(self._cache):
            metadata = self._cache.execute(
                "SELECT record_ordinal, length(canonical_record) FROM raw_records "
                "WHERE projection_sha256=? ORDER BY record_ordinal", (key,),
            )
            for ordinal, size in metadata:
                budget.charge(size, "source_records")
                content, = self._cache.execute(
                    "SELECT canonical_record FROM raw_records WHERE projection_sha256=? AND record_ordinal=?",
                    (key, ordinal),
                ).fetchone()
                records.append(_typed(json.loads(content), SourcedRecord, decode=True))
        return tuple(records) if records else None

    def _store(self, projection: object, records: tuple[SourcedRecord, ...], budget: QueryBudget) -> None:
        key = hashlib.sha256(_canonical(projection)).hexdigest()
        with budget.sqlite(self._cache), self._cache:
            for ordinal, record in enumerate(records, start=1):
                content = _canonical(record)
                budget.charge(len(content), "source_records")
                self._cache.execute("INSERT INTO raw_records VALUES (?,?,?,?)",
                    (key,ordinal,content,hashlib.sha256(content).hexdigest()))

    def _indexed_records(self, projections, budget, adapter):
        """Materialize source facts from the authenticated bounded projection."""
        result = {}
        artifacts = {(Path(self._ot_contract.output_prefix) / item.relative_path).as_posix(): item
            for item in self._ot_artifacts}
        for projection in projections:
            key = hashlib.sha256(_canonical(projection)).hexdigest()
            if key in result:
                continue
            cached = self._cached(projection, budget)
            if cached is not None:
                result[key] = cached
                continue
            locator, raw = projection['source_locator'], projection['source_records']
            expected = projection['allele']
            identity = (expected['chromosome'], expected['position'], expected['reference'], expected['alternate'])
            if adapter == 'clinvar':
                parsed = raw['clinvar']
                # Reuse the strict parser to verify projected column/INFO agreement.
                if _canonical(links.parse_clinvar_vcf_row('\t'.join(parsed['columns']))) != _canonical(parsed):
                    raise WgsClinicalSessionError('indexed ClinVar source fields differ')
                columns = parsed['columns']
                if ((columns[0], int(columns[1]), columns[3], columns[4]) != identity
                        or columns[2] != projection['evidence']['variation_id']
                        or locator['path'] != self._clinvar_artifact.path
                        or locator['sha256'] != self._clinvar_artifact.sha256):
                    raise WgsClinicalSessionError('indexed ClinVar source identity differs')
                records = (links._source_record('ClinVar VCF', self._clinvar_contract.release,
                    locator['path'], locator['record_number'], locator['sha256'],
                    self.clinvar.manifest_sha256, parsed),)
            else:
                evidence, variant = raw['evidence'], raw['variant']
                target, condition = links.parse_open_targets_eva_row(evidence)
                variant_locator = locator['variant']
                if (links.parse_open_targets_variant_row(variant) != identity
                        or evidence['id'] != projection['evidence']['assertion_id']
                        or evidence['variantId'] != projection['evidence']['variant_id']
                        or (target, condition) != (projection['evidence']['target_id'],
                            projection['evidence']['disease_from_source_id'])):
                    raise WgsClinicalSessionError('indexed Open Targets source fields differ')
                records = []
                for dataset, item, row in (('Open Targets variant', variant_locator, variant),
                        ('Open Targets evidence_eva', locator, evidence)):
                    artifact = artifacts[item['path']]
                    if item['sha256'] != artifact.sha256:
                        raise WgsClinicalSessionError('indexed Open Targets source identity differs')
                    records.append(links._source_record(dataset, self._ot_contract.release,
                        artifact.relative_path, item['row_number'], artifact.sha256,
                        self.open_targets.manifest_sha256, row))
                records = tuple(records)
            self._store(projection, records, budget)
            result[key] = records
        return result

    def _clinvar_records(self, projections, budget):
        return self._indexed_records(projections, budget, 'clinvar')

    def _open_targets_records(self, projections, budget):
        return self._indexed_records(projections, budget, 'open_targets')

    @failure_boundary("clinical_replay", "clinical_response_bind")
    def use_case_cache(self, *, path, binding):
        from genosapiens.wgs_allele_query_cache import CaseAlleleQueryCache
        self._active()
        if self.query_cache is not None:
            raise WgsClinicalSessionError('case query cache is already bound')
        self.query_cache=CaseAlleleQueryCache(path,binding=binding,execution_config=self.execution_config)
        self._stack.callback(self.query_cache.close)
        return self.query_cache

    def _record_reference(self,record):
        content=_canonical(record)
        digest=hashlib.sha256(content).hexdigest()
        stored=self._cache.execute('SELECT canonical_record FROM raw_records WHERE record_sha256=? LIMIT 1',(digest,)).fetchone()
        if stored!=(content,):
            raise WgsClinicalSessionError('clinical raw record cache coverage is missing')
        return digest

    def _record_from_reference(self,digest,budget):
        if type(digest) is not str or len(digest)!=64 or any(c not in '0123456789abcdef' for c in digest):
            raise WgsClinicalSessionError('clinical raw record reference is invalid')
        metadata=self._cache.execute('SELECT projection_sha256,record_ordinal,length(canonical_record) '
            'FROM raw_records WHERE record_sha256=? LIMIT 1',(digest,)).fetchone()
        if metadata is None:
            raise WgsClinicalSessionError('clinical raw record reference is missing')
        projection,ordinal,size=metadata
        budget.charge(size,'source_records')
        content,=self._cache.execute('SELECT canonical_record FROM raw_records WHERE projection_sha256=? AND record_ordinal=?',
            (projection,ordinal)).fetchone()
        if hashlib.sha256(content).hexdigest()!=digest:
            raise WgsClinicalSessionError('clinical raw record cache changed')
        record=_typed(json.loads(content),SourcedRecord,decode=True)
        if _canonical(record)!=content:
            raise WgsClinicalSessionError('clinical raw record cache is invalid')
        return record

    def _query_fact_key(self,query):
        # Exact index scope is allele-wide; all gene/condition assertions remain
        # in the response. Population attachments do not change that query.
        return {'allele_id':query.allele.allele_id,'native':_plain(links._native(query.allele)),
            'chromosome':query.chromosome,'position':query.position,'reference':query.reference,
            'alternate':query.alternate,'reason':query.reason,'scope':'all-exact-gene-and-condition-assertions'}

    def _case_query_results(self,alleles,*,cache,budget=None):
        budget=budget or self._budget(len(alleles))
        queries=self._queries(alleles)
        keys=tuple(self._query_fact_key(q) for q in queries)
        facts=tuple(cache.get('clinical_association',key) for key in keys)
        missing=tuple(q.allele for q,fact in zip(queries,facts,strict=True) if fact is None)
        fresh={}
        if missing:
            batch=self._association_batch_uncached(missing,_budget=budget)
            for item in batch.proofs:fresh.setdefault(item.proof.allele_id,[]).append(item)
        result=[]
        for query,key,fact in zip(queries,keys,facts,strict=True):
            if fact is None:
                proofs=tuple(fresh.get(query.allele.allele_id,()))
                if not proofs:
                    raise WgsClinicalSessionError('clinical cache source coverage is incomplete')
                rows=[]
                for item in proofs:
                    data=_plain(item)
                    data['proof']['source_refs']=[self._record_reference(record) for record in item.proof.source_refs]
                    rows.append(data)
                status=('unsupported' if query.reason else 'uncertain' if any(
                    'applicable-public-source-record-unassessable' in item.proof.reasons for item in proofs) else
                    'positive' if any(item.proof.source_refs for item in proofs) else 'no_match')
                cache.put('clinical_association',key,status=status,value={'version':1,'proofs':rows})
            else:
                data=_plain(fact.value)
                if type(data) is not dict or set(data)!={'version','proofs'} or data['version']!=1 or type(data['proofs']) is not list:
                    raise WgsClinicalSessionError('clinical cached association is invalid')
                rows=[]
                for row in data['proofs']:
                    row['proof']['source_refs']=[_plain(self._record_from_reference(digest,budget))
                        for digest in row['proof']['source_refs']]
                    item=_typed(row,ShardedClinicalProof,decode=True)
                    if item.proof.allele_id!=query.allele.allele_id or item.proof.native!=links._native(query.allele):
                        raise WgsClinicalSessionError('clinical cached association binding changed')
                    budget.charge(len(_canonical(item)),'association_proofs')
                    rows.append(item)
                proofs=tuple(rows)
                status=('unsupported' if query.reason else 'uncertain' if any(
                    'applicable-public-source-record-unassessable' in item.proof.reasons for item in proofs) else
                    'positive' if any(item.proof.source_refs for item in proofs) else 'no_match')
                if fact.status!=status or {item.adapter for item in proofs}!={'clinvar','open_targets'}:
                    raise WgsClinicalSessionError('clinical cached association status changed')
            result.append(proofs)
        return tuple(result)

    def retention_batch(self, alleles):
        """Independent callers query and validate within a fresh operation."""
        with self._retention_window(alleles) as window:
            return tuple(window.take(self, allele) for allele in alleles)

    @contextmanager
    def _retention_window(self, alleles):
        """Own one bounded query and release typed reuse on success or failure."""
        if getattr(self, '_retention_owner', None) is not None:
            raise WgsClinicalSessionError('clinical result owner is already active')
        window = self._prepare_retention_window(alleles)
        self._retention_owner = window
        try:
            yield window
            window.finish()
        finally:
            window.close()
            self._retention_owner = None

    def _prepare_retention_window(self, alleles):
        """Batch unique exact facts without replaying raw public resources.

        This session cache includes completed empty queries. Its keys bind the
        index, reference and policy. It holds at most 256 entries and 16 MiB.
        """
        from genosapiens.wgs_admission_diagnostics import boundary, tag
        from genosapiens.variant_retention import ExactAllele, LookupCoverage, POLICY_SHA256
        budget = self._budget(len(alleles))
        by_identity = {q.allele.allele_id: q for q in self._queries(alleles)}
        queries = tuple(by_identity[a.allele_id] for a in alleles)
        with boundary('profile_clinical_bind'):
            if self._reference_sites is not None:
                self._reference_sites.assert_stable()
            coverage = self._index.retention_coverage()
            bindings = coverage['source_bindings']
            if (bindings['clinvar'] != self.clinvar.manifest_sha256
                    or bindings['open_targets'] != self.open_targets.manifest_sha256
                    or bindings['reference'] != self.reference.manifest_sha256):
                raise tag(WgsClinicalSessionError('selection index source binding differs'),
                        'profile_clinical_bind', 'source_binding_mismatch')
        reference_identity = next(key[1] for key in self._maps if key[0] == 'public_annotation_benchmark')
        complete = tuple(LookupCoverage(dataset, bindings[key], True, 'complete')
            for dataset, key in (('ClinVar VCF', 'clinvar'), ('Open Targets evidence_eva', 'open_targets'),
                                 ('gencc_submissions', 'curated'), ('evidence_clingen', 'curated')))
        requested, missing, records = [], {}, {}
        for query in queries:
            if query.reason is not None:
                requested.append(None)
                continue
            exact = ExactAllele(reference_identity, self.reference.manifest_sha256, query.chromosome,
                                query.position, query.reference, query.alternate)
            key = _canonical(dict(allele=exact, index=self.indexes.manifest_sha256, policy=POLICY_SHA256))
            requested.append((key, exact))
            if key in records or key in missing:
                continue
            cached = self._retention_cache.get(key)
            if cached is not None:
                self._retention_cache.move_to_end(key)
                budget.charge(len(cached), 'projection')
                records[key] = json.loads(cached)
                self.retention_metrics['cache_hits'] += 1
            else:
                missing[key] = PublicIndexQuery(query.chromosome, query.position, query.reference, query.alternate)
        if missing:
            with boundary('profile_clinical_query', blocking=True, reason='lookup_failed'):
                results = self._index.retention_many(tuple(missing.values()), budget=budget)
            self.retention_metrics['query_batches'] += 1
            self.retention_metrics['unique_queries'] += len(missing)
            for key, result in zip(missing, results, strict=True):
                resolved = None
                if result.status == 'retain_unknown':
                    from genosapiens.wgs_clinvar_reference_sites import ClinvarReferenceSites
                    with boundary('profile_clinical_query', reason='incomplete_required_evidence'):
                        if self._reference_sites is None:
                            self._reference_sites = ClinvarReferenceSites(self)
                        resolved = self._reference_sites.resolve(result.records, budget=budget, query=missing[key])
                if result.status not in {'matches', 'no_match'} and resolved is None:
                    raise tag(WgsClinicalSessionError('selection index has incomplete applicable source evidence'),
                        'profile_clinical_query', 'incomplete_required_evidence')
                records[key] = result.records if resolved is None else resolved
        return _RetentionWindow(self, queries, requested, records, missing, complete,
            bindings, reference_identity, budget)


    def association_batch(self,alleles,*,_budget=None):
        if self.query_cache is None:
            return self._association_batch_uncached(alleles,_budget=_budget)
        rows=self.query_cache.clinical_results(alleles,session=self,_budget=_budget)
        return ShardedClinicalAssociations(self.indexes.manifest_sha256,self.open_targets.manifest_sha256,
            self.clinvar.manifest_sha256,self.reference.manifest_sha256,
            tuple(proof for proofs in rows for proof in proofs))

    def _association_batch_uncached(
        self, alleles: tuple[AlleleEvidence, ...], *, _budget: QueryBudget | None = None,
    ) -> ShardedClinicalAssociations:
        """Return complete source-labelled reports for one bounded allele batch."""

        try:
            budget = _budget or self._budget(len(alleles))
            queries = self._queries(alleles)
            indexed = {}
            clinvar_projections = []
            open_targets_projections = []
            supported = tuple(query for query in queries if query.reason is None)
            with failure_boundary("clinical_index", "clinical_index_query"):
                results = self._index.query_many(tuple(PublicIndexQuery(
                    query.chromosome, query.position, query.reference, query.alternate,
                ) for query in supported), budget=budget) if supported else ()
            for query, result in zip(supported, results, strict=True):
                indexed[query.allele.allele_id] = result
                clinvar_projections.extend(
                    value for value in result.records
                    if value.get("source") == "clinvar"
                    and value.get("record_type") == "exact"
                )
                open_targets_projections.extend(
                    value for value in result.records
                    if value.get("source") == "open_targets"
                    and value.get("record_type") == "exact"
                )
            clinvar_rows = self._clinvar_records(tuple(clinvar_projections), budget)
            open_targets_rows = self._open_targets_records(
                tuple(open_targets_projections), budget,
            )
            proofs = []
            def append(proof):
                budget.charge(len(_canonical(proof)), "association_proofs")
                proofs.append(proof)
            for query in queries:
                result = indexed.get(query.allele.allele_id)
                records = () if result is None else result.records
                for adapter in ("open_targets", "clinvar"):
                    selected = tuple(
                        value for value in records
                        if value.get("source") == adapter
                        and value.get("record_type") == "exact"
                    )
                    provenance = links._provenance(
                        "bounded-" + adapter + "-public-index-session",
                        ({
                            "index_manifest_sha256": self.indexes.manifest_sha256,
                            "index_closure_sha256": self.indexes.closure_sha256,
                            "source_manifest_sha256": (
                                self.open_targets.manifest_sha256
                                if adapter == "open_targets"
                                else self.clinvar.manifest_sha256
                            ),
                            "raw_rows": "projected-source-facts-with-authenticated-locator",
                        },),
                        self._reference_proof,
                    )
                    if selected:
                        for projection in selected:
                            key = hashlib.sha256(
                                _canonical(projection)
                            ).hexdigest()
                            source_rows = (
                                open_targets_rows[key]
                                if adapter == "open_targets"
                                else clinvar_rows[key]
                            )
                            if adapter == "open_targets":
                                evidence = projection["evidence"]
                                target = evidence["target_id"]
                                condition = evidence["disease_from_source_id"]
                                outcome = (
                                    "admitted" if target or condition
                                    else "unresolved"
                                )
                                reasons = (
                                    "exact-variant-key-to-source-assertion"
                                    if target or condition
                                    else "source-target-and-condition-unresolved",
                                )
                                if condition is None:
                                    reasons += (
                                        "native-condition-unresolved; "
                                        "mapped-disease-not-equivalence",
                                    )
                                proof = links._proof(
                                    query.allele,
                                    "explicit_clinical_link",
                                    outcome,
                                    reasons,
                                    source_rows,
                                    provenance,
                                    target,
                                    condition,
                                )
                            else:
                                proof = links._proof(
                                    query.allele,
                                    "explicit_clinical_link",
                                    "unresolved",
                                    (
                                        "exact-allele-aggregate-retained; "
                                        "condition-assertion-link-unresolved",
                                    ),
                                    source_rows,
                                    provenance,
                                )
                            append(ShardedClinicalProof(adapter, proof))
                        if result.status == "retain_unknown":
                            append(ShardedClinicalProof(adapter, links._proof(
                                query.allele, "explicit_clinical_link", "unresolved",
                                ("applicable-public-source-record-unassessable",),
                                (), provenance,
                            )))
                    else:
                        reason = query.reason or (
                            "applicable-public-source-record-unassessable"
                            if result is not None
                            and result.status == "retain_unknown"
                            else "no-matched-source-assertion; "
                            "clinical-significance-unknown"
                        )
                        append(ShardedClinicalProof(
                            adapter,
                            links._proof(
                                query.allele,
                                "explicit_clinical_link",
                                "unresolved",
                                (reason,),
                                (),
                                provenance,
                            ),
                        ))
            return ShardedClinicalAssociations(
                self.indexes.manifest_sha256,
                self.open_targets.manifest_sha256,
                self.clinvar.manifest_sha256,
                self.reference.manifest_sha256,
                tuple(proofs),
            )
        except (WgsClinicalSessionError, WgsFailure, MemoryError):
            raise
        except Exception:
            raise WgsClinicalSessionError(
                "bounded clinical source query failed"
            ) from None

    def _clinical_rows(self,sources):
        if self.query_cache is None:
            return _adapt_clinical_sources(sources)
        by_hash={self._record_reference(record):record for record in sources}
        key={'source_records':tuple(by_hash),'scope':'complete-source-assertion-conflicts'}
        fact=self.query_cache.get('clinical_adapted',key)
        if fact is None:
            rows=_adapt_clinical_sources(sources)
            values=[]
            for row in rows:
                value=_plain(row)
                value['source']=self._record_reference(row.source)
                values.append(value)
            self.query_cache.put('clinical_adapted',key,status='positive' if rows else 'no_match',value=values)
            return rows
        if type(fact.value) is not tuple:
            raise WgsClinicalSessionError('cached clinical evidence is invalid')
        rows=[]
        for value in _plain(fact.value):
            if type(value) is not dict or value.get('source') not in by_hash:
                raise WgsClinicalSessionError('cached clinical source reference changed')
            value['source']=_plain(by_hash[value['source']])
            rows.append(_typed(value,ClinicalEvidence,decode=True))
        if fact.status!=('positive' if rows else 'no_match'):
            raise WgsClinicalSessionError('cached clinical evidence status changed')
        return tuple(rows)

    @failure_boundary("clinical_adaptation", "clinical_response_adapt")
    def annotate(
        self, alleles: tuple[AlleleEvidence, ...],
    ) -> clinical.ClinicalAnnotationResult:
        """Attach the same raw clinical rows without reopening whole sources."""

        try:
            budget = self._budget(len(alleles))
            alleles = links._inputs(alleles)
            batch = self.association_batch(alleles, _budget=budget)
            budget.origin = "clinical_adaptation"
            records = {value.allele_id: {} for value in alleles}
            gaps = {
                (value.allele_id, name): []
                for value in alleles for name in ("clinvar", "open_targets")
            }
            provenance = {}
            for item in batch.proofs:
                proof = item.proof
                provenance[item.adapter] = proof.provenance
                gaps[proof.allele_id, item.adapter].extend(proof.reasons)
                for record in proof.source_refs:
                    key = (
                        record.public_manifest_sha256,
                        record.source_record_id,
                    )
                    previous = records[proof.allele_id].setdefault(key, record)
                    if previous != record:
                        raise ValueError
            attached = []
            coverage = []
            stores = {}
            for allele in alleles:
                sources = sorted(
                    records[allele.allele_id].values(),
                    key=lambda value: (
                        value.source_dataset, value.source_record_id,
                    ),
                )
                clinical_rows = self._clinical_rows(sources)
                if any(value not in clinical_rows for value in allele.clinical or ()):
                    raise ValueError
                enriched = replace(
                    allele,
                    clinical=(
                        clinical_rows
                        or (None if allele.clinical is None else ())
                    ),
                )
                budget.charge(len(_canonical(enriched)), "annotation")
                attached.append(enriched)
                for name in ("clinvar", "open_targets"):
                    selected = tuple(
                        value for value in clinical_rows
                        if value.source.source_dataset.startswith(
                            "ClinVar" if name == "clinvar"
                            else "Open Targets"
                        )
                    )
                    assertions = tuple(
                        _reference(value.source) for value in selected
                        if value.assertion_id is not None
                    )
                    reasons = []
                    if assertions:
                        status = "assertions"
                    elif selected and name == "clinvar":
                        status = "aggregate_only"
                        reasons.append(
                            "aggregate-condition-assertion-link-unresolved"
                        )
                    else:
                        status = "unknown"
                        reasons.append(
                            "no-matched-source-assertion; "
                            "clinical-significance-unknown"
                        )
                    for value in selected:
                        if value.assertion_id is not None:
                            if not value.explicit_condition_link:
                                reasons.append(
                                    "native-condition-unresolved; "
                                    "mapped-disease-not-equivalence"
                                )
                            if value.target_id is None:
                                reasons.append(
                                    "source-target-unresolved; "
                                    "native-condition-retained"
                                )
                            if not value.significance or all(
                                child is None for child in value.significance
                            ):
                                reasons.append("source-significance-unknown")
                            if value.review_status is None:
                                reasons.append("source-review-status-unknown")
                            if value.evaluated_on is None:
                                reasons.append("source-evaluation-date-unknown")
                    if not selected:
                        reasons.extend(gaps[allele.allele_id, name])
                    conflict_refs = {
                        _reference(value.source) for value in selected
                        if value.conflict_refs
                    }
                    covered = clinical.ClinicalCoverage(
                        allele.allele_id,
                        allele.native_identity_id,
                        name,
                        status,
                        tuple(_reference(value.source) for value in selected),
                        assertions,
                        tuple(sorted(
                            conflict_refs,
                            key=lambda value: (
                                value.evidence_id, value.record_id,
                            ),
                        )),
                        tuple(sorted(set(reasons))),
                    )
                    budget.charge(len(_canonical(covered)), "annotation")
                    coverage.append(covered)
                for record in sources:
                    manifest = record.public_manifest_sha256
                    previous = stores.setdefault(manifest, {}).setdefault(
                        record.source_record_id, record,
                    )
                    if previous != record:
                        raise ValueError
                    budget.check()
            binding = {
                "schema": "genosapiens.wgs-bounded-clinical-session",
                "index_manifest_sha256": self.indexes.manifest_sha256,
                "index_closure_sha256": self.indexes.closure_sha256,
                "sources": {
                    "clinvar": self.clinvar.manifest_sha256,
                    "open_targets": self.open_targets.manifest_sha256,
                    "reference": self.reference.manifest_sha256,
                },
                "source_replay": _plain(provenance),
                "scope": "raw-rows-by-authenticated-index-locator",
            }
            digest = hashlib.sha256(_canonical(binding)).hexdigest()
            budget.charge(len(_canonical(binding)), "annotation")
            for values in stores.values():
                for record in values.values():
                    budget.charge(len(_canonical(record)), "annotation")
            return clinical.ClinicalAnnotationResult(
                tuple(attached),
                tuple(coverage),
                SourcedRecord(
                    "genosapiens.wgs-bounded-clinical-session",
                    "1",
                    digest,
                    digest,
                    None,
                    MappingProxyType(binding),
                ),
                tuple(
                    PublicSourcedEvidence(
                        "clinical-source-" + manifest,
                        manifest,
                        values,
                    )
                    for manifest, values in sorted(stores.items())
                ),
            )
        except (WgsClinicalSessionError, WgsFailure, MemoryError):
            raise
        except Exception:
            raise WgsClinicalSessionError(
                "bounded clinical annotation failed"
            ) from None

    def __exit__(self, *args: object) -> None:
        self._active()
        stack = self._stack
        self._stack = None
        try:
            stack.__exit__(*args)
        finally:
            self._reference_sites = None
            self._retention_cache.clear()
            self._retention_cache_bytes = 0
            self._files = None
            self._index = None
            self._maps = None
            self._fasta = None
            self._reference_proof = None
            self._ot_release = None
            self._ot_contract = None
            self._ot_artifacts = None
            self._clinvar_contract = None
            self._clinvar_artifact = None
            self._cache = None
