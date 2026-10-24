"""Context-owned protected population decisions before rich scientific work.

The original closed context stays authoritative. This owner publishes compact
decisions and a retained selection; it cannot qualify a release by itself.
"""
from contextlib import closing
from dataclasses import replace
import hashlib
import json
import os
import threading

from genosapiens.evidence_serialization import _canonical
from genosapiens.evidence_model import IdentityStage, AlleleReference, AlleleEvidence
from genosapiens.genotype_hypotheses import consuming_model_contexts, _requires_pairs
from genosapiens.population_annotation import PopulationRequest, PopulationResourceSession
from genosapiens.population_filter_candidate import (
    PopulationWorkPolicy, evaluate_source_candidate, evaluate_candidate,
    optional_frequency_support, context_rescue_vetoes,
)
from genosapiens.population_policy import ContextKey
from genosapiens.wgs_analysis_eligibility import AnalysisEligibilitySession
from genosapiens.wgs_contexts import ClosedContext
from genosapiens.wgs_evidence_store import EvidenceStoreError
from genosapiens.wgs_genotype_facts import ContextGenotypeFacts

MAX_DECISION_BYTES = 64 * 1024**2


def require(value, message='invalid population work-filter ownership'):
    if not value:
        raise EvidenceStoreError(message)


def population_records(allele):
    """Project the unchanged candidate predicate's source fields only."""
    fields = ('subset', 'stratum', 'ancestry', 'allele_count', 'allele_number',
        'homozygote_count', 'hemizygote_count', 'filtering_allele_frequency_99',
        'filters', 'callable', 'discordant')
    return tuple(dict(**{key: getattr(row, key) for key in fields},
        raw_context=dict((row.source.source_fields or {}).get('raw_context', {})))
        for row in allele.population or ())


class _ProtectionQueries:
    """Adapt the existing exact, bounded public query owner without new I/O."""
    def __init__(self, session):
        self.session = session

    def query_many(self, queries):
        return self.session.query_clinical(queries)


class ContextPopulationWorkFilter:
    def __init__(self, writer, closed, *, policy, resources, source_access,
                 population_session, clinical_annotator, query_cache, name):
        require(type(closed) is ClosedContext and type(policy) is PopulationWorkPolicy
            and type(source_access) is AnalysisEligibilitySession
            and type(population_session) is PopulationResourceSession
            and resources.writer is writer and resources.population.connection is writer.connection)
        require(source_access.owner_pid == os.getpid() and source_access.owner_thread == threading.get_ident()
            and not source_access.finished and source_access.source_sha256 == resources.source_manifest_sha256
            and source_access.routing.contract.reference_manifest_sha256 == resources.reference_binding.manifest_sha256)
        require(resources.parent_bindings.get('population_work_policy') == policy.sha256)
        from genosapiens.common_allele_filter import common_allele_filter_policy_as_json
        require(hashlib.sha256(common_allele_filter_policy_as_json(source_access.common_policy)).hexdigest()
            == json.loads(policy.content)['protection']['policy']['sha256'])
        if source_access.screening.merged_intake is None:
            source_access.public_session._assert_active()
            source_access.public_session._assert_clinical_stable()
        self.writer, self.connection, self.original = writer, writer.connection, closed
        self.policy, self.resources, self.source_access = policy, resources, source_access
        self.population_session, self.clinical_annotator, self.query_cache = population_session, clinical_annotator, query_cache
        self.name = name
        self.table = 'work_filter_' + hashlib.sha256(name.encode()).hexdigest()
        self.closed = False
        self.original_count = closed.member_count
        self.original_occurrences = self.retained_occurrences = 0
        self.retained_count = self.excluded_count = self.unresolved_count = 0
        self.input_bytes = 0
        try:
            self.connection.execute(f'CREATE TABLE {self.table} '
                '(id TEXT PRIMARY KEY, ordinal INTEGER UNIQUE, native TEXT, size INTEGER, disposition TEXT, content TEXT)')
            # Source genotype facts are compact shared interpretation, not a
            # constructed hypothesis universe or a second scoring pass.
            with ContextGenotypeFacts(writer, closed, sample_name=resources.sample_name,
                    source_manifest_sha256=resources.source_manifest_sha256,
                    reference_manifest_sha256=resources.reference_binding.manifest_sha256) as facts:
                self.original_genotype_inventory = facts.header.inventory_sha256
                for member, size in closed.iter_sized_members():
                    require(member.parent_source_manifest_sha256 == resources.source_manifest_sha256)
                    properties = facts.allele(member.allele.allele_id)
                    models = consuming_model_contexts(closed.seed.context, member.allele, properties)
                    value = self._assess(member, models)
                    value.update(allele_id=member.allele.allele_id, native_identity_id=member.allele.native_identity_id,
                        model_input=AlleleEvidence(member.allele.allele_id, member.allele.native_identity_id,
                            member.allele.identities), genotype_properties=properties,
                        identity_id=member.identity_id,
                        original_ordinal=member.retained_ordinal, source_allele_ids=member.source_allele_ids,
                        original_input_sha256=hashlib.sha256(member.as_bytes()).hexdigest(), original_input_bytes=size)
                    content = _canonical(value)
                    require(len(content) <= MAX_DECISION_BYTES, 'population work decision exceeds capacity')
                    self.connection.execute(f'INSERT INTO {self.table} VALUES (?,?,?,?,?,?)',
                        (member.allele.allele_id, member.retained_ordinal, member.allele.native_identity_id,
                         size, value['disposition'], content))
            # A mixed pair context needs the original potential partners to
            # preserve qualification/phase/completeness. Until a compact proof
            # can preserve those values, retain the implicated members.
            total, excluded = self.connection.execute(f"SELECT count(*),coalesce(sum(disposition='excluded'),0) FROM {self.table}").fetchone()
            require(total == closed.member_count)
            if _requires_pairs(closed.seed.context) and 0 < excluded < total:
                for key, content in self.connection.execute(f"SELECT id,content FROM {self.table} WHERE disposition='excluded'"):
                    value = json.loads(content)
                    value.update(disposition='retained_unresolved',
                        reasons=['excluded_partner_semantics_require_original_member'])
                    self.connection.execute(f'UPDATE {self.table} SET disposition=?,content=? WHERE id=?',
                        (value['disposition'], _canonical(value), key))
            refs = []
            for key, _, native, size, disposition, content in self.connection.execute(f'SELECT * FROM {self.table} ORDER BY ordinal'):
                count = len(json.loads(content)['source_allele_ids'])
                self.original_occurrences += count
                if disposition == 'excluded':
                    self.excluded_count += 1
                else:
                    self.retained_count += 1
                    self.retained_occurrences += count
                    self.unresolved_count += disposition == 'retained_unresolved'
                    self.input_bytes += size
                    refs.append(AlleleReference(key, native))
            require(self.retained_count + self.excluded_count == self.original_count)
            self.seed = replace(closed.seed, allele_refs=tuple(refs))
            digest = writer.add_sequence(name + '.members', 'population-work-member-decisions-v1', self.decisions())
            self.header = dict(schema='genosapiens.context-population-work-filter', version=1,
                context_ordinal=closed.ordinal, context_id=closed.context_id, original_seed=closed.seed,
                original_genotype_inventory_sha256=self.original_genotype_inventory,
                resources=resources.reference, assessment_domain_sha256=resources.assessment_domain_sha256,
                policy_id=policy.policy_id, policy_sha256=policy.sha256,
                protection_resource_bindings=source_access.resource_bindings,
                source_manifest_sha256=resources.source_manifest_sha256,
                reference_manifest_sha256=resources.reference_binding.manifest_sha256,
                original_member_count=self.original_count, retained_member_count=self.retained_count,
                excluded_member_count=self.excluded_count, unresolved_member_count=self.unresolved_count,
                original_member_occurrences=self.original_occurrences, retained_member_occurrences=self.retained_occurrences,
                excluded_member_occurrences=self.original_occurrences-self.retained_occurrences,
                retained_input_bytes=self.input_bytes, rich_excluded_hypothesis_count=None,
                skipped_scoring='not_evaluated', members=dict(name=name+'.members',
                    kind='population-work-member-decisions-v1', sha256=digest))
            self.sha256 = writer.add(name, 'context-population-work-filter-v1', self.header)
        except BaseException:
            writer.failed = True
            self.close()
            raise

    def _assess(self, member, models):
        from genosapiens.wgs_admission import iter_admission_source_rows, _occurrence
        from genosapiens.wgs_context_evidence import EnrichedContextEvidence
        access = self.source_access
        require(member.parent_working_set_manifest_sha256 == access.prepared_index.working_manifest_sha256
            and member.native_source_contract_sha256 == access.routing.contract.sha256,
            'population work source parent binding changed')
        keys = tuple(ContextKey.from_context(model) for model in models)
        base = dict(contexts=keys, occurrences=(), rescue_assessments=(), population=())
        if access.screening.merged_intake is not None:
            return dict(base, disposition='retained_unresolved', reasons=('merged_native_occurrence_scope_retained',))
        if not keys or any(key is None for key in keys):
            return dict(base, disposition='retained_unresolved', reasons=('model_context_unresolved',))
        with closing(iter_admission_source_rows(access.source_connection,
                ordinals=(member.retained_ordinal,), execution_config=access.config)) as rows:
            source = tuple(rows)
        require(len(source) == 1, 'population work source inventory incomplete')
        identity, payloads = source[0]
        require(identity[1] == member.identity_id and tuple(key for key, _ in payloads) == member.source_allele_ids,
            'population work source occurrence binding changed')
        normalized = tuple(i.identity_id for i in member.allele.identities if i.stage is IdentityStage.NORMALIZED)
        if len(normalized) != 1:
            return dict(base, disposition='retained_unresolved', reasons=('normalized_member_identity_unresolved',))
        requests = (PopulationRequest(member.allele, normalized[0]),)
        result = (self.population_session.query(requests) if self.query_cache is None else
            self.query_cache.population_results(requests, session=self.population_session))
        EnrichedContextEvidence._check_enrichment({member.allele.allele_id: member.allele}, result.alleles, 'population')
        enriched, = result.alleles
        records = population_records(enriched)
        base['population'] = records
        base['population_query_report'] = result.report
        if not optional_frequency_support(records).qualified:
            return dict(base, disposition='retained_unresolved', reasons=('frequency_exclusion_not_established',))
        proofs, retained = [], []
        bindings = dict(access.resource_bindings)
        for source_id, payload in payloads:
            facts = _occurrence(payload, source_allele_id=source_id,
                source_manifest=access.source_manifest, contract=access.routing.contract)
            decision, occurrence, protection = evaluate_source_candidate(value=facts, routing=access.routing,
                context=keys[0], records=records, population_allele_id=member.identity_id,
                population_reference_sha256=self.resources.reference_binding.manifest_sha256,
                protection_session=_ProtectionQueries(access.public_session), protection_policy=access.common_policy,
                clinical_manifest_sha256=bindings['clinical_manifest'], exception_manifest_sha256=bindings['exception_manifest'])
            decisions = [decision]
            decisions.extend(evaluate_candidate(context=key, occurrence=occurrence, records=records,
                protection=protection, protection_policy=access.common_policy) for key in keys[1:])
            proofs.append(dict(source_allele_id=source_id, payload_sha256=facts.payload_sha256,
                occurrence=occurrence, protection=protection, decisions=tuple(decisions)))
            retained.extend(reason for d in decisions if not d.excluded for reason in d.reasons)
        base['occurrences'] = tuple(proofs)
        if retained:
            return dict(base, disposition='retained_unresolved', reasons=tuple(sorted(set(retained))))
        if self.clinical_annotator is None:
            return dict(base, disposition='retained_unresolved', reasons=('context_clinical_rescue_unassessed',))
        clinical = self.clinical_annotator((enriched,))
        EnrichedContextEvidence._check_enrichment({enriched.allele_id: enriched}, clinical.alleles, 'clinical')
        enriched, = clinical.alleles
        # Retain only the existing rescue kernel's inputs, without transcripts,
        # call payloads, hypothesis objects or score ledgers.
        base['rescue_input'] = AlleleEvidence(enriched.allele_id, enriched.native_identity_id,
            enriched.identities, population=enriched.population, clinical=enriched.clinical)
        base['clinical_query'] = dict(coverage=clinical.coverage, bindings=clinical.bindings,
            public_sourced_evidence=clinical.public_sourced_evidence)
        assessed = []
        for key in keys:
            assessment, = self.resources.population.evaluate(key, (enriched,))
            assessed.append(dict(context=key, assessment=assessment))
            for proof in proofs:
                retained.extend(context_rescue_vetoes(assessment,
                    protection=proof['protection'], protection_policy=access.common_policy,
                    clinical_evidence=tuple(enriched.clinical or ())))
        base['rescue_assessments'] = tuple(assessed)
        return dict(base, disposition='retained_unresolved' if retained else 'excluded',
            reasons=tuple(sorted(set(retained))) if retained else ('complete_protected_population_exclusion',))

    @property
    def reference(self):
        self.require_owner(self.writer, self.original)
        return dict(name=self.name, kind='context-population-work-filter-v1', sha256=self.sha256)

    def require_owner(self, writer, closed):
        require(not self.closed and self.writer is writer and self.original is closed)

    def retains(self, member):
        self.require_owner(self.writer, self.original)
        row = self.connection.execute(f'SELECT ordinal,native,disposition FROM {self.table} WHERE id=?',
            (member.allele.allele_id,)).fetchone()
        require(row is not None and row[:2] == (member.retained_ordinal, member.allele.native_identity_id))
        return row[2] != 'excluded'

    def decisions(self):
        self.require_owner(self.writer, self.original)
        for content, in self.connection.execute(f'SELECT content FROM {self.table} ORDER BY ordinal'):
            yield json.loads(content)

    def population_assessment(self, key, allele, policy):
        """Reuse a live preconstruction result only for its exact inputs."""
        from genosapiens.evidence_serialization import _typed
        from genosapiens.population_policy import MemberAssessment
        self.require_owner(self.writer, self.original)
        require(policy is self.resources.population)
        row = self.connection.execute(f'SELECT content FROM {self.table} WHERE id=?', (allele.allele_id,)).fetchone()
        require(row is not None)
        value = json.loads(row[0])
        require(value['disposition'] != 'excluded')
        expected = AlleleEvidence(allele.allele_id, allele.native_identity_id, allele.identities,
            population=allele.population, clinical=allele.clinical)
        if _canonical(value.get('rescue_input')) != _canonical(expected):
            return None
        for saved in value['rescue_assessments']:
            if _canonical(saved['context']) == _canonical(key):
                return _typed(saved['assessment'], MemberAssessment, decode=True)
        return None

    def iter_public_records(self):
        """Keep query provenance needed by the original admission decisions."""
        from genosapiens.evidence_serialization import _typed
        from genosapiens.evidence_model import PublicSourcedEvidence
        from genosapiens.wgs_contexts import ClosedPublicRecord
        for row in self.decisions():
            for value in row.get('clinical_query', {}).get('public_sourced_evidence', ()):
                source = _typed(value, PublicSourcedEvidence, decode=True)
                if not source.records:
                    yield ClosedPublicRecord(source.evidence_id, source.public_manifest_sha256, None)
                for record in source.records.values():
                    yield ClosedPublicRecord(source.evidence_id, source.public_manifest_sha256, record)

    def close(self):
        if not self.closed:
            self.connection.execute(f'DROP TABLE IF EXISTS {self.table}')
            self.closed = True

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()


def load_work_selection(reader, reference, owner, resources):
    """Authenticate compact dispositions while restoring a finite diagnostic.

    Recompute the scheduling predicate from its saved inputs. Population rescue
    assessments retain their existing kernel values; completed source/query
    authority remains bound to the original producer and resource closure.
    """
    from genosapiens.candidate_admission import CandidateContextSeed
    from genosapiens.common_allele_filter import (
        PrefilterOccurrence, ClinicalProtectionFacts, load_common_allele_filter_policy,
    )
    from genosapiens.evidence_serialization import _typed
    from genosapiens.evidence_model import _sha
    from genosapiens.population_filter_candidate import WORK_POLICY_ID, WORK_POLICY_SHA256, load_population_work_policy
    from genosapiens.population_policy import MemberAssessment
    from genosapiens.genotype_properties import AlleleGenotypeProperties
    from genosapiens.wgs_pair_domain import MAX_CONTEXT_MEMBERS
    from genosapiens.wgs_case_resources import bound_reference
    from pathlib import Path
    value = reader.reconstruct(bound_reference(reader, reference, 'context-population-work-filter-v1'))
    require(resources is not None and value['resources'] == resources.reference
        and value['assessment_domain_sha256'] == resources.assessment_domain_sha256
        and value['source_manifest_sha256'] == resources.source_manifest_sha256
        and value['reference_manifest_sha256'] == resources.reference_binding.manifest_sha256
        and resources.parent_bindings.get('population_work_policy') == value['policy_sha256'])
    require(set(value) == {'schema', 'version', 'context_ordinal', 'context_id', 'original_seed',
        'original_genotype_inventory_sha256', 'resources', 'assessment_domain_sha256', 'policy_id',
        'policy_sha256', 'protection_resource_bindings', 'source_manifest_sha256', 'reference_manifest_sha256',
        'original_member_count', 'retained_member_count', 'excluded_member_count', 'unresolved_member_count',
        'original_member_occurrences', 'retained_member_occurrences', 'excluded_member_occurrences',
        'retained_input_bytes', 'rich_excluded_hypothesis_count', 'skipped_scoring', 'members'})
    require(value['schema'] == 'genosapiens.context-population-work-filter' and type(value['version']) is int
        and value['version'] == 1 and value['policy_id'] == WORK_POLICY_ID and value['policy_sha256'] == WORK_POLICY_SHA256)
    for key in ('source_manifest_sha256', 'reference_manifest_sha256', 'assessment_domain_sha256',
                'original_genotype_inventory_sha256'):
        _sha(value[key])
    original = _typed(value['original_seed'], CandidateContextSeed, decode=True)
    for key in ('context_ordinal', 'original_member_count', 'retained_member_count', 'excluded_member_count',
                'unresolved_member_count', 'original_member_occurrences', 'retained_member_occurrences',
                'excluded_member_occurrences', 'retained_input_bytes'):
        require(type(value[key]) is int and value[key] >= 0)
    require(value['original_member_count'] <= MAX_CONTEXT_MEMBERS)
    require(value['context_id'] == original.context.context_id == owner.context_id and value['context_ordinal'] == owner.ordinal)
    bindings = dict(value['protection_resource_bindings'])
    require(len(bindings) == len(value['protection_resource_bindings']) and set(bindings) == {
        'population_manifest', 'clinical_manifest', 'exception_manifest', 'precision_manifest',
        'all_sites_manifest', 'public_index_manifest'})
    for digest in bindings.values():
        _sha(digest)
    # The canonical policy's dependency is also bound in producer resources.
    content = (Path(__file__).resolve().parents[2] / 'config/track1/common-allele-prefilter-v1.json').read_bytes()
    policy = load_common_allele_filter_policy(content, hashlib.sha256(content).hexdigest())
    active = load_population_work_policy((Path(__file__).resolve().parents[2] /
        'config/track1/protected-population-work-filter-v1.json').read_bytes(), expected_sha256=WORK_POLICY_SHA256)
    expected_common = json.loads(active.content)['protection']['policy']['sha256']
    require(hashlib.sha256(content).hexdigest() == expected_common)
    table = owner.prefix + '_work_selection'
    owner.tables.append(table)
    reader.connection.execute(f'CREATE TABLE {table} (id TEXT PRIMARY KEY,ordinal INTEGER UNIQUE,native TEXT,input_sha TEXT,size INTEGER,retained INTEGER,content TEXT)')
    count = retained = excluded = unresolved = occurrences = retained_occurrences = retained_bytes = 0
    prior = 0
    refs, retained_refs = [], []
    for row in reader.iter_array(bound_reference(reader, value['members'], 'population-work-member-decisions-v1'),
            max_item_bytes=MAX_DECISION_BYTES):
        count += 1
        expected = {'contexts', 'occurrences', 'rescue_assessments', 'population', 'allele_id',
            'model_input', 'genotype_properties',
            'native_identity_id', 'identity_id', 'original_ordinal', 'source_allele_ids',
            'original_input_sha256', 'original_input_bytes', 'disposition', 'reasons'}
        require(set(row) in (expected, expected | {'population_query_report'},
            expected | {'population_query_report', 'rescue_input', 'clinical_query'}))
        require(type(row['original_ordinal']) is int and row['original_ordinal'] > prior
            and type(row['original_input_bytes']) is int and row['original_input_bytes'] > 0)
        prior = row['original_ordinal']
        _sha(row['original_input_sha256'])
        require(row['disposition'] in {'excluded', 'retained', 'retained_unresolved'}
            and type(row['source_allele_ids']) is list and row['source_allele_ids']
            and len(row['source_allele_ids']) == len(set(row['source_allele_ids'])))
        ref = AlleleReference(row['allele_id'], row['native_identity_id'])
        model_input = _typed(row['model_input'], AlleleEvidence, decode=True)
        properties = _typed(row['genotype_properties'], AlleleGenotypeProperties, decode=True)
        require(model_input.allele_id == properties.allele_id == ref.allele_id
            and model_input.native_identity_id == properties.native_identity_id == ref.identity_id
            and properties.source_allele_ids == tuple(row['source_allele_ids'])
            and not model_input.calls and not model_input.population and not model_input.clinical
            and not model_input.transcript_consequences and not model_input.mechanism_ids and not model_input.evidence_refs)
        expected_contexts = tuple(ContextKey.from_context(model) for model in
            consuming_model_contexts(original.context, model_input, properties))
        require(_canonical(expected_contexts) == _canonical(row['contexts']), 'population work consuming models changed')
        refs.append(ref)
        occurrences += len(row['source_allele_ids'])
        keep = row['disposition'] != 'excluded'
        reader.connection.execute(f'INSERT INTO {table} VALUES (?,?,?,?,?,?,?)',
            (row['allele_id'], prior, row['native_identity_id'], row['original_input_sha256'], row['original_input_bytes'], int(keep), _canonical(row)))
        if keep:
            retained += 1
            unresolved += row['disposition'] == 'retained_unresolved'
            retained_refs.append(ref)
            retained_occurrences += len(row['source_allele_ids'])
            retained_bytes += row['original_input_bytes']
            continue
        excluded += 1
        require(row['reasons'] == ['complete_protected_population_exclusion'])
        contexts = tuple(_typed(key, ContextKey, decode=True) for key in row['contexts'])
        require(contexts and len(contexts) == len(set(contexts))
            and tuple(proof['source_allele_id'] for proof in row['occurrences']) == tuple(row['source_allele_ids'])
            and len(row['rescue_assessments']) == len(contexts))
        rescue_input = _typed(row['rescue_input'], AlleleEvidence, decode=True)
        require(rescue_input.allele_id == row['allele_id'] and rescue_input.native_identity_id == row['native_identity_id']
            and rescue_input.identities == model_input.identities
            and rescue_input.calls is None and rescue_input.transcript_consequences is None
            and not rescue_input.mechanism_ids and not rescue_input.evidence_refs
            and _canonical(population_records(rescue_input)) == _canonical(row['population']))
        for proof in row['occurrences']:
            _sha(proof['payload_sha256'])
            occurrence = _typed(proof['occurrence'], PrefilterOccurrence, decode=True)
            protection = _typed(proof['protection'], ClinicalProtectionFacts, decode=True)
            require(occurrence.source_allele_id == proof['source_allele_id']
                and occurrence.normalized_allele_id == row['identity_id']
                and protection.clinical_manifest_sha256 == bindings['clinical_manifest']
                and protection.exception_manifest_sha256 == bindings['exception_manifest'])
            decisions = tuple(evaluate_candidate(context=key, occurrence=occurrence, records=row['population'],
                protection=protection, protection_policy=policy) for key in contexts)
            require(all(d.excluded for d in decisions) and _canonical(decisions) == _canonical(proof['decisions']))
            for key, saved in zip(contexts, row['rescue_assessments'], strict=True):
                require(_typed(saved['context'], ContextKey, decode=True) == key)
                assessment = _typed(saved['assessment'], MemberAssessment, decode=True)
                expected_assessment, = resources.population.evaluate(key, (rescue_input,))
                require(assessment.allele_id == row['allele_id'] and assessment.native_identity_id == row['native_identity_id']
                    and assessment == expected_assessment
                    and not context_rescue_vetoes(assessment, protection=protection, protection_policy=policy,
                        clinical_evidence=tuple(rescue_input.clinical or ())))
    require(tuple(refs) == original.allele_refs and replace(original, allele_refs=tuple(retained_refs)) == owner.seed
        and (count, retained, excluded, unresolved, occurrences, retained_occurrences, retained_bytes) == (
            value['original_member_count'], value['retained_member_count'], value['excluded_member_count'],
            value['unresolved_member_count'], value['original_member_occurrences'], value['retained_member_occurrences'],
            value['retained_input_bytes'])
        and occurrences-retained_occurrences == value['excluded_member_occurrences']
        and retained == owner.member_count and retained_bytes == owner.input_bytes
        and value['rich_excluded_hypothesis_count'] is None and value['skipped_scoring'] == 'not_evaluated')
    require(not _requires_pairs(original.context) or excluded in (0, count), 'mixed pair omission lacks exact partner scope')
    return value


def work_domain(members):
    """Compact evaluated/original-domain binding for downstream products.

    Counts are context-member occurrences, not distinct case-wide source
    records. No count or score is invented for skipped rich hypotheses.
    """
    if getattr(members, 'work_filter', None) is None:
        return {}
    return _work_domain(members.work_filter, members.work_filter_header)


def _work_domain(reference, value):
    keys = ('policy_id', 'policy_sha256', 'original_member_count', 'retained_member_count',
        'excluded_member_count', 'unresolved_member_count', 'original_member_occurrences',
        'retained_member_occurrences', 'excluded_member_occurrences', 'rich_excluded_hypothesis_count',
        'skipped_scoring')
    return {'population_work': {'evidence': reference, **{key: value[key] for key in keys}}}


def verify_saved_work_domain(reader, context, row):
    """Bind compact downstream counts without restoring excluded science."""
    from genosapiens.wgs_case_resources import bound_reference
    from genosapiens.population_filter_candidate import WORK_POLICY_SHA256
    members = reader.reconstruct(bound_reference(reader, context['members'], 'context-evidence-v1'))
    if 'work_filter' not in members:
        require('population_work' not in row)
        return
    reference = members['work_filter']
    value = reader.reconstruct(bound_reference(reader, reference, 'context-population-work-filter-v1'))
    require(value['policy_sha256'] == WORK_POLICY_SHA256
        and value['assessment_domain_sha256'] == context['assessment_domain_sha256']
        and value['resources'] == context['resources']
        and value['context_ordinal'] == context['ordinal']
        and value['context_id'] == context['context_id']
        and value['retained_member_count'] == context['member_count']
        and _work_domain(reference, value).get('population_work') == row.get('population_work'))
