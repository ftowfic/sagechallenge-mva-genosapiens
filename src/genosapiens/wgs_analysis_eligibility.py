"""Explicit assessment-domain policy; never a genotype/carriage certificate.

The original saved inputs, GT, source routing and read obligations stay intact.
This gate schedules assessment or preserves unresolved reference-call rescue.
"""
from collections.abc import Mapping
from dataclasses import asdict, dataclass, field
import hashlib
import json
import re
from functools import lru_cache
from enum import StrEnum
from types import MappingProxyType

from genosapiens.allele_identity import SequenceAllele, MosaicAllele
from genosapiens.evidence_model import IdentityStage, _typed, _typed_value, _validation_session
from genosapiens.vcf_qc import VcfQcStreamHeader
from genosapiens.wgs_admission import LedgerAdmissionInput, SourceTechnicalObligation
from genosapiens.common_allele_filter import PrefilterDecision
from genosapiens.evidence_serialization import _plain
from genosapiens.wgs_ledger import MAX_ROW_BYTES

RULE = 'deepvariant-r1.0-unconfirmed-homref-v1'
HEADER_FILTER = '##FILTER=<ID=RefCall,Description="Genotyping model thinks this site is reference.">'
DISPOSITIONS = frozenset(('assess', 'exclude_common', 'defer_refcall'))


class AnalysisEligibilityError(ValueError):
    """Static integrity diagnostic, with no original input values."""


def _sha(value):
    return type(value) is str and re.fullmatch('[0-9a-f]{64}', value) is not None


def _canonical(value):
    return (json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n').encode('ascii')


def _unique(pairs):
    result = {}
    for key, value in pairs:
        if key in result: raise AnalysisEligibilityError('duplicate eligibility field')
        result[key] = value
    return result


class OptionalAnalysisStatus(StrEnum):
    """Execution state, separate from whether a completed analysis found a match."""

    NOT_RUN = 'not_run'
    RUNNING = 'running'
    COMPLETE = 'complete'
    FAILED = 'failed'


# A single method identity binds screening, context closure and publication.
# This is descriptive provenance, never a selector for another implementation.
METHOD_BYTES = _canonical({
    'schema': 'genosapiens.wgs-analysis-method',
    'version': 1,
    'source_inventory': 'all_original_retained_identities_and_occurrences',
    'producer_validation': 'complete_inputs_annotations_and_source_proofs',
    'consumer_validation': 'complete_compact_coverage_and_requested_rich_evidence',
    'assessment_domain': 'supported_associations_of_assessed_inputs',
    'context_closure': 'global_nonempty_assessed_contexts',
    'pair_domain': 'all_relevant_distinct_source_events_in_closed_assessed_contexts',
    'outside_domain': 'unresolved_without_no_partner_certificate',
    'primary_readiness': 'primary_accounting_evidence_order_selection_selected_checks_and_packet',
    'optional_products': 'separate_immutable_publications_bound_to_primary',
    'optional_initial_status': 'not_run',
})
METHOD_SHA256 = hashlib.sha256(METHOD_BYTES).hexdigest()


def scientific_assessment_domain(component,parents,profile_component=None):
    """Current membership definition and its explicitly supplied provenance.

    Isolated scientific owners may have no source-screening publication. They
    retain an unbound original denominator, never a fabricated all-input count.
    Executable WGS builders always provide the authenticated component.
    """
    counts=None if component is None else dict(component['input_dispositions'])
    if counts is not None and (set(counts)!=DISPOSITIONS
            or any(type(v) is not int or v<0 for v in counts.values())
            or sum(counts.values())!=component['input_rows']):
        raise AnalysisEligibilityError('invalid scientific assessment domain')
    result=dict(schema='genosapiens.screened-scientific-domain',version=1,method_sha256=METHOD_SHA256,
        membership='supported_associations_of_assessed_inputs',
        pairs='all_relevant_distinct_source_events_in_closed_assessed_contexts',
        outside_domain='unresolved_without_no_partner_certificate',
        original_inputs=None if component is None else component['input_rows'],input_dispositions=counts,
        eligibility_manifest_sha256=None if component is None else component['manifest_sha256'],
        association_projection_sha256=parents.get('association_projection'),
        assessment_projection_sha256=parents.get('assessment_projection'))
    if profile_component is not None:
        from genosapiens.wgs_profile_projection import profile_component_parents, _cohort
        expected=profile_component_parents(profile_component)
        scope=profile_component['counts']; original=profile_component['identity']['parent']
        cohort=_cohort(profile_component['identity'])
        profile_domain=original if cohort is None else cohort
        if (component is None or any(parents.get(k)!=v for k,v in expected.items())
                or original['rows']!=component['input_rows']
                or original['occurrences']!=component['source_occurrences']
                or scope['rows']!=profile_domain['rows']
                or cohort is None and scope['occurrences']!=original['occurrences']
                or original['input_dispositions']!=counts
                or original['occurrence_dispositions']!=dict(component['occurrence_dispositions'])
                or sum(scope['dispositions'].values())!=scope['rows']
                or scope['advancing_rows']!=sum(scope['dispositions'][k] for k in ('selected','unresolved'))):
            raise AnalysisEligibilityError('profile scientific assessment domain changed')
        from genosapiens.track1_association_scope import POLICY_SHA256 as association_scope_sha256
        result.update(version=3,membership='supported_applicable_associations_of_profile_advancing_inputs',
            association_scope_sha256=association_scope_sha256,
            public_evidence='selected-association-source-closure-with-compact-exclusions',
            analysis_profile=profile_component['identity']['profile'],
            profile_manifest_sha256=profile_component['manifest_sha256'],
            profile_dispositions=dict(scope['dispositions']),advancing_inputs=scope['advancing_rows'],
            advancing_occurrences=scope['advancing_occurrences'],
            advancing_ordinals_sha256=scope['advancing_ordinals_sha256'])
        if cohort is not None:
            result['finite_profile_cohort']=dict(rows=cohort['rows'],occurrences=scope['occurrences'],
                ordered_ordinals_sha256=cohort['ordered_ordinals_sha256'],
                acceptance_sha256=cohort['acceptance_sha256'],full_parent_qualification=False)
    elif any(k.startswith('analysis_profile_') for k in parents):
        raise AnalysisEligibilityError('missing scientific profile domain')
    if 'population_work_policy' in parents:
        from genosapiens.population_filter_candidate import WORK_POLICY_SHA256
        if parents['population_work_policy'] != WORK_POLICY_SHA256:
            raise AnalysisEligibilityError('population work policy changed')
        result.update(version=4, population_work_policy_sha256=WORK_POLICY_SHA256,
            membership='retained_context_members_after_protected_population_work_filter',
            original_membership='supported_applicable_associations_of_profile_advancing_inputs' if profile_component
                else 'supported_associations_of_assessed_inputs',
            pairs='all_relevant_distinct_source_events_in_retained_contexts_with_conservative_partner_obligations',
            skipped_scientific_work='not_evaluated;compact_original_member_decisions_retained')
    return result


@dataclass(frozen=True)
class AnalysisEligibilityPolicy:
    common_filter_policy_sha256: str
    schema: str = 'genosapiens.wgs-analysis-eligibility-policy'
    version: int = 2
    method_sha256: str = METHOD_SHA256
    noncontribution_proof_version: int = 1
    refcall_rule_id: str = RULE
    phase_protection: bool = True
    read_rescue_mode: str = 'existing-interface-only'

    def __post_init__(self):
        if (not _sha(self.common_filter_policy_sha256)
                or self.schema != 'genosapiens.wgs-analysis-eligibility-policy'
                or type(self.version) is not int or self.version != 2
                or self.method_sha256 != METHOD_SHA256
                or type(self.noncontribution_proof_version) is not int or self.noncontribution_proof_version != 1
                or self.refcall_rule_id != RULE or self.phase_protection is not True
                or self.read_rescue_mode != 'existing-interface-only'):
            raise AnalysisEligibilityError('invalid analysis eligibility policy')

    def as_bytes(self): return _canonical(asdict(self))

    @property
    def sha256(self): return hashlib.sha256(self.as_bytes()).hexdigest()


def load_analysis_eligibility_policy(content, expected_sha256):
    try:
        if type(content) is not bytes or not 0 < len(content) <= 4096 or not _sha(expected_sha256):
            raise ValueError
        if hashlib.sha256(content).hexdigest() != expected_sha256: raise ValueError
        data = json.loads(content, object_pairs_hook=_unique)
        if type(data) is not dict or set(data) != set(AnalysisEligibilityPolicy.__dataclass_fields__): raise ValueError
        policy = AnalysisEligibilityPolicy(**data)
        if policy.as_bytes() != content: raise ValueError
        return policy
    except (TypeError, ValueError, UnicodeError):
        raise AnalysisEligibilityError('invalid analysis eligibility policy') from None


@dataclass(frozen=True, slots=True)
class AdmissionScreeningContext:
    header: VcfQcStreamHeader
    selected_sample: str
    policy: AnalysisEligibilityPolicy
    merged_intake: object = None
    _digest: str = field(init=False, repr=False, compare=False)

    def __post_init__(self):
        if (type(self.header) is not VcfQcStreamHeader
                or type(self.policy) is not AnalysisEligibilityPolicy
                or type(self.selected_sample) is not str or not self.selected_sample):
            raise AnalysisEligibilityError('invalid producer screening context')
        self.header.__post_init__(); self.policy.__post_init__()
        if self.merged_intake is None:
            if self.selected_sample not in self.header.sample_names:
                raise AnalysisEligibilityError('screening selected sample is absent')
        else:
            from genosapiens.wgs_read_rescue import MergedIntakeManifest
            if (type(self.merged_intake) is not MergedIntakeManifest
                    or self.selected_sample != self.merged_intake.sample_id
                    or self.header not in tuple(v.source_header for v in self.merged_intake.source_bindings)):
                raise AnalysisEligibilityError('invalid merged screening context')
            self.merged_intake.__post_init__()
        object.__setattr__(self, '_digest', hashlib.sha256(self.as_bytes()).hexdigest())

    def as_bytes(self):
        return _canonical(dict(header=_plain(self.header), selected_sample=self.selected_sample,
            policy=json.loads(self.policy.as_bytes()), merged_intake=_plain(self.merged_intake)))

    @property
    def sha256(self):
        return self._digest

    def derive(self, value):
        return derive_admission_screening_facts(value, header=self.header,
            selected_sample=self.selected_sample, merged_intake=self.merged_intake,
            policy=self.policy)


def load_admission_screening_context(content):
    try:
        if type(content) is not bytes or not 0 < len(content) <= 16*1024*1024:
            raise ValueError
        data = json.loads(content, object_pairs_hook=_unique)
        if type(data) is not dict or set(data) != {'header', 'selected_sample', 'policy', 'merged_intake'}:
            raise ValueError
        policy = _canonical(data['policy'])
        merged = data['merged_intake']
        if merged is not None:
            from genosapiens.wgs_read_rescue import MergedIntakeManifest
            merged = _typed(merged, MergedIntakeManifest, decode=True)
        value = AdmissionScreeningContext(_typed(data['header'], VcfQcStreamHeader, decode=True),
            data['selected_sample'], load_analysis_eligibility_policy(policy, hashlib.sha256(policy).hexdigest()), merged)
        if value.as_bytes() != content:
            raise ValueError
        return value
    except (ValueError, TypeError, KeyError, UnicodeError):
        raise AnalysisEligibilityError('invalid serialized producer screening context') from None


@dataclass(frozen=True)
class RefCallEligibility:
    can_defer: bool
    reasons: tuple[str, ...]
    source_allele_ids: tuple[str, ...]
    header_sha256: str


@dataclass(frozen=True, slots=True)
class AdmissionScreeningFacts:
    """Producer-derived facts, bound to canonical rich evidence by the store.

    Construction is structural validation only. Source/proof replay belongs to
    the producer, and neither this type nor its hashes confer producer authority.
    The exact source-call projection retains GT, phase, ALT and technical facts;
    transcript graphs remain in the independently bound rich payload.
    """
    retained_ordinal: int
    identity_id: str
    allele_id: str | None
    subject_id: str
    source_allele_ids: tuple[str, ...]
    required_rule_ids: tuple[str, ...]
    transport: str
    source_manifest_sha256: str
    working_manifest_sha256: str
    source_routing_sha256: str
    header_sha256: str
    selected_sample: str
    reference_manifest_sha256s: tuple[str, ...]
    normalized_identity_ids: tuple[str, ...]
    sequence_identities: tuple[tuple[str, SequenceAllele], ...]
    gene_ids: tuple[str, ...]
    technical_obligations: tuple[SourceTechnicalObligation, ...]
    source_call_facts: bytes
    refcall_can_defer: bool
    refcall_reasons: tuple[str, ...]
    merged_intake_sha256: str | None
    eligibility_policy_sha256: str
    method_sha256: str = METHOD_SHA256

    def __post_init__(self):
        try:
            calls = json.loads(self.source_call_facts, object_pairs_hook=_unique)
            if _canonical(calls) != self.source_call_facts: raise ValueError
        except (ValueError, TypeError, UnicodeError):
            raise AnalysisEligibilityError('invalid screening source projection') from None
        self._validate(calls)

    def _validate(self, calls):
        """Shared semantic checks for constructors and the owned JSON decoder."""
        if (type(self.retained_ordinal) is not int or self.retained_ordinal < 1
                or any(type(v) is not str or not v for v in (self.identity_id, self.subject_id, self.selected_sample))
                or self.allele_id is not None and (type(self.allele_id) is not str or not self.allele_id)
                or self.transport not in ('sequence_vep', 'native_vep', 'repeat_catalog', 'unresolved')
                or any(not _sha(v) for v in (self.source_manifest_sha256, self.working_manifest_sha256,
                    self.source_routing_sha256, self.header_sha256, self.method_sha256, self.eligibility_policy_sha256))
                or self.method_sha256 != METHOD_SHA256
                or self.merged_intake_sha256 is not None and not _sha(self.merged_intake_sha256)
                or type(self.refcall_can_defer) is not bool
                or self.merged_intake_sha256 is not None and self.refcall_can_defer
                or type(self.source_allele_ids) is not tuple or not self.source_allele_ids
                or len(set(self.source_allele_ids)) != len(self.source_allele_ids)
                or type(self.technical_obligations) is not tuple
                or any(type(v) is not SourceTechnicalObligation for v in self.technical_obligations)
                or tuple(v.source_allele_id for v in self.technical_obligations) != self.source_allele_ids
                or any(type(v) is not str or not v for rows in (self.source_allele_ids, self.required_rule_ids,
                    self.normalized_identity_ids, self.gene_ids, self.refcall_reasons) for v in rows)
                or not self.refcall_reasons
                or any(not _sha(v) for v in self.reference_manifest_sha256s)
                or any(type(rows) is not tuple or rows != tuple(sorted(set(rows))) for rows in (
                    self.reference_manifest_sha256s, self.normalized_identity_ids, self.gene_ids, self.refcall_reasons))
                or type(self.source_call_facts) is not bytes or not 0 < len(self.source_call_facts) <= MAX_ROW_BYTES):
            raise AnalysisEligibilityError('invalid admission screening facts')
        _typed(self.sequence_identities, tuple[tuple[str, SequenceAllele], ...])
        if (any(not key for key, _ in self.sequence_identities)
                or tuple(key for key, _ in self.sequence_identities)
                    != tuple(sorted({key for key, _ in self.sequence_identities}))):
            raise AnalysisEligibilityError('invalid screening sequence identities')
        try:
            if type(calls) is not list:
                raise ValueError
            if calls and {row['source_allele_id'] for row in calls} != set(self.source_allele_ids):
                raise ValueError
        except (ValueError, TypeError, KeyError):
            raise AnalysisEligibilityError('invalid screening source projection') from None

    def as_bytes(self):
        content = _canonical(_screening_facts_mapping(self, json.loads(self.source_call_facts)))
        if len(content) > MAX_ROW_BYTES:
            raise AnalysisEligibilityError('screening facts exceed row capacity')
        return content


def derive_admission_screening_facts(value, *, header, selected_sample, policy, merged_intake=None):
    """Derive facts while the producer owns the full strictly validated row.

    Ordinary inputs use the existing single RefCall predicate. Merged inputs
    retain assessment after their declared original/delta parents authenticate;
    the owning MergedIntakeIndex must have replayed the complete intake first.
    """
    value = _typed(value, LedgerAdmissionInput)
    if type(policy) is not AnalysisEligibilityPolicy:
        raise AnalysisEligibilityError('invalid screening policy')
    policy.__post_init__()
    if type(header) is not VcfQcStreamHeader:
        raise AnalysisEligibilityError('invalid screening header')
    _supported_refcall_header(header)
    if merged_intake is None:
        refcall = evaluate_refcall_eligibility(value, header=header, selected_sample=selected_sample)
    else:
        from genosapiens.wgs_read_rescue import MergedIntakeManifest, _bytes
        if (type(merged_intake) is not MergedIntakeManifest
                or value.parent_source_manifest_sha256 != merged_intake.sha256
                or value.parent_working_set_manifest_sha256 != merged_intake.sha256
                or value.native_source_contract_sha256 != hashlib.sha256(_bytes(merged_intake.source_bindings)).hexdigest()
                or selected_sample != merged_intake.sample_id):
            raise AnalysisEligibilityError('foreign merged screening parent')
        merged_intake.__post_init__()
        refcall = RefCallEligibility(False, ('merged_input_assessment_domain',),
            value.source_allele_ids, header.header_sha256)
    mappings = () if value.allele is None else value.allele.identities
    if value.unresolved_source is not None:
        source=value.unresolved_source
        allowed_inputs=({header.input_sha256} if merged_intake is None else
            {binding.input_sha256 for binding in merged_intake.source_bindings})
        if value.source_allele_ids != (source.source_allele_id,) or source.input_sha256 not in allowed_inputs:
            raise AnalysisEligibilityError('unresolved screening source binding mismatch')
    source_mappings = {m.identity_id: m for m in mappings if m.stage is IdentityStage.SOURCE}
    obligations = {o.source_allele_id: o for o in value.technical_obligations}
    calls = []
    for call in (() if value.allele is None else value.allele.calls or ()):
        fields = call.source.source_fields
        source_id = fields.get('source_allele_id')
        if source_id not in obligations or source_id not in source_mappings:
            raise AnalysisEligibilityError('foreign screening source occurrence')
        obligation = obligations[source_id]
        source_header, source_sample = header, selected_sample
        if merged_intake is not None:
            binding = next((v for v in merged_intake.source_bindings
                if v.input_sha256 == fields.get('source_input_sha256')), None)
            if binding is not None:
                source_header, source_sample = binding.source_header, binding.selected_sample_name
                if fields.get('source_manifest_sha256') != binding.source_manifest_sha256:
                    raise AnalysisEligibilityError('foreign merged screening source')
            elif fields.get('authorized_read_evidence') is None:
                raise AnalysisEligibilityError('undeclared merged screening source')
            mapping = fields.get('merged_intake_sample')
            if (not isinstance(mapping, Mapping) or mapping.get('manifest_sha256') != merged_intake.sha256
                    or mapping.get('sample_id') != selected_sample
                    or mapping.get('source_sample_name') != call.sample.sample_name):
                raise AnalysisEligibilityError('foreign merged screening sample')
        if fields.get('authorized_read_evidence') is None:
            _validate_screening_source_call(call, obligation, source_mappings[source_id],
                header=source_header, selected_sample=source_sample)
        else:
            if merged_intake is None:
                raise AnalysisEligibilityError('read addition requires its merged owner')
            from genosapiens.wgs_read_rescue import validate_screening_read_lineage
            validate_screening_read_lineage(call, obligation, merged_intake)
        calls.append(dict(source_allele_id=source_id, source_record_sha256=call.source.source_sha256,
            source_input_sha256=fields['source_input_sha256'], source_header_sha256=(
                source_header.header_sha256 if fields.get('authorized_read_evidence') is None else None),
            raw_line_sha256=fields.get('raw_line_sha256'), allele_number=fields.get('allele_number'),
            source_identity_sha256=hashlib.sha256(_canonical(_plain(source_mappings[source_id]))).hexdigest(),
            sample=None if call.sample is None else _plain(call.sample), filters=call.filters,
            technical_obligation=asdict(obligation),
            read_assessment=fields.get('authorized_read_assessment'), read_evidence=fields.get('authorized_read_evidence')))
    if calls and {c['source_allele_id'] for c in calls} != set(value.source_allele_ids):
        raise AnalysisEligibilityError('incomplete screening occurrence coverage')
    return AdmissionScreeningFacts(value.retained_ordinal, value.identity_id,
        None if value.allele is None else value.allele.allele_id, value.subject_id,
        value.source_allele_ids, value.required_rule_ids, value.transport,
        value.parent_source_manifest_sha256, value.parent_working_set_manifest_sha256,
        value.native_source_contract_sha256, header.header_sha256, selected_sample,
        tuple(sorted({m.reference_manifest_sha256 for m in mappings})),
        tuple(sorted(m.identity_id for m in mappings if m.stage is IdentityStage.NORMALIZED)),
        tuple(sorted((m.identity_id, m.identity) for m in mappings if type(m.identity) is SequenceAllele)),
        tuple(sorted({t.ensembl_gene_id for t in (() if value.allele is None else
            value.allele.transcript_consequences or ()) if t.ensembl_gene_id is not None})),
        value.technical_obligations, _canonical(_plain(calls)), refcall.can_defer, refcall.reasons,
        None if merged_intake is None else merged_intake.sha256, policy.sha256)


def _validate_screening_source_call(call, obligation, mapping, *, header, selected_sample):
    """Check declared raw bindings even when the scientific predicate assesses."""
    fields = call.source.source_fields
    columns = fields.get('columns')
    if (type(columns) is not tuple or len(columns) != 9+len(header.sample_names)
            or any(type(v) is not str for v in columns)
            or fields.get('source_input_sha256') != header.input_sha256
            or fields.get('source_allele_id') != call.source.source_record_id
            or fields.get('technical_flags') != obligation.technical_flags
            or fields.get('read_validation_required') is not obligation.read_validation_required
            or fields.get('filter_value') != columns[6]
            or call.filters != (None if columns[6] == '.' else () if columns[6] == 'PASS'
                else tuple(columns[6].split(';')))
            or obligation.source_filtered != (columns[6] not in ('PASS', '.'))):
        raise AnalysisEligibilityError('screening original source binding mismatch')
    line = '\t'.join(columns)
    if fields.get('raw_line_sha256') not in {hashlib.sha256((line+end).encode()).hexdigest()
            for end in ('', '\n', '\r\n')}:
        raise AnalysisEligibilityError('screening source line digest mismatch')
    number = fields.get('allele_number'); alternatives = columns[4].split(',')
    if type(number) is not int or not 1 <= number <= len(alternatives):
        raise AnalysisEligibilityError('screening source ALT binding missing')
    if type(mapping.identity) is SequenceAllele:
        source = mapping.identity
        if (source.contig, str(source.position), source.reference, source.alternate) != (
                columns[0], columns[1], columns[3], alternatives[number-1]):
            raise AnalysisEligibilityError('screening source ALT identity mismatch')
    if call.sample is not None:
        if call.sample.sample_name not in header.sample_names or call.sample.sample_name != obligation.sample_name:
            raise AnalysisEligibilityError('screening source sample missing')
        keys = columns[8].split(':')
        values = columns[9+header.sample_names.index(call.sample.sample_name)].split(':')
        if len(keys) != len(set(keys)) or len(keys) != len(values):
            raise AnalysisEligibilityError('screening source FORMAT conflict')
        raw = dict(zip(keys, values))
        if raw.get('GT') != call.sample.genotype:
            raise AnalysisEligibilityError('screening source genotype binding mismatch')
        for key, expected in (('PGT', call.sample.pgt), ('PID', call.sample.pid), ('PS', call.sample.phase_set)):
            if raw.get(key) not in (None, '.', expected) or expected is not None and raw.get(key) != expected:
                raise AnalysisEligibilityError('screening source phase binding mismatch')


def _screening_facts_mapping(facts, calls):
    # These fields have an exact, shallow schema. JSON already serializes the
    # immutable tuple containers and string enums; recursively copying every
    # scalar through dataclasses.asdict repeats the compact scan's work.
    value = {name: getattr(facts, name) for name in AdmissionScreeningFacts.__dataclass_fields__}
    value['sequence_identities'] = tuple((key, {name: getattr(sequence, name)
        for name in SequenceAllele.__dataclass_fields__}) for key, sequence in facts.sequence_identities)
    value['technical_obligations'] = tuple({name: getattr(obligation, name)
        for name in SourceTechnicalObligation.__dataclass_fields__} for obligation in facts.technical_obligations)
    value['source_call_facts'] = calls
    return value


def _screening_facts_from_mapping(value):
    """Build and fully validate once from an operation-owned parsed mapping."""
    if type(value) is not dict or set(value) != set(AdmissionScreeningFacts.__dataclass_fields__):
        raise ValueError
    value = dict(value)
    calls = value['source_call_facts']
    for name in ('source_allele_ids', 'required_rule_ids', 'reference_manifest_sha256s',
            'normalized_identity_ids', 'gene_ids', 'refcall_reasons'):
        if type(value[name]) is not list: raise ValueError
        value[name] = tuple(value[name])
    value['technical_obligations'] = _typed_value(value['technical_obligations'],
        tuple[SourceTechnicalObligation, ...], decode=True)
    value['sequence_identities'] = _screening_sequences_from_mapping(value['sequence_identities'])
    value['source_call_facts'] = _canonical(calls)
    # Populate only the exact dataclass fields, then run every semantic check.
    # The public constructor separately parses/validates its external byte field.
    result = object.__new__(AdmissionScreeningFacts)
    for name,item in value.items(): object.__setattr__(result,name,item)
    result._validate(calls)
    return result


def _screening_sequences_from_mapping(values):
    """Reuse identical immutable coordinates within this one parsed record.

    Original source/native/normalized IDs and their order all remain present.
    Only exact field values and scalar types share a strictly decoded child;
    differing fields/types still take the ordinary decoder and its full checks.
    The small table is discarded before this operation returns.
    """
    if type(values) not in (list,tuple):raise ValueError
    coordinates={};result=[]
    for pair in values:
        if type(pair) not in (list,tuple) or len(pair)!=2 or type(pair[1]) is not dict:raise ValueError
        key=_typed_value(pair[0],str)
        # The JSON-owned SequenceAllele mapping has only scalar fields. Type
        # tags distinguish e.g. true from 1 before either reaches validation.
        # An invalid container value remains unhashable and is rejected.
        identity=tuple((name,type(child),child) for name,child in sorted(pair[1].items()))
        sequence=coordinates.get(identity)
        if sequence is None:
            sequence=_typed_value(pair[1],SequenceAllele,decode=True)
            coordinates[identity]=sequence
        result.append((key,sequence))
    return _typed_value(tuple(result),tuple[tuple[str,SequenceAllele],...])


@_validation_session(isolated=True)
def load_admission_screening_facts(content):
    try:
        if type(content) is not bytes or not 0 < len(content) <= MAX_ROW_BYTES:
            raise ValueError
        value = json.loads(content, object_pairs_hook=_unique)
        result = _screening_facts_from_mapping(value)
        if _canonical(_screening_facts_mapping(result, value['source_call_facts'])) != content:
            raise ValueError
        return result
    except (ValueError, TypeError, KeyError, UnicodeError):
        raise AnalysisEligibilityError('invalid serialized admission screening facts') from None


def _simple(sequence):
    ref, alt = sequence.reference, sequence.alternate
    return (len(ref) == len(alt) == 1 or len(ref) == 1 < len(alt) and alt.startswith(ref)
            or len(alt) == 1 < len(ref) and ref.startswith(alt))


@lru_cache(maxsize=4)
def _supported_refcall_header(header):
    header.__post_init__()
    lines=tuple(line.rstrip('\r\n') for line in header.raw_lines)
    return (tuple(line for line in lines if line.startswith('##DeepVariant_version='))==('##DeepVariant_version=1.0.0',)
        and tuple(line for line in lines if line.startswith('##FILTER=<ID=RefCall,'))==(HEADER_FILTER,))


def evaluate_refcall_eligibility(value, *, header, selected_sample):
    """Require every original selected-sample occurrence; default to assessment.

    Input/header owners authenticate bytes before calling this pure predicate.
    It also checks all source/header/call/ALT links needed to certify deferral.
    Missing/unsupported scientific facts retain assessment. Broken declared
    source identities fail, rather than masquerading as source absence.
    """
    if type(value) is not LedgerAdmissionInput or type(header) is not VcfQcStreamHeader:
        raise AnalysisEligibilityError('invalid eligibility source types')
    value.__post_init__()
    def assess(reason):
        return RefCallEligibility(False, (reason,), value.source_allele_ids, header.header_sha256)
    # Declared source integrity precedes scientific selection, including for
    # unsupported headers and ALT-bearing calls that will remain assessable.
    if value.allele is not None and value.allele.calls:
        observed = set()
        original_obligations = {row.source_allele_id:row for row in value.technical_obligations}
        for call in value.allele.calls:
            fields = call.source.source_fields
            if not isinstance(fields, Mapping):
                raise AnalysisEligibilityError('eligibility source fields missing')
            source_id = fields.get('source_allele_id')
            if source_id not in value.source_allele_ids or source_id != call.source.source_record_id:
                raise AnalysisEligibilityError('foreign eligibility source occurrence')
            if (fields.get('source_manifest_sha256') != value.parent_source_manifest_sha256
                    or fields.get('source_input_sha256') != header.input_sha256):
                raise AnalysisEligibilityError('eligibility source/header binding mismatch')
            columns = fields.get('columns')
            if type(columns) is not tuple or len(columns) != 9+len(header.sample_names) or any(type(v) is not str for v in columns):
                raise AnalysisEligibilityError('eligibility source columns missing')
            line = '\t'.join(columns)
            if fields.get('raw_line_sha256') not in {hashlib.sha256((line+end).encode()).hexdigest() for end in ('','\n','\r\n')}:
                raise AnalysisEligibilityError('eligibility source line digest mismatch')
            obligation = original_obligations[source_id]
            if (fields.get('technical_flags') != obligation.technical_flags
                    or fields.get('read_validation_required') is not obligation.read_validation_required
                    or fields.get('filter_value') != columns[6]
                    or obligation.source_filtered != (columns[6] not in ('PASS','.'))):
                raise AnalysisEligibilityError('eligibility original technical obligation mismatch')
            observed.add(source_id)
        if observed != set(value.source_allele_ids):
            raise AnalysisEligibilityError('incomplete eligibility occurrence coverage')
    if not _supported_refcall_header(header):
        return assess('unsupported_refcall_header')
    if selected_sample not in header.sample_names:
        return assess('selected_sample_unavailable')
    if value.allele is None or not value.allele.calls:
        return assess('source_calls_unavailable')
    native = next(m.identity for m in value.allele.identities if m.identity_id == value.allele.native_identity_id)
    if type(native) not in (SequenceAllele, MosaicAllele):
        return assess('protected_native_class')
    sequence = native.sequence if type(native) is MosaicAllele else native
    if type(native) is MosaicAllele and (native.source_genotype not in ('0/0', '0|0') or native.ploidy != 2):
        return assess('mosaic_source_genotype_not_homref')
    if re.fullmatch(r'(?:chr)?(?:[1-9]|1[0-9]|2[0-2])', sequence.contig) is None or not _simple(sequence):
        return assess('protected_sequence_class')
    obligations = {row.source_allele_id:row for row in value.technical_obligations}
    source_mappings = {m.identity_id:m for m in value.allele.identities if m.stage is IdentityStage.SOURCE}
    calls = {}
    for call in value.allele.calls:
        fields = call.source.source_fields
        if not isinstance(fields, Mapping): return assess('source_fields_unavailable')
        source_id = fields.get('source_allele_id')
        if source_id not in obligations or source_id != call.source.source_record_id:
            raise AnalysisEligibilityError('foreign eligibility source occurrence')
        if (fields.get('source_manifest_sha256') != value.parent_source_manifest_sha256
                or fields.get('source_input_sha256') != header.input_sha256):
            raise AnalysisEligibilityError('eligibility source/header binding mismatch')
        if fields.get('authorized_read_assessment') is not None:
            # Preserve independent evidence for the existing read/genotype
            # owners; mere presence here never certifies its carriage.
            return assess('independent_read_assessment_present')
        calls.setdefault(source_id, []).append(call)
    if set(calls) != set(value.source_allele_ids):
        raise AnalysisEligibilityError('incomplete eligibility occurrence coverage')
    sample_column = 9 + header.sample_names.index(selected_sample)
    for source_id in value.source_allele_ids:
        obligation = obligations[source_id]
        if obligation.sample_name != selected_sample: return assess('selected_sample_mismatch')
        if type(native) is MosaicAllele and not obligation.read_validation_required:
            return assess('mosaic_not_unconfirmed_screen')
        for call in calls[source_id]:
            sample, fields = call.sample, call.source.source_fields
            if sample is None or sample.sample_name != selected_sample: return assess('selected_sample_call_unavailable')
            if sample.genotype not in ('0/0','0|0') or sample.ploidy != 2:
                return assess('genotype_not_complete_homref')
            if sample.pgt is not None and (re.fullmatch(r'0[|/]0', sample.pgt) is None):
                return assess('phase_genotype_conflict')
            columns = fields.get('columns')
            if type(columns) is not tuple or len(columns) != 9+len(header.sample_names):
                raise AnalysisEligibilityError('eligibility source columns missing')
            line = '\t'.join(columns)
            if fields.get('raw_line_sha256') not in {hashlib.sha256((line+end).encode()).hexdigest() for end in ('','\n','\r\n')}:
                raise AnalysisEligibilityError('eligibility source line digest mismatch')
            if columns[6] != 'RefCall' or fields.get('filter_value') != 'RefCall' or call.filters != ('RefCall',):
                return assess('filter_not_supported_refcall')
            keys, values = columns[8].split(':'), columns[sample_column].split(':')
            if len(keys) != len(set(keys)) or len(keys) != len(values):
                raise AnalysisEligibilityError('eligibility source FORMAT conflict')
            raw = dict(zip(keys, values))
            if raw.get('GT') != sample.genotype:
                raise AnalysisEligibilityError('eligibility source genotype binding mismatch')
            for key, expected in (('PGT',sample.pgt), ('PID',sample.pid), ('PS',sample.phase_set)):
                if raw.get(key) not in (None,'.',expected) or expected is not None and raw.get(key) != expected:
                    raise AnalysisEligibilityError('eligibility source phase binding mismatch')
            mapping = source_mappings.get(source_id)
            number = fields.get('allele_number'); alternatives = columns[4].split(',')
            if (mapping is None or type(mapping.identity) is not SequenceAllele or type(number) is not int
                    or not 1 <= number <= len(alternatives)):
                raise AnalysisEligibilityError('eligibility source ALT binding missing')
            original = mapping.identity
            if (original.contig, str(original.position), original.reference, original.alternate) != (
                    columns[0], columns[1], columns[3], alternatives[number-1]):
                raise AnalysisEligibilityError('eligibility source ALT identity mismatch')
            if not _simple(original): return assess('protected_source_representation')
    return RefCallEligibility(True, ('supported_unconfirmed_refcall_homref',), value.source_allele_ids, header.header_sha256)


@dataclass(frozen=True)
class AnalysisEligibilityDecision:
    retained_ordinal: int
    identity_id: str
    allele_id: str | None
    source_allele_ids: tuple[str, ...]
    disposition: str
    reasons: tuple[str, ...]
    admission_binding_sha256: str
    input_sha256: str
    proof_sha256: str
    eligibility_policy_sha256: str
    source_manifest_sha256: str
    reference_manifest_sha256: str
    resource_bindings: tuple[tuple[str, str], ...]
    proof_references: tuple[tuple[str, str, str], ...]
    schema: str = 'genosapiens.wgs-analysis-eligibility-decision'
    version: int = 2

    def __post_init__(self):
        if (type(self.retained_ordinal) is not int or self.retained_ordinal < 1
                or type(self.identity_id) is not str or not self.identity_id
                or self.allele_id is not None and (type(self.allele_id) is not str or not self.allele_id)
                or type(self.disposition) is not str or self.disposition not in DISPOSITIONS
                or type(self.version) is not int or self.version != 2
                or self.schema != 'genosapiens.wgs-analysis-eligibility-decision'
                or type(self.source_allele_ids) is not tuple or not self.source_allele_ids
                or any(type(v) is not str or not v for v in self.source_allele_ids)
                or len(set(self.source_allele_ids)) != len(self.source_allele_ids)
                or type(self.reasons) is not tuple or not self.reasons
                or any(type(v) is not str or not v for v in self.reasons)
                or self.reasons != tuple(sorted(set(self.reasons)))
                or any(not _sha(v) for v in (self.admission_binding_sha256, self.input_sha256, self.proof_sha256,
                    self.eligibility_policy_sha256, self.source_manifest_sha256, self.reference_manifest_sha256))
                or type(self.resource_bindings) is not tuple
                or any(type(v) is not tuple or len(v) != 2 for v in self.resource_bindings)
                or self.resource_bindings != tuple(sorted(self.resource_bindings))
                or len(dict(self.resource_bindings)) != len(self.resource_bindings)
                or any(type(key) is not str or not key or not _sha(value) for key,value in self.resource_bindings)
                or type(self.proof_references) is not tuple
                or any(type(v) is not tuple or len(v) != 3 or type(v[0]) is not str or not v[0]
                    or v[1] != 'witnesses' or not _sha(v[2])
                    or v[0] not in ('refcall_header','common_exclusion','population_noncontribution')
                    for v in self.proof_references)
                or self.proof_references != tuple(sorted(set(self.proof_references)))):
            raise AnalysisEligibilityError('invalid analysis eligibility decision')
        required = {'defer_refcall': 'refcall_header', 'exclude_common': 'common_exclusion'}.get(self.disposition)
        if required is not None and not any(role == required for role,path,digest in self.proof_references):
            raise AnalysisEligibilityError('eligibility disposition lacks its rule witness')

    def as_bytes(self): return _canonical(_eligibility_decision_mapping(self))


def _eligibility_decision_mapping(decision):
    # All nested fields are tuples of JSON scalars. Preserve their exact order
    # and multiplicity without a second recursive object-copy traversal.
    return {name: getattr(decision, name) for name in AnalysisEligibilityDecision.__dataclass_fields__}


def _eligibility_decision_from_mapping(data):
    if type(data) is not dict or set(data) != set(AnalysisEligibilityDecision.__dataclass_fields__): raise ValueError
    data=dict(data)
    for field in ('source_allele_ids', 'reasons', 'resource_bindings', 'proof_references'):
        if type(data[field]) is not list: raise ValueError
        if field in ('resource_bindings', 'proof_references'):
            if any(type(v) is not list for v in data[field]): raise ValueError
            data[field] = tuple(tuple(v) for v in data[field])
        else: data[field] = tuple(data[field])
    return AnalysisEligibilityDecision(**data)


def load_analysis_eligibility_decision(content):
    try:
        if type(content) is not bytes or not 0 < len(content) <= 1024*1024: raise ValueError
        data = json.loads(content, object_pairs_hook=_unique)
        decision = _eligibility_decision_from_mapping(data)
        if decision.as_bytes() != content: raise ValueError
        return decision
    except (TypeError, ValueError, UnicodeError):
        raise AnalysisEligibilityError('invalid serialized eligibility decision') from None


@dataclass(frozen=True)
class EligibilityBatch:
    decisions: tuple[AnalysisEligibilityDecision, ...]

    def __post_init__(self):
        if (type(self.decisions) is not tuple or any(type(v) is not AnalysisEligibilityDecision for v in self.decisions)
                or any(b.retained_ordinal <= a.retained_ordinal for a,b in zip(self.decisions,self.decisions[1:]))
                or len({v.identity_id for v in self.decisions}) != len(self.decisions)):
            raise AnalysisEligibilityError('invalid ordered eligibility batch')

    @property
    def assessed_subjects(self):
        return tuple(row.allele_id or row.identity_id for row in self.decisions if row.disposition == 'assess')

    def require_inputs(self,inputs,*,input_rows_sha256=None):
        if type(inputs) is not tuple or len(inputs)!=len(self.decisions):
            raise AnalysisEligibilityError('eligibility batch input coverage mismatch')
        for value,decision in zip(inputs,self.decisions,strict=True):
            if (type(value) is not LedgerAdmissionInput or value.retained_ordinal!=decision.retained_ordinal
                    or value.identity_id!=decision.identity_id or value.subject_id!=(decision.allele_id or decision.identity_id)
                    or value.source_allele_ids!=decision.source_allele_ids
                    or value.parent_source_manifest_sha256!=decision.source_manifest_sha256):
                raise AnalysisEligibilityError('eligibility batch original occurrence mismatch')
        if self.decisions and len({row.eligibility_policy_sha256 for row in self.decisions})!=1:
            raise AnalysisEligibilityError('eligibility batch policy conflict')
        if input_rows_sha256 is not None and hashlib.sha256(_canonical([row.input_sha256 for row in self.decisions])).hexdigest()!=input_rows_sha256:
            raise AnalysisEligibilityError('eligibility batch original digest mismatch')
        selected=set(self.assessed_subjects)
        return tuple(value for value in inputs if value.subject_id in selected)


def evaluate_analysis_eligibility(facts, *, policy, common_decisions=()):
    """One decision predicate over producer facts and fresh common evidence."""
    if (type(policy) is not AnalysisEligibilityPolicy or type(facts) is not AdmissionScreeningFacts
            or facts.eligibility_policy_sha256 != policy.sha256):
        raise AnalysisEligibilityError('invalid eligibility policy owner')
    if facts.merged_intake_sha256 is not None:
        if common_decisions:
            raise AnalysisEligibilityError('ordinary common filtering does not apply to merged intake')
        return 'assess', ('merged_input_assessment_domain',)
    if facts.refcall_can_defer: return 'defer_refcall', facts.refcall_reasons
    if common_decisions:
        if (type(common_decisions) is not tuple or any(type(v) is not PrefilterDecision for v in common_decisions)
                or tuple(v.source_allele_id for v in common_decisions) != facts.source_allele_ids
                or any(v.policy_sha256 != policy.common_filter_policy_sha256
                    or v.normalized_allele_id not in facts.normalized_identity_ids for v in common_decisions)):
            raise AnalysisEligibilityError('common eligibility occurrence/policy mismatch')
        if all(v.disposition == 'exclude_common' for v in common_decisions):
            return 'exclude_common', ('all_occurrences_proved_common',)
    return 'assess', tuple(sorted(set((*facts.refcall_reasons, 'conservative_assessment'))))


def evaluate_profile_intake(value, *, upstream_decision, molecular_sources,
        selection, catalog, route_sources, clinical_session, regulatory_index, reference, _clinical_window=None,
        _record_owner=None, _route_owner=None):
    """Apply the shared scope policy after validated upstream intake screening.

    The original input digest is independent of admission proof payloads.
    Upstream common/RefCall states and their witnesses remain separate from scope.
    """
    from genosapiens.wgs_admission import LedgerAdmissionInput
    from genosapiens.track1_analysis_profiles import (AnalysisProfileSelection,
        ProfileCatalog, ProfileCandidateFacts, ProfileRouteFacts, ROUTES, evaluate_profile_candidate)
    from genosapiens.wgs_profile_projection import derive_profile_candidate_facts, profile_route_sources
    if _record_owner is not None:
        from genosapiens.wgs_admission_projection import _ProfileRecord
        if type(_record_owner) is not _ProfileRecord:
            raise AnalysisEligibilityError('profile input owner differs')
        input_sha256 = _record_owner.input_digest(value)
    else:
        input_sha256 = hashlib.sha256(value.as_bytes()).hexdigest()
    if (type(value) is not LedgerAdmissionInput or type(upstream_decision) is not AnalysisEligibilityDecision
            or type(selection) is not AnalysisProfileSelection or type(catalog) is not ProfileCatalog
            or (upstream_decision.retained_ordinal,upstream_decision.identity_id,
                upstream_decision.allele_id,upstream_decision.source_allele_ids,
                upstream_decision.input_sha256,upstream_decision.source_manifest_sha256) !=
               (value.retained_ordinal,value.identity_id,None if value.allele is None else value.allele.allele_id,
                value.source_allele_ids,input_sha256,
                value.parent_source_manifest_sha256)):
        raise AnalysisEligibilityError('profile intake original/upstream binding differs')
    if _route_owner is None:
        expected = profile_route_sources(annotation_manifest_sha256=route_sources['molecular'],
            native_source_sha256=value.native_source_contract_sha256,
            clinical_session=clinical_session,regulatory_index=regulatory_index)
        if expected != route_sources or reference != clinical_session.reference:
            raise AnalysisEligibilityError('profile intake route sources differ')
    else:
        from genosapiens.wgs_profile_projection import _ProfileRouteContext
        if type(_route_owner) is not _ProfileRouteContext:
            raise AnalysisEligibilityError('profile route owner differs')
        _route_owner.require(route_sources, clinical_session, regulatory_index, reference,
            value.native_source_contract_sha256)
    if upstream_decision.disposition == 'assess':
        facts = derive_profile_candidate_facts(value,original_input_sha256=upstream_decision.input_sha256,
            upstream_disposition='assess',molecular_sources=molecular_sources,route_sources=route_sources,
            clinical_session=clinical_session,regulatory_index=regulatory_index,reference=reference,
            _clinical_window=_clinical_window, _route_owner=_route_owner)
    else:
        facts = ProfileCandidateFacts(value.identity_id,upstream_decision.input_sha256,
            upstream_decision.disposition,tuple(ProfileRouteFacts(r,'not_applicable',route_sources[r]) for r in ROUTES),
            source_occurrences=value.source_allele_ids)
    return facts,evaluate_profile_candidate(facts,policy=selection.policy,catalog=catalog,profile=selection.profile)


def _fail_closed_operation(function):
    from functools import wraps
    import time
    @wraps(function)
    def operation(self,*args,**kwargs):
        if getattr(self,'failed',False): raise AnalysisEligibilityError('eligibility owner previously failed')
        started=time.monotonic_ns()
        try: return function(self,*args,**kwargs)
        except BaseException:
            self.failed=True
            raise
        finally:
            if hasattr(self,'metrics'):
                key='seal_nanoseconds' if function.__name__=='finish' else 'write_nanoseconds'
                self.metrics[key]+=time.monotonic_ns()-started
    return operation


class AnalysisEligibilitySession:
    """One ordered, bounded coordinator over authenticated original indexes.

    Connections and the public resource session are opened and authenticated by
    their existing owners. This object neither rebuilds inputs nor publishes a
    completion marker. Witnesses are handed to the durable stream per window.
    """
    def __init__(self, *, policy, common_policy, prepared_index, source_manifest,
                 source_connection, filter_state, routing, header, selected_sample,
                 execution_config, public_session, precision, resource_bindings, case_cache=None, merged_intake=None):
        screening=AdmissionScreeningContext(header,selected_sample,policy,merged_intake)
        from genosapiens.bounded_wgs_config import BoundedWgsExecutionConfig, bounded_wgs_execution_config_as_json
        from genosapiens.wgs_source_routing import SourceRoutingContext
        if (type(execution_config) is not BoundedWgsExecutionConfig
                or type(routing) is not SourceRoutingContext
                or type(resource_bindings) is not tuple or resource_bindings!=tuple(sorted(resource_bindings))
                or len(dict(resource_bindings))!=len(resource_bindings)
                or any(not _sha(v) for k,v in resource_bindings)):
            raise AnalysisEligibilityError('invalid screening session binding')
        self.screening=screening
        self.source_sha256=source_manifest.sha256 if merged_intake is None else merged_intake.sha256
        if merged_intake is None:
            self._validate_original(policy=policy,common_policy=common_policy,prepared_index=prepared_index,
                source_manifest=source_manifest,source_connection=source_connection,filter_state=filter_state,
                routing=routing,header=header,selected_sample=selected_sample,execution_config=execution_config,
                public_session=public_session,precision=precision,resource_bindings=resource_bindings)
        elif (not any(source.source_manifest_sha256==source_manifest.sha256
                    and source.native_contract_sha256==routing.contract.sha256
                    for source in merged_intake.source_bindings)
                or merged_intake.execution_config_sha256!=hashlib.sha256(
                    bounded_wgs_execution_config_as_json(execution_config)).hexdigest()
                or any(value is not None for value in (source_connection,filter_state,public_session))):
            raise AnalysisEligibilityError('merged screening source ownership is invalid')
        self.policy,self.common_policy,self.prepared_index=policy,common_policy,prepared_index
        import os
        import threading
        self.owner_pid=os.getpid()
        self.owner_thread=threading.get_ident()
        self.source_manifest,self.source_connection,self.filter_connection=source_manifest,source_connection,None if filter_state is None else filter_state._connection
        self.routing,self.header,self.selected_sample=routing,header,selected_sample
        self.config,self.public_session,self.precision=execution_config,public_session,precision
        self.resource_bindings,self.case_cache=resource_bindings,case_cache
        self.completed_inputs=self.completed_occurrences=0
        self.counts={key:0 for key in DISPOSITIONS}; self.occurrences={key:0 for key in DISPOSITIONS}
        self.witnesses={}; self.finished=False; self.reconsidered_occurrences=0; self.public_query_batches=0

    @staticmethod
    def _validate_original(*,policy,common_policy,prepared_index,source_manifest,source_connection,
            filter_state,routing,header,selected_sample,execution_config,public_session,precision,resource_bindings):
        from genosapiens.common_allele_filter import CommonAlleleFilterPolicy, common_allele_filter_policy_as_json
        from genosapiens.wgs_prepared_case import PreparedIndex
        from genosapiens.wgs_source_routing import SourceRoutingContext
        from genosapiens.population_clinical_session import PopulationClinicalResourceSession
        from genosapiens.bounded_wgs_config import BoundedWgsExecutionConfig
        from genosapiens.wgs_filter import _FilterState
        from genosapiens.wgs_ledger import LedgerManifest
        from genosapiens.wgs_annotation import _projection_binding
        from genosapiens.gnomad_faf95_adapter import GnomadFaf95PrecisionQualification
        if (type(policy) is not AnalysisEligibilityPolicy or type(common_policy) is not CommonAlleleFilterPolicy
                or hashlib.sha256(common_allele_filter_policy_as_json(common_policy)).hexdigest() != policy.common_filter_policy_sha256
                or type(prepared_index) is not PreparedIndex or type(routing) is not SourceRoutingContext
                or type(execution_config) is not BoundedWgsExecutionConfig
                or type(public_session) is not PopulationClinicalResourceSession
                or type(precision) is not GnomadFaf95PrecisionQualification
                or type(source_manifest) is not LedgerManifest
                or type(filter_state) is not _FilterState or filter_state._closed or not filter_state._sealed
                or prepared_index.source_manifest_sha256 != source_manifest.sha256
                or routing.contract.parent_source_manifest_sha256 != source_manifest.sha256
                or type(header) is not VcfQcStreamHeader or header.input_sha256 != source_manifest.source.input_sha256
                or header.input_size_bytes != source_manifest.source.input_size_bytes
                or header.reference_concordance_sha256 != source_manifest.source.concordance_sha256
                or type(selected_sample) is not str or selected_sample not in header.sample_names
                or type(resource_bindings) is not tuple or resource_bindings != tuple(sorted(resource_bindings))
                or len(dict(resource_bindings)) != len(resource_bindings)
                or any(not _sha(v) for k,v in resource_bindings)):
            raise AnalysisEligibilityError('invalid eligibility session binding')
        required={'population_manifest','precision_manifest','clinical_manifest','exception_manifest',
                  'all_sites_manifest','public_index_manifest'}
        if not required <= dict(resource_bindings).keys():
            raise AnalysisEligibilityError('incomplete eligibility resource bindings')
        resources=dict(resource_bindings)
        if (public_session._population_arguments['joint'].manifest_sha256!=resources['population_manifest']
                or public_session._population_arguments['all_sites'].bundle_manifest_sha256!=resources['all_sites_manifest']
                or public_session._population_arguments['reference'].manifest_sha256!=routing.contract.reference_manifest_sha256
                or precision.publication_manifest_sha256!=resources['population_manifest']
                or precision.report_sha256!=resources['precision_manifest']
                or public_session._verified.manifest_sha256!=resources['public_index_manifest']):
            raise AnalysisEligibilityError('foreign eligibility public resource binding')
        public_session._assert_active(); public_session._assert_clinical_stable()
        envelope=_canonical(dict(schema='genosapiens.wgs-header',version=1,classification='controlled',
            llm_policy='no-llm',selected_sample_name=selected_sample,header=asdict(header)))
        if (len(envelope)!=source_manifest.header.size_bytes
                or hashlib.sha256(envelope).hexdigest()!=source_manifest.header.sha256):
            raise AnalysisEligibilityError('foreign eligibility original header artifact')
        if (filter_state.source_manifest_sha256!=source_manifest.sha256
                or filter_state.policy_sha256!=policy.common_filter_policy_sha256
                or filter_state.resource_sha256!=tuple(resources[k] for k in ('population_manifest','clinical_manifest','exception_manifest','precision_manifest'))
                or filter_state.working_binding()[0]!=prepared_index.working_manifest_sha256
                or filter_state._connection.execute('SELECT filter_manifest_sha256 FROM filter_meta WHERE singleton=1').fetchone()!=(prepared_index.filter_manifest_sha256,)
                or source_connection.execute('PRAGMA query_only').fetchone()!=(1,)):
            raise AnalysisEligibilityError('foreign eligibility original index binding')
        meta=source_connection.execute('SELECT binding_sha256,identities,occurrences,complete FROM projection_meta WHERE singleton=1').fetchone()
        if meta is None or meta[0]!=_projection_binding(prepared_index,routing,execution_config) or meta[3]!=1:
            raise AnalysisEligibilityError('incomplete eligibility retained source projection')

    def _witness(self, role, content):
        digest=hashlib.sha256(content).hexdigest()
        existing=self.witnesses.setdefault(digest,content)
        if existing != content: raise AnalysisEligibilityError('eligibility witness digest conflict')
        return role,'witnesses',digest

    @_fail_closed_operation
    def evaluate_window(self, records):
        from genosapiens.wgs_admission_projection import EncodedAdmissionRow
        from genosapiens.wgs_admission import iter_admission_source_rows, _occurrence
        from genosapiens.wgs_source_routing import classify_source_occurrence
        from genosapiens.wgs_public_filter import (_prefilter_occurrence,_population_request,_clinical_requests,
                                                  _population_filter_records,clinical_protection_facts)
        from genosapiens.evidence_serialization import to_json
        from genosapiens.population_annotation import _PopulationScreeningResult, _PopulationScreeningRow
        from genosapiens.common_allele_filter import prefilter_bypass_reasons,evaluate_common_allele_filter
        from dataclasses import replace
        import os,threading
        if (self.owner_pid!=os.getpid() or self.owner_thread!=threading.get_ident() or self.finished or type(records) is not tuple or not records
                or len(records)>min(self.config.public_queries.batch_rows,self.config.sqlite.cursor_batch_rows)
                or any(not isinstance(row, EncodedAdmissionRow) for row in records)):
            raise AnalysisEligibilityError('invalid eligibility window')
        self.witnesses={}
        compact_facts = tuple(row.authenticated_facts() for row in records)
        for offset,value in enumerate(compact_facts,1):
            if value.retained_ordinal != self.completed_inputs+offset:
                raise AnalysisEligibilityError('eligibility window ordinal gap')
            if (value.source_manifest_sha256 != self.source_sha256
                    or value.header_sha256 != self.header.header_sha256
                    or value.selected_sample != self.selected_sample
                    or value.eligibility_policy_sha256 != self.policy.sha256
                    or value.merged_intake_sha256 != (None if self.screening.merged_intake is None else self.screening.merged_intake.sha256)):
                raise AnalysisEligibilityError('foreign eligibility window source')
            if self.screening.merged_intake is None and any(o.sample_name!=self.selected_sample for o in value.technical_obligations):
                raise AnalysisEligibilityError('foreign eligibility selected sample')
        preliminary={}; candidate_ids=set()
        source_ids=tuple(source_id for value in compact_facts for source_id in value.source_allele_ids)
        if len(set(source_ids)) != len(source_ids): raise AnalysisEligibilityError('duplicate eligibility occurrence')
        reconsider = tuple(value for value in compact_facts if not value.refcall_can_defer
            and value.merged_intake_sha256 is None
            and value.refcall_reasons != ('independent_read_assessment_present',))
        lookup_ids = tuple(key for value in reconsider for key in value.source_allele_ids)
        reasons={key:set() for key in lookup_ids}
        # Existing primary-key lookups, bounded per window, never shard scans.
        for start in range(0,len(lookup_ids),512):
            keys=lookup_ids[start:start+512]
            cursor=self.filter_connection.execute('SELECT source_allele_id,reason FROM reasons WHERE source_allele_id IN ('+
                ','.join('?' for _ in keys)+')',keys)
            try:
                for source_id,reason in cursor: reasons[source_id].add(reason)
            finally: cursor.close()
        for value in compact_facts:
            disposition, reason = evaluate_analysis_eligibility(value, policy=self.policy)
            if disposition == 'defer_refcall':
                refs=(self._witness('refcall_header',_canonical(dict(schema='genosapiens.wgs-refcall-header-proof',version=1,
                    rule=RULE,selected_sample=self.selected_sample,header=asdict(self.header)))),)
                preliminary[value.retained_ordinal]=(disposition,reason,refs)
            else:
                preliminary[value.retained_ordinal]=(disposition,reason,())
                # A historic reason selects fresh work only. It never proves a
                # new exclusion; all occurrences are reconstructed below.
                if any('exomes_counts_missing' in reasons.get(key, ()) for key in value.source_allele_ids):
                    candidate_ids.add(value.retained_ordinal)
        original={}
        if candidate_ids:
            for identity,payloads in iter_admission_source_rows(self.source_connection,
                    ordinals=tuple(sorted(candidate_ids)),
                    execution_config=self.config):
                original[identity[0]]=(identity,payloads)
            if set(original) != candidate_ids: raise AnalysisEligibilityError('eligibility source index coverage mismatch')
        pending=[]; by_ordinal={}
        for value in compact_facts:
            ordinal=value.retained_ordinal
            if ordinal not in candidate_ids: continue
            identity,payloads=original[ordinal]
            if identity[1] != value.identity_id or tuple(k for k,p in payloads) != value.source_allele_ids:
                raise AnalysisEligibilityError('eligibility source identity coverage mismatch')
            items=[]; protected=False
            for source_id,payload in payloads:
                facts=_occurrence(payload,source_allele_id=source_id,source_manifest=self.source_manifest,contract=self.routing.contract)
                route=classify_source_occurrence(facts,contract=self.routing.contract,renderer=self.routing.renderer)
                occurrence=_prefilter_occurrence(facts,self.routing.renderer,self.routing.contract,route=route)
                protected |= bool(prefilter_bypass_reasons(occurrence))
                items.append((facts,occurrence))
            if protected:
                preliminary[ordinal]=('assess',('original_occurrence_protected',),())
            else:
                by_ordinal[ordinal]=items; pending.extend((ordinal,facts,occurrence) for facts,occurrence in items)
        decisions_by_source={}; proofs_by_source={}; bindings=dict(self.resource_bindings)
        for start in range(0,len(pending),self.config.public_queries.batch_rows):
            part=pending[start:start+self.config.public_queries.batch_rows]
            requests=tuple(_population_request(facts,reference_manifest_sha256=self.routing.contract.reference_manifest_sha256)
                for ordinal,facts,occurrence in part)
            from genosapiens.wgs_failures import WgsFailure
            self.public_query_batches+=1
            try:
                population=self.public_session._population_screening_facts(requests,_case_cache=self.case_cache)
                clinical=self.public_session.query_clinical(_clinical_requests(tuple((facts,occurrence) for ordinal,facts,occurrence in part),self.routing.renderer))
            except WgsFailure as error:
                if error.category != 'timeout': raise
                for ordinal,facts,occurrence in part:
                    by_ordinal.pop(ordinal,None)
                    preliminary[ordinal]=('assess',('fresh_public_evidence_timed_out',),())
                continue
            if type(population) is not _PopulationScreeningResult:
                raise AnalysisEligibilityError('eligibility population response mismatch')
            proofs=population.proofs
            expected={request.allele.allele_id:request for request in requests}
            responses={}
            for fact in population.rows:
                if (type(fact) is not _PopulationScreeningRow or fact.allele_id not in expected
                        or fact.allele_id in responses
                        or fact.normalized_identity_id != expected[fact.allele_id].normalized_identity_id
                        or fact.input_sha256 != hashlib.sha256(to_json(expected[fact.allele_id].allele)).hexdigest()):
                    raise AnalysisEligibilityError('eligibility population response mismatch')
                responses[fact.allele_id]=fact
            if set(responses)!=set(expected) or len(clinical)!=len(part):
                raise AnalysisEligibilityError('eligibility public response coverage mismatch')
            for (ordinal,facts,occurrence),protected in zip(part,clinical,strict=True):
                selected=tuple(p for p in proofs if p.normalized_allele_id==occurrence.normalized_allele_id)
                # Identical normalized requests may produce identical witnesses.
                selected=tuple({p.as_bytes():p for p in selected}.values())
                response=responses[occurrence.source_allele_id]
                pop=_population_filter_records(response.population,expected[occurrence.source_allele_id].allele.identities,precision=self.precision,
                    population_manifest_sha256=bindings['population_manifest'],precision_manifest_sha256=bindings['precision_manifest'],noncontribution_proofs=selected)
                protect=clinical_protection_facts(protected.status,protected.records,clinical_manifest_sha256=bindings['clinical_manifest'],exception_manifest_sha256=bindings['exception_manifest'])
                decision=evaluate_common_allele_filter(occurrence,pop,protect,policy=self.common_policy)
                decisions_by_source[occurrence.source_allele_id]=decision
                refs=[self._witness('common_exclusion',_canonical(dict(schema='genosapiens.wgs-analysis-common-witness',version=1,
                    occurrence=asdict(occurrence),population=asdict(pop),clinical=asdict(protect),decision=asdict(decision),
                    population_query_receipts=tuple(json.loads(receipt) for receipt in response.query_receipts))))]
                refs.extend(self._witness('population_noncontribution',p.as_bytes()) for p in selected)
                proofs_by_source[occurrence.source_allele_id]=refs
            self.reconsidered_occurrences+=len(part)
        result=[]
        for row,value in zip(records,compact_facts,strict=True):
            disposition,reason,refs=preliminary[value.retained_ordinal]
            if value.retained_ordinal in by_ordinal:
                common=tuple(decisions_by_source[k] for k in value.source_allele_ids)
                disposition,reason=evaluate_analysis_eligibility(value,policy=self.policy,common_decisions=common)
                refs=tuple(sorted(set(ref for key in value.source_allele_ids for ref in proofs_by_source[key])))
            result.append(AnalysisEligibilityDecision(value.retained_ordinal,value.identity_id,row.allele_id,
                value.source_allele_ids,disposition,reason,row.binding_sha256,row.input_sha256,
                row.proof_sha256,self.policy.sha256,value.source_manifest_sha256,
                self.routing.contract.reference_manifest_sha256,self.resource_bindings,refs))
        batch=EligibilityBatch(tuple(result))
        required_witnesses={ref[2] for row in batch.decisions for ref in row.proof_references}
        # Each completed window owns a fresh immutable map of immutable bytes.
        # A later window may replace session.witnesses while this one is queued.
        self.witnesses=MappingProxyType({key:value for key,value in self.witnesses.items()
            if key in required_witnesses})
        for decision in batch.decisions:
            self.counts[decision.disposition]+=1; self.occurrences[decision.disposition]+=len(decision.source_allele_ids)
        self.completed_inputs+=len(result); self.completed_occurrences+=len(source_ids)
        return batch

    @_fail_closed_operation
    def finish(self, *, expected_inputs, expected_occurrences):
        import os,threading
        if (self.owner_pid!=os.getpid() or self.owner_thread!=threading.get_ident() or self.finished or self.completed_inputs != expected_inputs or self.completed_occurrences != expected_occurrences
                or sum(self.counts.values()) != expected_inputs or sum(self.occurrences.values()) != expected_occurrences):
            raise AnalysisEligibilityError('incomplete eligibility session')
        self.finished=True
        return dict(inputs=self.completed_inputs,occurrences=self.completed_occurrences,
                    input_dispositions=dict(self.counts),occurrence_dispositions=dict(self.occurrences))


@dataclass(frozen=True)
class EligibilityShard:
    relative_path: str
    sha256: str
    size_bytes: int
    uncompressed_sha256: str
    uncompressed_size_bytes: int
    row_count: int
    first_ordinal: int
    last_ordinal: int

    def __post_init__(self):
        if (type(self.relative_path) is not str or not re.fullmatch(r'(?:decisions|witnesses)-[0-9]{8}\.jsonl\.gz',self.relative_path)
                or not _sha(self.sha256) or not _sha(self.uncompressed_sha256)
                or any(type(v) is not int or v<1 for v in (self.size_bytes,self.uncompressed_size_bytes,self.row_count,self.first_ordinal,self.last_ordinal))
                or self.size_bytes>128*1024*1024 or self.uncompressed_size_bytes>64*1024*1024
                or self.last_ordinal-self.first_ordinal+1!=self.row_count):
            raise AnalysisEligibilityError('invalid eligibility shard')


@dataclass(frozen=True)
class AnalysisEligibilityManifest:
    admission_binding_sha256: str
    eligibility_policy_sha256: str
    source_manifest_sha256: str
    reference_manifest_sha256: str
    input_rows: int
    source_occurrences: int
    input_dispositions: tuple[tuple[str,int], ...]
    occurrence_dispositions: tuple[tuple[str,int], ...]
    ordered_decisions_sha256: str
    ordered_original_inputs_sha256: str
    ordered_original_proofs_sha256: str
    ordered_screening_sha256: str
    shards: tuple[EligibilityShard, ...]
    witness_shards: tuple[EligibilityShard, ...]
    witness_count: int
    dependencies: tuple[tuple[str,str,str,int], ...]
    resource_bindings: tuple[tuple[str,str], ...]
    schema: str = 'genosapiens.wgs-analysis-eligibility-manifest'
    version: int = 3
    status: str = 'complete'
    analysis_domain_status: str = 'complete_for_policy'
    rescue_assessment_status: str = 'not_performed'
    index_status: str = 'derived_rebuildable'

    def __post_init__(self):
        if (self.schema!='genosapiens.wgs-analysis-eligibility-manifest' or type(self.version) is not int or self.version!=3
                or self.status!='complete' or self.analysis_domain_status!='complete_for_policy'
                or self.rescue_assessment_status!='not_performed' or self.index_status!='derived_rebuildable'
                or any(not _sha(v) for v in (self.admission_binding_sha256,self.eligibility_policy_sha256,
                    self.source_manifest_sha256,self.reference_manifest_sha256,self.ordered_decisions_sha256,
                    self.ordered_original_inputs_sha256,self.ordered_original_proofs_sha256,self.ordered_screening_sha256))
                or type(self.input_rows) is not int or self.input_rows<0
                or type(self.source_occurrences) is not int or self.source_occurrences<self.input_rows):
            raise AnalysisEligibilityError('invalid eligibility manifest binding')
        for counts,total in ((self.input_dispositions,self.input_rows),(self.occurrence_dispositions,self.source_occurrences)):
            if (type(counts) is not tuple or any(type(v) is not tuple or len(v)!=2 or type(v[1]) is not int or v[1]<0 for v in counts)
                    or tuple(k for k,v in counts)!=tuple(sorted(DISPOSITIONS)) or sum(v for k,v in counts)!=total):
                raise AnalysisEligibilityError('invalid eligibility disposition reconciliation')
        previous=0
        if type(self.shards) is not tuple: raise AnalysisEligibilityError('invalid eligibility shard collection')
        for number,shard in enumerate(self.shards,1):
            if (type(shard) is not EligibilityShard or shard.first_ordinal!=previous+1
                    or shard.relative_path!=f'decisions-{number:08d}.jsonl.gz'):
                raise AnalysisEligibilityError('invalid ordered eligibility shard coverage')
            previous=shard.last_ordinal
        if previous!=self.input_rows: raise AnalysisEligibilityError('incomplete eligibility shard coverage')
        previous=0
        if type(self.witness_count) is not int or self.witness_count < 0 or type(self.witness_shards) is not tuple:
            raise AnalysisEligibilityError('invalid eligibility witness coverage')
        for number,shard in enumerate(self.witness_shards,1):
            if (type(shard) is not EligibilityShard or shard.first_ordinal != previous+1
                    or shard.relative_path != f'witnesses-{number:08d}.jsonl.gz'):
                raise AnalysisEligibilityError('invalid eligibility witness shard coverage')
            previous=shard.last_ordinal
        if previous != self.witness_count:
            raise AnalysisEligibilityError('incomplete eligibility witness shards')
        if (type(self.dependencies) is not tuple or not self.dependencies or self.dependencies!=tuple(sorted(self.dependencies))
                or any(type(v) is not tuple or len(v)!=4 or type(v[0]) is not str or not v[0] or type(v[1]) is not str
                    or not (v[1].startswith('/') if v[0]!='screening_context' else v[1]=='screening-context.json')
                    or '..' in v[1].split('/') or not _sha(v[2]) or type(v[3]) is not int or v[3]<1 for v in self.dependencies)
                or len({v[0] for v in self.dependencies})!=len(self.dependencies)
                or type(self.resource_bindings) is not tuple or self.resource_bindings!=tuple(sorted(self.resource_bindings))
                or any(type(v) is not tuple or len(v)!=2 or type(v[0]) is not str or not v[0] or not _sha(v[1]) for v in self.resource_bindings)
                or len(dict(self.resource_bindings))!=len(self.resource_bindings)):
            raise AnalysisEligibilityError('invalid retained eligibility dependencies')

    def as_bytes(self): return _canonical(asdict(self))

    @property
    def sha256(self): return hashlib.sha256(self.as_bytes()).hexdigest()


def load_analysis_eligibility_manifest(content, expected_sha256):
    try:
        if type(content) is not bytes or not 0<len(content)<=16*1024*1024 or not _sha(expected_sha256): raise ValueError
        if hashlib.sha256(content).hexdigest()!=expected_sha256: raise ValueError
        data=json.loads(content,object_pairs_hook=_unique)
        if type(data) is not dict or set(data)!=set(AnalysisEligibilityManifest.__dataclass_fields__): raise ValueError
        for field in ('input_dispositions','occurrence_dispositions','dependencies','resource_bindings'):
            if type(data[field]) is not list or any(type(v) is not list for v in data[field]): raise ValueError
            data[field]=tuple(tuple(v) for v in data[field])
        for name in ('shards','witness_shards'):
            if type(data[name]) is not list: raise ValueError
            data[name]=tuple(EligibilityShard(**v) for v in data[name])
        result=AnalysisEligibilityManifest(**data)
        if result.as_bytes()!=content: raise ValueError
        return result
    except (TypeError,ValueError,UnicodeError):
        raise AnalysisEligibilityError('invalid serialized eligibility manifest') from None


def _eligibility_shard_writer(root,config,*,prefix='decisions'):
    """Reuse scientific writer row/byte bounds and shard-boundary persistence."""
    import gzip,os,time
    from genosapiens.wgs_ledger import _HashingWriter
    from genosapiens.wgs_shards import _ShardWriter,_bytes
    class CompressedWriter(_HashingWriter):
        def fileno(self):
            return self.stream.fileno()

    class Writer(_ShardWriter):
        def __init__(self,*args,**kwargs):
            super().__init__(*args,**kwargs)
            self.metrics=dict(shards=0,shard_close_nanoseconds=0,fsync_nanoseconds=0)
        def _open(self):
            relative=f'{prefix}-{len(self.descriptors)+1:08d}.jsonl.gz'
            self._relative=relative; self._path=self.root/relative
            self.raw=os.fdopen(os.open(self._path,os.O_WRONLY|os.O_CREAT|os.O_EXCL,0o600),'wb')
            self.compressed=CompressedWriter(self.raw)
            self._handle=gzip.GzipFile(filename='',mode='wb',fileobj=self.compressed,mtime=0,compresslevel=1)
            self._digest=hashlib.sha256()
        def _close(self):
            if self._handle is None: return
            started=time.monotonic_ns()
            self._handle.close(); self.raw.flush()
            sync_started=time.monotonic_ns()
            try:os.fsync(self.raw.fileno())
            finally:self.metrics['fsync_nanoseconds']+=time.monotonic_ns()-sync_started
            self.raw.close()
            shard=EligibilityShard(self._relative,self.compressed.digest.hexdigest(),self.compressed.size,
                self._digest.hexdigest(),self._size,self._rows,int(self._first),int(self._last))
            self.descriptor_bytes+=len(_bytes(shard))
            if self.descriptor_bytes>8*1024*1024: raise AnalysisEligibilityError('eligibility shard metadata exceeds capacity')
            self.descriptors.append(shard)
            self.metrics['shards']+=1
            self.metrics['shard_close_nanoseconds']+=time.monotonic_ns()-started
            self._handle=self._digest=self._path=self._relative=None
            self._rows=self._size=0; self._first=self._last=None
    return Writer(root,'analysis_eligibility',filename_prefix=prefix,
        max_rows=config.source_shards.max_uncompressed_rows,max_bytes=config.source_shards.max_uncompressed_bytes)


def _open_eligibility_index(path,config):
    from genosapiens.wgs_contexts import _open_index
    db=_open_index(path,config)
    db.executescript('''
        CREATE TABLE binding(manifest_sha256 TEXT NOT NULL);
        CREATE TABLE witnesses(sha256 TEXT PRIMARY KEY);
        CREATE TABLE decisions(ordinal INTEGER PRIMARY KEY,identity_id TEXT NOT NULL UNIQUE,
            disposition TEXT NOT NULL,decision_json BLOB NOT NULL);
        CREATE TABLE sources(source_allele_id TEXT PRIMARY KEY,ordinal INTEGER NOT NULL REFERENCES decisions(ordinal));
        CREATE TABLE genes(ordinal INTEGER NOT NULL REFERENCES decisions(ordinal),gene_id TEXT);
        CREATE INDEX genes_by_id ON genes(gene_id,ordinal);
        CREATE TABLE sequences(ordinal INTEGER NOT NULL REFERENCES decisions(ordinal),identity_id TEXT NOT NULL,
            contig TEXT NOT NULL,position INTEGER NOT NULL,reference TEXT NOT NULL,alternate TEXT NOT NULL);
        CREATE INDEX sequences_exact ON sequences(contig,position,reference,alternate,ordinal);
    ''')
    db.commit()
    return db


class AnalysisEligibilityWriter:
    """Durable complete decisions and an independently rebuildable index."""
    def __init__(self,root,*,execution_config,admission_binding_sha256,eligibility_policy_sha256,
                 source_manifest_sha256,reference_manifest_sha256,expected_inputs,expected_occurrences,
                 dependencies,resource_bindings):
        from pathlib import Path
        self.root=Path(root)
        if self.root.exists(): raise AnalysisEligibilityError('eligibility scratch already exists')
        self.root.mkdir(mode=0o700)
        self.config=execution_config; self.expected_inputs=expected_inputs; self.expected_occurrences=expected_occurrences
        self.bindings=(admission_binding_sha256,eligibility_policy_sha256,source_manifest_sha256,reference_manifest_sha256)
        self.dependencies=dependencies; self.resources=resource_bindings
        if any(not _sha(v) for v in self.bindings) or any(type(v) is not int or v<0 for v in (expected_inputs,expected_occurrences)):
            raise AnalysisEligibilityError('invalid eligibility writer binding')
        self.writer=_eligibility_shard_writer(self.root,execution_config)
        self.witness_writer=_eligibility_shard_writer(self.root,execution_config,prefix='witnesses')
        self.witness_count=0
        self.db=_open_eligibility_index(self.root/'index.sqlite',execution_config)
        self.input_digest=hashlib.sha256(); self.proof_digest=hashlib.sha256()
        self.screening_digest=hashlib.sha256()
        self.rows=self.occurrences=0; self.counts={k:0 for k in DISPOSITIONS}; self.occurrence_counts={k:0 for k in DISPOSITIONS}
        self.finished=False; self.manifest=None; self.windows=0
        self.metrics=dict(commits=0,commit_nanoseconds=0,write_nanoseconds=0,
            seal_nanoseconds=0,committed_inputs=0,committed_occurrences=0)

    def _commit(self):
        import time
        started=time.monotonic_ns()
        try:self.db.commit()
        finally:self.metrics['commit_nanoseconds']+=time.monotonic_ns()-started
        self.metrics['commits']+=1
        self.metrics['committed_inputs']=self.rows
        self.metrics['committed_occurrences']=self.occurrences

    @_fail_closed_operation
    def add_window(self,records,batch,*,witnesses):
        from genosapiens.wgs_admission_projection import AdmissionScreeningRecord, AdmissionProjectionError
        self._begin_window(records,batch,witnesses=witnesses)
        for encoded,decision in zip(records,batch.decisions,strict=True):
            try:
                compact=AdmissionScreeningRecord.from_encoded(encoded,decision)
            except AdmissionProjectionError:
                raise AnalysisEligibilityError('eligibility original input/proof coverage mismatch') from None
            self._add_authenticated(encoded,compact,compact.as_bytes())
        self._commit_window()

    def _begin_window(self,records,batch,*,witnesses):
        if (self.finished or type(records) is not tuple or type(batch) is not EligibilityBatch
                or len(records)!=len(batch.decisions) or not records
                or len(records)>self.config.public_queries.batch_rows
                or type(witnesses) not in (dict,MappingProxyType)):
            raise AnalysisEligibilityError('invalid persisted eligibility window')
        required={ref[2] for row in batch.decisions for ref in row.proof_references}
        if set(witnesses)!=required: raise AnalysisEligibilityError('eligibility witness coverage mismatch')
        for digest,content in witnesses.items():
            if (not _sha(digest) or type(content) is not bytes or not 0<len(content)<=MAX_ROW_BYTES
                    or hashlib.sha256(content).hexdigest()!=digest):
                raise AnalysisEligibilityError('eligibility witness digest mismatch')
            if not self.db.execute('SELECT 1 FROM witnesses WHERE sha256=?',(digest,)).fetchone():
                payload=json.loads(content,object_pairs_hook=_unique)
                if _canonical(payload)!=content:
                    raise AnalysisEligibilityError('noncanonical eligibility witness')
                self.witness_count+=1
                self.witness_writer.write(dict(sha256=digest,content=payload),key=str(self.witness_count))
                self.db.execute('INSERT INTO witnesses VALUES (?)',(digest,))

    def _add_authenticated(self,encoded,compact,compact_bytes):
        # Called only with a fresh private record from this operation, including
        # joint publication. Every store-specific binding and order check remains.
        from genosapiens.wgs_shards import encode_row
        decision=compact.decision
        if (decision.retained_ordinal!=self.rows+1 or decision.retained_ordinal!=encoded.retained_ordinal
                or decision.identity_id!=encoded.identity_id or decision.allele_id!=encoded.allele_id
                or decision.source_allele_ids!=compact.facts.source_allele_ids
                or decision.input_sha256!=encoded.input_sha256 or decision.proof_sha256!=encoded.proof_sha256
                or (decision.admission_binding_sha256,decision.eligibility_policy_sha256,
                    decision.source_manifest_sha256,decision.reference_manifest_sha256)!=self.bindings
                or decision.resource_bindings!=self.resources):
            raise AnalysisEligibilityError('eligibility original input/proof coverage mismatch')
        canonical=encode_row(decision)
        self.writer.write(decision,key=str(decision.retained_ordinal),canonical_bytes=canonical)
        self.db.execute('INSERT INTO decisions VALUES (?,?,?,?)',(decision.retained_ordinal,decision.identity_id,decision.disposition,canonical.content))
        self.db.executemany('INSERT INTO sources VALUES (?,?)',((key,decision.retained_ordinal) for key in decision.source_allele_ids))
        self.db.executemany('INSERT INTO genes VALUES (?,?)',((decision.retained_ordinal,key) for key in (encoded.gene_ids or (None,))))
        self.db.executemany('INSERT INTO sequences VALUES (?,?,?,?,?,?)',((decision.retained_ordinal,key,
            sequence.contig,sequence.position,sequence.reference,sequence.alternate)
            for key,sequence in compact.facts.sequence_identities))
        self.rows+=1; self.occurrences+=len(decision.source_allele_ids)
        self.counts[decision.disposition]+=1; self.occurrence_counts[decision.disposition]+=len(decision.source_allele_ids)
        self.input_digest.update((encoded.input_sha256+'\n').encode()); self.proof_digest.update((encoded.proof_sha256+'\n').encode())
        self.screening_digest.update(compact_bytes)

    def _commit_window(self):
        self._commit(); self.windows+=1

    @_fail_closed_operation
    def finish(self):
        from genosapiens.wgs_ledger import publish_bytes
        if self.finished or self.rows!=self.expected_inputs or self.occurrences!=self.expected_occurrences:
            raise AnalysisEligibilityError('incomplete persisted eligibility coverage')
        shards=self.writer.close()
        witnesses=self.witness_writer.close()
        manifest=AnalysisEligibilityManifest(*self.bindings,self.rows,self.occurrences,tuple(sorted(self.counts.items())),
            tuple(sorted(self.occurrence_counts.items())),self.writer.stream.hexdigest(),self.input_digest.hexdigest(),
            self.proof_digest.hexdigest(),self.screening_digest.hexdigest(),shards,witnesses,self.witness_count,
            tuple((role,'screening-context.json' if role=='screening_context' else path,digest,size)
                for role,path,digest,size in self.dependencies),self.resources)
        # The index is rebuildable. SQLite enforces its keys during insertion.
        # Required coverage and digests close above, without a second index scan.
        from pathlib import Path
        from genosapiens.wgs_shards import read_bounded_metadata
        retained=set()
        for role,path,digest,size in self.dependencies:
            content=read_bounded_metadata(Path(path))
            if len(content)!=size or hashlib.sha256(content).hexdigest()!=digest:
                raise AnalysisEligibilityError('original eligibility dependency changed before sealing')
            if digest not in retained:
                publish_bytes(self.root,f'dependency-{digest}.json',content)
                retained.add(digest)
        self.db.execute('INSERT INTO binding VALUES (?)',(manifest.sha256,)); self._commit()
        publish_bytes(self.root,'manifest.json',manifest.as_bytes())
        rescue=_canonical(dict(schema='genosapiens.wgs-refcall-rescue-manifest',version=1,
            eligibility_manifest_path='manifest.json',eligibility_manifest_sha256=manifest.sha256,
            deferred_inputs=self.counts['defer_refcall'],deferred_occurrences=self.occurrence_counts['defer_refcall'],
            assessment_status='not_performed',index_path='index.sqlite',index_status='derived_rebuildable',
            query_contract='original-source-id/original-ordinal/exact-sequence/gene-or-unknown-v1',
            dependencies=manifest.dependencies))
        publish_bytes(self.root,'rescue-manifest.json',rescue)
        self.finished=True; self.manifest=manifest
        return manifest

    def close(self):
        # Closing an unfinished owner releases resources without sealing its
        # partial shards. Only finish() can publish their durable descriptors.
        for writer in (self.writer,self.witness_writer):
            if writer._handle is not None:
                writer._handle.close();writer.raw.close();writer._handle=None
        self.db.close()


def iter_eligibility_decisions(root,manifest,*,scratch_directory=None):
    """Bounded replay with independent cross-shard identity/source uniqueness."""
    from pathlib import Path
    import sqlite3,tempfile
    from contextlib import closing
    with tempfile.TemporaryDirectory(prefix='.eligibility-replay-',dir=Path(root).parent if scratch_directory is None else scratch_directory) as temporary:
        with closing(sqlite3.connect(Path(temporary)/'uniqueness.sqlite')) as db:
            db.execute('PRAGMA journal_mode=DELETE'); db.execute('PRAGMA synchronous=FULL')
            db.execute('PRAGMA cache_size=-65536'); db.execute('PRAGMA mmap_size=0')
            db.executescript('CREATE TABLE identities(id TEXT PRIMARY KEY);CREATE TABLE sources(id TEXT PRIMARY KEY);'
                'CREATE TABLE witnesses(sha256 TEXT PRIMARY KEY,used INTEGER NOT NULL);')
            _load_eligibility_witnesses(root,manifest,db)
            with closing(_iter_eligibility_decisions(root,manifest,witness_db=db)) as rows:
                for row in rows:
                    try:
                        db.execute('INSERT INTO identities VALUES (?)',(row.identity_id,))
                        db.executemany('INSERT INTO sources VALUES (?)',((key,) for key in row.source_allele_ids))
                    except sqlite3.IntegrityError:
                        raise AnalysisEligibilityError('duplicate eligibility identity or source occurrence') from None
                    if row.retained_ordinal%4096==0: db.commit()
                    yield row
            if db.execute('SELECT 1 FROM witnesses WHERE used=0 LIMIT 1').fetchone():
                raise AnalysisEligibilityError('unreferenced eligibility witness')


def _load_eligibility_witnesses(root,manifest,db):
    """Authenticate bounded witness shards into a disposable membership index."""
    import gzip,os,stat,sqlite3
    from pathlib import Path
    count=0
    for shard in manifest.witness_shards:
        path=Path(root)/shard.relative_path
        if any(p.is_symlink() for p in (path,*path.parents)):
            raise AnalysisEligibilityError('eligibility witness shard path changed')
        with os.fdopen(os.open(path,os.O_RDONLY|os.O_NOFOLLOW),'rb') as raw:
            before=os.fstat(raw.fileno()); compressed=hashlib.sha256()
            if not stat.S_ISREG(before.st_mode) or before.st_size!=shard.size_bytes:
                raise AnalysisEligibilityError('eligibility witness shard size mismatch')
            while part:=raw.read(1024*1024): compressed.update(part)
            if compressed.hexdigest()!=shard.sha256:
                raise AnalysisEligibilityError('eligibility witness shard digest mismatch')
            raw.seek(0); digest=hashlib.sha256();size=rows=0
            with gzip.GzipFile(fileobj=raw,mode='rb') as handle:
                while content:=handle.readline(MAX_ROW_BYTES+1):
                    size+=len(content);rows+=1
                    if len(content)>MAX_ROW_BYTES or size>shard.uncompressed_size_bytes or rows>shard.row_count:
                        raise AnalysisEligibilityError('eligibility witness shard exceeds capacity')
                    item=json.loads(content,object_pairs_hook=_unique)
                    if (type(item) is not dict or set(item)!={'sha256','content'}
                            or _canonical(item)!=content or not _sha(item['sha256'])
                            or hashlib.sha256(_canonical(item['content'])).hexdigest()!=item['sha256']):
                        raise AnalysisEligibilityError('eligibility witness content mismatch')
                    try: db.execute('INSERT INTO witnesses VALUES (?,0)',(item['sha256'],))
                    except sqlite3.IntegrityError:
                        raise AnalysisEligibilityError('duplicate eligibility witness') from None
                    count+=1;digest.update(content)
                    if count%4096==0: db.commit()
            after=os.fstat(raw.fileno())
            if (rows!=shard.row_count or size!=shard.uncompressed_size_bytes or digest.hexdigest()!=shard.uncompressed_sha256
                    or (before.st_dev,before.st_ino,before.st_size,before.st_mtime_ns,before.st_ctime_ns)!=
                       (after.st_dev,after.st_ino,after.st_size,after.st_mtime_ns,after.st_ctime_ns)):
                raise AnalysisEligibilityError('eligibility witness final reconciliation failed')
    if count!=manifest.witness_count:
        raise AnalysisEligibilityError('eligibility witness coverage is incomplete')
    db.commit()


def _iter_eligibility_decisions(root,manifest,*,witness_db):
    """Bounded strict replay; count/digest closure is required on exhaustion."""
    import gzip,os,stat
    from pathlib import Path
    from genosapiens.wgs_shards import read_bounded_metadata
    root=Path(root)
    if type(manifest) is not AnalysisEligibilityManifest: raise AnalysisEligibilityError('invalid eligibility replay manifest')
    stream=hashlib.sha256(); inputs=hashlib.sha256(); proofs=hashlib.sha256()
    count=occurrences=0; totals={k:0 for k in DISPOSITIONS}; occurrence_totals={k:0 for k in DISPOSITIONS}
    for shard in manifest.shards:
        path=root/shard.relative_path
        if any(p.is_symlink() for p in (path,*path.parents)): raise AnalysisEligibilityError('eligibility shard path changed')
        with os.fdopen(os.open(path,os.O_RDONLY|os.O_NOFOLLOW),'rb') as raw:
            before=os.fstat(raw.fileno()); compressed=hashlib.sha256()
            if not stat.S_ISREG(before.st_mode) or before.st_size!=shard.size_bytes:
                raise AnalysisEligibilityError('eligibility compressed shard size mismatch')
            while part:=raw.read(1024*1024): compressed.update(part)
            if compressed.hexdigest()!=shard.sha256: raise AnalysisEligibilityError('eligibility compressed shard digest mismatch')
            raw.seek(0); digest=hashlib.sha256(); size=rows=0
            with gzip.GzipFile(fileobj=raw,mode='rb') as handle:
                while content:=handle.readline(1024*1024+1):
                    size+=len(content); rows+=1
                    if len(content)>1024*1024 or size>shard.uncompressed_size_bytes or rows>shard.row_count:
                        raise AnalysisEligibilityError('eligibility shard exceeds declared capacity')
                    decision=load_analysis_eligibility_decision(content)
                    if (decision.retained_ordinal!=count+1 or
                            (decision.admission_binding_sha256,decision.eligibility_policy_sha256,
                             decision.source_manifest_sha256,decision.reference_manifest_sha256)!=
                            (manifest.admission_binding_sha256,manifest.eligibility_policy_sha256,
                             manifest.source_manifest_sha256,manifest.reference_manifest_sha256)
                            or decision.resource_bindings!=manifest.resource_bindings):
                        raise AnalysisEligibilityError('eligibility replay original binding mismatch')
                    for role,collection,expected in decision.proof_references:
                        if not witness_db.execute('SELECT 1 FROM witnesses WHERE sha256=?',(expected,)).fetchone():
                            raise AnalysisEligibilityError('eligibility witness is missing from sealed shards')
                        witness_db.execute('UPDATE witnesses SET used=1 WHERE sha256=? AND used=0',(expected,))
                    digest.update(content); stream.update(content); count+=1
                    occurrences+=len(decision.source_allele_ids); totals[decision.disposition]+=1
                    occurrence_totals[decision.disposition]+=len(decision.source_allele_ids)
                    inputs.update((decision.input_sha256+'\n').encode()); proofs.update((decision.proof_sha256+'\n').encode())
                    yield decision
            after=os.fstat(raw.fileno())
            if (rows!=shard.row_count or size!=shard.uncompressed_size_bytes or digest.hexdigest()!=shard.uncompressed_sha256
                    or (before.st_dev,before.st_ino,before.st_size,before.st_mtime_ns,before.st_ctime_ns)!=
                       (after.st_dev,after.st_ino,after.st_size,after.st_mtime_ns,after.st_ctime_ns)):
                raise AnalysisEligibilityError('eligibility shard final reconciliation failed')
    if (count!=manifest.input_rows or occurrences!=manifest.source_occurrences or tuple(sorted(totals.items()))!=manifest.input_dispositions
            or tuple(sorted(occurrence_totals.items()))!=manifest.occurrence_dispositions
            or (stream.hexdigest(),inputs.hexdigest(),proofs.hexdigest())!=(manifest.ordered_decisions_sha256,
                manifest.ordered_original_inputs_sha256,manifest.ordered_original_proofs_sha256)):
        raise AnalysisEligibilityError('eligibility replay final reconciliation failed')


def iter_analysis_eligibility(windows,*,session,writer):
    """Evaluate inside the original decode pass; publish only at full closure."""
    if type(session) is not AnalysisEligibilitySession or type(writer) is not AnalysisEligibilityWriter:
        raise AnalysisEligibilityError('invalid eligibility stream owners')
    for window in windows:
        batch=session.evaluate_window(window)
        writer.add_window(window,batch,witnesses=session.witnesses)
        yield batch
    session.finish(expected_inputs=writer.expected_inputs,expected_occurrences=writer.expected_occurrences)
    writer.finish()


def iter_refcall_rescue_candidates(root,*,manifest):
    """Original source IDs remain addressable; this view confers no carriage."""
    for decision in iter_eligibility_decisions(root,manifest):
        if decision.disposition=='defer_refcall': yield decision


def rebuild_eligibility_index(root,*,manifest,admission_reader,execution_config):
    """Recover only the derived index from sealed decisions and original rows."""
    from contextlib import closing
    from pathlib import Path
    from genosapiens.wgs_admission_projection import AdmissionProjectionReader
    if (type(admission_reader) is not AdmissionProjectionReader or admission_reader.binding_sha256!=manifest.admission_binding_sha256
            or admission_reader.row_count!=manifest.input_rows
            or admission_reader.completion.screening_sha256!=manifest.ordered_screening_sha256):
        raise AnalysisEligibilityError('foreign eligibility rebuild input')
    with closing(_open_eligibility_index(Path(root)/'index.sqlite',execution_config)) as db, \
            closing(iter_eligibility_decisions(root,manifest)) as decisions,closing(admission_reader.iter_screening_records()) as rows:
        for count,(decision,record) in enumerate(zip(decisions,rows,strict=True),1):
            if decision != record.decision:
                raise AnalysisEligibilityError('eligibility rebuild source coverage mismatch')
            db.execute('INSERT INTO decisions VALUES (?,?,?,?)',(decision.retained_ordinal,decision.identity_id,decision.disposition,decision.as_bytes()))
            db.executemany('INSERT INTO sources VALUES (?,?)',((key,decision.retained_ordinal) for key in decision.source_allele_ids))
            db.executemany('INSERT INTO genes VALUES (?,?)',((decision.retained_ordinal,key) for key in (record.facts.gene_ids or (None,))))
            db.executemany('INSERT INTO sequences VALUES (?,?,?,?,?,?)',((decision.retained_ordinal,key,
                sequence.contig,sequence.position,sequence.reference,sequence.alternate)
                for key,sequence in record.facts.sequence_identities))
            if count%execution_config.public_queries.batch_rows==0: db.commit()
        if db.execute('PRAGMA foreign_key_check').fetchall() or db.execute('PRAGMA integrity_check').fetchall()!=[('ok',)]:
            raise AnalysisEligibilityError('rebuilt eligibility index integrity failure')
        db.execute('INSERT INTO binding VALUES (?)',(manifest.sha256,)); db.commit()


def cleanup_eligibility_index(root,*,expected_manifest_sha256):
    """Remove only the rebuildable cache; the exact eligibility subtree lives on."""
    import os,stat
    from pathlib import Path
    from genosapiens.wgs_shards import read_bounded_metadata
    root=Path(root)
    manifest=load_analysis_eligibility_manifest(read_bounded_metadata(root/'manifest.json'),expected_manifest_sha256)
    rescue=json.loads(read_bounded_metadata(root/'rescue-manifest.json'),object_pairs_hook=_unique)
    if rescue.get('eligibility_manifest_sha256')!=manifest.sha256 or rescue.get('index_path')!='index.sqlite':
        raise AnalysisEligibilityError('foreign eligibility cleanup view')
    path=root/'index.sqlite'
    try: info=path.lstat()
    except FileNotFoundError: return
    if (any(p.is_symlink() for p in (path,*path.parents)) or not stat.S_ISREG(info.st_mode)
            or info.st_nlink!=1 or info.st_uid!=os.getuid() or stat.S_IMODE(info.st_mode)!=0o600
            or any((root/name).exists() for name in ('index.sqlite-journal','index.sqlite-wal','index.sqlite-shm'))):
        raise AnalysisEligibilityError('eligibility index is unsafe or active')
    path.unlink()


def analysis_eligibility_component(root, manifest):
    """Small scientific reference to retained evidence beside its science root."""
    from pathlib import Path
    from genosapiens.wgs_shards import read_bounded_metadata
    root=Path(root)
    if type(manifest) is not AnalysisEligibilityManifest or root.name!='eligibility':
        raise AnalysisEligibilityError('invalid scientific eligibility root')
    rescue=read_bounded_metadata(root/'rescue-manifest.json')
    value=dict(schema='genosapiens.wgs-analysis-eligibility-component',version=1,
        policy_sha256=manifest.eligibility_policy_sha256,
        manifest_path='eligibility/manifest.json',manifest_sha256=manifest.sha256,
        rescue_manifest_path='eligibility/rescue-manifest.json',rescue_manifest_sha256=hashlib.sha256(rescue).hexdigest(),
        input_rows=manifest.input_rows,source_occurrences=manifest.source_occurrences,
        input_dispositions=manifest.input_dispositions,occurrence_dispositions=manifest.occurrence_dispositions,
        analysis_domain_status=manifest.analysis_domain_status,rescue_assessment_status=manifest.rescue_assessment_status)
    verify_analysis_eligibility_component(root.parent,value,full=False)
    return value


def eligibility_component_parents(value):
    return dict(analysis_eligibility_policy=value['policy_sha256'],
        analysis_eligibility_manifest=value['manifest_sha256'],
        analysis_eligibility_rescue_manifest=value['rescue_manifest_sha256'])


def verify_analysis_eligibility_component(scope,value,*,full=True,scratch_directory=None,eligibility_root=None):
    """Authenticate bounded roots, original descriptor references and witnesses."""
    from pathlib import Path
    from genosapiens.wgs_shards import read_bounded_metadata
    scope=Path(scope)
    keys={'schema','version','policy_sha256','manifest_path','manifest_sha256','rescue_manifest_path',
        'rescue_manifest_sha256','input_rows','source_occurrences','input_dispositions','occurrence_dispositions',
        'analysis_domain_status','rescue_assessment_status'}
    if (type(value) is not dict or set(value)!=keys
            or value['schema']!='genosapiens.wgs-analysis-eligibility-component'
            or type(value['version']) is not int or value['version']!=1
            or value['manifest_path']!='eligibility/manifest.json'
            or value['rescue_manifest_path']!='eligibility/rescue-manifest.json'
            or any(not _sha(value[k]) for k in ('policy_sha256','manifest_sha256','rescue_manifest_sha256'))):
        raise AnalysisEligibilityError('invalid scientific eligibility component')
    root=scope/'eligibility' if eligibility_root is None else Path(eligibility_root)
    if (not root.is_absolute() or root.name!='eligibility'
            or any(p.is_symlink() for p in (root,*root.parents))):
        raise AnalysisEligibilityError('invalid retained eligibility root')
    manifest=load_analysis_eligibility_manifest(read_bounded_metadata(root/'manifest.json'),value['manifest_sha256'])
    if (manifest.eligibility_policy_sha256!=value['policy_sha256']
            or any(value[k]!=getattr(manifest,k) for k in ('input_rows','source_occurrences','analysis_domain_status','rescue_assessment_status'))
            or _canonical(value['input_dispositions'])!=_canonical(manifest.input_dispositions)
            or _canonical(value['occurrence_dispositions'])!=_canonical(manifest.occurrence_dispositions)):
        raise AnalysisEligibilityError('scientific eligibility coverage changed')
    rescue=read_bounded_metadata(root/'rescue-manifest.json')
    expected=_canonical(dict(schema='genosapiens.wgs-refcall-rescue-manifest',version=1,
        eligibility_manifest_path='manifest.json',eligibility_manifest_sha256=manifest.sha256,
        deferred_inputs=dict(manifest.input_dispositions)['defer_refcall'],
        deferred_occurrences=dict(manifest.occurrence_dispositions)['defer_refcall'],
        assessment_status='not_performed',index_path='index.sqlite',index_status='derived_rebuildable',
        query_contract='original-source-id/original-ordinal/exact-sequence/gene-or-unknown-v1',dependencies=manifest.dependencies))
    if rescue!=expected or hashlib.sha256(rescue).hexdigest()!=value['rescue_manifest_sha256']:
        raise AnalysisEligibilityError('scientific rescue coverage changed')
    for role,path,digest,size in manifest.dependencies:
        # Original paths remain provenance. The current publication retains
        # authenticated descriptor bytes so cleanup cannot remove its parents.
        content=read_bounded_metadata(root/f'dependency-{digest}.json')
        if len(content)!=size or hashlib.sha256(content).hexdigest()!=digest:
            raise AnalysisEligibilityError('original eligibility dependency changed')
    if full:
        for row in iter_eligibility_decisions(root,manifest,scratch_directory=scratch_directory): pass
    return manifest


def eligibility_retained_files(scope,value,*,eligibility_root=None):
    """Stream immutable leaves for the existing scientific publication inventory."""
    from pathlib import Path
    import os,re
    root=Path(scope)/'eligibility' if eligibility_root is None else Path(eligibility_root)
    manifest=verify_analysis_eligibility_component(scope,value,full=False,eligibility_root=root)
    for name,digest in (('manifest.json',manifest.sha256),('rescue-manifest.json',value['rescue_manifest_sha256'])):
        path=root/name
        yield path,'eligibility/'+name,digest,path.stat().st_size
    for shard in manifest.shards:
        yield root/shard.relative_path,'eligibility/'+shard.relative_path,shard.sha256,shard.size_bytes
    for shard in manifest.witness_shards:
        yield root/shard.relative_path,'eligibility/'+shard.relative_path,shard.sha256,shard.size_bytes
    for digest,size in sorted({(digest,size) for _,_,digest,size in manifest.dependencies}):
        name=f'dependency-{digest}.json'
        yield root/name,'eligibility/'+name,digest,size


def scientific_eligibility_component(artifact,directory,*,eligibility_root=None):
    """Reconstruct exact expected parents from their bounded retained roots.

    The existing scientific reader separately verifies the normalized component
    reference against these values; this does not open another science store.
    """
    from pathlib import Path
    from genosapiens.wgs_shards import read_bounded_metadata
    actual={k:v for k,v in artifact.parent_bindings if k.startswith('analysis_eligibility_')}
    if not actual: return None
    if set(actual)!={'analysis_eligibility_policy','analysis_eligibility_manifest','analysis_eligibility_rescue_manifest'}:
        raise AnalysisEligibilityError('incomplete scientific eligibility parents')
    root=Path(directory).parent/'eligibility' if eligibility_root is None else Path(eligibility_root)
    manifest=load_analysis_eligibility_manifest(read_bounded_metadata(root/'manifest.json'),actual['analysis_eligibility_manifest'])
    value=analysis_eligibility_component(root,manifest)
    if actual!=eligibility_component_parents(value):
        raise AnalysisEligibilityError('scientific eligibility parent reconstruction differs')
    return value
