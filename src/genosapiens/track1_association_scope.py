"""Focused association scope, distinct from biological support and scoring.

Original inputs and the complete regulatory index remain the source of record.
This owner selects rich evaluation/publication and records compact exclusions.
"""
from collections import Counter
from dataclasses import asdict, dataclass, replace
from decimal import Decimal, InvalidOperation
import hashlib

from genosapiens.evidence_model import SourcedRecord
from genosapiens.evidence_serialization import _owned_canonical_bytes
from genosapiens.wgs_ledger import canonical

SCHEMA = 'genosapiens.focused-association-scope'
POLICY = dict(schema=SCHEMA, version=2,
    molecular='exact-selected-transcript-source-and-target',
    expansion='positive-source-match-only; compact-other-source-accounting',
    grouping='same-subject-rule-adapter-state-target-condition-reasons-context-geometry-regulatory-support',
    regulatory='maximum-finite-score-per-element-raw-target-method',
    regulatory_ties='source-record-id-then-source-sha256',
    uncertainty='outside-competitive-scope; not-benign-or-no-partner',
    original_sources='complete-immutable-upstream-inputs-and-regulatory-index',
    curated='exact-functional-source-facts-from-sealed-retention-projection')
POLICY_SHA256 = hashlib.sha256(canonical(POLICY)).hexdigest()


def _digest(value):
    return hashlib.sha256(_owned_canonical_bytes(value)).hexdigest()


def source_locator(source):
    return canonical(dict(dataset=source.source_dataset, release=source.source_release,
        record_id=source.source_record_id, source_sha256=source.source_sha256,
        manifest_sha256=source.public_manifest_sha256)).decode().strip()


def _summary(kind, value):
    content = dict(schema=SCHEMA, version=1, policy_sha256=POLICY_SHA256, kind=kind, **value)
    digest = _digest(content)
    return SourcedRecord(SCHEMA, 'v1', kind+'-'+digest, digest, None, content)


def select_regulatory_sources(sources):
    """One representative per method/relationship; never compare methods."""
    from genosapiens.admission_links import ENCODE_RE2G
    retained, groups = [], {}
    commitments = hashlib.sha256()
    total = missing = 0
    prior = []
    for source in sources:
        if source.source_dataset == SCHEMA and source.source_fields.get("kind") == "regulatory-support-selection":
            prior.append(source)
        if source.source_dataset != ENCODE_RE2G:
            retained.append(source)
            continue
        row = source.source_fields['record']
        total += 1
        commitments.update(canonical((source.source_record_id, source.source_sha256,
            source.public_manifest_sha256, row['legacy_id'], row['source_target_id'],
            row['method'], row['score'])))
        try:
            score = Decimal(row['score'])
        except (InvalidOperation, ValueError, TypeError):
            missing += 1
            continue
        if not score.is_finite():
            missing += 1
            continue
        key = (row['legacy_id'], row['source_target_id'], row['method'])
        order = (source.source_record_id, source.source_sha256)
        previous = groups.get(key)
        if (previous is None or score > previous[0]
                or score == previous[0] and order < previous[1]):
            groups[key] = score, order, source
    if prior:
        if len(prior)!=1:
            raise ValueError('duplicate regulatory support selection')
        data=prior[0].source_fields
        counts={k:data[k] for k in ('total','selected','missing','omitted')}
        original=tuple(groups[key][2] for key in sorted(groups))
        selected_ids={s.source_record_id for s in original}
        if (data['policy_sha256']!=POLICY_SHA256 or prior[0].source_sha256!=_digest(data)
                or prior[0].source_record_id!='regulatory-support-selection-'+prior[0].source_sha256
                or prior[0].public_manifest_sha256 is not None
                or set(data['selected_source_records'])!=selected_ids
                or counts['selected']!=len(groups) or total!=len(groups) or missing
                or any(type(v) is not int or v<0 for v in counts.values())
                or counts['total']!=counts['selected']+counts['omitted']
                or counts['missing']>counts['omitted']):
            raise ValueError('preselected regulatory support differs')
        return (*retained,*original),counts
    if not total:
        return tuple(retained), dict(total=0, selected=0, missing=0, omitted=0)
    selected = tuple(groups[key][2] for key in sorted(groups))
    counts = dict(total=total, selected=len(selected), missing=missing, omitted=total-len(selected))
    summary = _summary('regulatory-support-selection', dict(**counts,
        source_locations_sha256=commitments.hexdigest(),
        selected_source_records=tuple(source.source_record_id for source in selected),
        interpretation='representative-support-only; no-tissue-or-causal-inference'))
    return (*retained, *selected, summary), counts


@dataclass(frozen=True, slots=True)
class _PendingRegulatorySource:
    """Private selection view of an already parsed, authenticated index row."""
    source_dataset: str
    source_release: str
    source_record_id: str
    source_sha256: str
    public_manifest_sha256: str
    source_fields: dict

    def materialize(self):
        return SourcedRecord(self.source_dataset,self.source_release,self.source_record_id,
            self.source_sha256,self.public_manifest_sha256,self.source_fields)


def select_parsed_regulatory_sources(rows, *, artifact_path, source_sha256, manifest_sha256):
    """Apply the same selector before rich construction of omitted index rows.

Only the regulatory index calls this with rows produced by the strict rE2G
parser after source-byte checks. No raw-input parser or source authentication
is skipped. Validate shared source metadata on the first row, even if all scores
are missing. All selected records use the ordinary SourcedRecord constructor.
"""
    from genosapiens.admission_links import ENCODE_RE2G
    def pending():
        first=True
        for ordinal,row in rows:
            source=_PendingRegulatorySource(ENCODE_RE2G,'Registry-V4',
                f'{artifact_path}#row={ordinal}',source_sha256,manifest_sha256,
                {'artifact_path':artifact_path,'row_ordinal':ordinal,'record':row})
            if first:
                first=False
                source.materialize()
            yield source
    selected,counts=select_regulatory_sources(pending())
    return tuple(source.materialize() if type(source) is _PendingRegulatorySource else source
        for source in selected),counts


@dataclass(frozen=True, slots=True)
class ScopedAssociations:
    associations: tuple
    raw_count: int
    surviving_proofs: int
    excluded_proofs: int
    grouped_count: int
    exclusion_summaries: int
    regulatory_rows: int
    regulatory_rows_selected: int


@dataclass(frozen=True, slots=True)
class AssociationScopeSummary:
    raw_count: int
    surviving_proofs: int
    excluded_proofs: int
    grouped_count: int
    exclusion_summaries: int
    regulatory_rows: int
    regulatory_rows_selected: int
    policy_sha256: str = POLICY_SHA256

    def __post_init__(self):
        counts = (self.raw_count,self.surviving_proofs,self.excluded_proofs,
            self.grouped_count,self.exclusion_summaries,self.regulatory_rows,self.regulatory_rows_selected)
        if (self.policy_sha256 != POLICY_SHA256 or any(type(v) is not int or v<0 for v in counts)
                or self.raw_count != self.surviving_proofs+self.excluded_proofs
                or self.grouped_count>self.surviving_proofs
                or self.exclusion_summaries>self.excluded_proofs
                or self.regulatory_rows_selected>self.regulatory_rows):
            raise ValueError('invalid focused association scope accounting')

    @classmethod
    def from_result(cls, result):
        return cls(*(getattr(result,name) for name in ('raw_count','surviving_proofs','excluded_proofs',
            'grouped_count','exclusion_summaries','regulatory_rows','regulatory_rows_selected')))


def validate_scoped_associations(profile, associations, summary):
    """Validate selected-domain meaning on independently decoded discovery."""
    from genosapiens.admission_links import ENCODE_RE2G
    from genosapiens.wgs_profile_projection import profile_association_applicability
    if type(summary) is not AssociationScopeSummary:
        raise ValueError('focused association scope summary is missing')
    summary.__post_init__()
    expected_curated = tuple(v for v in curated_associations(profile) if v.state == 'match')
    actual_curated = tuple(v for v in associations if v.adapter == 'curated' and v.state == 'match')
    def source_entries(values):
        return sorted(_digest(replace(value, source_refs=(source,)))
                      for value in values for source in value.source_refs)
    if source_entries(expected_curated) != source_entries(actual_curated):
        raise ValueError('curated association differs from retained exact source facts')
    selected = replace(profile, decision={**profile.decision, 'routes':tuple(
        row for row in profile.decision['routes'] if row['outcome']=='selected')})
    routes = selected.decision['routes']
    native = any(row['route']=='native' for row in routes)
    molecular = {(row['target_namespace'],row['target_id'],locator)
        for row in routes if row['route']=='molecular' for locator in row['source_refs']}
    positive = excluded = exclusion_count = 0
    for value, allowed in zip(associations,profile_association_applicability(selected,associations),strict=True):
        if value.state=='match':
            if not allowed:
                raise ValueError('focused association is outside selected target scope')
            if value.adapter=='annotation' and not native:
                target=value.target.split('.')[0] if value.namespace=='Ensembl' and value.target else value.target
                if not any((value.namespace,target,source_locator(source)) in molecular for source in value.source_refs):
                    raise ValueError('focused association lacks a qualifying transcript proof')
                for source in value.source_refs:
                    if (source.source_dataset=='Ensembl VEP' and (source.source_fields or {}).get('csq') is not None
                            and (value.namespace,target,source_locator(source)) not in molecular):
                        raise ValueError('focused association carries an excluded transcript proof')
            if value.adapter=='encode':
                seen=set(); selected_ids=set(); selections=[]
                for source in value.source_refs:
                    if source.source_dataset==SCHEMA:
                        selections.append(source)
                    if source.source_dataset!=ENCODE_RE2G:continue
                    row=source.source_fields['record'];key=(row['legacy_id'],row['source_target_id'],row['method'])
                    if key in seen or not Decimal(row['score']).is_finite():
                        raise ValueError('focused regulatory support selection differs')
                    seen.add(key); selected_ids.add(source.source_record_id)
                if not seen or len(selections)!=1:
                    raise ValueError('focused regulatory support is missing')
                selection=selections[0]; data=selection.source_fields
                if (data['kind']!='regulatory-support-selection' or data['policy_sha256']!=POLICY_SHA256
                        or selection.source_sha256!=_digest(data)
                        or selection.source_record_id!='regulatory-support-selection-'+selection.source_sha256
                        or selection.public_manifest_sha256 is not None
                        or set(data['selected_source_records'])!=selected_ids
                        or data['selected']!=len(seen)
                        or any(type(data[k]) is not int or data[k]<0 for k in ('total','selected','missing','omitted'))
                        or data['total']!=data['selected']+data['omitted'] or data['missing']>data['omitted']):
                    raise ValueError('focused regulatory selection accounting differs')
            positive+=1
        else:
            if (value.state!='unresolved' or value.target is not None or value.condition is not None
                    or value.reasons!=('focused-association-scope-excluded','outside-competitive-scope-not-benign')
                    or len(value.source_refs)!=1):
                raise ValueError('focused excluded route is not compact scope accounting')
            source=value.source_refs[0];data=source.source_fields
            if (source.source_dataset!=SCHEMA or source.public_manifest_sha256 is not None
                    or data['policy_sha256']!=POLICY_SHA256 or data['kind']!='source-scope-exclusion'
                    or data['subject']!=value.subject or data['rule']!=value.rule or data['adapter']!=value.adapter
                    or tuple(data['original_input_ref'])!=profile.original_ref
                    or data['screening_sha256']!=profile.screening_sha256
                    or data['profile_policy_sha256']!=profile.decision['policy_sha256']
                    or source.source_sha256!=_digest(data)
                    or source.source_record_id!='source-scope-exclusion-'+source.source_sha256
                    or sum(data['reasons'].values())!=data['excluded_associations']):
                raise ValueError('focused scope exclusion source binding differs')
            excluded+=1;exclusion_count+=data['excluded_associations']
    if (positive,excluded,exclusion_count)!=(summary.grouped_count,summary.exclusion_summaries,summary.excluded_proofs):
        raise ValueError('focused scope row coverage differs')


def curated_associations(profile):
    """Expose exact functional targets without inventing clinical classifications."""
    from genosapiens import candidate_admission as admission
    from genosapiens.wgs_profile_projection import _facts
    from genosapiens.variant_retention import RetentionFacts, select_variant
    facts = _facts(profile.facts)
    if profile.decision['retention_disposition'] != 'retained':
        raise ValueError('deferred input cannot construct curated associations')
    subject = profile.original_ref[3]
    if subject is None:
        subjects = profile.source_allele_ids
    else:
        subjects = (subject,)
    result, seen = [], set()
    for route in facts.routes:
        for evidence in route.evidence:
            fact = evidence.functional
            if fact is None:
                continue
            selected = select_variant(RetentionFacts(facts.identity_id, facts.exact_allele, 'assess',
                facts.source_occurrences, (), (), (fact,), (), facts.public_coverage),
                established_genes=frozenset(), catalog_sha256=profile.decision['catalog_sha256'])
            if not selected.retained:
                continue
            # Source disease identifiers remain native, without an inferred
            # ontology equivalence or a pathogenic clinical classification.
            document = asdict(fact)
            digest = _digest(document)
            if digest in seen:
                continue
            seen.add(digest)
            source = SourcedRecord(fact.source.dataset, fact.extraction_policy, fact.source.record_id,
                digest, fact.source.release_sha256, {'curated_functional_allele': document})
            for condition in fact.disease_ids:
                for subject in subjects:
                    result.append(admission._Association(subject, admission.inv.RULE_IDS[4], 'curated',
                        'match', 'Ensembl' if fact.gene_id else None, fact.gene_id, condition,
                        ('exact-curated-disease-relevant-damaging-functional-allele',), (source,)))
    return tuple(result) or tuple(admission._Association(subject, admission.inv.RULE_IDS[4], 'curated',
        'no_match', None, None, None, ('complete-curated-lookup-without-qualifying-exact-allele',))
        for subject in subjects)


def scope_associations(profile, associations):
    """Operate on one authenticated input; never hide missing adapter coverage."""
    from genosapiens import candidate_admission as admission
    from genosapiens.wgs_profile_projection import ProfileDecisionRow, profile_association_applicability
    if type(profile) is not ProfileDecisionRow or type(associations) is not tuple:
        raise ValueError('invalid focused association scope input')
    routes = tuple(row for row in profile.decision['routes'] if row['outcome']=='selected')
    selected = replace(profile, decision={**profile.decision, 'routes':routes})
    applicable = profile_association_applicability(selected, associations)
    native = any(row['route']=='native' for row in routes)
    molecular = {(row['target_namespace'], row['target_id'], locator)
        for row in routes if row['route']=='molecular' for locator in row['source_refs']}
    survivors, excluded = {}, {}
    retained = omitted = regulatory_total = regulatory_selected = 0
    for value, allowed in zip(associations, applicable, strict=True):
        reason = None
        if value.state != 'match':
            reason = 'source-'+value.state
        elif not allowed:
            reason = 'outside-selected-target-route'
        elif value.adapter == 'annotation' and not native:
            target = value.target.split('.')[0] if value.namespace=='Ensembl' and value.target else value.target
            if not any((value.namespace,target,source_locator(source)) in molecular for source in value.source_refs):
                reason = 'outside-selected-transcript-proof'
        if reason is None and value.adapter == 'encode':
            sources, counts = select_regulatory_sources(value.source_refs)
            regulatory_total += counts['total']; regulatory_selected += counts['selected']
            if counts['selected']==0:
                reason = 'regulatory-finite-score-unavailable'
            else:
                value = replace(value, source_refs=sources)
        if reason is not None:
            key = (value.subject,value.rule,value.adapter)
            entry = excluded.setdefault(key, Counter())
            entry[reason] += 1
            omitted += 1
            continue
        retained += 1
        # Geometry and regulatory interpretation must survive coalescing.
        interpretation = tuple(sorted(admission._record_bytes(source) for source in value.source_refs
            if source.source_dataset == 'genosapiens.admission-geometry' or value.adapter == 'encode'))
        key = (admission._bytes(replace(value,source_refs=(),evidence_refs=())), interpretation)
        if key not in survivors:
            survivors[key] = [value,{},{}]
        group = survivors[key]
        for source in value.source_refs:
            group[1][admission._record_bytes(source)] = source
        for reference in value.evidence_refs:
            group[2][admission._record_bytes(reference)] = reference
    result = []
    for key in sorted(survivors):
        value, sources, references = survivors[key]
        result.append(replace(value,source_refs=tuple(sources[key] for key in sorted(sources)),
            evidence_refs=tuple(references[key] for key in sorted(references))))
    for (subject,rule,adapter), reasons in sorted(excluded.items()):
        summary = _summary('source-scope-exclusion', dict(subject=subject, rule=rule, adapter=adapter,
            original_input_ref=profile.original_ref, screening_sha256=profile.screening_sha256,
            profile_policy_sha256=profile.decision['policy_sha256'],
            excluded_associations=sum(reasons.values()), reasons=dict(sorted(reasons.items()))))
        result.append(admission._Association(subject,rule,adapter,'unresolved',None,None,None,
            ('focused-association-scope-excluded','outside-competitive-scope-not-benign'),(summary,)))
    if {(a.subject,a.rule,a.adapter) for a in result}!={(a.subject,a.rule,a.adapter) for a in associations}:
        raise ValueError('focused association required coverage differs')
    return ScopedAssociations(tuple(result),len(associations),retained,omitted,
        len(survivors),len(excluded),regulatory_total,regulatory_selected)
