"""Strict codecs for completed scientific projections, without scientific replay.

These helpers confer no producer authority. Their caller must own a compatible
completed-science capability and the live authenticated evidence reader.
"""
from genosapiens.wgs_evidence_store import EvidenceStoreError


def require(value, message="invalid completed scientific projection"):
    if not value:
        raise EvidenceStoreError(message)


def integer(value, *, minimum=0, maximum=2**63-1):
    require(type(value) is int and minimum <= value <= maximum)
    return value


def fields(value, names):
    require(type(value) is dict and set(value) == set(names))
    return value


def child(reader, value, kind):
    fields(value, ("name", "sha256"))
    reference = reader.reference(value["name"], kind=kind)
    require(reference.sha256 == value["sha256"], "completed scientific child binding changed")
    return reference


def header(reader, reference, *, kind, schema, keys, version=1):
    require(reference.kind == kind)
    value = fields(reader.reconstruct(reference), ("schema", "version", *keys))
    require(value["schema"] == schema and type(value["version"]) is int and value["version"] == version)
    return value


def observed_rows(rows, progress, subphase, total):
    from genosapiens.wgs_runtime import analysis_progress
    analysis_progress(progress, 'publication', subphase, unit='records', completed=0, total=total)
    count = 0
    for count, row in enumerate(rows, 1):
        yield row
        if count % 64 == 0:
            analysis_progress(progress, 'publication', subphase, unit='records', completed=count, total=total)
    analysis_progress(progress, 'publication', subphase, unit='records', completed=count, total=total)


def scalar_value(kind, value):
    """Validate the saved scalar wire shape without recomputing decisions."""
    from genosapiens.evidence_model import _exact
    from genosapiens.population_policy import PopulationOutcome
    if kind == 'members':
        require(type(value) is list and len(value) == 3 and type(value[2]) is bool)
        if value[0] is not None:
            _exact(value[0])
        PopulationOutcome(value[1])
        return
    keys = ('model_state', 'completeness', 'member_states', 'clinical_support', 'alternate_present',
        'conditional_singleton_search_satisfied', 'read_validation_required')
    extra = ('constructible',) if kind == 'classes' else ('population_penalty', 'population_contribution',
        'population_outcome', 'unknown_rescue_members', 'population_excluded_members')
    fields(value, (*keys, *extra))
    require(all(type(value[k]) is str and value[k] for k in ('model_state', 'completeness')))
    require(all(type(value[k]) is bool for k in ('conditional_singleton_search_satisfied', 'read_validation_required')))
    width = 2 if kind == 'classes' else 1
    for key in ('member_states', 'clinical_support', 'alternate_present'):
        require(type(value[key]) is list and len(value[key]) == width)
    require(all(type(v) is str and v for v in value['member_states']))
    require(all(type(v) is bool for v in value['clinical_support']))
    # Saved genotype properties retain unknown presence; do not turn it into
    # false or reject a scientifically valid unresolved genotype.
    require(all(v is None or type(v) is bool for v in value['alternate_present']))
    if kind == 'classes':
        require(type(value['constructible']) is bool)
    else:
        require((value['population_penalty'] is None) == (value['population_contribution'] is None))
        if value['population_penalty'] is not None:
            _exact(value['population_penalty'])
            _exact(value['population_contribution'])
        PopulationOutcome(value['population_outcome'])
        integer(value['unknown_rescue_members'], maximum=width)
        require(type(value['population_excluded_members']) is list and len(value['population_excluded_members']) == width
            and all(type(v) is bool for v in value['population_excluded_members']))


def outcome_chunks(reader, tape_name, tape_value, *, catalog_sha256, catalog_count,
                   pair_count, domain_sha256=None, mapping=None):
    """Authenticate and decode a saved compact tape once, optionally remap codes."""
    from array import array
    import hashlib
    import sys
    import zlib
    from genosapiens.wgs_pair_outcomes import PairOutcomeTape, _decode, MAX_COMPRESSED_BYTES
    from genosapiens.evidence_serialization import _canonical
    tape = PairOutcomeTape(**tape_value)
    reference = reader.reference(tape_name, kind='pair-outcome-tape-v1')
    require(_canonical(reader.reconstruct(reference)) == _canonical(tape_value))
    require((tape.catalog_sha256, tape.catalog_count, tape.pair_count) ==
        (catalog_sha256, catalog_count, pair_count))
    require(domain_sha256 is None or tape.domain_sha256 == domain_sha256)
    root = child(reader, dict(name=tape.chunks_name, sha256=tape.chunks_sha256), 'pair-outcome-chunks-v1')
    counts = array('Q', [0]) * integer(catalog_count, maximum=1_000_000)
    count, digest = 0, hashlib.sha256()
    for ordinal, chunk in enumerate(reader.iter_array(root, max_item_bytes=2*MAX_COMPRESSED_BYTES), 1):
        raw, values = _decode(chunk, first=count+1, remaining=pair_count-count, catalog_count=catalog_count)
        digest.update(raw)
        count += len(values)
        for value in values:
            counts[value] += 1
        if mapping is not None:
            require(len(mapping) == catalog_count)
            mapped = array('I', (mapping[value] for value in values))
            if sys.byteorder != 'little':
                mapped.byteswap()
            raw = mapped.tobytes()
        yield ordinal, len(values), zlib.compress(raw, 6)
    require(count == pair_count and digest.hexdigest() == tape.stream_sha256 and all(counts))


def publish_selected_index(writer, context):
    """Stream only missing lookup ordinals from the live computed owners.

    Singleton construction and score arrays already share hypothesis-id order;
    pair outcomes/classes already have exact ordinal tapes. Neither is copied.
    Member enrichment/representation use source order, and population uses
    context/id order, so only their compact lookup keys need an extra index.
    """
    name = context.name + '.selected_index'
    members = context.members
    pair_table = None if context.universe.pairs is None else context.universe.pairs.prefix + '_members'
    source_order = f'(SELECT id,row_number() OVER (ORDER BY ordinal) AS ordinal FROM {members.prefix}_members)'
    member_sql = f'SELECT id,ordinal FROM {source_order} ORDER BY id'
    if pair_table is not None:
        member_sql = (f'SELECT m.id,m.ordinal,p.position FROM {source_order} m LEFT JOIN '
            f'(SELECT id,row_number() OVER (ORDER BY id)-1 AS position FROM {pair_table}) p '
            'ON p.id=m.id ORDER BY m.id')
    def member_rows():
        count = 0
        for position, row in enumerate(writer.connection.execute(member_sql)):
            key, ordinal = row[:2]
            require(type(ordinal) is int and 1 <= ordinal <= members.member_count)
            query, = writer.connection.execute(
                f'SELECT query_ordinal FROM {context.representations.table} WHERE id=?', (key,)).fetchone()
            yield [key, ordinal - 1, ordinal - 1, position,
                   None if pair_table is None else row[2], query]
            count += 1
        require(count == members.member_count)
    def population_rows():
        count = 0
        for position, (key, allele) in enumerate(writer.connection.execute(
                f'SELECT context,id FROM {context.population.prefix}_members ORDER BY context,id')):
            yield [key, allele, position]
            count += 1
        require(count == context.population.assessment_count)
    arrays = {}
    for key, rows in (('members', member_rows()), ('population', population_rows())):
        kind = 'selected-reconstruction-' + key + '-v1'
        digest = writer.add_sequence(name + '.' + key, kind, rows)
        arrays[key] = dict(name=name + '.' + key, kind=kind, sha256=digest)
    value = dict(schema='genosapiens.selected-reconstruction-index', version=1,
        context_ordinal=context.ordinal, context_id=context.context_id,
        resources=context.resources.reference, assessment_domain_sha256=context.assessment_domain_sha256,
        dependencies={key: context.header[key] for key in (
            'members', 'facts', 'universe', 'population', 'phenotype', 'scores',
            'representations', 'eligibility', 'order')},
        member_count=members.member_count, population_count=context.population.assessment_count,
        singleton_count=context.universe.singleton_count, pair_count=context.pair_count,
        hypothesis_count=context.hypothesis_count, arrays=arrays)
    digest = writer.add(name, 'selected-reconstruction-index-v1', value)
    return dict(name=name, kind='selected-reconstruction-index-v1', sha256=digest)


class SelectedReconstructionIndex:
    """Exact-owner bounded seeks; no rich context or scientific restoration."""
    def __init__(self, reader, context, resources):
        from genosapiens.wgs_case_resources import bound_reference
        require(type(context.get('version')) is int and context['version'] == 2)
        root = bound_reference(reader, context['selected_reconstruction'], 'selected-reconstruction-index-v1')
        value = header(reader, root, kind='selected-reconstruction-index-v1',
            schema='genosapiens.selected-reconstruction-index', keys=(
                'context_ordinal', 'context_id', 'resources', 'assessment_domain_sha256',
                'dependencies', 'member_count', 'population_count', 'singleton_count',
                'pair_count', 'hypothesis_count', 'arrays'))
        require(value['resources'] == context['resources'] == resources.reference
            and value['assessment_domain_sha256'] == context['assessment_domain_sha256'] == resources.assessment_domain_sha256
            and value['context_ordinal'] == context['ordinal'] and value['context_id'] == context['context_id'])
        fields(value['dependencies'], ('members', 'facts', 'universe', 'population', 'phenotype',
            'scores', 'representations', 'eligibility', 'order'))
        require(all(value['dependencies'][key] == context[key] for key in value['dependencies']))
        for key in ('member_count', 'pair_count', 'hypothesis_count'):
            require(integer(value[key]) == context[key])
        integer(value['member_count'], maximum=100_000)
        integer(value['population_count'])
        integer(value['singleton_count'])
        fields(value['arrays'], ('members', 'population'))
        self.reader, self.value = reader, value
        self.roots = {key: bound_reference(reader, ref, 'selected-reconstruction-' + key + '-v1')
                      for key, ref in value['arrays'].items()}
        for key, count in (('members', value['member_count']), ('population', value['population_count'])):
            metadata = reader._lookup(self.roots[key].sha256)
            require(metadata[0] == 'array' and metadata[2] == count)

    def _find(self, kind, key):
        self.reader._check_files()
        lower, upper = 0, self.value['member_count' if kind == 'members' else 'population_count']
        width = 1 if kind == 'members' else 2
        while lower < upper:
            middle = (lower + upper) // 2
            row = self.reader.array_item(self.roots[kind], middle, max_item_bytes=1024**2)
            require(type(row) is list and len(row) == (6 if kind == 'members' else 3))
            require(all(type(item) is str and item for item in row[:width]))
            found = tuple(row[:width])
            if found < key:
                lower = middle + 1
            elif found > key:
                upper = middle
            else:
                return row
        raise EvidenceStoreError('selected reconstruction locator is missing')

    def member(self, allele_id):
        row = self._find('members', (allele_id,))
        for position in row[1:4]:
            integer(position, maximum=self.value['member_count'] - 1)
        require(row[1] == row[2])
        if row[4] is not None:
            integer(row[4], maximum=self.value['member_count'] - 1)
        if row[5] is not None:
            integer(row[5], maximum=self.value['member_count'] - 1)
        return row

    def population(self, context_key, allele_id):
        from genosapiens.evidence_serialization import _canonical
        row = self._find('population', (_canonical(context_key), allele_id))
        integer(row[2], maximum=self.value['population_count'] - 1)
        return row[2]


def _saved_pair_code(reader, value, position, *, domain, catalog, count):
    from genosapiens.wgs_pair_outcomes import PairOutcomeTape, CHUNK_PAIRS, _decode
    tape = PairOutcomeTape(**value['tape'])
    root = reader.reference(value['tape_name'], kind='pair-outcome-tape-v1')
    require(reader.reconstruct(root) == value['tape'])
    require(tape.domain_sha256 == domain['membership_sha256']
        and tape.pair_count == domain['pair_count'] and tape.catalog_sha256 == catalog['sha256']
        and tape.catalog_count == count)
    integer(position, minimum=1, maximum=tape.pair_count)
    chunks = reader.reference(tape.chunks_name, kind='pair-outcome-chunks-v1')
    require(chunks.sha256 == tape.chunks_sha256)
    offset = ((position - 1) // CHUNK_PAIRS) * CHUNK_PAIRS
    chunk = reader.array_item(chunks, offset // CHUNK_PAIRS, max_item_bytes=1024**2)
    _, values = _decode(chunk, first=offset + 1, remaining=tape.pair_count - offset,
        catalog_count=tape.catalog_count)
    return values[(position - 1) % CHUNK_PAIRS]


def selected_records(reader, case, resources, ordinal, ordered_rows, *, reference_names):
    """Read selected singleton/pair detail using saved, exact producer facts.

    Initial sealing has no saved groups yet. It derives them from these checked
    representations; later consumers additionally check the sealed groups.
    """
    from genosapiens.evidence_serialization import _canonical, _typed
    from genosapiens.candidate_admission import _hash
    from genosapiens.evidence_model import (
        ScoreLedger, HypothesisContext, EvidenceReference, EnumerationStatus,
    )
    from genosapiens.evidence_fusion import (
        HypothesisEvidence, EvidenceCharge, _score, _inheritance_charges, population_assessment_charges,
    )
    from genosapiens.genotype_hypotheses import ConstructedHypothesis, construct_evaluated_hypothesis, _record
    from genosapiens.genotype_model import ModelMemberFacts
    from genosapiens.partner_qualification import PartnerQualificationAssessment
    from genosapiens.population_policy import ContextKey, MemberAssessment, summarize_population_members
    from genosapiens.phenotype_scoring import PhenotypeScore, _freeze
    from genosapiens.wgs_context_population import PopulationAssessmentValues
    from genosapiens.wgs_admission import LedgerAdmissionInput
    from genosapiens.wgs_case_resources import bound_reference
    from genosapiens.wgs_genotype_facts import GenotypeFactHeader
    from genosapiens.wgs_pair_decisions import PairDecision
    from genosapiens.wgs_pair_evidence import assemble_evaluated_pair
    from genosapiens.wgs_phase_patterns import PairPhasePattern
    from genosapiens.wgs_score_decisions import ScoreDecision
    from genosapiens.track1_submission import SubmissionRepresentationAssessment, _allele_trace_key
    if not ordered_rows:
        return
    require(len(ordered_rows) <= 10 and len({row['rank'] for row in ordered_rows}) == len(ordered_rows)
        and all(row['locator'][0] == ordinal for row in ordered_rows))
    directory = case.context(ordinal)
    context = reader.reconstruct(bound_reference(reader, directory['context'], 'scientific-case-context-v1'))
    require(all(context[key] == directory[key] for key in (
        'ordinal', 'context_id', 'member_count', 'pair_count', 'hypothesis_count')))
    index = SelectedReconstructionIndex(reader, context, resources)
    def component(key, kind):
        return reader.reconstruct(bound_reference(reader, context[key], kind))
    def ref(value):
        return {key: value[key] for key in ('name', 'sha256')}
    def item(value, kind, position, maximum=128*1024**2):
        return reader.array_item(child(reader, value, kind), position, max_item_bytes=maximum)
    universe = component('universe', 'complete-hypothesis-universe-v1')
    scores = component('scores', 'complete-context-scores-v1')
    population = component('population', 'context-population-evidence-v1')
    members = component('members', 'context-evidence-v1')
    representations = component('representations', 'complete-submission-representations-v1')
    require(scores['hypotheses'] == population['hypotheses'] == ref(context['universe'])
        and scores['population'] == ref(context['population'])
        and scores['policy_sha256'] == resources.fusion_policy_sha256
        and _typed(scores['policy'], type(resources.fusion_policy), decode=True) == resources.fusion_policy
        and population['policy'] == resources.population_reference
        and population['context_evidence'] == representations['members'] == ref(context['members'])
        and representations['genotype_facts'] == context['facts']
        and representations['reference'] == reader.reconstruct(
            bound_reference(reader, context['resources'], 'scientific-case-resources-v1'))['reference']
        and population['assessment_count'] == index.value['population_count']
        and scores['singleton_count'] == universe['singleton_count'] == index.value['singleton_count']
        and members['member_count'] == representations['count'] == index.value['member_count'])
    selected_groups = (None if not hasattr(case, 'selected_groups') else
        {group.source_rank: group for group in case.selected_groups})
    for ordered in ordered_rows:
        rank = integer(ordered['rank'], minimum=1, maximum=case.hypothesis_count)
        saved_locator, saved_total, saved_ties = case.orders['primary'].at_rank(rank)
        require(ordered['locator'] == saved_locator and ordered['total'] == saved_total
            and ordered['ties'] == saved_ties and ordered['mode'] == 'primary'
            and ordered['hypothesis_id'] == case.identities.identity(saved_locator))
        _, left, right = ordered['locator']
        integer(left, minimum=1)
        integer(right)
        pair_position = None
        if right == 0:
            integer(left, maximum=universe['singleton_count'])
            constructed = _typed(item(universe['singletons'], 'constructed-singletons-v1', left-1),
                ConstructedHypothesis, decode=True)
            saved_score = item(scores['singletons'], 'singleton-score-evidence-v1', left-1)
            fields(saved_score, ('hypothesis_id', 'evidence', 'ledger'))
            require(saved_score['hypothesis_id'] == constructed.hypothesis.hypothesis_id)
            evidence = _typed(saved_score['evidence'], HypothesisEvidence, decode=True)
            ledger = _typed(saved_score['ledger'], ScoreLedger, decode=True)
        else:
            positions = case.inputs._context(ordinal).positions
            require(positions is not None)
            pair_position = positions.ordinal(left, right)
            pairs = reader.reconstruct(child(reader, universe['pair_evidence'], 'complete-pair-evidence-v1'))
            domain = pairs['domain_summary']
            require(domain['member_count'] == directory['member_count']
                and domain['pair_count'] == pairs['pair_count'] == directory['pair_count']
                and pairs['genotype_facts'] == context['facts'])
            keys, templates = [], []
            for position in (left, right):
                row = item(pairs['domain'], 'pair-domain-members-v1', position-1, 1024**2)
                require(type(row) is list and len(row) == 5 and row[0] == position and row[4] is False)
                key = row[1]
                locator = index.member(key)
                require(locator[3] == position-1 and locator[4] is not None)
                template = item(pairs['members'], 'pair-member-evidence-v1', locator[4])
                fields(template, ('allele_id', 'decision_class', 'model_member', 'qualification'))
                require(template['allele_id'] == key and template['model_member']['allele']['identity_id'] == row[2])
                keys.append(key)
                templates.append(template)
            code = _saved_pair_code(reader, pairs, pair_position, domain=domain,
                catalog=pairs['outcomes'], count=pairs['outcome_count'])
            outcome = item(pairs['outcomes'], 'pair-computed-outcomes-v1', code, 1024**2)
            fields(outcome, ('ordinal', 'decision', 'pattern'))
            require(outcome['ordinal'] == code)
            catalog = reader.reconstruct(child(reader, pairs['decisions'], 'pair-decision-catalog-v1'))
            require(catalog['context'] == pairs['context'] and catalog['outcome_count'] == pairs['decision_count'])
            decision_row = item(dict(name=catalog['outcomes_name'], sha256=catalog['outcomes_sha256']),
                'pair-decision-results-v1', integer(outcome['decision'], maximum=catalog['outcome_count']-1), 1024**2)
            require(decision_row['ordinal'] == outcome['decision']
                and (decision_row['left'], decision_row['right']) == tuple(row['decision_class'] for row in templates))
            decision = _typed(decision_row['result'], PairDecision, decode=True)
            pattern_row = item(pairs['patterns'], 'pair-phase-patterns-v1',
                integer(outcome['pattern'], maximum=pairs['pattern_count']-1), 16*1024**2)
            import hashlib
            require(pattern_row['ordinal'] == outcome['pattern'] and pattern_row['sha256'] ==
                hashlib.sha256(_canonical(pattern_row['pattern']).encode('ascii')).hexdigest())
            pattern = _typed(pattern_row['pattern'], PairPhasePattern, decode=True)
            require(pattern.phase.value == decision_row['phase'] and decision.constructible)
            facts = component('facts', 'genotype-facts-v1')
            genotype = _typed(facts['header'], GenotypeFactHeader, decode=True)
            require(facts['input_evidence'] == ref(context['members'])
                and genotype.input_sha256 == resources.source_manifest_sha256
                and genotype.reference_manifest_sha256 == resources.reference_binding.manifest_sha256
                and genotype.sample_name == resources.sample_name)
            pair_context = _typed(pairs['context'], HypothesisContext, decode=True)
            _, _, _, model, qualified = assemble_evaluated_pair(pair_position, tuple(keys), decision, pattern,
                tuple(_typed(row['model_member'], ModelMemberFacts, decode=True) for row in templates),
                tuple(_typed(row['qualification'], PartnerQualificationAssessment, decode=True) for row in templates),
                context=pair_context, header=genotype, facts_reference=context['facts'],
                pair_reference=universe['pair_evidence'],
                locus_policy_sha256=_hash(tuple(sorted(resources.hypothesis_policy.loci, key=_canonical))))
            search = reader.reconstruct(child(reader, universe['search'], 'complete-partner-search-v1'))
            require(search['context_evidence'] == ref(context['members'])
                and search['policy_sha256'] == resources.hypothesis_policy_sha256)
            constructed = construct_evaluated_hypothesis(context=pair_context, admission_context_id=directory['context_id'],
                references=tuple(m.allele for m in model.members), properties=tuple(m.genotype for m in model.members),
                assessments=qualified, phase=model.phase, phase_evidence=model.phase_evidence, model=model,
                enumeration=EnumerationStatus.COMPLETE,
                enumeration_evidence=(_record('partner-search', dict(context_id=directory['context_id'], search_evidence=universe['search'])),),
                sample_name=universe['sample_name'], observation_ids=tuple(universe['observation_ids']),
                incompleteness_penalty=universe['incompleteness_penalty'],
                outside_context_partner_uncertainty=search['outside_context_partner_uncertainty'])
            require(constructed is not None)
        hypothesis = constructed.hypothesis
        require(hypothesis.hypothesis_id == ordered['hypothesis_id']
            and constructed.admission_context_id == directory['context_id'])
        key = ContextKey.from_context(hypothesis.context)
        assessed, selected_values = [], set()
        alleles, locators, traces = {}, {}, {}
        for subject in sorted(hypothesis.evidence_subject_refs, key=lambda v: v.allele_id):
            location = index.member(subject.allele_id)
            locators[subject.allele_id] = location
            saved = item(members['members'], 'context-enriched-members-v1', location[1])
            fields(saved, ('admitted', 'enriched', 'batch'))
            member = _typed(saved['enriched'], LedgerAdmissionInput, decode=True)
            require(member.allele.allele_id == subject.allele_id and member.allele.native_identity_id == subject.identity_id)
            alleles[subject.allele_id] = member.allele
            row = item(population['members'], 'context-population-members-v1', index.population(key, subject.allele_id))
            fields(row, ('context', 'selected', 'assessment'))
            require(_canonical(row['context']) == _canonical(key) and type(row['selected']) is bool)
            assessment = _typed(row['assessment'], MemberAssessment, decode=True)
            require(assessment.allele_id == subject.allele_id and assessment.native_identity_id == subject.identity_id)
            assessed.append(assessment)
            selected_values.add(row['selected'])
        require(len(selected_values) == 1)
        selected, = selected_values
        penalty, contribution, outcome = summarize_population_members(tuple(assessed))
        population_assessment = PopulationAssessmentValues(key if selected else None, tuple(assessed), penalty, contribution, outcome,
            'population_only;not_all_checks_or_scientific_acceptance' if selected else 'no_exact_context_policy')
        for subject in hypothesis.allele_refs:
            location = locators[subject.allele_id]
            result = _typed(item(representations['representations'], 'member-submission-representations-v1', location[2]),
                SubmissionRepresentationAssessment, decode=True)
            require(result.allele_id == subject.allele_id and result.native_identity_id == subject.identity_id
                and result.status == 'eligible' and result.trace is not None and location[5] is not None)
            query = representations['queries']
            integer(location[5], maximum=query['count']-1)
            saved_query = item(ref(query), 'submission-reference-queries-v1', location[5], 1024**2)
            require(saved_query['contig'] == reference_names.get(result.trace.rendered.contig.removeprefix('chr')))
            require(saved_query['bases'] == result.trace.rendered.reference,
                'selected reference observation does not replay')
            require(saved_query['position'] == result.trace.rendered.position and saved_query['length'] == len(result.trace.rendered.reference))
            traces[subject.allele_id] = result.trace
        if right:
            pair_scores = reader.reconstruct(child(reader, scores['pairs'], 'complete-pair-scores-v1'))
            require(pair_scores['hypotheses'] == ref(context['universe']) and pair_scores['population'] == ref(context['population']))
            sources = reader.reconstruct(child(reader, pair_scores['sources'], 'pair-score-source-evidence-v1'))
            require(sources['hypotheses'] == ref(context['universe']) and sources['policy_sha256'] == resources.fusion_policy_sha256)
            shared = reader.reconstruct(child(reader, sources['shared'], 'pair-score-shared-sources-v1'))
            charges = list(_typed(shared['charges'], tuple[EvidenceCharge, ...], decode=True))
            for member_id in hypothesis.evidence_subject_ids:
                row = item(sources['members'], 'pair-score-member-sources-v1', locators[member_id][4])
                require(row['allele_id'] == member_id)
                charges.extend(_typed(row['charges'], tuple[EvidenceCharge, ...], decode=True))
            phenotype = component('phenotype', 'context-phenotype-evidence-v1')
            score = PhenotypeScore(_freeze(phenotype['score']))
            refs = _typed(phenotype['profile_refs'], tuple[EvidenceReference, ...], decode=True)
            evidence = HypothesisEvidence(hypothesis, score, (*charges, *_inheritance_charges(constructed, resources.fusion_policy),
                *population_assessment_charges(population_assessment, alleles, resources.population_reference['policy_sha256'])), refs, None)
            ledger = _score(evidence, resources.fusion_policy, resources.fusion_policy_sha256)
            catalog = reader.reconstruct(child(reader, pair_scores['decisions'], 'score-decision-catalog-v1'))
            require(catalog['policy_sha256'] == resources.fusion_policy_sha256)
            score_code = _saved_pair_code(reader, pair_scores, pair_position, domain=domain,
                catalog=pair_scores['outcomes'], count=catalog['outcome_count'] + int(pair_scores['no_hypothesis_code'] is not None))
            row = item(pair_scores['outcomes'], 'pair-score-outcomes-v1', score_code, 1024**2)
            require(row['ordinal'] == score_code and type(row['score']) is int)
            saved = item(catalog['outcomes'], 'score-computed-outcomes-v1', row['score'], 1024**2)
            require(saved['ordinal'] == row['score'])
            decision = _typed(saved['result'], ScoreDecision, decode=True)
            require(ledger.total == decision.total and tuple((c.name, c.raw, c.normalized) for c in ledger.components)
                == tuple((c.name, c.raw, c.normalized) for c in decision.components))
            eligibility_code = case.inputs._pair_code(ordinal, pair_position-1) - case.inputs._context(ordinal).offset
        else:
            row = reader.connection.execute(f'SELECT code FROM {case.inputs.prefix}_singletons WHERE context=? AND ordinal=?',
                (ordinal, left)).fetchone()
            require(row is not None)
            eligibility_code = row[0]
        reasons, _, _ = case.inputs._class(ordinal, eligibility_code)
        rendered = tuple(sorted(traces.values(), key=_allele_trace_key))
        require(not reasons and evidence.hypothesis == hypothesis and ledger.hypothesis_id == hypothesis.hypothesis_id
            and ledger.total == ordered['total'])
        if selected_groups is not None:
            require(ordered['rank'] in selected_groups and rendered == selected_groups[ordered['rank']].alleles)
        yield dict(context_ordinal=ordinal, context_id=directory['context_id'], constructed=constructed,
            evidence=evidence, ledger=ledger, population_assessment=population_assessment,
            population_policy=resources.population_reference, eligibility_reasons=reasons,
            alleles=tuple(alleles[subject.allele_id] for subject in hypothesis.allele_refs),
            allele_traces=rendered, ordering=ordered)


def selected_singletons(reader, case, resources, ordinal, ordered_rows):
    """Decode saved singleton detail; only aggregate already assessed members.

    No resource corpus, hypothesis constructor, phenotype scorer or score fusion
    runs here. Pair-rich detail remains the explicit selected-context fallback.
    """
    from genosapiens.evidence_serialization import _canonical, _typed
    from genosapiens.genotype_hypotheses import ConstructedHypothesis
    from genosapiens.evidence_fusion import HypothesisEvidence
    from genosapiens.evidence_model import ScoreLedger
    from genosapiens.population_policy import ContextKey, MemberAssessment, summarize_population_members
    from genosapiens.wgs_context_population import PopulationAssessmentValues
    from genosapiens.wgs_admission import LedgerAdmissionInput
    from genosapiens.wgs_case_resources import bound_reference
    from genosapiens.track1_submission import SubmissionRepresentationAssessment, _allele_trace_key
    require(all(row['locator'][0] == ordinal and row['locator'][2] == 0 for row in ordered_rows))
    wanted = {row['locator'][1]: row for row in ordered_rows}
    require(len(wanted) == len(ordered_rows) <= 10)
    directory = case.context(ordinal)
    context = reader.reconstruct(bound_reference(reader, directory['context'], 'scientific-case-context-v1'))
    universe = reader.reconstruct(bound_reference(reader, context['universe'], 'complete-hypothesis-universe-v1'))
    scores = reader.reconstruct(bound_reference(reader, context['scores'], 'complete-context-scores-v1'))
    population = reader.reconstruct(bound_reference(reader, context['population'], 'context-population-evidence-v1'))
    def ref(value):
        return {k: value[k] for k in ('name', 'sha256')}
    require(scores['hypotheses'] == population['hypotheses'] == ref(context['universe'])
        and scores['population'] == ref(context['population']) and scores['policy_sha256'] == resources.fusion_policy_sha256
        and population['policy'] == resources.population_reference and population['context_evidence'] == ref(context['members']))
    constructed, count = {}, 0
    for count, row in enumerate(reader.iter_array(child(reader, universe['singletons'], 'constructed-singletons-v1'),
            max_item_bytes=128*1024**2), 1):
        if count in wanted:
            item = _typed(row, ConstructedHypothesis, decode=True)
            require(item.hypothesis.hypothesis_id == wanted[count]['hypothesis_id']
                and item.admission_context_id == directory['context_id'])
            constructed[item.hypothesis.hypothesis_id] = item
    require(count == universe['singleton_count'] and len(constructed) == len(wanted))
    score_rows, count = {}, 0
    for row in reader.iter_array(child(reader, scores['singletons'], 'singleton-score-evidence-v1'), max_item_bytes=128*1024**2):
        count += 1
        fields(row, ('hypothesis_id', 'evidence', 'ledger'))
        key = row['hypothesis_id']
        if key in constructed:
            require(key not in score_rows)
            evidence, ledger = _typed(row['evidence'], HypothesisEvidence, decode=True), _typed(row['ledger'], ScoreLedger, decode=True)
            require(evidence.hypothesis == constructed[key].hypothesis)
            score_rows[key] = evidence, ledger
    require(count == scores['singleton_count'] == universe['singleton_count'] and set(score_rows) == set(constructed))
    population_needed = {(_canonical(ContextKey.from_context(item.hypothesis.context)), subject.allele_id): subject.identity_id
        for item in constructed.values() for subject in item.hypothesis.evidence_subject_refs}
    population_rows, count = {}, 0
    for row in reader.iter_array(child(reader, population['members'], 'context-population-members-v1'), max_item_bytes=64*1024**2):
        count += 1
        fields(row, ('context', 'selected', 'assessment'))
        require(type(row['selected']) is bool and type(row['assessment']) is dict)
        key = _canonical(row['context']), row['assessment'].get('allele_id')
        if key in population_needed:
            require(key not in population_rows)
            assessment = _typed(row['assessment'], MemberAssessment, decode=True)
            require(assessment.native_identity_id == population_needed[key])
            population_rows[key] = row['selected'], assessment
    require(count == population['assessment_count'] and set(population_rows) == set(population_needed))
    allele_ids = {ref.allele_id for item in constructed.values() for ref in item.hypothesis.allele_refs}
    members = reader.reconstruct(bound_reference(reader, context['members'], 'context-evidence-v1'))
    alleles, count = {}, 0
    for row in reader.iter_array(child(reader, members['members'], 'context-enriched-members-v1'), max_item_bytes=69*1024**2):
        count += 1
        fields(row, ('admitted', 'enriched', 'batch'))
        key = row['enriched']['allele']['allele_id']
        if key in allele_ids:
            require(key not in alleles)
            item = _typed(row['enriched'], LedgerAdmissionInput, decode=True)
            alleles[key] = item.allele
    require(count == directory['member_count'] == members['member_count'] and set(alleles) == allele_ids)
    representations = reader.reconstruct(bound_reference(reader, context['representations'], 'complete-submission-representations-v1'))
    traces, count = {}, 0
    for row in reader.iter_array(child(reader, representations['representations'], 'member-submission-representations-v1'), max_item_bytes=128*1024**2):
        count += 1
        key = row.get('allele_id')
        if key in allele_ids:
            assessment = _typed(row, SubmissionRepresentationAssessment, decode=True)
            require(key not in traces and assessment.status == 'eligible' and assessment.trace is not None
                and assessment.native_identity_id == alleles[key].native_identity_id)
            traces[key] = assessment.trace
    require(count == directory['member_count'] and set(traces) == allele_ids)
    groups = {group.source_rank: group for group in case.selected_groups}
    for ordered in ordered_rows:
        item = constructed[ordered['hypothesis_id']]
        evidence, ledger = score_rows[ordered['hypothesis_id']]
        key = ContextKey.from_context(item.hypothesis.context)
        assessed, selected_values = [], set()
        for subject in sorted(item.hypothesis.evidence_subject_refs, key=lambda value: value.allele_id):
            selected, member = population_rows[_canonical(key), subject.allele_id]
            selected_values.add(selected)
            assessed.append(member)
        require(len(selected_values) == 1)
        selected, = selected_values
        penalty, contribution, outcome = summarize_population_members(tuple(assessed))
        population_assessment = PopulationAssessmentValues(key if selected else None, tuple(assessed), penalty, contribution, outcome,
            'population_only;not_all_checks_or_scientific_acceptance' if selected else 'no_exact_context_policy')
        singleton = reader.connection.execute(f'SELECT code FROM {case.inputs.prefix}_singletons WHERE context=? AND ordinal=?',
            (ordinal, ordered['locator'][1])).fetchone()
        require(singleton is not None)
        reasons, _, _ = case.inputs._class(ordinal, singleton[0])
        rendered = tuple(sorted((traces[subject.allele_id] for subject in item.hypothesis.allele_refs), key=_allele_trace_key))
        require(not reasons and ledger.total == ordered['total'] and rendered == groups[ordered['rank']].alleles)
        yield dict(context_ordinal=ordinal, context_id=directory['context_id'], constructed=item, evidence=evidence, ledger=ledger,
            population_assessment=population_assessment, population_policy=resources.population_reference, eligibility_reasons=reasons,
            alleles=tuple(alleles[subject.allele_id] for subject in item.hypothesis.allele_refs), allele_traces=rendered, ordering=ordered)
