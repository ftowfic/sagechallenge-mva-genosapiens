"""Build a compact scientific roster before merging rich context seeds."""
from dataclasses import asdict
from collections import OrderedDict
from contextlib import closing
import hashlib
import json
import time

from genosapiens.track1_science_protection import iter_profile_protection
from genosapiens.track1_science_roster import ContextProjection, ScienceRoster, strict_association_category
from genosapiens.track1_science_scope import ScienceScopePolicy
from genosapiens.wgs_ledger import canonical


def _equivalence_index(evidence):
    """Index exactly the distance-zero portion of the ordinary disease graph.

    Only equivalent relations have zero cost, and the ordinary graph inserts
    them in both directions. Selection never uses nonzero distances or path
    proofs. Those remain available to the ordinary scientific kernels.
    """
    from genosapiens.candidate_admission import _relation_edges
    parents, sizes = {}, {}
    def root(value):
        parents.setdefault(value, value)
        while parents[value] != value:
            parents[value] = parents[parents[value]]
            value = parents[value]
        return value
    for left, edges in _relation_edges(evidence).items():
        for right, distance, _ in edges:
            if distance != 0:
                continue
            a, b = root(left), root(right)
            if a != b:
                if sizes.get(a, 1) < sizes.get(b, 1):
                    a, b = b, a
                parents[b] = a
                sizes[a] = sizes.get(a, 1) + sizes.get(b, 1)
    return {value: root(value) for value in parents}


class ScientificContextSelector:
    def __init__(self, *, policy, profile_reader, phenotype_scores, gene_evidence, roster_connection=None):
        from genosapiens.phenotype_scoring import PreparedPhenotypeScores
        from genosapiens.evidence_model import DiseaseGeneEvidence
        if (type(policy) is not ScienceScopePolicy or type(phenotype_scores) is not PreparedPhenotypeScores
                or type(gene_evidence) is not tuple
                or any(type(value) is not DiseaseGeneEvidence for value in gene_evidence)):
            raise ValueError('invalid scientific selector dependencies')
        self.policy, self.profile_reader, self.phenotype_scores = policy, profile_reader, phenotype_scores
        self.gene_evidence = gene_evidence
        self.roster, self.metrics = None, {}
        self.roster_connection = roster_connection
        self._profile_collection = None
        self.assertions = {}
        for public in gene_evidence:
            for assertion in public.store.assertions:
                record = public.store.records[assertion.record_id]
                key = (public.evidence_id, assertion.record_id, public.public_manifest_sha256)
                value = (record.source_fields.get('source_dataset', record.table),
                    assertion.validity, assertion.evidence_kind, assertion.polarity)
                previous = self.assertions.setdefault(key, value)
                if previous != value:
                    raise ValueError('scientific association assertion conflicts')

    def collect_profile_rows(self):
        """Produce protection during the ordinary association reader's replay.

        The private index is usable only after this complete live stream closes.
        It is never a persisted checkpoint or a replacement for policy replay.
        """
        from genosapiens.track1_science_protection import iter_validated_profile_facts, extract_protection
        if self._profile_collection is not None or self.roster_connection is None or self.profile_reader is None:
            raise ValueError('scientific profile collection owner is invalid')
        state = dict(reader=self.profile_reader, connection=self.roster_connection,
            complete=False, metrics={'profile_rows': 0}, seconds=0, count=0)
        digest = hashlib.sha256()
        self._profile_collection = state
        state['connection'].execute('CREATE TABLE science_profile_protection('
            'retained_ordinal INTEGER PRIMARY KEY,original_ref BLOB NOT NULL,content BLOB NOT NULL)')
        with closing(iter_validated_profile_facts(self.profile_reader, metrics=state['metrics'])) as rows:
            while True:
                started = time.monotonic()
                try:
                    row, facts = next(rows)
                except StopIteration:
                    state['seconds'] += time.monotonic() - started
                    break
                selected = row.decision['disposition'] in {'selected', 'unresolved'}
                if selected and row.original_ref[3] is not None:
                    reference = canonical(row.original_ref)
                    content = canonical(asdict(extract_protection(facts)))
                    state['connection'].execute('INSERT INTO science_profile_protection VALUES (?,?,?)',
                        (row.original_ref[1], reference, content))
                    digest.update(reference + b'\n' + content + b'\n')
                    state['count'] += 1
                state['seconds'] += time.monotonic() - started
                if selected:
                    yield row
        self.profile_reader.require_stable()
        state['sha256'] = digest.hexdigest()
        state['complete'] = True

    def _support(self, seed):
        states, supported, uncertain = set(), False, False
        for entry in seed.support_ledger:
            ref = entry.source_ref
            source = self.assertions.get((ref.evidence_id, ref.record_id, ref.public_manifest_sha256))
            if source is None:
                raise ValueError('scientific association source is missing')
            dataset, validity, kind, polarity = source
            supported |= entry.positive and entry.requirement_state == 'exact_alternative'
            uncertain |= bool(entry.conflict_refs)
            states.add(strict_association_category(dataset, validity,
                direct=entry.scope == 'direct', positive=entry.positive and polarity == 'positive',
                causal=kind == 'curated', exact=entry.applicability in {'exact', 'equivalent'}
                    and entry.requirement_state == 'exact_alternative', conflict=bool(entry.conflict_refs)))
        state = next((key for key in ('strong', 'weak', 'nonqualifying') if key in states), 'unknown')
        return state, supported, uncertain

    def build(self, connection, payloads, *, page_rows=256):
        from genosapiens.candidate_admission import _normalise_curie
        from genosapiens.wgs_context_payloads import expand, header
        from genosapiens.wgs_contexts import _context_keys, _load_seed
        if self.roster is not None:
            raise ValueError('scientific selector already consumed its domain')
        started = time.monotonic()
        connection.execute('CREATE TABLE science_member_protection(allele_id TEXT PRIMARY KEY,content BLOB NOT NULL)')
        connection.execute('CREATE TABLE science_disease_scores(disease_id TEXT PRIMARY KEY,score INTEGER,available INTEGER)')
        profile_metrics = {'profile_rows': 0}
        if self.profile_reader is not None:
            def included(ref):
                row = connection.execute('SELECT projection_binding_sha256,retained_ordinal,identity_id,allele_id,input_sha256,input_size '
                    'FROM inputs WHERE retained_ordinal=?', (ref[1],)).fetchone()
                if row is None:
                    return False
                if tuple(row) != ref:
                    raise ValueError('scientific protection and context input differ')
                return ref[3] is not None
            state = self._profile_collection
            if state is None:
                for ref, protection in iter_profile_protection(self.profile_reader, include=included, metrics=profile_metrics):
                    connection.execute('INSERT INTO science_member_protection VALUES (?,?)',
                        (ref[3], canonical(asdict(protection))))
            else:
                if (not state['complete'] or state['reader'] is not self.profile_reader
                        or state['connection'] is not self.roster_connection):
                    raise ValueError('scientific profile collection is incomplete or foreign')
                self.profile_reader.require_stable()
                profile_metrics = dict(state['metrics'])
                digest, count = hashlib.sha256(), 0
                for reference, content in state['connection'].execute(
                        'SELECT original_ref,content FROM science_profile_protection ORDER BY retained_ordinal'):
                    digest.update(reference + b'\n' + content + b'\n')
                    count += 1
                    ref = tuple(json.loads(reference))
                    if included(ref):
                        connection.execute('INSERT INTO science_member_protection VALUES (?,?)', (ref[3], content))
                if count != state['count'] or digest.hexdigest() != state['sha256']:
                    raise ValueError('scientific profile collection content changed')
                self.profile_reader.require_stable()
                state['connection'].execute('DROP TABLE science_profile_protection')
            missing, = connection.execute('SELECT count(*) FROM inputs i WHERE i.allele_id IS NOT NULL '
                'AND NOT EXISTS(SELECT 1 FROM science_member_protection p WHERE p.allele_id=i.allele_id)').fetchone()
            if missing:
                raise ValueError('scientific protection domain is incomplete')
        protection_seconds = time.monotonic() - started
        roster = ScienceRoster(connection if self.roster_connection is None else self.roster_connection, self.policy)
        equivalents = _equivalence_index(self.gene_evidence)
        support_cache = OrderedDict()
        for context_id in _context_keys(connection, page_rows):
            grouping, = connection.execute('SELECT grouping_bytes FROM contexts WHERE context_id=?', (context_id,)).fetchone()
            metadata = json.loads(grouping)
            disease = metadata['disease_id']
            scored = connection.execute('SELECT score,available FROM science_disease_scores WHERE disease_id=?', (disease,)).fetchone()
            if scored is None:
                score = self.phenotype_scores.score(disease)
                scored = (score.ledger['total_units'], int(score.ledger['profile_status'] == 'available'))
                connection.execute('INSERT INTO science_disease_scores VALUES (?,?,?)', (disease, *scored))
            context_sha = hashlib.sha256(grouping).hexdigest()
            native_disease = metadata['native_disease_id']
            disease_group = equivalents.get(native_disease, native_disease)
            protections = []
            for content, in connection.execute('SELECT p.content FROM members m JOIN science_member_protection p USING(allele_id) '
                    'WHERE m.context_id=? ORDER BY m.allele_id', (context_id,)):
                protections.append(json.loads(content))
            # Missing protection inputs do not authorize a negative claim.
            uncertain = self.profile_reader is None or any(p['uncertain'] for p in protections)
            clinical = False
            native = False
            for protection in protections:
                for fact in protection['facts']:
                    target = fact['gene_id']
                    if target is not None and target.split('.')[0] != (metadata['gene_id'] or '').split('.')[0]:
                        continue
                    if fact['kind'] == 'native':
                        native = True
                        uncertain |= fact['uncertain']
                        continue
                    conditions = (_normalise_curie(condition) for condition in fact['disease_ids'])
                    if any(equivalents.get(condition, condition) == disease_group for condition in conditions):
                        clinical = True
                    else:
                        # A positive allele assertion with unresolved condition
                        # applicability cannot support symptom-only deferral.
                        uncertain = True
                    uncertain |= fact['uncertain']
            seeds = connection.execute('SELECT p.canonical_seed FROM ('
                'SELECT seed_sha256 FROM baseline_seeds WHERE context_id=? UNION '
                'SELECT seed_sha256 FROM seeds WHERE context_id=?) s '
                'JOIN seed_payloads p ON p.sha256=s.seed_sha256 ORDER BY s.seed_sha256', (context_id, context_id))
            found = False
            for content, in seeds:
                seed = _load_seed(expand(content, payloads.get))
                if seed.context.context_id != context_id:
                    raise ValueError('scientific context seed identity differs')
                support_key = header(content)['shared'].get('support_ledger')
                support = support_cache.get(support_key) if support_key is not None else None
                if support is None:
                    support = self._support(seed)
                    if support_key is not None:
                        support_cache[support_key] = support
                        if len(support_cache) > 1024:
                            support_cache.popitem(last=False)
                elif support_key is not None:
                    support_cache.move_to_end(support_key)
                association, supported, conflicts = support
                members = connection.execute('SELECT i.allele_id,i.input_sha256 FROM members m JOIN inputs i USING(allele_id) '
                    'WHERE m.context_id=? ORDER BY i.allele_id', (context_id,)) if not found else ()
                roster.add(ContextProjection(context_id, context_sha, scored[0], bool(scored[1]),
                    clinical or native and supported, uncertain or conflicts, association), members)
                found = True
            if not found:
                raise ValueError('scientific context has no seed')
        self.roster = roster.close()
        self.metrics = dict(**profile_metrics, protection_seconds=protection_seconds,
            selection_seconds=time.monotonic() - started,
            profile_preparation_seconds=0 if self._profile_collection is None else self._profile_collection['seconds'],
            profile_protection_reused=self._profile_collection is not None, **roster.counts)
        return roster
