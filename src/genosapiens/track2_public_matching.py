"""Exact offline reuse of reviewed public results, without case-readiness transfer.

The public installation owns the source reports and reviewed identity mappings.
Case packages contain only matching decisions and small public result summaries.
No public review is rebound to a case and no drug analysis is rerun here.
"""
import copy
from dataclasses import asdict
from datetime import datetime
import json
from pathlib import Path
import re
import shutil

from genosapiens.track2_operator_inputs import require, small_json, write_new
from genosapiens.track2_published_genotype import canonical, sha256, _source_file, load_published_genotype

SCHEMA = 'genosapiens.track2-public-matching/v1'
MATCH_SCHEMA = 'genosapiens.track2-public-match-decisions/v1'
PREFIX = 'public-match:'
IDENTITY = ('reference_identity_id', 'reference_role', 'contig', 'position', 'reference', 'alternate')


def _accepted(review):
    require(review['decision'] == 'accepted' and review['reviewer_type'] == 'human'
        and review['reviewer'] and review['rationale']
        and review['independence'] in {'independent', 'self_review_disclosed'}
        and datetime.fromisoformat(review['reviewed_at'].replace('Z', '+00:00')).tzinfo is not None,
        'Public evidence requires its actual accepted source review.')


def _mapping(mapping, review, data, source_hashes):
    require(set(mapping) == {'schema', 'classification', 'public_input_sha256', 'disease_id',
        'gene_id', 'reference_manifest_sha256', 'alleles', 'phase', 'sources'}
        and mapping['schema'] == 'genosapiens.track2-public-genotype-mapping/v1'
        and mapping['classification'] == 'public', 'Unsupported public identity mapping.')
    require(mapping['public_input_sha256'] == data['input_sha256']
        and re.fullmatch(r'(?:OMIM|ORPHA|MONDO):[0-9]+', mapping['disease_id'])
        and re.fullmatch(r'ENSG[0-9]+', mapping['gene_id'])
        and re.fullmatch(r'[0-9a-f]{64}', mapping['reference_manifest_sha256']),
        'Public mapping requires exact disease, gene, reference and source bindings.')
    alleles = {a['allele_id']: a for a in data['genotype']['alleles']}
    mapped = {a['public_allele_id']: a for a in mapping['alleles']}
    require(len(mapped) == len(mapping['alleles']) and set(mapped) == set(alleles),
        'Public mapping must cover the complete published genotype.')
    keys = []
    for key, item in mapped.items():
        require(set(item) == {*IDENTITY, 'public_allele_id', 'copies'}
            and item['copies'] == alleles[key]['copies'] and type(item['copies']) is int
            and type(item['position']) is int and item['position'] > 0
            and all(isinstance(item[k], str) and item[k] for k in IDENTITY if k != 'position')
            and item['reference_role'] == 'challenge_alignment'
            and re.fullmatch('[ACGT]+', item['reference']) and re.fullmatch('[ACGT]+', item['alternate'])
            and item['reference'] != item['alternate'], 'Public allele mapping is incomplete.')
        keys.append(tuple(item[k] for k in IDENTITY))
    require(len(set(keys)) == len(keys) and mapping['phase'] == (
        'not_applicable' if data['genotype']['zygosity'] == 'homozygous' else 'in_trans'),
        'Public mapping copies or phase disagree with the published genotype.')
    require(set(review) == {'schema', 'mapping_sha256', 'source_hashes', 'reviewer', 'reviewer_type',
        'independence', 'decision', 'rationale', 'limitations', 'reviewed_at'}
        and review['schema'] == 'genosapiens.track2-public-mapping-review/v1'
        and review['mapping_sha256'] == sha256(canonical(mapping))
        and review['source_hashes'] == source_hashes and source_hashes and review['limitations'],
        'Public mapping review is incomplete or stale.')
    _accepted(review)


def export_public_evidence(bundle, plan, envelope, report, output, *, mapping=None,
                           mapping_review=None, allow_synthetic=False):
    """Retain completed public results and their bindings, without scientific replay."""
    from genosapiens.track2_synthesis import SCHEMA as PLAN_SCHEMA
    from genosapiens.track2_freeze import load_frozen_public_envelope
    root, target, report = Path(bundle), Path(output), Path(report)
    public = load_published_genotype(root, allow_synthetic=allow_synthetic)
    plan_path = Path(plan)
    p = small_json(plan_path)
    review_path = plan_path.with_name('plan-review.json')
    review = small_json(review_path)
    _accepted(review)
    manifest = small_json(report/'manifest.json')
    frozen = load_frozen_public_envelope(envelope)
    require(p['schema'] == PLAN_SCHEMA and p['classification'] == 'public'
        and p['input_origin'] == 'published_patient_genotype'
        and manifest['classification'] == 'public'
        and manifest['analysis_origin'] == asdict(public.admitted_origin())
        and review['plan_sha256'] == sha256(plan_path.read_bytes())
        and manifest['bindings']['plan'] == review['plan_sha256']
        and manifest['bindings']['plan_review'] == sha256(review_path.read_bytes())
        and all(manifest['bindings'].get(k) == v for k, v in review['source_hashes'].items())
        and frozen['manifest_checksums'] == manifest['bindings']
        and manifest['frozen_envelope_sha256'] == sha256(Path(envelope).read_bytes())
        and p['scope_sha256'] == frozen['public_scope']['content_sha256'],
        'Public report, plan, reviews and frozen envelope disagree.')
    # Authenticate the retained result under its recorded runtime. Regenerating
    # figure records with today's renderer would reinterpret an accepted output.
    for name, digest in manifest['outputs'].items():
        require(sha256(_source_file(report,name).read_bytes()) == digest, 'Public report output changed.')
    require(small_json(report/'readiness.json') == manifest['readiness']
        and manifest['readiness']['research_report_ready'],
        'Public research report is not complete.')
    require(not target.exists() and target.resolve() == target.absolute(),
        'Public evidence export requires a new real directory.')
    target.mkdir(parents=True, mode=0o700)
    files = {}
    def retain(path, name):
        raw = path.read_bytes()
        files[name] = write_new(target/name, raw)
    for name in ('genotype.json', 'review.json'):
        retain(_source_file(root, name), 'input/'+name)
    for source in public.data['sources']:
        retain(_source_file(root, source['path']), 'input/'+source['path'])
    for path, name in ((plan_path, 'plan.json'), (review_path, 'plan-review.json'),
                       (Path(envelope), 'envelope.json'), (report/'manifest.json', 'report/manifest.json')):
        retain(path, name)
    for name in manifest['outputs']:
        retain(_source_file(report, name), 'report/'+name)
    mapping_hashes = {}
    require((mapping is None) == (mapping_review is None), 'Supply both mapping and its actual review.')
    if mapping is not None:
        m = small_json(mapping)
        for source in m['sources']:
            require(set(source) == {'path', 'sha256', 'citation', 'url', 'locator'}
                and all(isinstance(v, str) and v for v in source.values())
                and source['url'].startswith('https://'), 'Mapping source provenance is incomplete.')
            path = _source_file(Path(mapping).parent, source['path'])
            require(sha256(path.read_bytes()) == source['sha256'], 'Mapping source changed.')
            require(source['path'] not in mapping_hashes, 'Duplicate mapping source.')
            mapping_hashes[source['path']] = source['sha256']
            retain(path, 'mapping-sources/'+source['path'])
        _mapping(m, small_json(mapping_review), {**public.data, 'input_sha256':public.input_sha256}, mapping_hashes)
        files['mapping.json'] = write_new(target/'mapping.json', m)
        retain(Path(mapping_review), 'mapping-review.json')
    return write_new(target/'evidence.json', {'schema':SCHEMA, 'classification':'public',
        'evidence_origin':public.data['evidence_origin'], 'files':files,
        'mapping_status':'reviewed' if mapping is not None else 'missing_reviewed_mapping'})


def empty_public_evidence(output):
    return write_new(Path(output)/'evidence.json', {'schema':SCHEMA, 'classification':'public',
        'evidence_origin':None, 'files':{}, 'mapping_status':'no_public_evidence_installed'})


def load_public_evidence(root, *, allow_synthetic=False):
    root = Path(root)
    raw = (root/'evidence.json').read_bytes()
    data = small_json(root/'evidence.json')
    require(set(data) == {'schema','classification','evidence_origin','files','mapping_status'}
        and data['schema'] == SCHEMA and data['classification'] == 'public', 'Invalid public evidence installation.')
    require(data['evidence_origin'] in (None, 'primary_literature') or
        allow_synthetic and data['evidence_origin'] == 'independently_synthetic',
        'Synthetic public evidence cannot support an actual case.')
    actual = {p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file()}
    require(actual == set(data['files']) | {'evidence.json'}, 'Public evidence file closure changed.')
    for name, digest in data['files'].items():
        require(sha256(_source_file(root, name).read_bytes()) == digest, 'Installed public evidence changed.')
    if not data['files']:
        require(data['evidence_origin'] is None and data['mapping_status'] == 'no_public_evidence_installed',
            'Empty public evidence has an invalid status.')
        return {'sha256':sha256(raw), 'mapping_status':data['mapping_status'], 'summary':None, 'mapping':None}
    public = load_published_genotype(root/'input', allow_synthetic=allow_synthetic)
    require(public.data['evidence_origin'] == data['evidence_origin'], 'Public evidence origin changed.')
    plan = small_json(root/'plan.json')
    review = small_json(root/'plan-review.json')
    _accepted(review)
    manifest = small_json(root/'report/manifest.json')
    from genosapiens.track2_freeze import load_frozen_public_envelope
    from genosapiens.track2_input import source_mechanism_inputs
    frozen = load_frozen_public_envelope(root/'envelope.json')
    source_mechanism_inputs(public, plan['mechanism_inputs'])
    require(manifest['analysis_origin'] == asdict(public.admitted_origin())
        and manifest['bindings']['plan'] == data['files']['plan.json'] == review['plan_sha256']
        and manifest['bindings']['plan_review'] == data['files']['plan-review.json']
        and all(manifest['bindings'].get(k) == v for k,v in review['source_hashes'].items())
        and frozen['manifest_checksums'] == manifest['bindings']
        and manifest['frozen_envelope_sha256'] == data['files']['envelope.json']
        and manifest['readiness'] == small_json(root/'report/readiness.json')
        and all(data['files'].get('report/'+k) == v for k,v in manifest['outputs'].items()),
        'Installed public result bindings changed.')
    mapping = None
    if data['mapping_status'] == 'reviewed':
        mapping = small_json(root/'mapping.json')
        hashes = {s['path']:s['sha256'] for s in mapping['sources']}
        require(all(data['files'].get('mapping-sources/'+k) == v for k,v in hashes.items()),
            'Installed mapping source bindings changed.')
        _mapping(mapping, small_json(root/'mapping-review.json'),
            {**public.data, 'input_sha256':public.input_sha256}, hashes)
    else:
        require(data['mapping_status'] == 'missing_reviewed_mapping'
            and 'mapping.json' not in data['files'] and 'mapping-review.json' not in data['files'],
            'Public evidence mapping status changed.')
    # Small summaries reference the original report's complete gates and figures.
    tables = small_json(root/'report/tables.json')
    lincs_path = root/'report/lincs/results.json'
    lincs = small_json(lincs_path) if lincs_path.exists() else None
    outcomes = small_json(root/'report/outcomes.json')
    outcome_fields = ('objective_id','primary_candidate_ids','research_candidate_ids',
        'scientific_outcome','track2_objective_status')
    names = {} if lincs is None else {r['medication_entity_key']:
        next((v['name'] for v in r['cells'].values() if v.get('name')),r['medication_entity_key'])
        for r in lincs['cross_cell']}
    summary = {'public_scenario_id':public.data['scenario_id'], 'input_sha256':public.input_sha256,
        'input_review':public.review, 'plan_review':review,
        'report_manifest_sha256':data['files']['report/manifest.json'],
        'mechanism_inputs':plan['mechanism_inputs'],
        'outcomes':[{k:row[k] for k in outcome_fields} for row in outcomes],
        'source_readiness':small_json(root/'report/readiness.json'),
        'source_catalog':tables['evidence_catalog'],
        'lincs_research_eligibility':None if lincs is None else {
            'policy':lincs['research_eligibility_policy'],
            'eligible_medications':[{**row, 'name':names[row['medication_entity_key']]}
                for row in lincs['candidate_decisions'] if row['research_eligible']]},
        'limitations':public.data['limitations']}
    return {'sha256':sha256(raw), 'mapping_status':data['mapping_status'],
        'mapping':mapping, 'summary':summary}


def _identity(components, scenario):
    """Use actual normalized identity and genotype-property records from the builder."""
    evidence = {e['evidence_id']:e for e in components['evidence']}
    proofs = [evidence[k] for k in scenario['evidence_ids'] if k.startswith('scoped-model:')]
    require(len(proofs) == 1, 'Public matching requires the actual constructed genotype proof.')
    row = proofs[0]
    require(sha256(row['assertion'].encode()) == row['source_identity'], 'Genotype proof changed.')
    proof = json.loads(row['assertion'])
    genotype = proof['genotype']; hypothesis = genotype['hypothesis']
    require(hypothesis['hypothesis_id'] == scenario['genotype_id'], 'Foreign genotype matching proof.')
    alleles = {}
    for key in scenario['evidence_ids']:
        if not key.startswith('scoped-allele:'):
            continue
        row = evidence[key]
        require(sha256(row['assertion'].encode()) == row['source_identity'], 'Allele proof changed.')
        allele = json.loads(row['assertion']); alleles[allele['allele_id']] = allele
    require(set(alleles) == set(proof['normalized_to_package_allele_ids'])
        and set(proof['normalized_to_package_allele_ids'].values()) == set(scenario['allele_ids']),
        'Matching allele domain differs from the constructed genotype.')
    records=[]
    for member in genotype['members']:
        prop=member['genotype']; allele=alleles[member['allele_id']]
        if (prop['conflicting_source_genotypes'] or prop['minimum_alt_copies'] is None
            or prop['minimum_alt_copies'] != prop['maximum_alt_copies']
            or prop['reported_ploidies'] != [2]):
            return None, 'unresolved_copy_number'
        identities=[i for i in allele['identities'] if i['stage']=='normalized']
        require(len(identities)==1 and identities[0]['reference_manifest_sha256'] == components['scope']['reference_sha256'],
            'Normalized matching identity or reference is ambiguous.')
        records.append({**{k:identities[0]['identity'][k] for k in IDENTITY}, 'copies':prop['minimum_alt_copies']})
    require(len(records)==len(alleles), 'Genotype member coverage is incomplete.')
    identity={'disease_id':scenario['disease_id'], 'gene_id':scenario['gene_id'],
        'reference_manifest_sha256':components['scope']['reference_sha256'],
        'phase':hypothesis['phase'], 'alleles':sorted(records,key=canonical)}
    if hypothesis['completeness'] != 'complete':
        return identity, 'incomplete_genotype'
    if hypothesis['phase'] not in {'not_applicable','in_trans'}:
        return identity, 'unresolved_or_incompatible_phase'
    return identity, None


def reviewed_functional_check(packet, expected, context, allele, transcript_ids):
    """Inspect a frozen reviewed import without transferring genotype effects.

    A mapped unrelated allele yields a completed empty lookup. An overlapping
    allele without a transcript or a structured effect stays unresolved.
    """
    from genosapiens.track2_scientific_planning import validate_packet
    from genosapiens.evidence_serialization import _plain
    if packet is None:
        return {'state':'unavailable','sources':[],'matches':[],'contrary':False}
    plan,_=validate_packet(packet,expected)
    files=packet['files']; export=json.loads(files['evidence.json'])
    if not {'mapping.json','mapping-review.json','input/genotype.json'} <= set(files):
        return {'state':'unavailable','sources':[expected],'matches':[],'contrary':False}
    mapping=json.loads(files['mapping.json']); data=json.loads(files['input/genotype.json'])
    _mapping(mapping,json.loads(files['mapping-review.json']),
        {**data,'input_sha256':sha256(files['input/genotype.json'].encode())},
        {s['path']:s['sha256'] for s in mapping['sources']})
    native=next(i for i in allele.identities if i.identity_id==allele.native_identity_id)
    identity=_plain(native.identity)
    matches=[]; incomplete=False
    if ((context.gene_id,context.native_disease_id,native.reference_manifest_sha256)==
            (mapping['gene_id'],mapping['disease_id'],mapping['reference_manifest_sha256'])):
        for row in mapping['alleles']:
            if any(identity.get(k)!=row[k] for k in IDENTITY): continue
            public=next(a for a in data['genotype']['alleles'] if a['allele_id']==row['public_allele_id'])
            transcript=public['transcript_id']
            if transcript and public['transcript_version'] is not None and '.' not in transcript:
                transcript+='.'+str(public['transcript_version'])
            if not transcript:
                incomplete=True
                continue
            if transcript not in transcript_ids: continue
            facts=[a for r in plan['mechanism_inputs'].values() for a in r['allele_evidence']
                   if a['allele_id']==public['allele_id'] and a['curation_status']=='curated']
            if not facts or any(a['curated_consequence'] not in {
                    'loss_of_function','gain_of_function','dominant_negative','hypomorphic','other'} for a in facts):
                incomplete=True
            for fact in facts:
                matches.append({'allele_id':public['allele_id'],'transcript_id':transcript,
                    'value':fact['curated_consequence'],'evidence_ids':fact['evidence_ids'],
                    'source_sha256':expected,'lineage':'original_reviewed_import',
                    'experimental_context':[e for e in data['evidence'] if e['evidence_id'] in fact['evidence_ids']]})
    return {'state':'incomplete' if incomplete else 'complete','sources':[expected],
        'matches':matches,'contrary':any(r['value'] in {'gain_of_function','dominant_negative','hypomorphic','other'}
                                        for r in matches),
        'rule':'reviewed-native-allele-gene-subtype-transcript; no-joint-to-allele-transfer/v1'}


def matching_decisions(components, public):
    mapping=public['mapping']; decisions=[]
    for scenario in components['scenarios']:
        identity, reason = _identity(components, scenario)
        reasons=[]
        if public['mapping_status'] != 'reviewed': reasons.append(public['mapping_status'])
        if reason: reasons.append(reason)
        if mapping is not None and identity is not None:
            expected={k:mapping[k] for k in ('disease_id','gene_id','reference_manifest_sha256','phase')}
            expected['alleles']=sorted([{k:v for k,v in a.items() if k!='public_allele_id'}
                for a in mapping['alleles']],key=canonical)
            reasons.extend(k+'_differs' for k in expected if expected[k]!=identity[k])
        decisions.append({'scenario_id':scenario['scenario_id'], 'genotype_id':scenario['genotype_id'],
            'status':'not_matched' if reasons else 'matched', 'reasons':reasons,
            'identity_sha256':None if identity is None else sha256(canonical(identity))})
    return {'schema':MATCH_SCHEMA, 'classification':components['selection']['classification'],
        'public_evidence_sha256':public['sha256'], 'coverage':public['mapping_status'],
        'decisions':decisions, 'matched_public_result':public['summary']
            if any(d['status']=='matched' for d in decisions) else None,
        'case_applicability':'not_assessed', 'case_candidate_eligibility':'not_assessed'}


def attach_matching(components, public):
    require(not any(e['evidence_id'].startswith(PREFIX) for e in components['evidence']),
        'Public matching has already been attached.')
    result=matching_decisions(components, public); raw=canonical(result)
    parts=copy.deepcopy(components)
    parts['evidence'].append({'evidence_id':PREFIX+sha256(raw)[:24], 'source_identity':sha256(raw),
        'lineage_key':public['sha256'], 'assertion':raw.decode(), 'evidence_role':'source_observation'})
    parts['evidence'].sort(key=lambda e:e['evidence_id'])
    return parts


def retained_matching(components):
    rows=[e for e in components['evidence'] if e['evidence_id'].startswith(PREFIX)]
    require(len(rows)<=1, 'Duplicate public matching decisions.')
    if not rows:return None
    row=rows[0]; raw=row['assertion'].encode(); result=json.loads(raw)
    require(result['schema']==MATCH_SCHEMA and row['source_identity']==sha256(raw)
        and row['evidence_id']==PREFIX+sha256(raw)[:24]
        and row['lineage_key']==result['public_evidence_sha256'], 'Public matching decisions changed.')
    return result


def verify_matching(components, public):
    require(retained_matching(components)==matching_decisions(components,public),
        'Public match no longer agrees with its source or constructed genotype.')


def render_matching(components, *, report_link=False):
    result=retained_matching(components)
    if result is None:return ''
    lines=['## Public evidence matching','',
        'An exact identity match links source evidence. It does not establish case applicability or medication benefit.', '',
        'Installed mapping coverage: '+result['coverage'].replace('_',' ')+'.','',
        '| Scenario | Result | Reason |','|---|---|---|']
    for row in result['decisions']:
        lines.append('| '+row['scenario_id']+' | '+row['status'].replace('_',' ')+' | '+
            (', '.join(r.replace('_',' ') for r in row['reasons']) or 'Exact reviewed genotype and disease identity')+' |')
    source=result['matched_public_result']
    if source:
        lines+=['','### Matched public-study results','',
            'The source report retains medication eligibility, scientific gates, hazards, limitations and LINCS figures.',
            'Case technical assessment and mechanism readiness remain separate.',
            'Public scenario: `'+source['public_scenario_id']+'`.','',
            'Public input reviewer: '+source['input_review']['reviewer']+'.',
            'Public plan reviewer: '+source['plan_review']['reviewer']+'.','']
        if report_link:
            lines+=['[Read the retained public report](matched-public-report/report.md).','']
        for facts in source['mechanism_inputs'].values():
            lines+=['Public mechanism direction: '+facts['direction'].replace('_',' ')+'.']
        for outcome in source['outcomes']:
            lines+=['Public research candidates: '+(', '.join(outcome['research_candidate_ids']) or 'none')+'.']
        if source['lincs_research_eligibility']:
            lines+=['','| Public LINCS research candidate | Scientific assessment |','|---|---|']
            for medication in source['lincs_research_eligibility']['eligible_medications']:
                lines+=['| '+medication['name'].replace('|','\\|')+' | '+
                    medication['scientific_assessment'].replace('_',' ')+' |']
        lines+=['','These are public research candidates. Source matching did not assess case candidate eligibility.','']
    return '\n'.join(lines)+'\n'


def copy_matched_report(public_root, output):
    """Copy only authenticated public report products into the controlled output."""
    root=Path(public_root)/'report'; target=Path(output)/'matched-public-report'
    manifest=small_json(root/'manifest.json'); target.mkdir(mode=0o700)
    hashes={}
    for name in (*manifest['outputs'], 'manifest.json'):
        source=_source_file(root,name); destination=target/name
        destination.parent.mkdir(mode=0o700,parents=True,exist_ok=True)
        shutil.copyfile(source,destination)
        hashes['matched-public-report/'+name]=sha256(destination.read_bytes())
    return hashes
