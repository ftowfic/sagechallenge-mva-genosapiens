"""Exploratory reports from completed scientific records without matrix access."""
import json
from pathlib import Path
import shutil

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_published_genotype import canonical,sha256,_source_file
from genosapiens.track2_run_state import files,atomic_record
from genosapiens.track2_scientific_storage import publish_json,publish_chunks,json_identity


def _text(value):
    if isinstance(value,(list,dict)):value=json.dumps(value,sort_keys=True)
    return str(value if value is not None else 'Unavailable').replace('|','&#124;').replace('\n',' ')


def candidate_table(rows):
    header=['Medication','Gene','Scenario','Mechanism route','LINCS route','Both routes','Direction agreement','LINCS method','Conflicts']
    yield '| '+' | '.join(header)+' |\n'+'|'+'---|'*len(header)+'\n'
    for r in rows:
        values=[r['name'],r['gene_symbol'],r['scenario_id'],r['mechanism_nominated'],r['lincs_nominated'],r['both_routes'],
            r['direction_agreement'],r['lincs_methods'],r['applicability_conflicts']+r['target_group_warnings']+r['condition_conflicts']]
        yield '| '+' | '.join(_text(v) for v in values)+' |\n'
    if not rows:yield '| None | Unavailable | Unavailable | False | False | False | Unavailable | Unavailable | See coverage |\n'


def coverage_table(plan,result):
    objectives={o['objective_id']:o for o in plan['objectives']}
    yield '| Gene | Objective | Route | Completion | Source limits |\n|---|---|---|---|---|\n'
    for r in result['route_outcomes']:
        values=[objectives[r['objective_id']]['gene_symbol'],r['objective_id'],r['route'],r['state'],r['reason_codes_result']]
        yield '| '+' | '.join(_text(v) for v in values)+' |\n'


def native_composition(plan,result):
    from genosapiens.track2_assessment import assess_native_plan
    from genosapiens.track2_report_tables import exploratory_tables
    from genosapiens.track2_no_candidate import native_objective_outcomes
    assessment=assess_native_plan(plan['exploratory_discovery'],result)
    return dict(assessment=assessment.as_dict(),tables=exploratory_tables(plan,result),
        outcomes=native_objective_outcomes(plan['exploratory_discovery'],result,assessment))


def discovery_section(plan,result,assessment):
    if result['schema']=='genosapiens.track2-discovery-results/v2':
        from genosapiens.track2_tissue_report_text import RULES
        yield '## Tissue expression evidence\n\n'
        yield 'Native discovery uses genetic tissue queries and separate drug conditions. These results remain exploratory research hypotheses.\n\n'
        if plan['input_origin']=='published_patient_genotype':
            from genosapiens.track2_public_discovery import source_section
            yield source_section(plan)
        yield from candidate_table(result['concordance'])
        yield '\n[Every evaluated medication](discovery/tissue-scores.json) includes nonnominees and negative scores.\n\n'
        yield '[All conditions](discovery/tissue-conditions.csv) retain raw correlations and exact experimental metadata.\n\n'
        yield '[Genetic contributors and source values](discovery/tissue-queries.json) retain the frozen query evidence.\n\n'
        for rule in RULES: yield rule+'\n\n'
        yield '## Route coverage\n\n'
        yield from coverage_table(plan['exploratory_discovery'],result)
        yield '\nTissue support does not satisfy the legacy same-cell eligibility rule. Clinical gates remain unchanged.\n\n'
        yield 'Patient tissue relevance, functional protein availability, and rescue efficacy require separate evidence.\n\n'
        return
    from genosapiens.track2_discovery_recurrence import recurrence_views,recurrence_chunks
    from genosapiens.track2_report_summary import exploratory_summary
    nested=plan['exploratory_discovery']
    yield exploratory_summary(plan,result,assessment)+'\n\n'
    if plan['input_origin']=='published_patient_genotype':
        from genosapiens.track2_public_discovery import source_section
        yield source_section(plan)
    yield '## Scenario scope\n\n'
    yield '| Scenario | Gene | Mechanism acceptance | Native eligibility | Objective or exclusion |\n|---|---|---|---|---|\n'
    from genosapiens.track2_mechanism import assess_mechanism_readiness
    for row in plan['mechanism_assessment']['scenarios']:
        eligibility=row['discovery_eligibility'];decision=assess_mechanism_readiness(row['mechanism_inputs'])
        values=[row['scenario_id'],row['gene_id'],decision.readiness_state,eligibility['outcome'],
            eligibility['objective_ids'] or eligibility['reasons']]
        yield '| '+' | '.join(_text(v) for v in values)+' |\n'
    yield '\n## Medication nominations\n\n'
    yield from candidate_table(result['concordance'])
    yield '\n\n'.join([
        '[Complete candidate export](discovery/candidates.csv). [All mechanism annotations](discovery/mechanisms.csv).','',
        'Assumed direction remains an assumption. Target-group membership does not establish gene-specific activity.',
        'Same-gene inhibitors retain a direction conflict. Other-gene inhibition does not establish indirect compensation.',''])
    yield from recurrence_chunks(recurrence_views(nested,result['concordance']))
    yield '\n## Route coverage\n\n'
    yield from coverage_table(nested,result)
    yield '\n\n'.join(['',
        '[Complete route coverage](discovery/route-coverage.csv). [All local results](discovery/results.json).','',
        '## Evidence limits','',
        'The source catalog defines search coverage. A missing result does not establish that no medication or mechanism exists.',
        'Cell relevance to patient tissue remains unknown. Residual protein availability and rescue efficacy require separate evidence.',
        'Expression matching does not establish production of absent protein. Genetic and drug exposure times can differ.',
        'Paired analysis uses both genetic directions in one cell. Each fallback uses one available direction.',
        'Opposite correlation and joint score remain unavailable for fallbacks. Repeated fallbacks do not establish paired research eligibility.',
        'Paired research eligibility requires two passing paired cells. Any observed failed cell or condition conflict prevents this eligibility.',
        'G0–G7 and clinical tiers remain unassessed for native nominations. These results do not establish submission readiness.','',
        ('Published genotype evidence retains its original review and source limits.' if plan['input_origin']=='published_patient_genotype' else
         'Optional study evidence: '+('retained with separate applicability limits.' if plan['mechanism_assessment']['automatic_inputs']['public_packet'] else 'unavailable.')),'',
        '## Follow-up','',
        '- Review each mechanism assertion and its source.\n- Assess tissue relevance and functional protein availability.\n'
        '- Test the stated hypothesis in an appropriate disease model.'])


def _csv(root,name,rows,fields):
    import csv,io
    def chunks():
        buffer=io.StringIO(newline='')
        writer=csv.DictWriter(buffer,fieldnames=fields,extrasaction='ignore')
        writer.writeheader()
        yield buffer.getvalue().encode()
        for row in rows:
            buffer.seek(0);buffer.truncate()
            writer.writerow({key:json.dumps(value,sort_keys=True) if isinstance(value,(list,dict)) else value
                for key in fields for value in (row.get(key),)})
            yield buffer.getvalue().encode()
    publish_chunks(root/name,chunks())


def synthesize_native(bundle,plan,review,bindings,plan_path,output,*,input_snapshot_path,expected_snapshot_sha256,
                      controlled_store,discovery_path,reference_root,figures_root,_evidence_only=False):
    from genosapiens.track2_discovery import discovery_bindings,load_completed_discovery,_read
    from genosapiens.track2_exploratory_planning import load_plan_resources
    from genosapiens.track2_discovery_union import reconcile_results
    from genosapiens.track2_discovery_resources import file_digest
    from genosapiens.track2_storage import output_storage
    from genosapiens.track2_synthesis import policies
    from genosapiens.track2_assessment import assess_native_plan
    resources=load_plan_resources(plan_path,plan['exploratory_discovery'],reference_root)
    try:
        if plan['input_origin']=='published_patient_genotype':
            from genosapiens.track2_public_discovery import check_public_plan
            check_public_plan(plan,resources,bundle=bundle)
        expected=discovery_bindings(plan,resources,plan['exploratory_discovery']['runtime_images']['discovery'])
        from genosapiens.track2_report_discovery_compatibility import load_report_discovery
        result=load_report_discovery(discovery_path,expected,allow_partial=True)
        mechanisms=_read(Path(discovery_path)/'mechanism-results.json')
        jobs={r['job_id']:_read(Path(discovery_path)/r['path']) if r['completion_sha256'] else
            dict(job_id=r['job_id'],medications=[]) for r in result['jobs']}
        queries=None
        if result['schema']=='genosapiens.track2-discovery-results/v2':
            from genosapiens.track2_tissue_discovery import retained_tissue_jobs
            jobs,queries=retained_tissue_jobs(discovery_path,result)
        union=reconcile_results(plan['exploratory_discovery'],mechanisms,jobs,resources['connection'],resources['policy'],result['route_outcomes'],queries=queries)
        require(all(json_identity(result[k])==json_identity(v) for k,v in union.items()),'Candidate reconciliation differs from retained scientific results.')
    finally:resources['connection'].close()
    bindings={**bindings,'policy':sha256(canonical(policies()[0]))}
    snapshot_raw=Path(input_snapshot_path).read_bytes();snapshot=json.loads(snapshot_raw)
    require(sha256(snapshot_raw)==expected_snapshot_sha256 and snapshot['bindings']==bindings
        and snapshot['schema']=='genosapiens.track2-private-input-snapshot/v2','Native snapshot bindings differ.')
    figures_root=Path(figures_root)
    from genosapiens.track2_lincs_report import native_render_binding
    from genosapiens.track2_scientific_completion import check_completion
    check_completion(figures_root,stage='figures',binding=native_render_binding(discovery_path))
    assessment=assess_native_plan(plan['exploratory_discovery'],result)
    composition=native_composition(plan,result)
    readiness=dict(assessment_status=assessment.status,input_report_ready=bundle.input_report_ready,
        research_report_ready=bundle.input_report_ready and assessment.complete,
        submission_ready=False,classification=plan['classification'])
    with output_storage(output,classification=plan['classification'],controlled_store=controlled_store,kind='track2-research-report') as target:
        target.mkdir(mode=0o700,parents=True,exist_ok=False);discovery=target/'discovery';discovery.mkdir(mode=0o700)
        for name,data in [('plan-source.json',plan),('input-snapshot.json',snapshot),('assessment.json',composition['assessment']),
                ('tables.json',composition['tables']),('outcomes.json',composition['outcomes']),('readiness.json',readiness),
                ('discovery/results.json',result),('discovery/concordance.json',result['concordance']),('discovery/mechanism-results.json',mechanisms)]:
            publish_json(target/name,data)
        if review is not None:atomic_record(target/'plan-review-source.json',review)
        fields=['name','canonical_active_ingredient_id','gene_id','gene_symbol','scenario_id','objective_id',
            'mechanism_nominated','lincs_nominated','both_routes','direction_agreement','lincs_methods','paired_research_eligible',
            'fallback_recurrence','applicability_conflicts','target_group_warnings','condition_conflicts']
        _csv(discovery,'candidates.csv',result['concordance'],fields)
        if 'evaluated_tissue_scores' in result:
            from genosapiens.track2_tissue_exports import publish_tissue_exports
            publish_tissue_exports(discovery,plan,result,discovery_path)
            if (figures_root/'tissues').is_dir():
                shutil.copytree(figures_root/'tissues',target/'lincs/tissues')
        from genosapiens.track2_discovery_recurrence import recurrence_views,VIEWS
        recurrence=recurrence_views(plan['exploratory_discovery'],result['concordance'])
        publish_json(discovery/'lincs-recurrence.json',recurrence)
        for view in ([] if 'evaluated_tissue_scores' in result else VIEWS):
            rows=(row for scope in recurrence['scopes'] for row in scope['views'][view])
            _csv(discovery,view+'-recurrence.csv',rows,['objective_id','gene_id','gene_symbol','scenario_id',
                'review_rank','name','ingredient_key','canonical_active_ingredient_id','agreement_group','passing_cells',
                'tested_cell_count','selected_cell_count','paired_research_eligible','fallback_recurrence',
                'objective_failed_cells','objective_condition_conflicts','cells'])
        _csv(discovery,'route-coverage.csv',result['route_outcomes'],['objective_id','route','planned_state','state','executed','reason_codes_result','result_references'])
        unique={}
        for records,field in ((result['concordance'],'mechanisms'),(mechanisms,'source_records')):
            for row in records:
                for annotation in row[field]:
                    digest=json_identity(dict(objective_id=row['objective_id'],**annotation))['sha256']
                    unique[digest]=(row['objective_id'],annotation)
        annotations=(dict(objective_id=unique[k][0],**unique[k][1]) for k in sorted(unique))
        _csv(discovery,'mechanisms.csv',annotations,['objective_id','mechanism_id','action','direction_class','description',
            'target_name','target_type','targets','references','group_warning','assumption_warning','source'])
        def lincs_links():
            for job in result['jobs']:
                if not job['completion_sha256']:continue
                local=jobs[job['job_id']];relative=Path('lincs')/local['gene_id']/job['job_id'];folder=target/relative
                folder.mkdir(parents=True,mode=0o700)
                shutil.copytree(Path(discovery_path)/'jobs'/job['job_id'],folder/'science')
                # Reports own one copy of each shared block, across all jobs.
                for block in local.get('shared_blocks',[]):
                    from genosapiens.track2_shared_blocks import check_shared_reference
                    source=check_shared_reference(discovery_path,block).parent
                    destination=discovery/'extractions'/block['block_id']
                    if not destination.exists():shutil.copytree(source,destination)
                shutil.copytree(figures_root/local['gene_id']/job['job_id'],folder/'figures')
                query=_read(folder/'science/query.json')
                from genosapiens.track2_figure_stream import condition_rows
                conditions=condition_rows(folder/'science/conditions.jsonl')
                if query['schema']=='genosapiens.track2-tissue-query/v1':
                    _csv(folder,'query-features.csv',query['genes'],['gene_id','gene_symbol','oe','kd','desired',
                        'selected','reason','source_values','missing_contributors'])
                    _csv(folder,'medications.csv',local['medications'],['name','medication_entity_key','cell_id','tissue_id',
                        'method','n_conditions','n_valid','n_undefined','witnesses','undefined_reasons','weaker_nomination'])
                    _csv(folder,'conditions.csv',conditions,
                        ['sig_id','pert_id','pert_iname','cell_id','pert_time','pert_time_unit','pert_dose','pert_dose_unit',
                         'method','r_oe','r_kd','r_desired','state','reason'])
                else:
                    _csv(folder,'query-features.csv',query['genes'],['gene_id','gene_symbol','oe','kd','desired','opposite','selected','reason'])
                    _csv(folder,'medications.csv',local['medications'],['name','medication_entity_key','cell_id','method','rank','state','r_desired','r_opposite','joint','times_hours','n_conditions','conflict','n_undefined'])
                    _csv(folder,'conditions.csv',conditions,['sig_id','pert_id','pert_iname','cell_id','pert_time','pert_time_unit','pert_dose','pert_dose_unit','method','r_desired','r_opposite','joint','state','reason'])
                yield '\n\n'.join(['### '+local['gene_id']+' / '+local['cell_id']+' / '+(local['method'] or 'unavailable'),
                    '[Selected signatures and QC]('+str(relative/'science/query.json')+'). [All conditions]('+str(relative/'conditions.csv')+').',
                    '[All local medications]('+str(relative/'medications.csv')+').'])+'\n\n'
                from genosapiens.track2_figure_stream import condition_rows
                for figure in condition_rows(folder/'figures/figure-data.jsonl'):
                    yield '!['+figure['alt']+']('+str(relative/'figures'/(figure['id']+'.png'))+')\n\n'+figure['caption']+'\n\n'
        from genosapiens.track2_report_language import report_terms,report_acknowledgement
        def report_text():
            yield from discovery_section(plan,result,assessment)
            yield '\n\n## LINCS details\n\n'
            yield from lincs_links()
            yield '\n\n## Source versions\n\n'
            for row in resources['descriptor']['source_manifests']:
                yield '`'+_text(row)+'`\n\n'
            yield report_terms(lincs=True,diagnosis=True,exploratory=True,
                tissue=result['schema']=='genosapiens.track2-discovery-results/v2')+'\n\n'+report_acknowledgement()+'\n'
        publish_chunks(target/'scientific-details.md',(part.encode() for part in report_text()))
        from genosapiens.track2_report_paper import compose_paper
        from genosapiens.track2_report_ranking import selected_native_cells
        from genosapiens.track2_report_gene_metadata import load_gene_metadata
        from genosapiens.track2_report_context import retained_gene_ids
        objectives=plan['exploratory_discovery']['objectives']
        metadata=load_gene_metadata(Path(__file__).resolve().parents[2], retained_gene_ids(bundle, plan, objectives))
        if _evidence_only:
            from genosapiens.track2_report_paper import export_scientific_supplements
            from genosapiens.track2_report_evidence import publish_native_evidence
            model=export_scientific_supplements(target,bundle=bundle,plan=plan,objectives=objectives,
                nominations=result['concordance'],coverage=result['route_outcomes'],
                selected_cells=selected_native_cells(plan['exploratory_discovery']),gene_metadata=metadata)
            return publish_native_evidence(target,model=model,bindings=bindings,readiness=readiness,
                snapshot_sha256=expected_snapshot_sha256,discovery_bindings=result['bindings'],
                source_versions=resources['descriptor']['source_manifests'])
        paper=compose_paper(target,bundle=bundle,plan=plan,objectives=objectives,
            nominations=result['concordance'],coverage=result['route_outcomes'],
            selected_cells=selected_native_cells(plan['exploratory_discovery']),readiness=readiness,
            title='# Exploratory repurposing research',gene_metadata=metadata)
        publish_chunks(target/'report.md',[paper])
        from genosapiens.track2_scientific_completion import publish_completion
        publish_completion(target,stage='native-report',binding=bindings,filename='manifest.json',
            result=dict(input_snapshot_sha256=expected_snapshot_sha256,discovery_bindings=result['bindings'],readiness=readiness))
        validate_native_report(target)
        return dict(**readiness,manifest_sha256=file_digest(target/'manifest.json'))


def export_native_evidence(*args, **kwargs):
    """Produce native scientific exports without invoking a paper composer."""
    require('_evidence_only' not in kwargs, 'The evidence export mode is fixed.')
    return synthesize_native(*args, **kwargs, _evidence_only=True)


def validate_native_report(output, *, page_context=None):
    if page_context is None:
        from genosapiens.track2_checked_pages import CheckedPageStreams
        owned=CheckedPageStreams()
        try:return validate_native_report(output,page_context=owned)
        finally:owned.close()
    from genosapiens.track2_discovery import _read,validate_job
    from genosapiens.track2_scientific_completion import check_completion,CONTROL_BYTES
    from genosapiens.track2_operator_inputs import small_json
    root=Path(output);header=small_json(root/'manifest.json',maximum=CONTROL_BYTES)
    checked=check_completion(root,stage='native-report',binding=header['binding'],filename='manifest.json')
    manifest=dict(bindings=checked['binding'],**checked['result'])
    readiness=validate_native_science(root,manifest,page_context=page_context)
    with (root/'report.md').open() as stream:
        require(stream.readline(128)=='# Exploratory repurposing research\n','Native report headline changed.')
    if (root/'presentation.json').is_file():
        from genosapiens.track2_report_checks import check_paper
        check_paper(root)
    return readiness


def validate_native_science(root, manifest, *, page_context=None):
    """Replay scientific meaning independently of a report headline or layout."""
    if page_context is None:
        from genosapiens.track2_checked_pages import CheckedPageStreams
        owned=CheckedPageStreams()
        try:return validate_native_science(root,manifest,page_context=owned)
        finally:owned.close()
    from genosapiens.track2_discovery import _read,validate_job
    root=Path(root)
    plan=_read(root/'plan-source.json');result=_read(root/'discovery/results.json');snapshot=_read(root/'input-snapshot.json')
    if plan['input_origin']=='published_patient_genotype':
        from genosapiens.track2_public_discovery import checked_public_provenance
        public=checked_public_provenance(plan)
        require(manifest['bindings']['input']==public.input_sha256
            and manifest['bindings']['input_review']==public.review_sha256,
            'The public report input identities differ.')
    require(snapshot['bindings']==manifest['bindings'] and json_identity(snapshot)['sha256']==manifest['input_snapshot_sha256']
        and manifest['bindings']['plan']==json_identity(plan)['sha256'] and result['bindings']==manifest['discovery_bindings'],
        'Native report scientific bindings differ.')
    require(snapshot['schema']=='genosapiens.track2-private-input-snapshot/v2'
        and snapshot['classification']==plan['classification']==manifest['readiness']['classification'],
        'Native report classification differs from its retained evidence.')
    composition=native_composition(plan,result)
    from genosapiens.track2_discovery_recurrence import recurrence_views
    require(json_identity(_read(root/'discovery/lincs-recurrence.json'))==
        json_identity(recurrence_views(plan['exploratory_discovery'],result['concordance'])),
        'Native recurrence differs from retained evidence.')
    for name in ['assessment','tables','outcomes']:
        require(json_identity(_read(root/(name+'.json')))==json_identity(composition[name]),'Native report decisions differ from retained evidence.')
    readiness=_read(root/'readiness.json')
    require(readiness==manifest['readiness'] and readiness['assessment_status']==composition['assessment']['status']
        and readiness['research_report_ready']==(readiness['input_report_ready'] and readiness['assessment_status']=='complete')
        and not readiness['submission_ready'],'Native report readiness changed.')
    shared_checked=set()
    for job in result['jobs']:
        if not job['completion_sha256']:continue
        matches=list((root/'lincs').glob('*/'+job['job_id']+'/science'))
        require(len(matches)==1,'Native report is missing a scientific job.')
        validate_job(matches[0],revision_root=root/'discovery',shared_checked=shared_checked)
        from genosapiens.track2_figure_stream import check_figure_pages
        check_figure_pages(matches[0].parent/'figures',_read(matches[0]/'query.json'),
            _read(matches[0]/'result.json'),matches[0]/'conditions.jsonl',context=page_context)
    from genosapiens.track2_report_tissue_checks import check_tissue_scores
    check_tissue_scores(root)
    return readiness
