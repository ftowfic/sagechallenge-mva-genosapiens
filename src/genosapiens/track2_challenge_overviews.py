"""Research overviews from the same selections and retained claims as the paper."""
import json
from pathlib import Path

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_report_model import relative_destination
from genosapiens.track2_scientific_reader import read_scientific


def _record(root, reference):
    path,_,pointer=relative_destination(reference).partition('#')
    file=Path(root)/path
    require(file.resolve()==file.absolute(), 'An overview source contains a symbolic link.')
    value=read_scientific(file)
    if pointer:
        require(pointer.startswith('/'), 'An overview source requires a JSON pointer.')
        for part in pointer[1:].split('/'):
            part=part.replace('~1','/').replace('~0','~')
            value=value[part] if isinstance(value,dict) else value[int(part)]
    return value


def _interpretations(value):
    """Read existing interpretation records. Do not infer an assay from prose."""
    if isinstance(value,dict):
        if value.get('policy')=='track2-functional-observation-interpretation/1.0.0':
            yield value
        else:
            for child in value.values():yield from _interpretations(child)
    elif isinstance(value,list):
        for child in value:yield from _interpretations(child)


def mechanism_records(context, selected_scenarios, context_path):
    observations={o['observation_id']:o for o in context['observations']}
    rows=[]
    for index,hypothesis in enumerate(context['hypotheses']):
        if selected_scenarios and hypothesis['scenario_id'] not in selected_scenarios:continue
        source=context_path+'#/hypotheses/'+str(index)
        variants=[]
        for oid in hypothesis['observation_ids']:
            observation=observations[oid]
            descriptions=[observation.get('reported_protein'),observation.get('reported_coding')]
            for transcript in observation.get('transcripts',[]):
                descriptions.extend((transcript.get('protein'),transcript.get('coding')))
            variants.append(dict(label=next((v for v in descriptions if v),observation['label']),
                origin='Observed variant' if context['input_origin']=='diagnosis_guided_vcf' else 'Reported genotype',
                source_reference=source,observation_id=oid,applicability='Variant identity does not establish its functional effect.'))
        facts=hypothesis['mechanism_inputs'];effects=[]
        for allele in facts.get('allele_evidence',[]):
            status=allele['curation_status']
            consequence=allele.get('curated_consequence') if status=='curated' else allele.get('predicted_consequence')
            effects.append(dict(label=consequence or 'Functional effect unresolved',
                origin='Curated interpretation' if status=='curated' else 'Prediction' if status=='predicted' else 'Unresolved',
                source_reference=source,evidence_ids=allele['evidence_ids'],
                applicability='An interpretation does not establish rescue in this genotype.'))
        proof=(facts.get('prediction_acceptance') or {}).get('evidence',{}).get('scientific_proof')
        interpreted=list(_interpretations(json.loads(proof))) if proof else []
        for proof_record in hypothesis.get('construction_proofs',[]):interpreted.extend(_interpretations(proof_record))
        functional=[];seen=set()
        for record in interpreted:
            for observation in record['observations']:
                key=observation['statement'],tuple(record['citations'])
                if key in seen:continue
                seen.add(key)
                functional.append(dict(label=observation['statement'],origin='Source-reported observation with agent interpretation',
                    source_reference=source,citations=record['citations'],authority=record['authority'],
                    human_review=record['human_review'],scope=observation['scope'],
                    applicability='Transfer to the analyzed genotype and tissue remains unconfirmed.',
                    measured_molecular_direction=record['measured_molecular_direction'],
                    tissue_generalization=record['tissue_generalization'],treatment_response=record['treatment_response']))
        # Explicit functional records can exist without a generated interpretation.
        # Quote these records with their original scope, without classifying prose.
        evidence_ids={key for allele in facts.get('allele_evidence',[]) for key in allele['evidence_ids']}
        for evidence_index,evidence in enumerate(context.get('source_evidence',[])):
            if evidence.get('kind')!='functional' or evidence['evidence_id'] not in evidence_ids:continue
            functional.append(dict(label=evidence['statement'],origin='Source-reported functional evidence (quotation)',
                source_reference=context_path+'#/source_evidence/'+str(evidence_index),
                evidence_ids=[evidence['evidence_id']],locators=evidence.get('locators',[]),
                patient_labels=evidence.get('patient_labels',[]),model_ids=evidence.get('model_ids',[]),
                limitations=evidence.get('limitations',[]),
                applicability='The source scope does not establish transfer to the analyzed genotype or tissue.'))
        if not functional:
            functional=[dict(label='Applicable process evidence unresolved',origin='Unresolved',source_reference=source,
                applicability='The retained record does not establish this biological step.')]
        objectives=[dict(label=o.get('desired_modulation','Unresolved').replace('_',' ')+' function',
            origin='Proposed intervention',objective_id=o['objective_id'],source_reference=source,
            applicability='The proposed direction does not establish drug action or functional rescue.') for o in hypothesis['objectives']]
        rows.append(dict(hypothesis=hypothesis['label'],scenario_id=hypothesis['scenario_id'],
            gene_labels=[g.get('symbol') or g['gene_id'] for g in hypothesis['genes']],
            variants=variants,effects=effects,functional=functional,objectives=objectives,source_reference=source))
    return rows


def overview_records(root, model):
    context=_record(root,model['artifacts']['context']['path'])
    selected={r['scenario_id'] for r in model['selections']}
    mechanisms=mechanism_records(context,selected,model['artifacts']['context']['path'])
    labels={oid:h['label'] for h in model['summary']['hypotheses'] for oid in h['objective_ids']}
    coverage=[dict(hypothesis=labels[o['objective_id']],objective_id=o['objective_id'],
        routes=o['routes'],source_reference=model['artifacts']['coverage']['path']) for o in model['summary']['objectives']]
    records=[dict(id='challenge-coverage',number=1,kind='coverage',rows=coverage,
        historical_ties=model['summary']['historical_ties'],derived_methods=model['summary']['derived_methods'],
        caption='Figure 1 identifies completed, negative, unavailable, and unfinished searches. Labels retain route scope and historical ties. '
            'Coverage does not rank medication efficacy.')]
    if mechanisms:
        records.append(dict(id='challenge-mechanism',number=2,kind='mechanism',rows=mechanisms,
            caption='Figure 2 separates variant observations, interpretations, reported functional evidence, and proposed intervention steps. '
                'Arrows connect the hypothesis steps. They do not demonstrate rescue. Source data retain applicability limits and evidence references.'))
    evidence=[]
    for selected in model['selections']:
        if selected['group']=='mechanism':continue
        source=_record(root,selected['source_reference'])
        item=dict(name=selected['name'],gene=selected.get('gene_symbol') or selected['gene_id'],
            hypothesis=labels.get(selected['objective_id'],'Objective'),objective_id=selected['objective_id'],
            method=selected['group'],expression=selected['expression'],source_reference=selected['source_reference'])
        if selected['expression']:
            item['unit']='tissue';item['states']=[]
            for row in source['tissue_score']['coverage']:
                if row['method']!=selected['group']:continue
                state=row['state']
                if row['method']!='paired_oe_kd' and not row['incomplete_jobs'] and row['n_valid']:
                    state='single_direction_nomination' if row['weaker_nomination'] else 'single_direction_no_nomination'
                item['states'].append(dict(label=row['tissue_id'],state=state,source_state=row['state'],
                    n_valid=row['n_valid'],n_undefined=row['n_undefined'],conflict=row['conflict'],
                    query_sha256=row['query_sha256']))
        else:
            item['unit']='cell';item['states']=[dict(label=r['cell']['cell_id'],method=r['cell']['method'],
                state=r['cell']['state'],conflict=bool(r['cell'].get('conflict')),
                source_reference=r['source_reference']) for r in source['cells']]
        evidence.append(item)
    if evidence:
        records.append(dict(id='challenge-evidence',number=3,kind='evidence',rows=evidence,
            caption='Figure 3 uses the main table selection. Each row retains its expression method, score, calibration identity, and tissue or cell states. '
                'Single-direction evidence does not supply paired support. Complete supplements retain every condition and detailed figure.'))
    return records


def pinned_runtime():
    return read_scientific(Path(__file__).resolve().parents[2]/'config/track2/discovery-runtime.json')['image']


def render_overviews(root, model, destination, *, runtime_image, run_renderer=None):
    """Render a new private figure component in the pinned R runtime."""
    import subprocess
    from genosapiens.track2_scientific_storage import publish_json, file_identity, json_identity
    from genosapiens.track2_challenge_revision import inventory
    script=Path(__file__).resolve().parents[2]/'scripts/render_track2_challenge_overviews.R'
    require(runtime_image==pinned_runtime(), 'Challenge rendering requires the configured pinned runtime image.')
    destination=Path(destination).absolute()
    require(destination.resolve()==destination, 'The overview destination contains a symbolic link.')
    destination.mkdir(mode=0o700,parents=True,exist_ok=False)
    records=overview_records(root,model)
    publish_json(destination/'overview-data.json',records)
    if run_renderer is None:
        subprocess.run(['Rscript',str(script),str(destination/'overview-data.json'),str(destination)],
            check=True,capture_output=True,text=True)
    else:
        run_renderer(script,destination/'overview-data.json',destination)
    outputs=inventory(destination)
    expected={'overview-data.json','runtime.json'} | {r['id']+suffix for r in records for suffix in ('.png','.svg')}
    require(set(outputs)==expected, 'The challenge overview inventory is incomplete.')
    manifest=dict(schema='genosapiens.track2-challenge-figures/v1',
        binding=dict(data_sha256=json_identity(records)['sha256'],policy_sha256=model['policy_sha256'],
            renderer_sha256=file_identity(script)['sha256'],generator_sha256=file_identity(Path(__file__))['sha256'],
            runtime_image=runtime_image),figures=records,outputs=outputs)
    publish_json(destination/'manifest.json',manifest)
    return manifest


def check_overviews(root, model, figures):
    from genosapiens.track2_scientific_storage import file_identity, json_identity
    from genosapiens.track2_challenge_revision import inventory
    figures=Path(figures)
    manifest=read_scientific(figures/'manifest.json')
    script=Path(__file__).resolve().parents[2]/'scripts/render_track2_challenge_overviews.R'
    require(manifest['schema']=='genosapiens.track2-challenge-figures/v1', 'The challenge figure contract is unsupported.')
    records=overview_records(root,model)
    expected=dict(data_sha256=json_identity(records)['sha256'],policy_sha256=model['policy_sha256'],
        renderer_sha256=file_identity(script)['sha256'],generator_sha256=file_identity(Path(__file__))['sha256'],
        runtime_image=pinned_runtime())
    require(manifest['binding']==expected,
            'The challenge figure dependencies changed.')
    require(json_identity(manifest['figures'])==json_identity(records)
            and json_identity(read_scientific(figures/'overview-data.json'))==json_identity(records),
            'Challenge figure data differ from the paper evidence.')
    require(manifest['outputs']==inventory(figures,omit={'manifest.json'}), 'A challenge figure asset changed.')
    require(set(manifest['outputs'])=={'overview-data.json','runtime.json'} |
        {r['id']+suffix for r in records for suffix in ('.png','.svg')}, 'A challenge figure asset is missing.')
    return manifest


def figure_references(root, model, component):
    """Bind checked assets to the same presentation selection used by the table."""
    from genosapiens.track2_challenge_revision import reference
    manifest=check_overviews(root,model,component)
    prefix=Path(component).relative_to(Path(root)).as_posix()+'/'
    return [dict(label='Figure '+str(row['number']),caption=row['caption'],
        png=reference(root,prefix+row['id']+'.png'),svg=reference(root,prefix+row['id']+'.svg'),
        data=reference(root,prefix+'overview-data.json'),
        provenance=dict(kind=row['kind'],record_id=row['id'],
            manifest=reference(root,prefix+'manifest.json'))) for row in manifest['figures']]
