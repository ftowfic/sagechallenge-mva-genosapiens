"""Column meanings for retained scientific tables and their supplementary copies."""
import re
from pathlib import Path

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_report_model import table_caption

MEANINGS = {
    'Exclusion': 'Retained reason that excludes this candidate from the stated eligibility scope.',
'Reported coding change':'Coding change reported by the source publication. Missing transcript mapping remains unresolved.',
'Reported protein change':'Protein change reported by the source publication. It does not establish measured protein abundance.',
'Copies':'Allele-copy count reported for the published genotype.',
'Review rank':'Original LINCS review rank under its retained source policy.',
'Local agreement':'Original summary of passing cells, failures, and condition conflicts.',
'Reason / next step':'Retained decision reason and required scientific follow-up.',
'Query genes':'Number of selected query features for this cell and method.',
'Tested medications':'Number of medication entities with retained local results.',
'Local nominations':'Number of medications that passed the local method.',
'Failed':'Number of failed comparisons under the source method.',
'Unassessed':'Number of comparisons without an assessed result.',
'Agreement':'Original summary of cell-specific agreement and conflicts.',
'Pass / tested':'Counts of passing and tested cells in the original view.',
'Direction +':'Number of cells with positive desired-direction correlation.',
'Review score':'Original LINCS review priority. It is separate from the main presentation score.',
'LINCS result':'Original category of cell-specific expression evidence.',
'Medications':'Number of medication entities in the stated source category.',
'What happens next':'Original research follow-up requirement for this category.',
'Check':'Scientific gate or evidence check under the original policy.',
'LINCS discovery candidates without drug-specific review':'Original gate coverage for expression nominations without independent medication review.',
'Coverage State':'Original coverage state for the stated scientific gate.',
'Gate Reason Codes':'Exact source reason codes for the gate decision.',
'Gate Review State':'Original review status for the gate evidence.',
'G0':'Original input and identity gate decision.',
'G1':'Original current-approval gate decision.',
'G2':'Original lesion and mechanism qualification gate decision.',
'G3':'Original drug-action and contradiction gate decision.',
'G4':'Original evidence-independence gate decision.',
'G5':'Original exposure feasibility gate decision.',
'G6':'Original safety and hazard gate decision.',
'G7':'Original experimental validation-plan gate decision.',
'Review state':'Original scientific review status and its authority.',
'Research eligibility':'Original research eligibility decision. It does not establish treatment benefit.',
'Scientific tier':'Original scientific evidence tier.',
'Scope':'Declared evidence, selection, or experimental scope of this row.',
'Coverage':'Original evidence availability and completeness within the declared bounds.',
'Rank':'Original table rank under the source ordering policy. It is not a treatment recommendation.',
'Hypothesis':'Retained proposed mechanism or intervention, separate from a proven effect.',
'Gene':'Gene identifier or symbol from the retained source mapping.',
'Method':'Actual comparison or interpretation method used for this row.',
'Type':'Source record or observation category.',
'Provenance':'Source identities, versions, and methods that produced this record.',
'Active Ingredient Id':'Canonical active ingredient identifier after source identity resolution.',
'Approval Jurisdictions':'Jurisdictions identified by the retained approval evidence.',
'Causal Paths':'Retained causal chains supporting the candidate rationales.',
'Child And General Applicability':'Separate source assessments for children and the broader disease context.',
'Chromosome Instability And Other Safety':'Retained chromosome-instability and other hazard assessments.',
'Current Management':'Source-described management context. It is not a treatment recommendation.',
'Evidence Dimensions':'Original evidence assessments across the declared scientific dimensions.',
'G0 G7':'Original decisions for gates G0 through G7. Presentation scores do not change them.',
'Generation Rationales':'Source explanations for the original candidate nominations.',
'Generation Routes':'Scientific routes that nominated this ingredient.',
'Independent Rationale Count':'Original count of independent causal rationales. It is distinct from E.',
'Mechanism Alternative Robustness':'Source assessment of support across competing mechanisms.',
'Negative Evidence':'Retained findings that oppose or limit the proposed hypothesis.',
'Patient Context':'Available patient-specific context and its source limitations.',
'Priority Tier':'Original scientific evidence tier, independent of presentation rank.',
'Proposed Actions':'Source-bound research actions or follow-up requirements.',
'Reachable Polypharmacology':'Other target actions at potentially achievable exposure.',
'Research Eligibility':'Original decision on inclusion as a research candidate.',
'Research Submission Review':'Separate retained review of research submission eligibility.',
'Resolved Drug Identity':'Source-bound mapping to the canonical active ingredient.',
'Route Coverage':'Actual route states and declared search bounds.',
'Selected Endpoint Id':'Exact identifier for the selected biological endpoint.',
'Target Hypothesis Ids':'Exact identifiers for the retained target hypotheses.',
'Target Liabilities':'Source-supported hazards from modulating the target.',
'Target Site Exposure':'Source assessment of medication exposure at the relevant biological site.',
'Uncertainties':'Retained evidence gaps and interpretation limits.',
'Gate Id':'Exact gate identifier from G0 through G7.',
'Gate Outcome':'Original pass, fail, or unresolved gate decision.',
'Gate Reasons':'Source reasons for the gate decision.',
'Gate Evidence Ids':'Exact source evidence identifiers supporting the gate decision.',
'Evidence Gaps':'Missing evidence required by the original scientific assessment.',
'Resolution Experiments':'Source-defined experiments that could resolve the evidence gaps.',
'AD':'Read counts for reference and alternate alleles in source order.',
'AD/DP state':'Availability and consistency of allele counts and total read depth.',
'ALT copies':'Observed or possible alternate-allele copy count.',
'Allele':'Retained allele identifier or source allele description.',
'Allele Evidence':'Source evidence for the individual allele effects.',
'Alternate':'Alternate allele sequence at the stated location.',
'Alternatives':'Retained competing interpretations or transcript mappings.',
'Any result?':'Whether the declared route returned a result within its bounds.',
'Both routes':'Whether both evidence routes nominated this ingredient for this objective.',
'Causal Proximity':'How directly the measured endpoint relates to the disease mechanism.',
'Causal Role':'Whether the target acts directly or through a supported indirect path.',
'Child Applicability':'Applicability of the source evidence to children.',
'Completion':'Actual assessment or route completion state.',
'Confidence Inputs':'Retained evidence inputs for confidence assessment.',
'Conflicts':'Contradictory source observations or conditions that remain unresolved.',
'Consequence':'Retained transcript consequence term. It does not establish measured function.',
'Contract SHA-256':'Checksum of the exact contract used by the source operation.',
'Copy number':'Retained copy-number estimate and its source interpretation.',
'Coverage outcome':'Actual search outcome within the declared bounds.',
'DP':'Total read depth at this observation.',
'Decision Reasons':'Source reasons for the recorded assessment decision.',
'Derivation':'Source method or retained mapping that produced this field.',
'Direction':'Retained direction of the effect or required intervention.',
'Direction / result':'Reported action direction or source result for this target.',
'Direction agreement':'Whether the retained actions support the objective direction.',
'Effect Magnitude':'Retained effect estimate. Source records define its units and uncertainty.',
'Endpoint Class':'Category of the measured or proposed outcome.',
'Evidence IDs':'Exact identifiers for the supporting source records.',
'Evidence type':'Source evidence category. It is distinct from the two-route score count.',
'Executed?':'Whether the source operation actually ran.',
'Expected Direction':'Direction required for the selected endpoint.',
'Expected ploidy':'Expected chromosome-copy count under the retained context.',
'Exposure':'Source assessment of activity at achievable medication exposure.',
'FORMAT key':'Exact sample-field identifier from the VCF FORMAT declaration.',
'GT':'Source sample genotype.',
'General MVA Applicability':'Applicability to mosaic variegated aneuploidy beyond the exact subtype.',
'Genotype Completeness':'Whether the retained genotype satisfies the stated disease model.',
'Geometry':'Retained structural observation shape and interval representation.',
'Harm Conditions':'Conditions under which the proposed intervention could cause harm.',
'Inheritance Consistency':'Agreement with the retained disease inheritance model.',
'Input Origin':'Published genotype or diagnosis-guided sample origin.',
'Intervention Objective':'Exact retained objective for the proposed target modulation.',
'Joint Evidence':'Evidence for the combined genotype effect.',
'LINCS method':'Actual paired, knockdown-only, or overexpression-only comparison.',
'LINCS route':'Nomination status from the local expression route.',
'Local results':'Retained cell-specific scores and result states.',
'Meaning':'Definition of the term within this report.',
'Meaningful Interpretation':'Supported meaning of the endpoint result within its source scope.',
'Measurement':'Assay or readout used to assess the endpoint.',
'Mechanism Confidence':'Retained confidence assessment for the proposed mechanism.',
'Mechanism acceptance':'Retained mechanism readiness or prediction acceptance decision.',
'Mechanism route':'Nomination status from the target-drug linkage route.',
'Missing Direction':'Genetic perturbation direction without usable evidence.',
'Modifiability':'Evidence that intervention can change the selected endpoint.',
'Name':'Retained preferred source name.',
'Native eligibility':'Original discovery eligibility, independent of presentation scores.',
'Nominal source interval':'Original structural interval before any retained normalization.',
'Objective':'Exact intervention objective identifier.',
'Objective condition conflicts':'Cells with condition conflicts within this objective.',
'Objective failed cells':'Cells with failed comparisons within this objective.',
'Objective or exclusion':'Retained intervention objective or reason for excluding that hypothesis.',
'Observations':'Source measurements or reported findings.',
'Observed ploidy':'Observed chromosome-copy count in the source genotype.',
'Outcome':'Retained assessment result within its declared scope.',
'PGT':'Physical phasing genotype from the source FORMAT field.',
'PID':'Physical phasing identifier from the source FORMAT field.',
'PS':'Phase-set identifier from the source FORMAT field.',
'Paired research eligibility':'Original research eligibility under the paired-expression policy.',
'Pass / tested / selected':'Counts of passing, tested, and selected cells within the stated view.',
'Phase':'Retained relationship of allele copies across chromosomes.',
'Phase state':'Availability and interpretation of the retained phase evidence.',
'Policy Version':'Exact source policy version used for the decision.',
'PolyPhen':'Retained PolyPhen prediction. It does not establish molecular mechanism.',
'Predicted direction':'Predicted functional direction, distinct from a measured effect.',
'Prediction Acceptance':'Recorded decision on whether the prediction permits hypothesis generation.',
'Predictor disagreement':'Retained disagreement among consequence or effect predictors.',
'Quality':'Source quality measurements and their original acceptance state.',
'Raw INFO':'Unchanged source VCF INFO value.',
'Raw QUAL':'Unchanged source VCF quality value.',
'Raw value':'Unchanged source value before typed interpretation.',
'Readiness State':'Original mechanism or assessment readiness state.',
'Ready For Drug Retrieval':'Whether the source target assessment permits medication retrieval.',
'Reason':'Source explanation for the recorded state.',
'Reasons':'Source explanations for the recorded state.',
'Reference declaration':'Exact reference identity or assembly declaration.',
'Release':'Version of the retained source resource.',
'Relevance':'Evidence that the endpoint applies to this disease mechanism.',
'Required Modulation':'Required increase or decrease in target activity.',
'Required copies':'Affected-allele copy count required by the disease model.',
'Research status':'Original research assessment or eligibility state.',
'Resolution':'Status of the source identity or interpretation mapping.',
'Result':'Retained result, including unavailable and incomplete states.',
'Route':'Scientific search route under the retained plan.',
'Row ID':'Exact stable identifier for the source evidence row.',
'SIFT':'Retained SIFT prediction. It does not establish molecular mechanism.',
'Safety':'Retained hazard assessment, including unresolved risks.',
'Scenario':'Exact identifier for the separate genotype or mechanism hypothesis.',
'Scope hits':'Retained observations that meet the declared scope rules.',
'Selection':'Original inclusion or exclusion decision under the source policy.',
'Signed Rescue Path':'Retained sequence of causal effects and intervention directions.',
'Source ALT occurrences':'Count of alternate-allele occurrences in the source records.',
'Source SHA-256':'Checksum of the exact retained source bytes.',
'Source allele':'Original allele identity before retained normalization.',
'Source genotype':'Original genotype representation and its source attribution.',
'Source limits':'Declared source bounds and interpretation limitations.',
'Source observation':'Exact identifier or description of the original observation.',
'Source selection outcome':'Original source inclusion or exclusion decision.',
'Stage':'Biological stage or original workflow stage, as specified by the source table.',
'State':'Original result state under the applicable source policy.',
'Subtype Applicability':'Applicability to the exact disease subtype or a broader context.',
'Support':'Retained evidence supporting the stated interpretation.',
'Target Id':'Exact intervention target identifier.',
'Technical Support':'Source assessment of the observation quality and technical validity.',
'Term':'Expression defined in this report.',
'Tissue':'Tissue or cell context of the source evidence.',
'Transcript':'Retained transcript identifier and version when available.',
'Transcript Relevance':'Evidence that the retained transcript applies to the disease context.',
'Translational Maturity':'Original stage of experimental support for the proposed intervention.',
'Typed QUAL':'Source quality value after typed parsing. Missing values remain missing.',
'What we found':'Retained route findings within the declared search bounds.',
'Balance':'Alternate-allele depth fraction under the retained quality calculation.',
'Agreement group':'Original summary of passing cells, failed cells, and condition conflicts.',
}
START = '<!-- report-column-caption -->'
END = '<!-- /report-column-caption -->'


def scientific_text(text):
    """Recover legacy scientific wording for its unchanged semantic checks."""
    return re.sub(re.escape(START)+r'.*?'+re.escape(END)+r'\n\n', '', text, flags=re.S)


def annotate_tables(text, filename, cell_ids=()):
    from genosapiens.track2_report_paper import DEFINITIONS
    from genosapiens.track2_report_supplements import FIGURE_DEFINITIONS
    definitions = {key.casefold():value for source in (DEFINITIONS,FIGURE_DEFINITIONS,MEANINGS) for key,value in source.items()}
    for cell in cell_ids:
        definitions[cell.casefold()]='Local joint correlation for '+cell+'. Missing values remain unavailable.'
    rows = scientific_text(text).splitlines(keepends=True)
    result, captions = [], {}
    for index, line in enumerate(rows):
        if (line.startswith('|') and index+1<len(rows) and rows[index+1].startswith('|')
                and set(rows[index+1].replace('|','').strip()) <= {'-',':',' '}):
            columns = [v.strip() for v in line.strip().strip('|').split('|')]
            for column in columns:
                match=re.fullmatch(r'(.+) \((\d+) genes\)',column)
                if match and match[1] in cell_ids:
                    definitions[column.casefold()]='Local joint correlation for '+match[1]+'. The query contains '+match[2]+' selected genes.'
            require(all(key.casefold() in definitions for key in columns),
                    'A supplementary table header has no report definition: '+', '.join(k for k in columns if k.casefold() not in definitions))
            label = 'Table '+filename.replace('.md','').replace('/','-')+'-'+str(len(captions)+1)
            caption = table_caption('Retained scientific evidence',columns,{key:definitions[key.casefold()] for key in columns},
                scope='This table retains the original source scope and scientific decisions.',
                rules=['A missing value does not establish a negative result.'],data='supplements/column-dictionaries.json')
            captions[label] = caption
            result += [START+'\n\n', '**'+label+'. Retained scientific evidence.**\n\n']
            for key in columns: result += ['**'+key+':** '+caption['definitions'][key]+'\n\n']
            prefix='../' if '/' in filename else ''
            result += ['[Column dictionary]('+prefix+'supplements/column-dictionaries.json).\n\n',END+'\n\n']
        result.append(line)
    return ''.join(result), captions


def annotate_supplements(root):
    root = Path(root)
    captions = {}
    from genosapiens.track2_scientific_reader import read_scientific
    cells=list(read_scientific(root/'lincs/results.json')['cells']) if (root/'lincs/results.json').is_file() else []
    # Imported, accepted public reports retain their original bytes and reviews.
    for path in [*root.glob('*.md'), *(root/'lincs').glob('*.md')]:
        if path.name == 'report.md':continue
        relative=str(path.relative_to(root))
        text, local = annotate_tables(path.read_text(), relative, cells)
        path.write_text(text)
        captions.update(local)
        if path.parent == root/'lincs':
            import json
            from genosapiens.track2_scientific_storage import file_identity
            from genosapiens.track2_published_genotype import canonical
            manifest_path=root/'lincs/manifest.json'
            manifest=json.loads(manifest_path.read_bytes())
            manifest['outputs'][relative]=file_identity(path)['sha256']
            manifest_path.write_bytes(canonical(manifest))
    return captions

CSV_MEANINGS = {
'target_group_warnings':'Unresolved applicability of grouped target annotations.',
'canonical_active_ingredient_id':'Canonical active ingredient identifier.',
'ingredient_key':'Stable ingredient identity used for deterministic ordering.',
'medication_entity_key':'Source medication identity after the original resolution rules.',
'gene_id':'Gene identifier in the source namespace. Local query features can use a different namespace from the intervention target.', 'gene_symbol':'Retained gene symbol.',
'scenario_id':'Exact identifier for the separate hypothesis.', 'objective_id':'Exact intervention objective identifier.',
'cell_id':'Cell line used for this experiment.',
'lincs_nominated':'Original nomination flag from the expression route.',
'mechanism_nominated':'Original nomination flag from the target-drug linkage route.',
'lincs_methods':'Actual expression comparison methods for this objective.',
'paired_research_eligible':'Original research eligibility under the paired-expression policy.',
'fallback_recurrence':'Original recurrence decision for the single-direction methods.',
'applicability_conflicts':'Source conflicts that limit applicability to this objective.',
'condition_conflicts':'Retained conflicts among admitted medication conditions.',
'conflict':'Whether retained conditions disagree under the source policy.',
'passing_cells':'Distinct passing cells in the stated source view.',
'tested_cell_count':'Number of cells with assessed local results in this view.',
'selected_cell_count':'Number of cells selected by the original plan.',
'cells':'Complete cell-specific results for this source row.',
'planned_state':'Route state declared by the retained plan.',
'executed':'Whether the route actually ran.',
'reason_codes_result':'Exact source reason codes for the actual route outcome.',
'result_references':'References to retained route result records.',
'mechanism_id':'Exact source identifier for the reported mechanism record.',
'action':'Reported medication action on the stated target.',
'direction_class':'Original classification of reported, assumed, or unresolved action direction.',
'description':'Unchanged source mechanism description.',
'target_name':'Reported target name.', 'target_type':'Source target category.',
'references':'Original source citations and record identifiers.',
'group_warning':'Unresolved applicability of a multi-target annotation.',
'assumption_warning':'Source limitation for an assumed action direction.',
'source':'Retained source identity and provenance.',
'targets':'All target identifiers in the original mechanism annotation.',
'query_sha256':'Checksum of the exact cell-specific query.',
'feature_class':'Original category of the measured gene feature.',
'oe':'Overexpression response value in the source expression units.',
'kd':'Knockdown response value in the source expression units.',
'desired':'Desired query response in the source expression units.',
'opposite':'Opposite query response in the source expression units.',
'desired_direction':'Categorical desired direction for this selected feature.',
'selected':'Whether the original query-selection rule retained this feature.',
'r_desired':'Spearman correlation with the desired query. It has no physical unit.',
'r_opposite':'Spearman correlation with the opposite query. Missing directions remain unavailable.',
'joint':'Minimum of desired correlation and negative opposite correlation. It is unavailable without both comparisons.',
'n_conditions':'Number of admitted medication conditions for this local result.',
'n_undefined':'Number of conditions with undefined comparisons.',
'n_agree':'Number of conditions that pass the actual local method.',
'agreement_fraction':'Fraction of admitted conditions that pass the local method.',
'gene_count':'Number of selected features in this cell-specific query.',
'times_hours':'Admitted measurement times in hours.',
'sig_id':'Exact signature identifier for this condition.',
'pert_id':'Source perturbation identifier.', 'pert_iname':'Source perturbation name.',
'pert_time':'Perturbation duration in the accompanying time unit.',
'pert_time_unit':'Unit of the perturbation duration.',
'pert_dose':'Perturbation dose in the accompanying dose unit.',
'pert_dose_unit':'Unit of the perturbation dose.',
'drug':'Retained medication name for this condition.',
'dose':'Source medication dose in the accompanying dose unit.',
'dose_unit':'Unit of the source medication dose.',
'time_hours':'Measurement time in hours.',
'modz':'Source standardized expression value. It is not fold change or clinical benefit.',
'finite_score_cells':'Number of cells with finite local scores.',
'desired_direction_cells':'Number of cells with positive desired-direction correlation.',
'opposed_direction_cells':'Number of cells with correlation opposed to the desired direction.',
'joint_score_sum':'Sum of finite local joint scores in the original view.',
'joint_priority':'Original paired-cell count multiplied by the sum of finite local joint scores.',
'desired_priority':'Original positive-correlation cell count multiplied by the sum of finite desired correlations.',
'failed_cells':'Cells with failed local comparisons.',
'conflicting_cells':'Cells with retained condition conflicts.',
'missing_cells':'Selected cells without an assessed local result.',
'research_eligible':'Original research eligibility decision, independent of presentation rank.',
'eligibility_bases':'Original scientific criteria supporting research eligibility.',
'tissue_id':'Normalized LINCS primary_site tissue identity.',
'P':'Number of supporting evaluable tissues.', 'C':'Number of conflicting evaluable tissues.',
'N':'Number of evaluable paired tissues.', 'numerator':'Exact supporting count minus conflicting count.',
'denominator':'Exact evaluable tissue count. Zero means the normalized score is unavailable.',
'nominated':'Whether a paired or weaker expression route nominates this medication.',
'provisional':'Whether the planned scientific scope includes unfinished work.',
'declared_tissues':'Complete planned tissue inventory for this objective.',
'supporting_tissues':'Distinct tissues with both supporting directional observations.',
'conflicting_tissues':'Distinct tissues with both conflicting directional observations.',
'evaluable_tissues':'Distinct completed tissues with a ready paired query and a valid condition.',
'method_labels':'Actual query methods. Single-direction methods remain weaker.',
'coverage':'Complete tissue states, counts, witnesses, and unavailable or incomplete reasons.',
'evaluable':'Whether the tissue belongs in the paired denominator.',
'support':'Whether both supporting directions exist within this medication and tissue query.',
'contribution':'Support minus conflict for an evaluable tissue. Missing means unavailable.',
'witnesses':'Exact condition identifiers for each supporting and conflicting direction.',
'condition_ids':'All admitted condition identifiers retained in this summary.',
'n_valid':'Number of conditions with defined comparisons under the selected method.',
'coverage_reasons':'Reasons for unavailable, incomplete, or partly undefined coverage.',
'incomplete_jobs':'Planned scientific jobs without a valid completion.',
'weaker_nomination':'A single-direction nomination outside the paired tissue denominator.',
'legacy_research_eligibility':'The legacy same-cell rule is not applicable to tissue flags.',
'job_id':'Exact dependent scientific job identifier.',
'quality':'Retained source quality metrics.', 'gene_ids':'Frozen measured feature identities in vector order.',
'values':'Individual condition values on the complete frozen feature set.',
'r_oe':'Raw signed Spearman correlation with pooled OE.',
'r_kd':'Raw signed Spearman correlation with pooled KD.',
'query_reference':'Package source for the pooled query and its contributors.',
'condition_reference':'Package source for the exact condition records.',
'source_values':'Selected genetic contributor values before separate directional medians.',
'missing_contributors':'Contributor values that are missing or nonfinite for this feature.',
'undefined_reasons':'Condition identifiers and reasons for undefined or unavailable comparisons.',
}


def export_index(root, captions):
    """Sidecars link immutable scientific exports to the shared dictionaries."""
    import csv
    from genosapiens.track2_scientific_storage import publish_json
    root=Path(root)
    definitions={k.casefold():v for k,v in MEANINGS.items()}
    index=[]
    paths=[*root.glob('*.json'), *(root/'discovery').glob('*.json'), *(root/'discovery').glob('*.csv'),
           *(root/'lincs').rglob('*.csv'), *(root/'supplements').glob('*.json')]
    for path in sorted(paths):
        if path.name in {'manifest.json','column-dictionaries.json','export-index.json'}:continue
        relative=str(path.relative_to(root))
        label='Export '+relative
        if path.suffix=='.csv':
            with path.open() as stream:columns=next(csv.reader(stream),[])
            meanings={key:CSV_MEANINGS.get(key,definitions.get(key.replace('_',' ').casefold())) for key in columns}
            if 'tissue' in path.name or any(part.startswith('tissue-drug-') for part in path.parts):
                meanings.update({k:v for k,v in {
                    'query_sha256':'Checksum of the frozen pooled tissue query.',
                    'oe':'Median OE value from the fixed tissue contributors.',
                    'kd':'Median KD value from the fixed tissue contributors.',
                    'conflict':'Whether both opposite directions exist within this medication and tissue query.',
                    'n_conditions':'Number of admitted conditions retained within this tissue scope.'}.items() if k in columns})
            if (path.parent.name in {'genetic','drugs'} and path.parent.parent.name=='science'
                    and path.name.startswith('block-') and columns[0]=='gene_id'):
                meanings.update({key:'Standardized expression (MODZ) for source signature '+key+
                    '. These values do not establish clinical benefit.' for key in columns[1:]})
            require(all(meanings.values()),'A CSV export column lacks a definition: '+', '.join(k for k,v in meanings.items() if not v))
            captions[label]=table_caption('Complete scientific export',columns,meanings,
                scope='Original scientific records retain their own query and objective scope.',
                rules=['Missing values remain unavailable. These values do not establish treatment benefit.'],data=relative)
        linked=[key for key,value in captions.items() if value['data']==relative]
        if relative=='tables.json':linked=[key for key in captions if key.startswith('Table tables-')]
        if relative=='supplements/objective-scores.json':linked=['combined','lincs','mechanism']
        index.append(dict(path=relative,column_dictionary='supplements/column-dictionaries.json',captions=linked))
    publish_json(root/'supplements/export-index.json',index)
    return captions
