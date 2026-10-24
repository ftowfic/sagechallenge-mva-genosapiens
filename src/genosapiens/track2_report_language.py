"""Shared explanations for the terms used in Track 2 reports."""


def report_acknowledgement():
    """Exact required quotation from the local official challenge rules."""
    return ('## Acknowledgement\n\n'
        'The challenge requires the following wording. It remains verbatim.\n\n'
        '> This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership with the MVA Society, '
        'Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment Consortium for Science), with prize sponsorship from AWS and Anthropic. '
        'We are deeply grateful to the child and their family who generously contributed their data and their story to advance research into this rare disease. '
        'We acknowledge their trust in making this Hackathon possible.')


def report_terms(*, lincs=False, diagnosis=False, exploratory=False, tissue=False):
    terms = [
        ('Assessment status', 'Completion of medication research within the declared scope.'),
        ('Not performed', 'No supported objective reached a medication search or a valid frozen assessment.'),
        ('Incomplete', 'Assessment started, but required work or a selected scenario remains incomplete.'),
        ('Complete', 'Every required objective has complete assessment evidence within the declared scope.'),
        ('Research eligibility', 'Permission to include a medication as a research candidate. It does not establish treatment benefit.'),
        ('Mechanism', 'The biological process that explains a disease effect or a medication effect.'),
        ('Genotype', 'The genetic variants considered together for the reported condition.'),
        ('Functional product', 'A protein that can perform the biological function under study.'),
        ('Expression', 'The measured RNA response. It does not establish functional protein availability.'),
        ('Rescue', 'An improvement in the specified biological defect in the tested model.'),
        ('Exposure', 'The amount of medication that reaches the relevant tissue over time.'),
        ('Endpoint', 'The measured result used to assess an intervention.'),
        ('G0–G7', 'The checks for input, medication identity, mechanism, action, independent evidence, exposure, safety and validation.'),
        ('SHA-256', 'A file hash used to identify the exact file contents.'),
    ]
    if lincs:
        terms += [
            ('LINCS', 'Library of Integrated Network-based Cellular Signatures. This report uses its measured cell responses.'),
            ('Cell line', 'Cells grown under controlled laboratory conditions. Each reported comparison uses one cell line.'),
            ('KD and OE', 'Knockdown (reduced target expression) and overexpression (extra target expression).'),
            ('shRNA', 'Short hairpin RNA, the reagent used for knockdown.'),
            ('MODZ', 'The standardized expression score in the retained LINCS profiles. It is not a fold change or a benefit measure.'),
            ('Local', 'Within one cell line.'),
            ('Condition', 'A recorded combination of cell line, medication, dose and exposure time.'),
            ('Joint score', 'The weaker of the two directional checks for a condition. Medication summaries combine these scores across conditions and times.'),
            ('Spearman correlation', 'A measure of similarity between the order of responses.'),
            ('GTEx and TPM', 'The reference expression resource and its RNA unit, transcripts per million.'),
            ('QTL', 'Quantitative trait locus. A genomic location associated with variation in a measured trait.'),
            ('FDR', 'False discovery rate, a statistical measure for multiple tests.'),
            ('FDA and QC', 'The US Food and Drug Administration and quality control.'),
        ]
    if diagnosis:
        terms += [
            ('Automatic acceptance', 'A machine decision under the named policy. It does not record human review.'),
            ('Prediction accepted', 'A compatible research direction assumption from accepted genotype, quality, transcript, and functional-source checks.'),
            ('Functional catalog', 'A frozen local source checked for applicable statements about protein function.'),
            ('Complete functional lookup', 'The declared source check finished. This does not establish complete literature coverage or a measured effect.'),
            ('Public study import', 'Retained evidence from a public study with its original source and mapping limits.'),
            ('Source limits', 'Recorded gaps in source coverage, identity mapping, or interpretation.'),
            ('Exploratory applicability', 'Permission to evaluate a hypothesis in the stated source context. It does not establish an exact genotype match.'),
            ('Affected copies', 'The number of allele copies required by the applicable disease model.'),
            ('Ploidy', 'The chromosome-copy count reported or expected at a locus.'),
            ('GT, GQ, and DP', 'The sample genotype, genotype quality, and read depth.'),
            ('ALT and allele balance', 'An alternate allele and its depth fraction under the existing sample-quality calculation.'),
            ('MANE', 'Reference transcripts shared by NCBI and Ensembl.'),
            ('SIFT and PolyPhen', 'Existing protein-damage predictors. Their results do not establish molecular mechanism.'),
            ('ExcessHet', 'A retained source heterozygosity diagnostic. This policy adds no threshold.'),
            ('HPO', 'Human Phenotype Ontology. It supplies standard terms for patient features.'),
            ('MVA', 'Mosaic variegated aneuploidy, the disease label used for the synthetic fixture.'),
            ('SNV', 'Single nucleotide variant. A change at one DNA base.'),
            ('Frozen input', 'The retained input with fixed contents for this analysis.'),
            ('Fixture', 'Invented test data used to check software behavior.'),
        ]
    if exploratory:
        terms=[row for row in terms if row[0]!='MVA']
        terms += [
            ('Hypothesis', 'A proposed explanation that requires further testing.'),
            ('Nomination', 'Selection for research review under the stated search rules.'),
            ('Paired analysis', 'A comparison that uses both overexpression and knockdown in the same cell line.'),
            ('KD-only reversal', 'A comparison with the reversed knockdown pattern. Overexpression evidence is missing.'),
            ('OE-only matching', 'A comparison with the overexpression pattern. Knockdown evidence is missing.'),
            ('Route overlap', 'Both searches nominate the same ingredient for the same gene, scenario, and objective.'),
            ('Recurrence', 'Passing results in multiple cell lines for the same ingredient and objective. Each view retains its method.'),
            ('Agreement group', 'A summary of passing cells, failed cells, and condition conflicts within the stated view.'),
            ('Joint priority', 'Passing paired cell count multiplied by the sum of finite local joint scores. Negative scores remain in the sum.'),
            ('Desired priority', 'Positive desired correlation cell count multiplied by the sum of finite local desired correlations.'),
            ('Ingredient key', 'The stable identifier used to order equal results for the same active ingredient.'),
            ('Unavailable score', 'The dash means that a score is unavailable. It does not mean zero or failure.'),
            ('Assumed direction', 'A policy assumption about drug action. The source does not establish that direction.'),
        ]
    if tissue:
        replacements={'Cell line':'Cells grown under controlled laboratory conditions.',
            'Paired analysis':'A comparison with separately pooled OE and KD from one tissue group.',
            'Unavailable score':'NA means no evaluable paired tissue evidence. It does not mean measured zero.'}
        omitted={'Joint score','Local','Recurrence','Agreement group','Joint priority','Desired priority'}
        terms=[(key,replacements.get(key,value)) for key,value in terms if key not in omitted]
    return '\n'.join(['## Terms used in this report', '', '| Term | Meaning |', '|---|---|']
                     + [f'| {term} | {meaning} |' for term, meaning in terms])
