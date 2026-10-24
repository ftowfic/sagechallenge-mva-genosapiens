"""Shared tissue report definitions and interpretation limits."""
COLUMNS=['Rank','Medication','Hypothesis','S','S_L','S_M','E','P','C','N','M','Method','Coverage']
DEFINITIONS={
    'Rank':'Competition rank uses exact fractions. Equal scores share rank.',
    'Medication':'Retained medication name and canonical active ingredient identity.',
    'Hypothesis':'Gene and objective representative, linked to its context.',
    'S':'Combined descriptive score. S = S_L + S_M.',
    'S_L':'Normalized expression score. S_L = (P - C) / N. NA means no evaluable tissue.',
    'S_M':'Mechanism score. S_M = I_M + 2M.',
    'E':'Number of nominating evidence routes, from zero to two.',
    'P':'Number of evaluable tissues with supporting directional evidence.',
    'C':'Number of evaluable tissues with conflicting directional evidence.',
    'N':'Number of tissues with a ready paired query and a valid comparison.',
    'M':'Number of applicable reported medication-target-action links.',
    'Method':'Actual genetic query methods. Single-direction methods remain weaker.',
    'Coverage':'Evaluable and declared tissue counts. Provisional marks unfinished scientific work.',
}
RULES=[
    'Support requires positive OE correlation and negative KD correlation within one medication and tissue query.',
    'Conflict requires negative OE correlation and positive KD correlation within the same scope.',
    'The directional witnesses can use different cells, doses, or times.',
    'Each witness condition needs both finite correlations. Exact zero supplies neither sign.',
    'Each tissue contributes at most one count to P, C, and N.',
    'Mixed evidence increments P and C. Its net contribution is zero.',
    'A valid tissue with neither flag also contributes zero. It retains a different evidence state.',
    'Unavailable LINCS evidence contributes zero only to the combined score. The LINCS component still displays NA.',
    'I_M equals one when M exceeds zero. Assumed actions and unresolved target groups earn no mechanism points.',
    'A supporting tissue can nominate a medication even when the net score is zero or negative.',
    'Each table selects at most five ingredients. Exact ties use name, ingredient, and objective for display order.',
    'Scores do not change scientific eligibility. Expression evidence does not establish treatment benefit.',
]


def glossary():
    return ('\n| Tissue group | Cells with the same resolved LINCS primary_site. |\n'
        '| Pooled signature | Per-feature median of selected genetic profiles within one tissue and direction. |\n'
        '| Contributor | One selected genetic signature from one cell and direction. |\n'
        '| Witness | A retained condition that supplies one directional observation. |\n'
        '| Evaluable tissue | A completed tissue with a ready paired query and at least one valid drug comparison. |\n'
        '| Mixed evidence | Supporting and conflicting directional evidence within one tissue. |\n'
        '| Provisional score | A score whose declared execution scope includes unfinished work. |\n'
        '| rOE | Spearman correlation with pooled OE. |\n'
        '| rKD | Raw Spearman correlation with pooled KD. |\n'
        '| -rKD | Negated KD correlation. Positive values indicate KD reversal. |\n')
