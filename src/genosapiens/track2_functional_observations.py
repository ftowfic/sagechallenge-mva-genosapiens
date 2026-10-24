"""Declared observation categories, without measured molecular-direction claims."""
import re

POLICY = 'track2-functional-observation-interpretation/1.0.0'
STATEMENTS = {
    'loss of protein expression': ('protein_expression_loss', 'source_reported_protein_expression'),
    'impairment of spindle assembly checkpoint': ('assay_impairment', 'spindle_assembly_checkpoint'),
}
COMPATIBLE_WITH_LOSS = frozenset({'loss_of_function', 'protein_expression_loss',
                                'assay_impairment', 'compatible_functional_observations'})


def observations(description):
    text = description.lower().replace('–', '-')
    if not any(phrase in text for phrase in STATEMENTS):
        return [], None
    if re.search(r'\b(?:not|no|without|may|might|possibly|putative|potential|uncertain|when|'
                 r'together|combination|combined|unless|if|cis|trans)\b', text):
        return [], 'effect_scope_or_negation_unresolved'
    result = []
    for clause in text.split(';'):
        clause = clause.strip().rstrip('.')
        if clause in STATEMENTS:
            category, scope = STATEMENTS[clause]
            result.append(dict(category=category, scope=scope, statement=clause))
        elif re.fullmatch(r'in [a-z0-9 -]+|dbsnp:rs\d+', clause):
            continue
        else:
            return [], 'functional_assertion_requires_interpretation'
    return result, None


def interpretation_record(description, citations):
    values, issue = observations(description)
    if issue or not values:
        return None
    return dict(policy=POLICY, authority='agent_interpretation', human_review=False,
        original_statement=description, citations=list(citations), observations=values,
        compatible_predicted_direction='loss_of_function', measured_molecular_direction='not_established',
        tissue_generalization='not_established', treatment_response='not_established')
