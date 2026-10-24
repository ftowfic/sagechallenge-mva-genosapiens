"""One selection statement for input dossiers and research reports."""


def render_selection(selection, *, deferred_path):
    coverage=selection['coverage']
    retention=coverage['retention']
    lines=['## Allele selection','',
        f"The selector assessed {retention['counts']['scoped']} normalized sequence alleles within the frozen scope.",
        f"Retained normalized alleles: {retention['counts']['retained']}. Deferred normalized alleles: {retention['counts']['deferred']}.",
        f"Source ALT occurrences for retained alleles: {retention['occurrences']['retained']}.",
        f"Source ALT occurrences for deferred alleles: {retention['occurrences']['deferred']}.",
        'The [deferred inventory]('+deferred_path+') contains identities, annotations, reasons, and source versions.',
        'Selector policy (exact identifier): `'+retention['policy_version']+'`.',
        'Selector digest (exact identifier): `'+retention['policy_sha256']+'`.',
        *retention['limitations'],
        'The source genome can contain a second disease allele outside the retained domain.',
        'Assessment completion and submission readiness require their separate evidence checks.','',
        'A normalized allele identifies one sequence change after reference normalization.',
        'A source ALT occurrence identifies one alternate entry in one input record.',
        'Deferral means that the allele has no supported retention route under the current policy.','',
        '### Source accounting','',
        'Each count below uses source ALT occurrences. The categories do not overlap.','',
        '| Source selection outcome | Source ALT occurrences |','|---|---|']
    labels={'quality_excluded':'FILTER quality exclusion','not_carried':'Explicit non-carrier',
        'outside_scope':'Outside the frozen scope','structural_review':'Structural review',
        'unsupported_class':'Unsupported class','unresolved_structural':'Unresolved structural scope',
        'no_sequence_change':'No sequence change'}
    lines += [f"| {label} | {coverage['source_accounting'][key]} |" for key,label in labels.items()]
    lines += [f"| Retained sequence allele occurrences | {retention['occurrences']['retained']} |",
        f"| Deferred sequence allele occurrences | {retention['occurrences']['deferred']} |",
        f"| Total source ALT occurrences | {sum(coverage['source_alternate_class_counts'].values())} |",'',
        'Structural review does not include structural genotype, mechanism, or medication assessment.',
        'Source class counts describe the input. They do not add independent evidence.']
    return '\n'.join(lines)+'\n'
