"""Read-only discovery admission for an exact tested reporting transition."""
from pathlib import Path

from genosapiens.track2_operator_inputs import require

# The reviewed transition adds package inventory, figure-attempt preservation,
# an immutable query snapshot, and a versioned narrative schema. Scoring is unchanged.
# Extend these identities only after fresh transition and snapshot tests.
PRODUCER_AGGREGATION = '4f3d284204b4668026d4ebad3f19b3b144cc499b3b70e20531bb137eded6fa91'
CONSUMER_AGGREGATION = '8493d9671e90f941baabc4df2261a8c9127a8645b7dd72596bbc07305b7bb00d'
# Baseline e2c932ee to the integrated reporter. Only package file inventories
# and the two new report schemas change in the aggregation dependency closure.
WORKFLOW_CONSUMER_AGGREGATION = 'a3b526e5cc2cb9af1ca0aecefd5dcfd866ac35a414b3f2c6ffede7e4056cbebe'
AGGREGATION_TRANSITIONS = frozenset({
    (PRODUCER_AGGREGATION, CONSUMER_AGGREGATION),
    (CONSUMER_AGGREGATION, WORKFLOW_CONSUMER_AGGREGATION),
})
SCIENTIFIC_CODE = 'bf8df3426d1c94eceffae6c26a686411ba1515948c22c619044a5508ad6380ba'
RULE = 'retained-discovery-for-challenge-report-v1'


def admitted_binding(saved, expected):
    """Never relabel retained products with consumer dependency identities."""
    if saved==expected:return saved
    fields={'plan_sha256','resources','policy_sha256','runtime','scientific_code_sha256','aggregation_code_sha256'}
    require(isinstance(saved,dict) and isinstance(expected,dict)
        and set(saved)==set(expected)==fields
        and (saved['aggregation_code_sha256'],expected['aggregation_code_sha256']) in AGGREGATION_TRANSITIONS
        and saved['scientific_code_sha256']==expected['scientific_code_sha256']==SCIENTIFIC_CODE
        and all(saved[key]==expected[key] for key in fields-{'aggregation_code_sha256'}),
        'Unsupported discovery transition for reporting. Preserve the original checkpoint and use a separately tested transition.')
    return saved


def load_report_discovery(path, expected, *, allow_partial=False):
    """Admit and replay old completion evidence without a cleanup or resume path."""
    from genosapiens.track2_discovery import _read,load_completed_discovery,scientific_identity
    saved=_read(Path(path)/'request.json')['bindings']
    binding=admitted_binding(saved,expected)
    if saved!=expected:
        require((Path(path)/'completion.json').is_file(),
            'This reporting transition requires complete discovery. Preserve partial evidence for a separate transition.')
        require(scientific_identity()==SCIENTIFIC_CODE,
            'The current scientific implementation lacks this discovery transition.')
    return load_completed_discovery(path,binding,allow_partial=allow_partial)
