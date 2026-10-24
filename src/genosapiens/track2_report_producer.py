"""Explicit report-only producer admission, independent of run resume rules."""
import hashlib
import json
from pathlib import Path

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_scientific_storage import file_identity

PREDECESSOR_NATIVE = '57429626bc1bb4ba3d344273e62ebe346f3925731712aca40af4f1697c66f37e'
# Baseline e2c932ee. The transition tests use reports made by that source tree.
PRECHANGE_NATIVE = '2f94e6d7dcadd5d9925d16cd6160b6e0e2836fefac53eabb45c731cdaec4f88f'
PRECHANGE_STUDY = '26f039ab35843658f3ac18769f51cd09526c016066be9a2962fb0053aa95e1fa'
TRANSITIONS = {
    (True, PREDECESSOR_NATIVE): 'retained-native-producer-57429626-report-only-v1',
    (True, PRECHANGE_NATIVE): 'retained-native-producer-e2c932ee-report-only-v1',
    (False, PRECHANGE_STUDY): 'retained-study-producer-e2c932ee-report-only-v1',
}


def code_identity(code):
    require(isinstance(code,dict) and bool(code),'The report producer lacks dependency evidence.')
    return hashlib.sha256(json.dumps(code,sort_keys=True,separators=(',',':'),ensure_ascii=True).encode()).hexdigest()


def current_code(*, native):
    from genosapiens.track2_construction_dependencies import python_dependencies
    from genosapiens.track2_report_finalization import VALIDATORS
    package=Path(__file__).resolve().parents[2]
    return python_dependencies(package,[
        ('genosapiens.track2_discovery_report','export_native_evidence') if native else
        ('genosapiens.track2_synthesis','export_study_evidence'),
        ('genosapiens.track2_report_fresh','run_fresh_report')],ignored=VALIDATORS)


def admit_producer(record, *, expected_production_sha256, expected_manifest_sha256, native):
    """Identity alone is insufficient: publication must also replay source semantics."""
    record=Path(record)
    require(file_identity(record)['sha256']==expected_production_sha256,
        'The production record differs from its independently expected identity.')
    value=json.loads(record.read_bytes())
    require(set(value)=={'schema','binding','manifest_sha256','result_sha256'}
        and value['schema']=='genosapiens.track2-report-production/v1'
        and value['manifest_sha256']==expected_manifest_sha256,
        'The production record does not bind the expected source report.')
    identity=code_identity(value['binding'].get('code'))
    if (native,identity) in TRANSITIONS:
        rule=TRANSITIONS[native,identity]
    else:
        require(identity==code_identity(current_code(native=native)),
            'Unsupported report producer. Preserve the original package, report, checkpoints, and failure evidence.')
        rule='current-report-producer-v1'
    return dict(schema='genosapiens.track2-report-producer-admission/v1',rule=rule,
        production_sha256=expected_production_sha256,source_manifest_sha256=expected_manifest_sha256,
        producer_code_binding_sha256=identity,scope='separate-report-revision-only',
        requires_current_scientific_admission=True,production=value)


def check_retained_producer(value, expected_manifest_sha256):
    require(value['schema']=='genosapiens.track2-report-producer-admission/v1'
        and value['rule'] in {*TRANSITIONS.values(),'current-report-producer-v1'}
        and value['source_manifest_sha256']==expected_manifest_sha256
        and value['production']['manifest_sha256']==expected_manifest_sha256
        and value['producer_code_binding_sha256']==code_identity(value['production']['binding']['code'])
        and value['scope']=='separate-report-revision-only'
        and value['requires_current_scientific_admission'] is True,
        'The retained producer admission differs from its source report.')
    if value['rule']!='current-report-producer-v1':
        require(any(rule==value['rule'] and identity==value['producer_code_binding_sha256']
            for (_,identity),rule in TRANSITIONS.items()),
            'The retained producer rule has another original code identity.')
