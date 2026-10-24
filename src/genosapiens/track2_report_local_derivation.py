"""Private current-derivation handoff, separate from serialized rankings."""
from contextlib import contextmanager
import copy

from genosapiens.track2_disk_memo import DiskMemo
from genosapiens.track2_operator_inputs import require
from genosapiens.track2_scientific_storage import json_identity

# These are exactly the aggregate fields checked by _check_observation.
_AGGREGATE_FIELDS=('status','strength','n_conditions','n_valid','n_undefined','n_cells','n_times',
                   'contributing_tissues','missing_tissues','tissue_strengths')
_PRODUCER=object()


class _CheckedDerivation:
    """Hold independent expected summaries from current checked derivation.

    Saved ranking rows and arbitrary derived_result arguments cannot create
    this authority. The producer closes the disk cache after its consumers.
    """
    def __init__(self, authority, evidence, calibration, result, *, validators=None, scientific_policy=None):
        require(authority is _PRODUCER, 'Summary seeds require the current derivation producer.')
        self.validators=validators
        self.scientific_policy=scientific_policy
        self.evidence=json_identity(evidence)
        self.policy=evidence['bindings']['policy']
        self.resource=calibration['expected_resource_sha256']
        self.result=json_identity(result)
        self.closed=False
        self.summaries=DiskMemo()
        try:
            from genosapiens.track2_report_calibration import aggregate_key
            for row in result['complete_export'] if result is not None else ():
                if not row.get('condition_scope_sha256'):continue
                key=aggregate_key(row['condition_scope_sha256'],row)
                summary=copy.deepcopy({field:row[field] for field in _AGGREGATE_FIELDS})
                require(key not in self.summaries or self.summaries[key]==summary,
                        'Checked derivation contains conflicting summaries.')
                self.summaries[key]=summary
        except BaseException:
            self.close()
            raise

    def seed(self, context, *, evidence, calibration, derived):
        """Seed only matching evidence, policies, validators, resources, and results.

        Each disk value owns its aggregate fields. Observation edits cannot
        change expected summaries. Rows without scope digests use normal replay.
        """
        require(not self.closed, 'The checked derivation owner has ended.')
        context._check_active()
        from genosapiens.track2_derived_strength import policy
        require(self.validators==context._validators and self.scientific_policy==json_identity(policy()),
                'The checked derivation scientific dependencies changed.')
        require(self.evidence==json_identity(evidence) and self.policy==evidence['bindings']['policy']
            and self.resource==calibration['expected_resource_sha256'] and self.result==json_identity(derived),
            'The checked derivation handoff differs from current report evidence.')
        from genosapiens.track2_report_calibration import aggregate_key
        for row in derived['complete_export'] if derived is not None else ():
            if not row.get('condition_scope_sha256'):continue
            key=aggregate_key(row['condition_scope_sha256'],row)
            summary=self.summaries[key]
            existing=context._calibration_aggregates.get(key)
            require(existing is None or all(existing[field]==summary[field] for field in _AGGREGATE_FIELDS),
                    'The calibration cache contains conflicting summaries.')
            context._calibration_aggregates[key]=summary

    def close(self):
        if not self.closed:
            self.closed=True
            self.summaries.close()


@contextmanager
def _checked_derivation(root, evidence, calibration):
    """Produce current scores and a temporary checked handoff for composition.

    This function runs derivation once and checks retained source identities.
    A resumed ranking checkpoint uses derived_result without this handoff.
    Standalone readers keep their ordinary semantic replay and disk cache.
    """
    from genosapiens.track2_report_package import derive_evidence
    from genosapiens.track2_report_contracts import check_evidence_contract
    from genosapiens.track2_checked_evidence import validator_identity
    from genosapiens.track2_derived_strength import policy
    validators=validator_identity()
    scientific_policy=json_identity(policy())
    result=derive_evidence(root,evidence,calibration)
    try:
        if result is None:
            yield None,None
            return
        check_evidence_contract(evidence,root=root)
        require(validator_identity()==validators and json_identity(policy())==scientific_policy,
                'Scientific dependencies changed during derivation.')
        handoff=_CheckedDerivation(_PRODUCER,evidence,calibration,result,
            validators=validators,scientific_policy=scientific_policy)
        try:yield result,handoff
        finally:handoff.close()
    finally:
        if hasattr(result,'document'):result.document.close()
