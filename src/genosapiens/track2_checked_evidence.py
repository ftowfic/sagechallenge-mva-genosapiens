"""Invocation-owned semantic results. No persistent success receipt authorizes reuse."""
import copy
from pathlib import Path

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_scientific_storage import file_identity, json_identity


def validator_identity():
    """Bind package-local validators, policies, schemas, and rendering checks."""
    root=Path(__file__).resolve().parents[2]
    files={}
    for directory,patterns in (('src/genosapiens',('*.py',)),('config',('*.json',)),
                               ('schemas',('*.json',)),('scripts',('*.py','*.R'))):
        for pattern in patterns:
            for path in (root/directory).rglob(pattern):
                files[path.relative_to(root).as_posix()]=file_identity(path)['sha256']
    return json_identity(files)['sha256']


class CheckedEvidenceContext:
    """Own current integrity checks and semantic results for one invocation."""
    def __init__(self, *, counters=None):
        from genosapiens.track2_checked_pages import CheckedPageStreams
        self.pages=CheckedPageStreams()
        self._validators=validator_identity()
        from genosapiens.track2_disk_memo import DiskMemo
        self._reports={};self._derived={};self._calibration={};self._calibration_aggregates=DiskMemo();self._closed=False
        self.counters={} if counters is None else counters

    def __enter__(self):
        self._check_active()
        return self

    def __exit__(self, *error):
        self.close()

    def close(self):
        self._closed=True;self._reports.clear();self._derived.clear();self.pages.close()
        self._calibration.clear();self._calibration_aggregates.close()

    def _check_active(self):
        require(not self._closed, 'The checked evidence invocation has ended.')
        require(validator_identity()==self._validators, 'Validator dependencies changed during the invocation.')

    def checked_report(self, root):
        from genosapiens.track2_challenge_revision import inventory
        from genosapiens.track2_report_metadata import checked_report_metadata
        self._check_active()
        before=inventory(root);key=json_identity(before)['sha256']
        if key not in self._reports:
            value=checked_report_metadata(root,context=self)
            require(inventory(root)==before, 'Scientific evidence changed during admission.')
            self.counters['report_semantic_passes']=self.counters.get('report_semantic_passes',0)+1
            self._reports[key]=copy.deepcopy(value)
        return copy.deepcopy(self._reports[key])

    def derived_report(self, root, legacy):
        from genosapiens.track2_challenge_revision import inventory, _derived
        from genosapiens.track2_scientific_reader import read_scientific
        self._check_active()
        before=inventory(root)
        if json_identity(before)['sha256'] not in self._reports:
            self.checked_report(root)
        require(json_identity(legacy)==json_identity(read_scientific(Path(root)/'presentation.json')),
                'Derived scoring received another source presentation.')
        key=json_identity(dict(files=before,presentation=json_identity(legacy)['sha256']))['sha256']
        if key not in self._derived:
            value=_derived(root,legacy)
            require(inventory(root)==before, 'Derived evidence changed during calculation.')
            self.counters['derived_semantic_passes']=self.counters.get('derived_semantic_passes',0)+1
            self._derived[key]=value
            if value is not None:
                from genosapiens.track2_report_calibration import aggregate_key
                for row in value['complete_export']:
                    if row.get('condition_scope_sha256'):
                        self._calibration_aggregates[aggregate_key(row['condition_scope_sha256'],row)]=copy.deepcopy(row)
        value=self._derived[key]
        return value if hasattr(value,'document') else copy.deepcopy(value)

    def checked_calibration(self, root, *, expected_manifest_sha256):
        from genosapiens.track2_challenge_revision import inventory
        from genosapiens.track2_report_calibration import _check_calibration
        self._check_active()
        before=inventory(root)
        require(before.get('manifest.json')==expected_manifest_sha256, 'The calibration manifest identity differs.')
        key=json_identity(before)['sha256']
        if key not in self._calibration:
            value=_check_calibration(root,expected_manifest_sha256=expected_manifest_sha256,
                aggregates=self._calibration_aggregates)
            require(inventory(root)==before, 'Calibration evidence changed during admission.')
            self.counters['calibration_semantic_passes']=self.counters.get('calibration_semantic_passes',0)+1
            self._calibration[key]=copy.deepcopy(value)
        return copy.deepcopy(self._calibration[key])
