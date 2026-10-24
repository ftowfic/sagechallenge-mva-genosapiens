"""An explicit maintenance supplement for the pinned operator package.

The baseline release and scientific/checkpoint identities stay intact. The new
package binding authenticates the actual code and the supplement records the
exact predecessor and changed modules. This is not a new science qualification.
"""

import hashlib
import json
import os
from functools import lru_cache
from pathlib import Path

from genosapiens.bounded_json import read_bounded_bytes
from genosapiens.wgs_ledger import canonical


PREDECESSOR_SHA256 = '346619d914bbf5841be2aef84af9b675c63bb5ab8011b1403f13d72abf724b24'
FIRST_REPAIR_SHA256 = 'ffaadf6119cec25a8ac50c88834237de6979994bb348883f89d1345862cfaa57'
SECOND_REPAIR_SHA256 = 'b11e5ea553dc2953f906b8a5768c494f84ba8365822fe133d1b85ba474ba0ec8'
THIRD_REPAIR_SHA256 = '80c6e7856d6a9b4a792083c48a8a8dfd3349efcd01b6b39438e1fa5cfbf45b2f'
FOURTH_REPAIR_SHA256 = 'e8adf9b8c287e48d2fe79059a5d6f8506458323552c6b56985ac889cc84c6c46'
FIFTH_REPAIR_SHA256 = '57742c22dbebb03bc8ac471bb52f7d86fca4f26787b181b06fffdc767bacde6d'
SIXTH_REPAIR_SHA256 = '264eb41305d3c1069fa8a9e269af75382c2f61d6525b7ff0a695fc5727fccb5c'
SEVENTH_REPAIR_SHA256 = '82dbafd7f8940acc71e398e665abbddd61168e9bbeb212913587a5e7f39ba520'
EIGHTH_REPAIR_SHA256 = 'a6fd891435be6d2c8ccd1b473db83df2c8b75a9dfa75b3664a86967b4f6df89c'
NINTH_REPAIR_SHA256 = 'd65a242a2bc9b353a30d44a2d507bf1d4c5c48fd0dde5987b739d3b3d611b4ce'
TENTH_REPAIR_SHA256 = 'a475b9e8fe12a902191890b8effffc04e552debe5d0c9052ca5ab384da095563'
ELEVENTH_REPAIR_SHA256 = '50d3a2cada716170020aa03cc1e4829cb13e806a203975039ad75304c4affb72'
TWELFTH_REPAIR_SHA256 = '2f194b7a8db01a09f819e8adcec3b485863a868c5ab5800d90d7181e803c4cab'
THIRTEENTH_REPAIR_SHA256 = '43adff367efc25fe8e4999e910bb6d75f557988f4d25ca145a1bba4a37f3a427'
COMPATIBLE_PREDECESSORS = (PREDECESSOR_SHA256, FIRST_REPAIR_SHA256, SECOND_REPAIR_SHA256,
                         THIRD_REPAIR_SHA256, FOURTH_REPAIR_SHA256, FIFTH_REPAIR_SHA256,
                         SIXTH_REPAIR_SHA256, SEVENTH_REPAIR_SHA256, EIGHTH_REPAIR_SHA256,
                         NINTH_REPAIR_SHA256, TENTH_REPAIR_SHA256, ELEVENTH_REPAIR_SHA256,
                         TWELFTH_REPAIR_SHA256, THIRTEENTH_REPAIR_SHA256)
REPAIR_NAME = 'nonprimary-alias-resume-v14'
MAINTENANCE_NAME = 'operator-maintenance.json'
ORIGINAL_MANIFEST = 'maintenance-original/package-bindings.json'
PRIOR_CHANGED_MODULES = (
    'phenotype_ingestion.py', 'reference_concordance.py',
    'wgs_operator_invocation.py', 'wgs_preparation.py', 'wgs_resume.py',
)
V3_CHANGED_MODULES = (*PRIOR_CHANGED_MODULES, 'wgs_failures.py')
V4_CHANGED_MODULES = (*V3_CHANGED_MODULES, 'submission_rendering.py',
                     'wgs_source_routing.py', 'wgs_filter.py')
V5_CHANGED_MODULES = V4_CHANGED_MODULES
V6_CHANGED_MODULES = (*V5_CHANGED_MODULES, 'wgs_public_analysis.py')
V7_CHANGED_MODULES = (*V6_CHANGED_MODULES, 'wgs_runtime.py', 'bounded_wgs_config.py',
                   'wgs_profile_projection.py', 'wgs_assessment_checkpoint.py',
                   'track1_public_validation_coordinator.py', 'track1_roster.py')
V8_CHANGED_MODULES = V7_CHANGED_MODULES
V9_CHANGED_MODULES = (*V8_CHANGED_MODULES, 'wgs_annotation.py')
V12_CHANGED_MODULES = (*V9_CHANGED_MODULES, 'wgs_ledger.py')
CHANGED_MODULES = (*V12_CHANGED_MODULES, 'wgs_admission.py', 'wgs_admission_projection.py', 'wgs_reporting.py')
CHANGED_SCRIPTS = ('check_bounded_host_capacity.py',)
PRIOR_ADDED_MODULES = ('operator_warnings.py', 'wgs_operator_maintenance.py')
V12_ADDED_MODULES = (*PRIOR_ADDED_MODULES, 'wgs_migration.py')
ADDED_MODULES = (*V12_ADDED_MODULES, 'wgs_annotation_coverage.py')


def sha(content):
    return hashlib.sha256(content).hexdigest()


def maintenance_baseline(files, read):
    """Verify the whole bounded code supplement, without reference/case reads."""
    if MAINTENANCE_NAME not in files:
        return None
    content = read(MAINTENANCE_NAME)
    maintenance = json.loads(content)
    original_content = read(ORIGINAL_MANIFEST)
    original = json.loads(original_content)
    repair = maintenance.get('repair')
    changed = {'python/genosapiens/' + name for name in (
        CHANGED_MODULES if repair in (REPAIR_NAME, 'unassessed-contig-resume-v13') else
        V12_CHANGED_MODULES if repair in ('failed-annotation-resume-v12', 'journal-head-migration-v11', 'output-migration-v10') else
        V9_CHANGED_MODULES if repair == 'interrupted-annotation-resume-v9' else
        V8_CHANGED_MODULES if repair == 'session-query-cache-resume-v8' else
        V7_CHANGED_MODULES if repair == 'remove-disk-capacity-gates-v7' else
        V6_CHANGED_MODULES if repair == 'analysis-checkpoint-package-path-v6' else
        V5_CHANGED_MODULES if repair == 'prepared-manifest-handoff-v5' else
        V4_CHANGED_MODULES if repair == 'source-handoff-and-contig-resume-v4' else
        V3_CHANGED_MODULES if repair == 'source-resume-and-input-notices-v3'
        else PRIOR_CHANGED_MODULES)}
    if repair in (REPAIR_NAME, 'unassessed-contig-resume-v13', 'failed-annotation-resume-v12', 'journal-head-migration-v11', 'output-migration-v10', 'interrupted-annotation-resume-v9', 'session-query-cache-resume-v8', 'remove-disk-capacity-gates-v7'):
        changed.update(CHANGED_SCRIPTS)
    added = {'python/genosapiens/' + name for name in (
        ADDED_MODULES if repair in (REPAIR_NAME, 'unassessed-contig-resume-v13') else
        V12_ADDED_MODULES if repair in ('failed-annotation-resume-v12', 'journal-head-migration-v11', 'output-migration-v10') else PRIOR_ADDED_MODULES)}
    backups = {'maintenance-original/' + name for name in changed}
    baseline = {row['relative_path']: row for row in original['files']}
    expected = dict(schema='genosapiens.operator-maintenance', version=1,
        repair=maintenance.get('repair'), predecessor_package_sha256=PREDECESSOR_SHA256,
        baseline_executable_stack_sha256=sha(read('executable-stack.json')),
        changed_modules=sorted(changed), added_modules=sorted(added))
    if expected['repair'] == REPAIR_NAME:
        expected['compatible_predecessor_package_sha256'] = list(COMPATIBLE_PREDECESSORS)
    elif expected['repair'] == 'unassessed-contig-resume-v13':
        expected['compatible_predecessor_package_sha256'] = list(COMPATIBLE_PREDECESSORS[:13])
    elif expected['repair'] == 'failed-annotation-resume-v12':
        expected['compatible_predecessor_package_sha256'] = list(COMPATIBLE_PREDECESSORS[:12])
    elif expected['repair'] == 'journal-head-migration-v11':
        expected['compatible_predecessor_package_sha256'] = list(COMPATIBLE_PREDECESSORS[:11])
    elif expected['repair'] == 'output-migration-v10':
        expected['compatible_predecessor_package_sha256'] = list(COMPATIBLE_PREDECESSORS[:10])
    elif expected['repair'] == 'interrupted-annotation-resume-v9':
        expected['compatible_predecessor_package_sha256'] = list(COMPATIBLE_PREDECESSORS[:9])
    elif expected['repair'] == 'session-query-cache-resume-v8':
        expected['compatible_predecessor_package_sha256'] = list(COMPATIBLE_PREDECESSORS[:8])
    elif expected['repair'] == 'remove-disk-capacity-gates-v7':
        expected['compatible_predecessor_package_sha256'] = list(COMPATIBLE_PREDECESSORS[:7])
    elif expected['repair'] == 'analysis-checkpoint-package-path-v6':
        expected['compatible_predecessor_package_sha256'] = list(COMPATIBLE_PREDECESSORS[:6])
    elif expected['repair'] == 'prepared-manifest-handoff-v5':
        expected['compatible_predecessor_package_sha256'] = list(COMPATIBLE_PREDECESSORS[:5])
    elif expected['repair'] == 'source-handoff-and-contig-resume-v4':
        expected['compatible_predecessor_package_sha256'] = list(COMPATIBLE_PREDECESSORS[:4])
    elif expected['repair'] == 'source-resume-and-input-notices-v3':
        expected['compatible_predecessor_package_sha256'] = list(COMPATIBLE_PREDECESSORS[:3])
    elif expected['repair'] == 'source-resume-and-input-notices-v2':
        expected['compatible_predecessor_package_sha256'] = list(COMPATIBLE_PREDECESSORS[:2])
    elif expected['repair'] != 'source-resume-and-input-notices-v1':
        raise ValueError('operator maintenance repair is unsupported')
    if (canonical(maintenance) != content or canonical(original) != original_content
            or maintenance != expected
            or sha(original_content) != PREDECESSOR_SHA256
            or set(files) != set(baseline) | added | backups | {MAINTENANCE_NAME, ORIGINAL_MANIFEST}):
        raise ValueError('operator maintenance binding is invalid')
    for name, row in baseline.items():
        if name in changed:
            old_content = read('maintenance-original/' + name)
            if sha(old_content) != row['sha256'] or len(old_content) != row['size_bytes']:
                raise ValueError('operator maintenance predecessor source changed')
            if files[name]['mode'] != row['mode']:
                raise ValueError('operator maintenance source mode changed')
        elif files[name] != row:
            raise ValueError('operator maintenance changed an unsupported file')
    return baseline


def maintenance_profile_adapter_sha256(root):
    """Use the original profile identity only inside an authenticated installation."""
    return _verified_profile_adapter_sha256(Path(root), os.environ.get('GENOSAPIENS_PACKAGE_SHA256'))


@lru_cache(maxsize=8)
def _verified_profile_adapter_sha256(root, installed):
    # Installed code is mounted read-only. Authenticate once per package/process.
    # This preserves checkpoint science identity while the package binds new code.
    from genosapiens.wgs_operator_invocation import _package
    _, _, _, files, read = _package(root, installed)
    baseline = maintenance_baseline(files, read)
    if baseline is None:
        raise ValueError('operator maintenance binding is missing')
    return baseline['python/genosapiens/wgs_profile_projection.py']['sha256']


def compatible_package_binding(root, recorded, installed):
    """Accept only an installed binding or the exact declared predecessor."""
    if recorded == installed:
        return True
    if recorded not in COMPATIBLE_PREDECESSORS:
        return False
    content = read_bounded_bytes(Path(root) / 'package-bindings.json')
    if sha(content) != installed:
        return False
    files = {row['relative_path']: row for row in json.loads(content)['files']}
    def read(name):
        data = read_bounded_bytes(Path(root) / name, reference_root=Path(root),
                                  maximum_bytes=16 * 1024**2)
        if sha(data) != files[name]['sha256'] or len(data) != files[name]['size_bytes']:
            raise ValueError('operator maintenance file changed')
        return data
    if maintenance_baseline(files, read) is None:
        return False
    maintenance = json.loads(read(MAINTENANCE_NAME))
    return recorded in maintenance.get('compatible_predecessor_package_sha256', [PREDECESSOR_SHA256])


def build_operator_maintenance(*, predecessor, destination, repository_root):
    """Create a distinct package. Never modify a predecessor or private output."""
    import stat
    from genosapiens.wgs_operator_invocation import _package
    from genosapiens.wgs_operator_package import _write

    predecessor, destination, repository = map(Path, (predecessor, destination, repository_root))
    predecessor_manifest = read_bounded_bytes(predecessor / 'package-bindings.json')
    predecessor_sha = sha(predecessor_manifest)
    if predecessor_sha not in COMPATIBLE_PREDECESSORS:
        raise ValueError('maintenance predecessor is unsupported')
    _, _, _, installed_files, installed_read = _package(predecessor, predecessor_sha)
    baseline = maintenance_baseline(installed_files, installed_read)
    original_files = installed_files if baseline is None else baseline
    original_manifest = (predecessor_manifest if baseline is None
                         else installed_read(ORIGINAL_MANIFEST))

    def original_read(name):
        backup = 'maintenance-original/' + name
        return installed_read(backup if baseline is not None and backup in installed_files else name)
    if destination.exists():
        raise ValueError('maintenance destination already exists')
    destination.mkdir(mode=0o700)
    changed = {'python/genosapiens/' + name for name in CHANGED_MODULES} | set(CHANGED_SCRIPTS)
    for name, row in original_files.items():
        data = original_read(name)
        if name in changed:
            _write(destination / 'maintenance-original' / name, data, row['mode'])
            source = (repository / 'scripts/operator_compatibility_repair' / name
                      if name in CHANGED_SCRIPTS else repository / 'src/genosapiens' / Path(name).name)
            data = read_bounded_bytes(source,
                                      maximum_bytes=16 * 1024**2)
        _write(destination / name, data, row['mode'])
    for name in ADDED_MODULES:
        _write(destination / 'python/genosapiens' / name,
               read_bounded_bytes(repository / 'src/genosapiens' / name), 0o400)
    _write(destination / ORIGINAL_MANIFEST, original_manifest, 0o400)
    supplement = canonical(dict(schema='genosapiens.operator-maintenance', version=1,
        repair=REPAIR_NAME, predecessor_package_sha256=PREDECESSOR_SHA256,
        compatible_predecessor_package_sha256=list(COMPATIBLE_PREDECESSORS),
        baseline_executable_stack_sha256=sha(original_read('executable-stack.json')),
        changed_modules=sorted(changed),
        added_modules=sorted('python/genosapiens/' + name for name in ADDED_MODULES)))
    _write(destination / MAINTENANCE_NAME, supplement, 0o400)
    manifest = json.loads(original_manifest)
    manifest['files'] = []
    for path in sorted(destination.rglob('*')):
        if path.is_file():
            data = read_bounded_bytes(path, maximum_bytes=16 * 1024**2)
            manifest['files'].append(dict(relative_path=path.relative_to(destination).as_posix(),
                sha256=sha(data), size_bytes=len(data), mode=stat.S_IMODE(path.stat().st_mode)))
    content = canonical(manifest)
    _write(destination / 'package-bindings.json', content, 0o400)
    _package(destination, sha(content))
    return sha(content)
