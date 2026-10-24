"""Reuse exact annotation requests independently of downstream selection facts."""
from pathlib import Path

from genosapiens.bounded_json import read_bounded_bytes
from genosapiens.evidence_serialization import _plain
from genosapiens.track2_operator_inputs import OperatorError, require, small_json, write_new
from genosapiens.track2_published_genotype import canonical, sha256
from genosapiens.track2_annotation_artifacts import AlleleArtifact, artifact_digest

SCHEMA = 'genosapiens.track2-annotation-checkpoint/v1'


def dependency_binding(requests, config, runtime, repository, *, explicit_empty_ast=True):
    from genosapiens.track2_construction_dependencies import python_dependencies
    repository = Path(repository).resolve()
    # Copied package paths are locations, not scientific version identities.
    # Bind each small contract's content. Resource roots and manifest digests
    # remain explicit. No installed data file is rehashed here.
    def version(value, key=''):
        if isinstance(value, dict):
            return {k: version(v, k) for k, v in value.items()}
        if isinstance(value, (tuple, list)):
            return [version(v) for v in value]
        if key in {'contract_path', 'identity_path'}:
            return {'sha256': sha256(read_bounded_bytes(value, maximum_bytes=16*1024*1024))}
        return value
    return {'requests_sha256': sha256(canonical(_plain(requests))),
        'resources': version(config), 'image': runtime['image'],
        'code': python_dependencies(repository, (
            ('genosapiens.vep_annotation', 'run_vep_annotation'),
            ('genosapiens.track2_scoped_interpretation', 'selected_alleles')),
            explicit_empty_ast=explicit_empty_ast)}


def artifacts(root):
    root = Path(root)
    require(root.is_dir() and root.resolve() == root.absolute(),
        'Annotation checkpoint directory is missing or linked.')
    paths = list(root.iterdir())
    require(paths and all(p.is_file() and not p.is_symlink() for p in paths),
        'Annotation checkpoint contains an invalid artifact.')
    return {p.name: (AlleleArtifact.capture(p) if p.name == 'alleles.json'
        else read_bounded_bytes(p, maximum_bytes=64*1024*1024)) for p in paths}


def load_or_run(root, *, binding, produce, recover=None):
    """Caller owns the run lock. Publish completion after complete output bytes.

    The caller replays annotation admission with the exact current requests.
    A matching checkpoint avoids tool execution, never scientific validation.
    """
    root = Path(root)
    require(root.resolve() == root.absolute(), 'Annotation checkpoint path is linked.')
    completion = root.parent/'annotation-completion.json'
    if completion.exists():
        try:
            record = small_json(completion)
        except OperatorError:
            record = {}
        if (set(record) == {'schema', 'binding', 'artifacts'} and record['schema'] == SCHEMA
                and record['binding'] == binding and root.is_dir()):
            current = artifacts(root)
            if {name: artifact_digest(raw) for name, raw in current.items()} == record['artifacts']:
                print('Reuse completed annotation for identical supplied calls and resources.', flush=True)
                return current
    if not completion.exists() and root.exists():
        recovery = root.parent/'annotation-recovery.json'
        require(recover is not None and recovery.is_file(),
            'Uncompleted annotation was preserved. Resume with a verified recovery record.')
        saved = small_json(recovery)
        require(saved.get('schema') == 'genosapiens.track2-annotation-recovery/v1'
            and saved.get('producer_binding') == {k:v for k,v in binding.items() if k != 'requests_sha256'},
            'Annotation recovery dependencies changed. Original annotation was preserved.')
        current = artifacts(root)
        require(sha256(current['binding.json']) == saved['binding_sha256'],
            'Annotation recovery publication changed. Original annotation was preserved.')
        recover(current)
        write_new(completion, {'schema': SCHEMA, 'binding': binding,
            'artifacts': {name: artifact_digest(raw) for name, raw in current.items()}})
        print('Recovered completed annotation after dependency and admission checks.', flush=True)
        return current
    if root.exists():
        require(root.is_dir() and not any(p.is_symlink() for p in root.rglob('*')),
            'Annotation checkpoint contains a symbolic link.')
        from uuid import uuid4
        preserved = root.parent/('preserved-'+uuid4().hex)
        preserved.mkdir(mode=0o700)
        root.rename(preserved/'annotation')
        for name in ('annotation-completion.json', 'annotation-recovery.json'):
            path = root.parent/name
            if path.exists(): path.rename(preserved/name)
    completion.unlink(missing_ok=True)
    print('Run annotation: its completion or dependencies changed.', flush=True)
    produce(root)
    current = artifacts(root)
    require('binding.json' in current and 'alleles.json' in current,
        'Annotation output lacks required scientific records.')
    write_new(completion, {'schema': SCHEMA, 'binding': binding,
        'artifacts': {name: artifact_digest(raw) for name, raw in current.items()}})
    return current
