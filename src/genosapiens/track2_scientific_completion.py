"""Compact completions over complete scientific artifacts and streamed indexes."""
import itertools
import json
import os
from pathlib import Path
import stat
import shutil

from genosapiens.track2_operator_inputs import require, small_json
from genosapiens.track2_scientific_storage import file_identity, json_chunks, json_identity, publish_chunks, ScientificCapacityError

SCHEMA = 'genosapiens.track2-scientific-completion/v1'
REFERENCE_SCHEMA = 'genosapiens.track2-artifact-reference/v1'
INDEX_SCHEMA = 'genosapiens.track2-artifact-index/v1'
CONTROL_BYTES = 2 * 1024**2
INDEX = 'artifact-index.jsonl'
PARTS = 'artifact-index-parts'
COMPLETIONS = frozenset({'completion.json', 'partial-completion.json', 'manifest.json'})
EXCLUDED = COMPLETIONS | frozenset({INDEX, PARTS, 'current-diagnostics.json', 'current-outcome.json', 'current-execution.json', '.discovery.lock'})


def artifact_path(root, relative):
    root = Path(root).absolute()
    relative = Path(relative)
    require(not relative.is_absolute() and relative.parts and '..' not in relative.parts,
        'An artifact reference has an invalid path.')
    target = root/relative
    require(not any(p.is_symlink() for p in (target, *target.parents)),
        'An artifact reference contains a symbolic link.')
    return target


def reference(root, relative, *, role, schema, dependencies, identity=None):
    target = artifact_path(root, relative)
    return dict(reference_schema=REFERENCE_SCHEMA, path=Path(relative).as_posix(), role=role,
        schema=schema, dependencies=dependencies, **(identity or file_identity(target)))


def check_reference(root, record, *, role, schema, dependencies):
    require(set(record) == {'reference_schema', 'path', 'role', 'schema', 'dependencies', 'size', 'sha256'}
        and record['reference_schema'] == REFERENCE_SCHEMA and record['role'] == role
        and record['schema'] == schema and record['dependencies'] == dependencies,
        'A scientific artifact reference has incompatible bindings.')
    target = artifact_path(root, record['path'])
    require(type(record['size']) is int and record['size'] >= 0,
        'A scientific artifact length is invalid.')
    require(file_identity(target) == {k: record[k] for k in ('size', 'sha256')},
        'A referenced scientific artifact is missing or changed.')
    return target


def artifact_files(root, *, pending=None, stage=None):
    """Visit one directory at a time, retaining no complete artifact inventory."""
    root = Path(root).absolute()
    require(root.resolve() == root, 'A scientific output contains a symbolic link.')
    def visit(folder):
        with os.scandir(folder) as entries:
            names = sorted(entry.name for entry in entries)
        for name in names:
            path = folder/name
            if path == pending:
                continue
            info = path.lstat()
            require(not stat.S_ISLNK(info.st_mode), 'A scientific output contains a symbolic link.')
            if folder == root and name in EXCLUDED:
                continue
            if folder == root and stage == 'assessment' and name == 'plan-review.json':
                continue
            if name.startswith('.scientific-'):
                require(False, 'An incomplete scientific publication requires recovery.')
            if stat.S_ISDIR(info.st_mode):
                yield from visit(path)
            else:
                require(stat.S_ISREG(info.st_mode), 'A scientific output has an invalid file type.')
                yield path.relative_to(root).as_posix()
    yield from visit(root)


def artifact_contract(name, *, root=None, binding=None, stage=None):
    """The index declares each file's format. Scientific readers check its meaning."""
    filename = Path(name).name
    roles = {
        'selection.json': ('job-selection', 'genosapiens.track2-lincs-selection/v1'),
        'query.json': ('genetic-query', 'genosapiens.track2-lincs-query/v2'),
        'result.json': ('job-result', 'genosapiens.track2-lincs-result/v2'),
        'request.json': ('discovery-request', 'genosapiens.track2-discovery-request/v1'),
        'discovery-results.json': ('discovery-result', 'genosapiens.track2-discovery-results/v1'),
        'mechanism-results.json': ('mechanism-results', 'genosapiens.track2-mechanism-results/v1'),
        'conditions.jsonl': ('scored-conditions', 'genosapiens.track2-scored-conditions/v1'),
        'figure-data.jsonl': ('figure-pages', 'genosapiens.track2-native-figure-page/v1'),
    }
    tissue_query=stage=='tissue-query' or any(part.startswith('tissue-query-') for part in Path(name).parts)
    tissue_drug=stage=='tissue-drug' or any(part.startswith('tissue-drug-') for part in Path(name).parts)
    if tissue_query or tissue_drug:
        roles.update({'query.json':('genetic-query','genosapiens.track2-tissue-query/v1'),
            'result.json':('job-result','genosapiens.track2-tissue-query-job-result/v1' if tissue_query else 'genosapiens.track2-tissue-drug-result/v1'),
            'conditions.jsonl':('scored-conditions','genosapiens.track2-tissue-condition/v1'),
            'selection.json':('job-selection','genosapiens.track2-tissue-query-selection/v1' if tissue_query else 'genosapiens.track2-tissue-drug-job/v1')})
    resources=(binding or {}).get('resources',{})
    if isinstance(resources,dict) and resources.get('query_policy_sha256'):
        roles['discovery-results.json']=('discovery-result','genosapiens.track2-discovery-results/v2')
    return roles.get(filename, ('stage-output', 'file-format:' + Path(name).suffix.lstrip('.')))


def publish_completion(root, *, stage, binding, result, filename='completion.json'):
    """Publish the index first and its small completion last under the output lock."""
    from genosapiens.track2_run_state import atomic_record
    root = Path(root)
    require(filename in COMPLETIONS, 'Invalid scientific completion path.')
    dependencies = json_identity(binding)['sha256']
    count = 0
    # Every index node fits the control bound. Parent nodes reference child
    # segments, so even an arbitrarily large inventory has a compact root.
    parts = root/PARTS
    if parts.exists():
        list(artifact_files(parts))  # Reject links before removing owned segments.
        shutil.rmtree(parts)
    serial = 0
    def segment(chunks):
        nonlocal serial
        name = PARTS+'/'+str(serial).zfill(8)+'.jsonl';serial += 1
        identity = publish_chunks(root/name,chunks)
        return reference(root,name,role='artifact-index-part',schema=INDEX_SCHEMA,
            dependencies=dependencies,identity=identity)
    def rows():
        nonlocal count
        for name in artifact_files(root, stage=stage):
            role, schema = artifact_contract(name,root=root,binding=binding,stage=stage)
            record = reference(root, name, role=role, schema=schema, dependencies=dependencies)
            count += 1
            yield record
    def nodes(records):
        chunks=[];size=0
        for row in records:
            encoded=b''.join(json_chunks(row))
            if len(encoded)>CONTROL_BYTES:
                raise ScientificCapacityError('An artifact reference exceeds its control bound.')
            if size+len(encoded)>CONTROL_BYTES:
                yield segment(chunks)
                chunks=[];size=0
            chunks.append(encoded);size+=len(encoded)
        if chunks:yield segment(chunks)
    # One level at a time keeps memory bounded independently of inventory size.
    # The initial segment is promoted to the root when it has no siblings.
    records=rows()
    while True:
        level=iter(nodes(records))
        first=next(level,None);second=next(level,None)
        require(first is not None,'A scientific completion requires artifact evidence.')
        if second is None:
            source=root/first['path']
            def chunks():
                with source.open('rb') as stream:
                    while chunk:=stream.read(64*1024):yield chunk
            identity=publish_chunks(root/INDEX,chunks())
            source.unlink()
            break
        records=itertools.chain((first,second),level)
    require(count > 0, 'A scientific completion requires artifact evidence.')
    record = dict(schema=SCHEMA, stage=stage, binding=binding, result=result, artifact_count=count,
        inventory=reference(root, INDEX, role='artifact-index', schema=INDEX_SCHEMA,
            dependencies=dependencies, identity=identity))
    # Keep admission independent of scientific artifact size.
    if json_identity(record)['size'] > CONTROL_BYTES:
        raise ScientificCapacityError('Scientific completion exceeds its control bound.')
    atomic_record(root/filename, record)
    return record


def check_completion(root, *, stage, binding, filename='completion.json'):
    root = Path(root)
    require(filename in COMPLETIONS, 'Invalid scientific completion path.')
    record = small_json(root/filename, maximum=CONTROL_BYTES)
    require(set(record) == {'schema', 'stage', 'binding', 'result', 'artifact_count', 'inventory'}
        and record['schema'] == SCHEMA and record['stage'] == stage and record['binding'] == binding,
        'The scientific completion contract or dependencies changed.')
    dependencies = json_identity(binding)['sha256']
    index = check_reference(root, record['inventory'], role='artifact-index', schema=INDEX_SCHEMA,
        dependencies=dependencies)
    count = 0
    segments=set()
    def rows(path):
        require(path.stat().st_size<=CONTROL_BYTES,'An artifact index node exceeds its control bound.')
        with path.open('rb') as stream:
            while line := stream.readline(CONTROL_BYTES + 1):
                require(len(line) <= CONTROL_BYTES and line.endswith(b'\n'), 'An artifact index record is incomplete.')
                row=json.loads(line)
                if row.get('role')=='artifact-index-part':
                    relative=Path(row['path'])
                    require(relative.parent==Path(PARTS) and relative.as_posix() not in segments,
                        'An artifact index segment is duplicated or misplaced.')
                    segments.add(relative.as_posix())
                    child=check_reference(root,row,role='artifact-index-part',schema=INDEX_SCHEMA,dependencies=dependencies)
                    yield from rows(child)
                else:yield row
    for name, row in itertools.zip_longest(artifact_files(root,stage=stage), rows(index)):
        require(name is not None and row is not None and row['path'] == name,
            'The scientific artifact inventory has missing or extra files.')
        role, schema = artifact_contract(name,root=root,binding=binding,stage=stage)
        check_reference(root, row, role=role, schema=schema, dependencies=dependencies)
        count += 1
    actual={PARTS+'/'+name for name in artifact_files(root/PARTS)} if (root/PARTS).exists() else set()
    require(actual==segments,'The artifact index has missing or extra segments.')
    require(type(record['artifact_count']) is int and record['artifact_count'] == count and count > 0,
        'The scientific artifact count changed.')
    return record
