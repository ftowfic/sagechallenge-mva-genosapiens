"""Explicit public source preparation with truthful receipts and bounded buffers."""
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import tempfile
import sys
import urllib.request
import zlib

CHUNK = 1024*1024


def canonical(value):
    return (json.dumps(value, sort_keys=True, separators=(',', ':'))+'\n').encode()


def identity(path):
    value = Path(path).stat()
    return [value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns]


def checked_path(root, relative):
    relative = Path(relative)
    path = Path(root).absolute()/relative
    if relative.is_absolute() or '..' in relative.parts or path.resolve() != path:
        raise ValueError('Public preparation paths must stay under the selected root without symbolic links.')
    return path


def write(path, value, *, preserve=False):
    content = canonical(value)
    path = Path(path)
    if path.resolve() != path.absolute():
        raise ValueError('Public metadata must not contain symbolic links.')
    path.parent.mkdir(parents=True, exist_ok=True)
    if preserve and path.exists():
        old = path.read_bytes()
        if old == content: return
        archive = path.parent/'preparation-history'/(path.name+'.'+hashlib.sha256(old).hexdigest())
        if archive.resolve() != archive.absolute():
            raise ValueError('Producer history must not contain symbolic links.')
        archive.parent.mkdir(exist_ok=True)
        if not archive.exists():
            with archive.open('xb') as stream: stream.write(old)
        elif archive.read_bytes() != old:
            raise ValueError('Original producer evidence changed.')
    fd, temporary = tempfile.mkstemp(prefix='.preparing-', dir=path.parent)
    try:
        with os.fdopen(fd, 'wb') as stream:
            stream.write(content); stream.flush(); os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally: Path(temporary).unlink(missing_ok=True)


def preserve_record(path):
    """Retain exact producer bytes before issuing a different local receipt."""
    path = Path(path)
    content = path.read_bytes()
    archive = checked_path(path.parent, 'original-producer-evidence/'+path.name+'.'+hashlib.sha256(content).hexdigest())
    archive.parent.mkdir(exist_ok=True)
    if archive.exists():
        if archive.read_bytes() != content: raise ValueError('Original producer evidence changed.')
    else:
        with archive.open('xb') as stream:
            stream.write(content); stream.flush(); os.fsync(stream.fileno())
    return archive


def authenticate(path, expected, size, *, algorithm='sha256'):
    digest = hashlib.new(algorithm)
    count = 0
    before = identity(path)
    with Path(path).open('rb') as stream:
        while chunk := stream.read(CHUNK):
            digest.update(chunk); count += len(chunk)
            if count > size: raise ValueError('Public source exceeds its declared size.')
    if count != size or digest.hexdigest() != expected or identity(path) != before:
        raise ValueError('Public source differs from its frozen integrity binding. Restore the exact source.')
    return before


def acquire(path, row, *, opener=urllib.request.urlopen):
    """Acquire only the declared file. Never replace a different installed file."""
    expected_size = row.get('size', row.get('bytes'))
    print('Prepare public source: '+path.name, flush=True, file=sys.stderr)
    if path.exists():
        authenticate(path, row['sha256'], expected_size)
        return 'verified_local_source'
    path.parent.mkdir(parents=True, exist_ok=True)
    pending = path.with_name(path.name+'.partial')
    if pending.is_symlink(): raise ValueError('Incomplete source must not be a symbolic link.')
    digest = hashlib.sha256(); size = 0
    with opener(row['url'], timeout=45) as response, pending.open('wb') as output:
        while chunk := response.read(CHUNK):
            size += len(chunk)
            if size > expected_size: raise ValueError('The source URL now supplies a different frozen release.')
            output.write(chunk); digest.update(chunk)
    if size != expected_size or digest.hexdigest() != row['sha256']:
        raise ValueError('Frozen source is unavailable at its declared URL. The partial file remains incomplete.')
    os.replace(pending, path)
    return 'downloaded_and_verified'


def completion(root, role, bindings, outputs):
    if existing := reusable(root, role, bindings): return existing
    record = dict(schema='genosapiens.portable-preparation/v1', status='complete', role=role,
        producer='genosapiens.portable_resources', completed_at=datetime.now(timezone.utc).isoformat(),
        contract_sha256=hashlib.sha256(canonical(bindings)).hexdigest(),
        outputs={str(p.relative_to(root)):identity(p) for p in outputs})
    path = checked_path(root, 'manifests/installation/'+role+'.json')
    write(path, record, preserve=True)
    # Separate producer evidence binds the exact record and its local filesystem
    # identity. Editing output identities or copying this pair cannot renew it.
    write(path.with_suffix('.seal.json'), dict(schema='genosapiens.local-completion-seal/v1',
        record_sha256=hashlib.sha256(canonical(record)).hexdigest(), record_identity=identity(path),
        root_identity=identity(root)[:2]), preserve=True)
    return record


def reusable(root, role, bindings):
    path = checked_path(root, 'manifests/installation/'+role+'.json')
    if not path.is_file(): return None
    with path.open('rb') as stream: content = stream.read(CHUNK+1)
    if len(content) > CHUNK: return None
    seal_path = path.with_suffix('.seal.json')
    if not seal_path.is_file() or seal_path.resolve() != seal_path: return None
    try:
        with seal_path.open('rb') as stream: seal = json.loads(stream.read(CHUNK+1))
        if (seal.get('schema') != 'genosapiens.local-completion-seal/v1'
                or seal.get('record_sha256') != hashlib.sha256(content).hexdigest()
                or seal.get('record_identity') != identity(path)
                or seal.get('root_identity') != identity(root)[:2]): return None
    except (ValueError, AttributeError): return None
    try: record = json.loads(content)
    except (ValueError, UnicodeError): return None
    if (not isinstance(record, dict) or record.get('schema') != 'genosapiens.portable-preparation/v1'
            or record.get('role') != role or not isinstance(record.get('outputs'), dict)
            or record.get('status') != 'complete'
            or record.get('contract_sha256') != hashlib.sha256(canonical(bindings)).hexdigest()):
        return None
    for name, previous in record.get('outputs', {}).items():
        output = checked_path(root, name)
        if not output.is_file() or identity(output) != previous: return None
    return record if record.get('outputs') else None


def prepare_metadata(root, role, pins, *, opener=urllib.request.urlopen):
    root = Path(root).absolute()
    rows = pins['lincs_metadata' if role == 'lincs-metadata' else 'fda']
    if existing := reusable(root, role, rows): return existing
    relative = 'lincs/GSE92742' if role == 'lincs-metadata' else 'fda/2026-10-16'
    manifest = checked_path(root, relative+('/metadata-downloads.json' if role == 'lincs-metadata' else '/download-manifest.json'))
    progress = checked_path(root, 'manifests/installation/'+role+'-progress.json')
    write(progress, dict(status='incomplete', role=role))
    outputs, receipts = [], []
    for row in rows:
        filename = row['name' if role == 'lincs-metadata' else 'file']
        path = checked_path(root, relative+'/'+filename)
        origin = acquire(path, row, opener=opener)
        receipts.append(dict(row, acquisition_kind=origin, checked_at=datetime.now(timezone.utc).isoformat()))
        outputs.append(path)
        write(progress, dict(status='incomplete', role=role, completed_files=len(outputs), required_files=len(rows)))
    write(manifest, receipts, preserve=True)
    record = completion(root, role, rows, [*outputs, manifest])
    write(progress, dict(status='complete', role=role))
    return record


def prepare_matrix(root, archive, pin):
    """Authenticate and decompress one archive in one bounded pass."""
    root, archive = Path(root).absolute(), Path(archive).absolute()
    if archive.resolve() != archive: raise ValueError('The source archive must not contain symbolic links.')
    if existing := reusable(root, 'lincs-matrix', pin): return existing
    directory = checked_path(root, 'lincs/GSE92742'); directory.mkdir(parents=True, exist_ok=True)
    matrix = checked_path(root, 'lincs/GSE92742/'+pin['decompressed_file'])
    progress = checked_path(root, 'manifests/installation/lincs-matrix-progress.json')
    if matrix.exists():
        # A crash can occur after matrix publication but before its receipts.
        # Explicit preparation authenticates both files before recovering them.
        write(progress, dict(status='incomplete', operation='recover_matrix_receipts'))
        authenticate(matrix, pin['decompressed_sha256'], pin['decompressed_bytes'])
        authenticate(archive, pin['sha512'], pin['compressed_bytes'], algorithm='sha512')
        return matrix_receipts(root, pin, matrix, 'verified_local_archive_and_matrix')
    write(progress, dict(status='incomplete', compressed_bytes=0, decompressed_bytes=0))
    pending = matrix.with_name(matrix.name+'.partial')
    if pending.is_symlink(): raise ValueError('Incomplete matrix must not be a symbolic link.')
    compressed = hashlib.sha512(); expanded = hashlib.sha256()
    received = written = 0; decoder = zlib.decompressobj(16+zlib.MAX_WBITS)
    before = identity(archive)
    print('Authenticate and decompress the declared LINCS archive.', flush=True, file=sys.stderr)
    with archive.open('rb') as source, pending.open('wb') as output:
        while chunk := source.read(CHUNK):
            received += len(chunk); compressed.update(chunk)
            if received > pin['compressed_bytes']: raise ValueError('The matrix archive exceeds its declared size.')
            while chunk:
                raw = decoder.decompress(chunk, CHUNK)
                chunk = decoder.unconsumed_tail
                written += len(raw)
                if written > pin['decompressed_bytes']: raise ValueError('The matrix exceeds its declared size.')
                output.write(raw); expanded.update(raw)
            if decoder.unused_data: raise ValueError('The matrix archive has unexpected trailing data.')
        output.flush(); os.fsync(output.fileno())
    if (not decoder.eof or received != pin['compressed_bytes'] or written != pin['decompressed_bytes']
            or compressed.hexdigest() != pin['sha512'] or expanded.hexdigest() != pin['decompressed_sha256']
            or identity(archive) != before):
        raise ValueError('Matrix acquisition is incomplete or differs from its publisher and content bindings.')
    os.replace(pending, matrix)
    return matrix_receipts(root, pin, matrix, 'verified_local_archive')


def matrix_receipts(root, pin, matrix, acquisition_kind):
    directory = matrix.parent
    receipt = directory/'matrix-integrity.json'
    transfer = directory/'matrix-transfer-state.json'
    write(receipt, dict(pin, publisher_match=True, producer='genosapiens.portable_resources.prepare_matrix',
        checked_at=datetime.now(timezone.utc).isoformat(), acquisition_kind=acquisition_kind), preserve=True)
    write(transfer, dict(status='download_complete_verified_local_archive', total_bytes=pin['compressed_bytes'],
        acquisition_kind=acquisition_kind, source_archive=pin['compressed_file']), preserve=True)
    record = completion(root, 'lincs-matrix', pin, [matrix, receipt, transfer])
    write(checked_path(root, 'manifests/installation/lincs-matrix-progress.json'),
        dict(status='complete', compressed_bytes=pin['compressed_bytes'], decompressed_bytes=pin['decompressed_bytes']))
    return record


def admit_copied_matrix(root, pin):
    """Authenticate copied matrix bytes and preserve the original producer receipts."""
    root = Path(root).absolute()
    if existing := reusable(root, 'lincs-matrix', pin): return existing
    receipt = checked_path(root, 'lincs/GSE92742/matrix-integrity.json')
    transfer = checked_path(root, 'lincs/GSE92742/matrix-transfer-state.json')
    original = json.loads(receipt.read_bytes()); acquisition = json.loads(transfer.read_bytes())
    if (any(original.get(k) != v for k,v in pin.items()) or original.get('publisher_match') is not True
            or not acquisition.get('status','').startswith('download_complete')
            or acquisition.get('total_bytes') != pin['compressed_bytes']):
        raise ValueError('Copied resources need valid original producer evidence. Restore the declared archive and receipts.')
    matrix = checked_path(root, 'lincs/GSE92742/'+pin['decompressed_file'])
    progress = checked_path(root, 'manifests/installation/lincs-matrix-progress.json')
    write(progress, dict(status='incomplete', operation='authenticate_copied_matrix'))
    print('Authenticate the copied LINCS matrix.', flush=True, file=sys.stderr)
    authenticate(matrix, pin['decompressed_sha256'], pin['decompressed_bytes'])
    record = completion(root, 'lincs-matrix', pin, [matrix, receipt, transfer])
    write(progress, dict(status='complete', operation='authenticate_copied_matrix'))
    return record


def prepare_open_targets_metadata(root, contract, *, opener=urllib.request.urlopen):
    root = Path(root).absolute()
    if existing := reusable(root, 'open-targets-metadata', contract): return existing
    progress = checked_path(root, 'manifests/installation/open-targets-metadata-progress.json')
    write(progress, dict(status='incomplete'))
    outputs = []
    base_url = contract['base_url'].removesuffix('output/').rstrip('/')
    for row in contract['metadata']:
        path = checked_path(root, contract['output_prefix']+'/'+row['name'])
        acquire(path, dict(row, size=row['size_bytes'], url=base_url+'/'+row['name']), opener=opener)
        outputs.append(path)
    record = completion(root, 'open-targets-metadata', contract, outputs)
    write(progress, dict(status='complete'))
    return record


def admit_manifest_outputs(root, relative, expected, *, role):
    """Explicitly authenticate one manifest's copied outputs, preserving its bytes."""
    root = Path(root).absolute()
    path = checked_path(root, relative)
    if path.stat().st_size > 8*CHUNK:
        raise ValueError('Source manifest exceeds the bounded metadata contract.')
    authenticate(path, expected, path.stat().st_size)
    manifest = json.loads(path.read_bytes())
    outputs = manifest['outputs']
    if isinstance(outputs, dict): outputs = list(outputs.values())
    bindings = dict(manifest_path=relative, manifest_sha256=expected)
    if existing := reusable(root, role, bindings): return existing
    progress = checked_path(root, 'manifests/installation/'+role+'-progress.json')
    write(progress, dict(status='incomplete'))
    paths = [path]
    for row in outputs:
        output = checked_path(root, row['path'])
        authenticate(output, row['sha256'], row['size_bytes'])
        paths.append(output)
    record = completion(root, role, bindings, paths)
    write(progress, dict(status='complete'))
    return record
