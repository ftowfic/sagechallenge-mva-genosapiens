"""Independent stage expectations and destination-owned atomic publication."""
from contextlib import contextmanager
import fcntl
import os
from pathlib import Path
import time
import uuid

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_scientific_storage import publish_json,json_identity,file_identity
from genosapiens.track2_scientific_reader import read_scientific
from genosapiens.track2_challenge_revision import inventory


def private_directory(path):
    path=Path(path).absolute()
    require(path.resolve()==path, 'Report operation links are unsupported.')
    path.mkdir(mode=0o700,parents=True,exist_ok=True)
    require(path.stat().st_uid==os.geteuid() and not path.stat().st_mode & 0o077,
            'Report operations require a private owned directory.')
    return path


@contextmanager
def destination_lock(destination, operation, *, identity, resume):
    destination=Path(destination).absolute();operation=Path(operation).absolute()
    require(destination.resolve()==destination and operation.resolve()==operation
        and not destination.is_relative_to(operation) and not operation.is_relative_to(destination),
        'Report publication requires independent paths without links.')
    require(operation.parent==destination.parent, 'Report staging must share the destination parent filesystem.')
    existed=operation.exists()
    require(resume or not existed and not destination.exists(), 'Existing report work requires explicit resume.')
    private_directory(operation)
    fd=os.open(operation/'.lock',os.O_CREAT|os.O_RDWR|os.O_NOFOLLOW,0o600)
    try:
        try:fcntl.flock(fd,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:raise ValueError('Another invocation owns this report destination.') from None
        request=dict(schema='genosapiens.track2-report-operation/v1',destination=str(destination),identity=identity)
        path=operation/'operation.json'
        if path.exists():
            require(resume and read_scientific(path)==request, 'Report resume has another operation identity.')
        else:
            require(not destination.exists() and not existed, 'The destination lacks matching operation ownership.')
            publish_json(path,request)
        yield operation
    finally:os.close(fd)


def code_binding(*modules):
    root=Path(__file__).resolve().parent
    return {name:file_identity(root/(name+'.py'))['sha256'] for name in modules}


def stage(operation, name, binding, produce, *, on_stage=None):
    """A coordinator-owned expectation precedes any later output reuse."""
    operation=Path(operation);stages=private_directory(operation/'stages')
    expected=stages/(name+'.json')
    if expected.exists():
        record=read_scientific(expected)
        require(record['binding']==binding and record['status']=='complete', 'Report stage dependencies changed.')
        relative=Path(record['attempt'])
        require(relative.parts[:1]==('attempts',) and len(relative.parts)==2 and '..' not in relative.parts,
                'The report stage has an invalid attempt reference.')
        output=operation/relative/'output'
        require(inventory(output)==record['outputs'], 'A completed report stage changed.')
        return output,True
    attempts=private_directory(operation/'attempts');attempt=private_directory(attempts/uuid.uuid4().hex)
    publish_json(attempt/'started.json',dict(stage=name,binding=binding,status='incomplete'))
    start=time.monotonic()
    try:
        produce(attempt/'output')
        outputs=inventory(attempt/'output')
        record=dict(stage=name,binding=binding,status='complete',outputs=outputs,
            attempt=attempt.relative_to(operation).as_posix())
        publish_json(expected,record)
        if on_stage is not None:on_stage(name)
    except BaseException as error:
        publish_json(attempt/'outcome.json',dict(status='incomplete',stage=name,error_type=type(error).__name__,
            elapsed_seconds=time.monotonic()-start))
        raise
    publish_json(attempt/'outcome.json',dict(status='complete',stage=name,elapsed_seconds=time.monotonic()-start))
    return attempt/'output',False


def atomic_publish(staging, destination):
    staging=Path(staging).absolute();destination=Path(destination).absolute()
    require(staging.resolve()==staging and destination.resolve()==destination and not destination.exists(),
            'Atomic report publication requires a new destination without links.')
    require(staging.stat().st_dev==destination.parent.stat().st_dev,
            'Report staging is on another filesystem.')
    # The launcher mounts the common parent. Separate sibling bind mounts can
    # still cause EXDEV even when st_dev matches, so never fall back to copying.
    os.rename(staging,destination)
    fd=os.open(destination.parent,os.O_RDONLY|os.O_DIRECTORY)
    try:os.fsync(fd)
    finally:os.close(fd)
