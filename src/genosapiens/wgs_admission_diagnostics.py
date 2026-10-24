"""Private admission metadata. This module never grants scientific completion.

The record contains finite labels, private binding digests, monotonic times,
and bounded window counters. It excludes exception text and scientific data.
The invocation token identifies one invocation. The admission digest identifies
its contract. Terminal receipt and failure digests authenticate terminal detail.
A last observation cannot establish a failure. Shared summaries omit all digests
and counters. Writers replace one private file atomically, at most once per
second during row work. Blocking boundaries and terminal states bypass that rate.
"""
from __future__ import annotations

from contextlib import contextmanager
from contextvars import ContextVar
import hashlib
import json
import math
import os
from pathlib import Path
import secrets
import sqlite3
import stat
import threading
import time

SCHEMA = 'genosapiens.admission-diagnostic'
VERSION = 1
MAX_BYTES = 16 * 1024
MAX_SECONDARY = 4
OPERATIONS = frozenset({
    'admission_result_collect', 'admission_worker_execute', 'admission_screening_receive',
    'profile_window_validate', 'profile_window_decode', 'profile_clinical_query',
    'profile_clinical_bind', 'profile_clinical_adapt', 'profile_row_convert',
    'profile_row_persist', 'admission_window_persist', 'admission_window_retire',
    'admission_complete', 'progress_publish', 'journal_publish', 'writer_close',
    'screening_close', 'worker_close', 'runtime_finalize', 'diagnostic_publish', 'monitor_close',
})
CATEGORIES = frozenset({'memory_error', 'interrupted', 'timeout', 'sqlite_error', 'os_error',
    'value_error', 'type_error', 'key_error', 'runtime_error', 'attribute_error',
    'index_error', 'assertion_error', 'iteration_error', 'other'})
REASONS = frozenset({'unknown', 'incomplete_required_evidence', 'source_binding_mismatch',
    'exact_identity_mismatch', 'malformed_evidence', 'lookup_failed', 'invalid_worker_envelope'})
RECOVERY = frozenset({'inspect_local_diagnostic_then_resume', 'check_required_evidence_then_resume',
    'check_storage_then_resume', 'check_resources_then_resume'})
STATES = frozenset({'initializing', 'active', 'failed', 'complete'})
COUNTERS = frozenset({'screened', 'pending', 'committed', 'cancelled', 'sealed'})
FAILURE_FIELDS = {'operation', 'category', 'reason', 'recovery'}
FIELDS = {'schema', 'version', 'invocation', 'admission_sha256', 'receipt_sha256',
    'failure_sha256', 'state', 'operation', 'operation_started', 'observed', 'counters',
    'primary', 'secondary', 'secondary_truncated'}
_CURRENT = ContextVar('admission_diagnostic_owner', default=None)


def digest(content):
    return hashlib.sha256(content).hexdigest()


def _digest(value):
    return type(value) is str and len(value) == 64 and all(c in '0123456789abcdef' for c in value)


def category(error):
    from genosapiens.wgs_failures import WgsFailure
    if isinstance(error, WgsFailure):
        return error.exception_kind if error.exception_kind in CATEGORIES else 'other'
    for family, label in ((MemoryError, 'memory_error'),
            ((KeyboardInterrupt, SystemExit, GeneratorExit), 'interrupted'),
            (TimeoutError, 'timeout'), (sqlite3.Error, 'sqlite_error'), (OSError, 'os_error'),
            (KeyError, 'key_error'), (IndexError, 'index_error'), (AttributeError, 'attribute_error'),
            (AssertionError, 'assertion_error'), ((StopIteration, StopAsyncIteration), 'iteration_error'),
            (ValueError, 'value_error'), (TypeError, 'type_error'), (RuntimeError, 'runtime_error')):
        if isinstance(error, family):
            return label
    return 'other'


def failure(operation, error, reason='unknown'):
    cat = category(error)
    recovery = ('check_required_evidence_then_resume' if reason != 'unknown' else
        'check_storage_then_resume' if cat in {'sqlite_error', 'os_error'} else
        'check_resources_then_resume' if cat == 'memory_error' else
        'inspect_local_diagnostic_then_resume')
    value = dict(operation=operation, category=cat, reason=reason, recovery=recovery)
    _validate_failure(value)
    return value


def _validate_failure(value):
    if (type(value) is not dict or set(value) != FAILURE_FIELDS
            or type(value['operation']) is not str or value['operation'] not in OPERATIONS
            or type(value['category']) is not str or value['category'] not in CATEGORIES
            or type(value['reason']) is not str or value['reason'] not in REASONS
            or type(value['recovery']) is not str or value['recovery'] not in RECOVERY):
        raise ValueError('admission diagnostic is unavailable')


def detail(error):
    value = getattr(error, '_admission_failure', None)
    try:
        _validate_failure(value)
    except ValueError:
        return None
    return dict(value)


def tag(error, operation, reason='unknown'):
    value = detail(error) or failure(operation, error, reason)
    error._admission_failure = value
    owner = _CURRENT.get()
    if owner is not None:
        owner.observe_failure(value)
    return error


@contextmanager
def boundary(operation, *, blocking=False, reason='unknown', known_errors=()):
    if operation not in OPERATIONS or reason not in REASONS:
        raise ValueError('invalid admission boundary')
    owner = _CURRENT.get()
    if owner is not None:
        owner.enter(operation, force=blocking)
    try:
        yield
    except BaseException as error:
        observed_reason = reason if not known_errors or isinstance(error, known_errors) else 'unknown'
        tag(error, operation, observed_reason)
        raise


def encode(value):
    validate(value)
    data = json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode() + b'\n'
    if len(data) > MAX_BYTES:
        raise ValueError('admission diagnostic is unavailable')
    return data


def validate(value):
    if type(value) is not dict or set(value) != FIELDS:
        raise ValueError('admission diagnostic is unavailable')
    if (value['schema'] != SCHEMA or type(value['version']) is not int or value['version'] != VERSION
            or not _digest(value['invocation'])
            or any(v is not None and not _digest(v) for v in
                (value['admission_sha256'], value['receipt_sha256'], value['failure_sha256']))
            or type(value['state']) is not str or value['state'] not in STATES
            or type(value['operation']) is not str or value['operation'] not in OPERATIONS):
        raise ValueError('admission diagnostic is unavailable')
    for key in ('operation_started', 'observed'):
        v = value[key]
        if type(v) not in (int, float) or not math.isfinite(v) or not 0 <= v <= 1e15:
            raise ValueError('admission diagnostic is unavailable')
    if value['observed'] < value['operation_started']:
        raise ValueError('admission diagnostic is unavailable')
    counters = value['counters']
    if (type(counters) is not dict or set(counters) != COUNTERS
            or any(type(v) is not int or not 0 <= v <= 2**63-1 for v in counters.values())
            or counters['sealed'] > counters['committed']):
        raise ValueError('admission diagnostic is unavailable')
    if value['primary'] is not None:
        _validate_failure(value['primary'])
    if (type(value['secondary']) is not list or len(value['secondary']) > MAX_SECONDARY
            or type(value['secondary_truncated']) is not bool
            or value['secondary_truncated'] and len(value['secondary']) != MAX_SECONDARY):
        raise ValueError('admission diagnostic is unavailable')
    for item in value['secondary']:
        _validate_failure(item)
    if ((value['state'] == 'failed') != (value['primary'] is not None)
            or (value['receipt_sha256'] is None) != (value['failure_sha256'] is None)
            or value['receipt_sha256'] is not None and
                (value['state'] != 'failed' or value['admission_sha256'] is None)):
        raise ValueError('admission diagnostic is unavailable')


def _pairs(items):
    value = {}
    for key, item in items:
        if key in value:
            raise ValueError('admission diagnostic is unavailable')
        value[key] = item
    return value


def read(path, *, invocation=None, admission_sha256=None, receipt_sha256=None, failure_sha256=None):
    """Reject special files before reading. Invalid detail returns no observation."""
    fd = -1
    try:
        path = Path(path)
        if any(p.is_symlink() for p in path.parents):
            return None
        fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
        before = os.fstat(fd)
        if (not stat.S_ISREG(before.st_mode) or not 0 < before.st_size <= MAX_BYTES
                or stat.S_IMODE(before.st_mode) & 0o077):
            return None
        content = os.read(fd, MAX_BYTES + 1)
        after = os.fstat(fd)
        signature = lambda s: (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns)
        if (signature(before) != signature(after) or signature(after) != signature(path.lstat())
                or len(content) != before.st_size):
            return None
        value = json.loads(content, object_pairs_hook=_pairs,
            parse_constant=lambda _: (_ for _ in ()).throw(ValueError()))
        validate(value)
        expected = dict(invocation=invocation, admission_sha256=admission_sha256,
            receipt_sha256=receipt_sha256, failure_sha256=failure_sha256)
        if any(v is not None and value[k] != v for k, v in expected.items()):
            return None
        return value
    except (OSError, ValueError, TypeError, RecursionError):
        return None
    finally:
        if fd >= 0:
            os.close(fd)


def publish(path, value, *, durable=True):
    """Replace one complete observation atomically.

    Active observations can disappear after power failure. Terminal records
    synchronize both the file and its directory. Neither authorizes admission.
    Scientific stores and checkpoint synchronization remain independent.
    """
    path = Path(path)
    content = encode(value)
    if any(p.is_symlink() for p in (path, *path.parents)):
        raise ValueError('admission diagnostic is unavailable')
    path.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
    temporary = path.with_name('.admission-diagnostic-' + secrets.token_hex(16))
    fd = os.open(temporary, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600)
    try:
        with os.fdopen(fd, 'wb') as stream:
            stream.write(content)
            stream.flush()
            if durable:
                os.fsync(stream.fileno())
        os.replace(temporary, path)
        if durable:
            directory = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
            try:
                os.fsync(directory)
            finally:
                os.close(directory)
    finally:
        temporary.unlink(missing_ok=True)


class Owner:
    """One invocation owns one terminal write attempt. Metadata is best effort."""
    def __init__(self, path, *, invocation=None, clock=time.monotonic, writer=publish):
        self.path, self.clock, self.writer = Path(path), clock, writer
        self.terminal_attempted = False
        self.primary_error = None
        self.last_write = -math.inf
        now = clock()
        self.value = dict(schema=SCHEMA, version=VERSION, invocation=invocation or secrets.token_hex(32),
            admission_sha256=None, receipt_sha256=None, failure_sha256=None, state='initializing',
            operation='admission_result_collect', operation_started=now, observed=now,
            counters={key: 0 for key in COUNTERS}, primary=None, secondary=[], secondary_truncated=False)

    @contextmanager
    def activate(self):
        token = _CURRENT.set(self)
        try:
            yield self
        finally:
            _CURRENT.reset(token)

    def bind(self, binding):
        self.value['admission_sha256'] = digest(binding)
        self.value['state'] = 'active'
        self.snapshot(force=True)

    def enter(self, operation, *, force=False):
        if self.value['primary'] is not None:
            return
        now = self.clock()
        if operation != self.value['operation']:
            self.value.update(operation=operation, operation_started=now)
        self.snapshot(force=force)

    def snapshot(self, *, counters=None, force=False):
        if counters is not None:
            self.value['counters'] = dict(counters)
        now = self.clock()
        self.value['observed'] = now
        if self.terminal_attempted or not force and now - self.last_write < 1:
            return
        self.last_write = now
        try:
            self.writer(self.path, self.value, durable=False)
        except Exception as error:
            self.secondary(error, 'diagnostic_publish')

    def observe_failure(self, value):
        if self.value['primary'] is None:
            self.value.update(primary=dict(value), state='failed', operation=value['operation'])

    def secondary(self, error, operation):
        if error is self.primary_error:
            return
        item = failure(operation, error)
        if len(self.value['secondary']) < MAX_SECONDARY:
            self.value['secondary'].append(item)
        else:
            self.value['secondary_truncated'] = True

    def capture(self, error, operation='admission_result_collect'):
        if self.primary_error is None:
            self.primary_error = error
        self.observe_failure(detail(error) or failure(operation, error))

    def terminal(self, *, receipt=None, primary=None):
        if self.terminal_attempted:
            return
        self.terminal_attempted = True
        if receipt is not None and primary is not None and self.value['primary'] is not None:
            self.value.update(receipt_sha256=digest(receipt), failure_sha256=digest(primary))
        self.value['observed'] = self.clock()
        try:
            self.writer(self.path, self.value)
        except Exception as error:
            self.secondary(error, 'diagnostic_publish')

    @contextmanager
    def cleanup(self, operation):
        try:
            yield
        except BaseException as error:
            if self.value['primary'] is None:
                self.capture(error, operation)
                raise
            self.secondary(error, operation)


@contextmanager
def preserving_close(manager, operation):
    """Preserve the operation failure when a disposable owner's close fails."""
    owner = _CURRENT.get()
    primary = None
    try:
        yield manager
    except BaseException as error:
        primary = error
        if owner is not None:
            owner.capture(error)
        raise
    finally:
        try:
            manager.close()
        except BaseException as error:
            if primary is None:
                tag(error, operation)
                raise
            if owner is not None:
                owner.secondary(error, operation)
            if error is not primary:
                from genosapiens.wgs_failures import run_cleanup
                run_cleanup((('resource_session_close', lambda: (_ for _ in ()).throw(error)),), primary=primary)


def current():
    return _CURRENT.get()


def worker_envelope(error):
    value = detail(error) or failure('admission_worker_execute', error)
    return {key: value[key] for key in ('operation', 'category', 'reason')}


def worker_error(value):
    from genosapiens.wgs_failures import WgsFailure, EXCEPTION_KINDS
    if (type(value) is not dict or set(value) != {'operation', 'category', 'reason'}
            or type(value['operation']) is not str or value['operation'] != 'admission_worker_execute'
            or type(value['category']) is not str or value['category'] not in CATEGORIES
            or type(value['reason']) is not str or value['reason'] not in REASONS):
        error = ValueError('admission worker response is invalid')
        return tag(error, 'admission_result_collect', 'invalid_worker_envelope')
    kind = value['category']
    error = WgsFailure('bounded_products',
        'memory_exhaustion' if kind == 'memory_error' else
        'interrupted' if kind == 'interrupted' else 'timeout' if kind == 'timeout' else 'unexpected',
        operation='annotation_replay', exception_kind=kind if kind in EXCEPTION_KINDS else 'other')
    error._admission_failure = dict(value, recovery=failure(value['operation'], error, value['reason'])['recovery'])
    return error


def terminal_summary(value):
    try:
        validate(value)
    except (ValueError, TypeError):
        return 'Admission failure detail is unavailable.'
    if value is None or value['state'] != 'failed' or value['receipt_sha256'] is None:
        return 'Admission failure detail is unavailable.'
    item = value['primary']
    return ('Admission failure: operation=' + item['operation'] + ', category=' + item['category']
        + ', reason=' + item['reason'] + ', recovery=' + item['recovery'] + '.')


class Monitor:
    """Read one bounded metadata file. Print at most one notice per 30 seconds."""
    def __init__(self, path, invocation, *, clock=time.monotonic, emit=print, interval=30):
        self.path, self.invocation, self.clock, self.emit = path, invocation, clock, emit
        self.interval = max(30, interval)
        self.last_notice = -math.inf
        self.stopped = threading.Event()
        self.thread = None

    def poll(self):
        now = self.clock()
        if now - self.last_notice < self.interval:
            return
        value = read(self.path, invocation=self.invocation)
        if value is None or value['state'] != 'active':
            return
        self.last_notice = now
        elapsed = max(0, int(now - value['operation_started']))
        self.emit('Admission last observation: ' + value['operation'] + '. Operation elapsed seconds: '
            + str(elapsed) + '. This observation does not establish progress or a failure cause.')

    def start(self):
        def run():
            while not self.stopped.wait(1):
                try:
                    self.poll()
                except Exception:
                    pass
        self.thread = threading.Thread(target=run, name='admission-operator-monitor', daemon=True)
        self.thread.start()
        return self

    def close(self):
        self.stopped.set()
        if self.thread is not None:
            self.thread.join(timeout=2)
            if self.thread.is_alive():
                raise RuntimeError('admission monitor shutdown failed')


def receipt_binding(workspace):
    """Read only bounded terminal metadata. Never read journals or case products."""
    from genosapiens.bounded_json import read_bounded_bytes
    from genosapiens.wgs_analysis import AnalysisRunReceipt, ANALYSIS_RECEIPT_NAME
    from genosapiens.evidence_model import _typed
    from genosapiens.wgs_failures import load_failure_evidence
    try:
        receipt = read_bounded_bytes(Path(workspace)/ANALYSIS_RECEIPT_NAME, maximum_bytes=MAX_BYTES)
        model = _typed(json.loads(receipt, object_pairs_hook=_pairs), AnalysisRunReceipt, decode=True)
        if model.as_bytes() != receipt:
            return None
        if model.status == 'passed':
            return None
        content = read_bounded_bytes(Path(workspace)/'failure-analysis.json', maximum_bytes=MAX_BYTES)
        value = json.loads(content, object_pairs_hook=_pairs)
        value = load_failure_evidence(content, receipt_content=receipt,
            expected_failure_code=None, expected_journal_heads=value.get('journal_heads'))
        if value['stage'] != 'analysis':
            return None
        return receipt, content, value
    except (OSError, ValueError, TypeError, KeyError, RuntimeError):
        return None


def finish_terminal(workspace, owner, *, emit=print, failed=False):
    bound = receipt_binding(workspace)
    expected = getattr(owner, 'terminal_binding', None)
    if bound is not None and expected == (digest(bound[0]), digest(bound[1])):
        owner.terminal(receipt=bound[0], primary=bound[1])
        value = read(owner.path, invocation=owner.value['invocation'],
            admission_sha256=owner.value['admission_sha256'],
            receipt_sha256=digest(bound[0]), failure_sha256=digest(bound[1]))
    else:
        owner.terminal()
        value = None
    if failed or owner.value['primary'] is not None or bound is not None:
        try:
            emit(terminal_summary(value))
        except Exception as error:
            owner.secondary(error, 'diagnostic_publish')


@contextmanager
def analysis_observation(workspace):
    """Standalone runner calls also own current invocation metadata."""
    owned = current() is None
    owner = current() or Owner(Path(workspace)/'analysis-work/state/admission-diagnostic.json')
    try:
        with owner.activate():
            yield owner
    finally:
        if owned:
            finish_terminal(workspace, owner)


def close_preserving(action, operation, *, primary=None):
    """Attempt cleanup once. Preserve the body failure and bounded static detail."""
    owner = current()
    try:
        action()
    except BaseException as error:
        if primary is None:
            tag(error, operation)
            raise
        if error is primary:
            return
        if owner is not None:
            owner.capture(primary)
            owner.secondary(error, operation)
        from genosapiens.wgs_failures import run_cleanup
        shared = 'worker_shutdown' if operation == 'worker_close' else 'resource_session_close'
        run_cleanup(((shared, lambda: (_ for _ in ()).throw(error)),), primary=primary)
