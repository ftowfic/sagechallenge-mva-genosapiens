#!/usr/bin/env python3
"""Host plumbing for public-only Docker preparation. Uses Python's standard library."""
import argparse
import csv
import fcntl
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile
import signal
import uuid


def mount(source, target, readonly=True):
    stream = io.StringIO()
    csv.writer(stream, lineterminator='').writerow(
        ['type=bind', 'source=' + str(source), 'target=' + str(target)] + (['readonly'] if readonly else []))
    return stream.getvalue()


def private_directory(path):
    if path.is_symlink():
        raise ValueError('The package directory must not be a symlink.')
    path.mkdir(mode=0o700, exist_ok=True)
    info = path.stat()
    if info.st_uid != os.getuid() or info.st_mode & 0o777 != 0o700:
        raise ValueError('The package directory must be owned by you with mode 0700.')


def preparation_command(repository, reference, store, image, name):
    args = ['docker', 'run', '--rm', '--init', '--network', 'none', '--read-only',
        '--user', f'{os.getuid()}:{os.getgid()}', '--cpus', '12', '--memory', '24g',
        '--memory-swap', '24g', '--pids-limit', '512', '--cap-drop', 'ALL',
        '--security-opt', 'no-new-privileges', '--tmpfs', '/tmp:rw,exec,nosuid,nodev,size=1g',
        '--env', 'PYTHONPATH=/workspace/src', '--env', 'PYTHONDONTWRITEBYTECODE=1',
        '--workdir', '/workspace', '--name', name, '--entrypoint', 'python']
    for name in ('src', 'config', 'schemas', 'scripts', 'docs', 'requirements', 'track1'):
        args += ['--mount', mount(repository / name, Path('/workspace') / name)]
    for name in ('releases', 'manifests', 'derived', 'downloads'):
        args += ['--mount', mount(reference / name, Path('/mnt/reference') / name)]
    # Only the public producers' output directories are writable.
    for producer in ('common-allele-public-indexes', 'track1-early-scope'):
        for name in ('manifests', 'derived'):
            path = reference / name / producer
            if path.is_symlink():
                raise ValueError('Public preparation output must not be a symlink.')
            path.mkdir(exist_ok=True)
            args += ['--mount', mount(path, Path('/mnt/reference') / name / producer, False)]
    args += ['--mount', mount(store, '/packages', False), image,
        '/workspace/scripts/build_track1_installation.py', '--repository-root', '/workspace',
        '--reference-root', '/mnt/reference', '--package-store', '/packages']
    return args


def ensure(repository, reference):
    repository, reference = repository.resolve(strict=True), reference.resolve(strict=True)
    image = json.loads((repository / 'config/track1/runtime-image.json').read_text())['image']
    if not re.fullmatch(r'genosapiens/track1-pipeline(?::[A-Za-z0-9._-]+)?@sha256:[0-9a-f]{64}', image):
        raise ValueError('The runtime image pin is invalid.')
    subprocess.run(['docker', 'image', 'inspect', image], check=True, stdout=subprocess.DEVNULL)
    store = repository / '.track1-packages'
    private_directory(store)
    descriptor = os.open(store / '.setup.lock', os.O_CREAT | os.O_RDWR | os.O_NOFOLLOW, 0o600)
    with os.fdopen(descriptor, 'r+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        name = 'track1-public-prepare-' + uuid.uuid4().hex
        try:
            result = subprocess.run(preparation_command(repository, reference, store, image, name),
                check=True, stdout=subprocess.PIPE, text=True, pass_fds=(lock.fileno(),))
        finally:
            # A cancelled Docker client must not leave its public producer running.
            subprocess.run(['docker', 'rm', '-f', name], check=False,
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=30)
        record = json.loads(result.stdout)
        name, sha = record['package_name'], record['package_sha256']
        if (not re.fullmatch(r'[0-9a-f]{64}', sha) or name != 'package-' + sha
                or record['runtime_image_digest'] != image.split('@')[-1]):
            raise ValueError('Prepared package identity is invalid.')
        package = store / name
        if hashlib.sha256((package / 'package-bindings.json').read_bytes()).hexdigest() != sha:
            raise ValueError('Prepared package metadata differs.')
        content = f'{image}\n{package}\n{sha}\n'
        target = repository / '.track1-installation'
        if target.is_symlink():
            raise ValueError('Installation metadata must not be a symlink.')
        if target.exists():
            old = target.read_bytes()
            if old == content.encode():
                return
            backup = store / ('installation-' + hashlib.sha256(old).hexdigest() + '.txt')
            if not backup.exists():
                with backup.open('xb') as stream:
                    stream.write(old)
        fd, temporary = tempfile.mkstemp(prefix='.track1-installation-', dir=repository)
        try:
            with os.fdopen(fd, 'w') as stream:
                stream.write(content)
                stream.flush()
                os.fsync(stream.fileno())
            os.replace(temporary, target)
        finally:
            Path(temporary).unlink(missing_ok=True)


def main():
    os.umask(0o077)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repository-root', type=Path, required=True)
    parser.add_argument('--reference-root', type=Path, required=True)
    parser.add_argument('--max-hours', type=float)
    args = parser.parse_args()
    if args.max_hours is not None and (not math.isfinite(args.max_hours) or args.max_hours <= 0):
        parser.error('--max-hours must be positive and finite.')
    def cancelled(signum, frame):
        raise TimeoutError('Preparation reached the invocation limit or received a stop signal.')
    signal.signal(signal.SIGTERM, cancelled)
    if args.max_hours is not None:
        signal.signal(signal.SIGALRM, cancelled)
        signal.setitimer(signal.ITIMER_REAL, args.max_hours * 3600)
    try:
        ensure(args.repository_root, args.reference_root)
    except subprocess.CalledProcessError as error:
        parser.exit(2, f'Track 1 preparation stopped (Docker exit {error.returncode}). Case processing did not start.\n')
    except (OSError, ValueError, KeyError, subprocess.TimeoutExpired) as error:
        parser.exit(2, f'Track 1 preparation failed: {error}\n')
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)


if __name__ == '__main__':
    main()
