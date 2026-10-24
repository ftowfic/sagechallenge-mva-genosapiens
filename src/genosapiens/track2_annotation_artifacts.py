"""File-backed allele arrays with complete byte authentication and bounded parsing."""
from contextlib import contextmanager
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import stat
from types import MappingProxyType

CHUNK_BYTES = 64 * 1024
MAX_ALLELE_BYTES = 64 * 1024 * 1024


def signature(info):
    return (info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns, info.st_ctime_ns)


@contextmanager
def regular_file(path):
    path = Path(path).absolute()
    if path.resolve() != path:
        raise ValueError('Annotation artifact paths must not contain links.')
    fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    with os.fdopen(fd, 'rb') as stream:
        before = os.fstat(stream.fileno())
        if not stat.S_ISREG(before.st_mode):
            raise ValueError('An annotation artifact is not a regular file.')
        yield stream, signature(before)
        if signature(os.fstat(stream.fileno())) != signature(before) or signature(path.lstat()) != signature(before):
            raise ValueError('An annotation artifact changed during reading.')


@dataclass(frozen=True)
class AlleleArtifact:
    path: Path
    sha256: str
    identity: tuple

    @classmethod
    def capture(cls, path):
        with regular_file(path) as (stream, identity):
            digest = hashlib.sha256()
            for chunk in iter(lambda: stream.read(CHUNK_BYTES), b''):
                digest.update(chunk)
        return cls(Path(path).absolute(), digest.hexdigest(), identity)

    def records(self, expected_sha256):
        """Exhaust the iterator before using records as accepted evidence.

        Only one serialized allele is buffered. Typed alleles remain owned by
        the existing admission API, whose working set still scales with them.
        """
        from genosapiens.vep_native_annotation import _strict_json
        if expected_sha256 != self.sha256:
            raise ValueError('Annotation allele digest does not match its publication.')
        with regular_file(self.path) as (stream, identity):
            if identity != self.identity:
                raise ValueError('The annotation allele artifact changed.')
            digest = hashlib.sha256()
            buffer = b''
            eof = False

            def fill():
                nonlocal buffer, eof
                chunk = stream.read(CHUNK_BYTES)
                digest.update(chunk)
                buffer += chunk
                eof = not chunk

            def space():
                nonlocal buffer
                while True:
                    buffer = buffer.lstrip(b' \t\r\n')
                    if buffer or eof:
                        return
                    fill()

            space()
            if not buffer.startswith(b'['):
                raise ValueError('Annotation alleles must form a JSON array.')
            buffer = buffer[1:]
            space()
            first = True
            decoder = json.JSONDecoder()
            while not buffer.startswith(b']'):
                if not first:
                    if not buffer.startswith(b','):
                        raise ValueError('Annotation alleles have an invalid separator.')
                    buffer = buffer[1:]
                    space()
                if not buffer.startswith(b'{'):
                    raise ValueError('An annotation allele must be an object.')
                while True:
                    try:
                        text = buffer.decode('utf-8')
                        _, end = decoder.raw_decode(text)
                        raw = text[:end].encode('utf-8')
                        if len(raw) > MAX_ALLELE_BYTES:
                            raise ValueError('One annotation allele exceeds the record limit.')
                        break
                    except (json.JSONDecodeError, UnicodeDecodeError):
                        if eof or len(buffer) > MAX_ALLELE_BYTES:
                            raise ValueError('An annotation allele is incomplete or exceeds the record limit.') from None
                        fill()
                value = _strict_json(raw.decode('utf-8'))
                buffer = buffer[len(raw):]
                yield value
                first = False
                space()
            buffer = buffer[1:]
            space()
            if buffer or not eof or digest.hexdigest() != expected_sha256:
                raise ValueError('Annotation allele bytes are incomplete or changed.')


@dataclass(frozen=True, repr=False)
class FileAnnotationInput:
    artifacts: object
    requests: tuple
    expected_binding_sha256: str

    def __post_init__(self):
        from genosapiens.candidate_admission import AnnotationInput
        values = dict(self.artifacts)
        if type(values.get('alleles.json')) is not AlleleArtifact:
            raise ValueError('A file annotation input requires an allele artifact.')
        # Preserve the established validation of every other field.
        AnnotationInput({**values, 'alleles.json': b'[]'}, self.requests, self.expected_binding_sha256)
        object.__setattr__(self, 'artifacts', MappingProxyType(values))


def artifact_digest(value):
    return value.sha256 if type(value) is AlleleArtifact else hashlib.sha256(value).hexdigest()
