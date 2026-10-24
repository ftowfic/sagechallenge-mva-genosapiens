"""Retrying, validated, atomic public-resource downloads."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import os
from pathlib import Path
import re
import time
from typing import Any, Callable
from urllib.error import URLError
from urllib.request import Request, urlopen
import uuid

from genosapiens.storage import CapacityPlan, require_capacity


SHA256 = re.compile(r"^[0-9a-f]{64}$")
SHA1 = re.compile(r"^[0-9a-f]{40}$")
MD5 = re.compile(r"^[0-9a-f]{32}$")


class DownloadError(RuntimeError):
    """Raised when a resource cannot be safely downloaded and published."""


class DownloadValidationError(DownloadError):
    """Raised when downloaded bytes do not match declared source identity."""


class CompletedResourceMismatch(DownloadError):
    """Raised when an existing published resource has unexpected contents."""


@dataclass(frozen=True, slots=True)
class DownloadRequest:
    url: str
    destination: Path
    expected_size_bytes: int
    expected_sha256: str | None = None
    expected_sha1: str | None = None
    expected_md5: str | None = None
    expected_generation: str | None = None
    attempts: int = 3
    timeout_seconds: int = 60

    def __post_init__(self) -> None:
        if not self.url.startswith("https://"):
            raise DownloadError("downloads require an HTTPS source")
        if not self.destination.is_absolute():
            raise DownloadError("download destination must be absolute")
        if not isinstance(self.expected_size_bytes, int) or self.expected_size_bytes <= 0:
            raise DownloadError("expected_size_bytes must be a positive integer")
        if (
            self.expected_sha256 is None
            and self.expected_sha1 is None
            and self.expected_md5 is None
            and self.expected_generation is None
        ):
            raise DownloadError(
                "a source SHA-256, SHA-1, MD5, or object generation is required"
            )
        if self.expected_sha256 is not None and not SHA256.fullmatch(
            self.expected_sha256
        ):
            raise DownloadError("expected_sha256 must be 64 lowercase hex characters")
        if self.expected_sha1 is not None and not SHA1.fullmatch(self.expected_sha1):
            raise DownloadError("expected_sha1 must be 40 lowercase hex characters")
        if self.expected_md5 is not None and not MD5.fullmatch(self.expected_md5):
            raise DownloadError("expected_md5 must be 32 lowercase hex characters")
        if self.expected_generation is not None and not self.expected_generation:
            raise DownloadError("expected_generation cannot be empty")
        if not isinstance(self.attempts, int) or self.attempts < 1:
            raise DownloadError("attempts must be at least one")
        if not isinstance(self.timeout_seconds, int) or self.timeout_seconds < 1:
            raise DownloadError("timeout_seconds must be at least one")


@dataclass(frozen=True, slots=True)
class DownloadResult:
    destination: Path
    sha256: str
    size_bytes: int
    source_generation: str | None
    attempts_used: int
    reused: bool
    md5: str | None = None
    sha1: str | None = None


def _checksums(path: Path) -> tuple[str, str, str, int]:
    sha256 = hashlib.sha256()
    sha1 = hashlib.sha1(usedforsecurity=False)
    md5 = hashlib.md5(usedforsecurity=False)
    size = 0
    with path.open("rb") as source:
        while chunk := source.read(1024 * 1024):
            sha256.update(chunk)
            sha1.update(chunk)
            md5.update(chunk)
            size += len(chunk)
    return sha256.hexdigest(), sha1.hexdigest(), md5.hexdigest(), size


def _validate_existing(request: DownloadRequest) -> DownloadResult | None:
    if not request.destination.exists():
        return None
    if not request.destination.is_file():
        raise CompletedResourceMismatch(
            f"published destination is not a file: {request.destination}"
        )
    digest, sha1, md5, size = _checksums(request.destination)
    if size != request.expected_size_bytes:
        raise CompletedResourceMismatch(
            f"published destination size mismatch: {request.destination}"
        )
    if (
        request.expected_sha256 is None
        and request.expected_sha1 is None
        and request.expected_md5 is None
    ):
        raise CompletedResourceMismatch(
            "an existing resource can only be reused with a declared checksum"
        )
    if request.expected_sha256 is not None and digest != request.expected_sha256:
        raise CompletedResourceMismatch(
            f"published destination checksum mismatch: {request.destination}"
        )
    if request.expected_sha1 is not None and sha1 != request.expected_sha1:
        raise CompletedResourceMismatch(
            f"published destination SHA-1 mismatch: {request.destination}"
        )
    if request.expected_md5 is not None and md5 != request.expected_md5:
        raise CompletedResourceMismatch(
            f"published destination MD5 mismatch: {request.destination}"
        )
    return DownloadResult(
        destination=request.destination,
        sha256=digest,
        size_bytes=size,
        source_generation=request.expected_generation,
        attempts_used=0,
        reused=True,
        md5=md5,
        sha1=sha1,
    )


def download_atomic(
    request: DownloadRequest,
    reference_root: str | Path,
    *,
    opener: Callable[..., Any] = urlopen,
    sleeper: Callable[[float], None] = time.sleep,
) -> DownloadResult:
    root = Path(reference_root).resolve()
    destination = request.destination
    destination.parent.mkdir(parents=True, exist_ok=True)
    if not destination.parent.resolve().is_relative_to(root):
        raise DownloadError("download destination must remain under reference_root")

    existing = _validate_existing(request)
    if existing is not None:
        return existing

    require_capacity(
        root,
        CapacityPlan(
            operation=f"download {destination.name}",
            source_bytes=request.expected_size_bytes,
        ),
    )
    partial_root = destination.parent / ".partial"
    partial_root.mkdir(mode=0o700, exist_ok=True)
    os.chmod(partial_root, 0o700)

    last_error: Exception | None = None
    for attempt in range(1, request.attempts + 1):
        partial = partial_root / f"{destination.name}.{uuid.uuid4().hex}.partial"
        generation: str | None = None
        try:
            http_request = Request(
                request.url,
                headers={"User-Agent": "genosapiens-reference-builder/0.1.0"},
            )
            digest = hashlib.sha256()
            sha1 = hashlib.sha1(usedforsecurity=False)
            md5 = hashlib.md5(usedforsecurity=False)
            size = 0
            with opener(http_request, timeout=request.timeout_seconds) as response:
                generation = response.headers.get("x-goog-generation")
                if (
                    request.expected_generation is not None
                    and generation != request.expected_generation
                ):
                    raise DownloadValidationError("source object generation mismatch")
                with partial.open("xb") as output:
                    os.chmod(partial, 0o600)
                    while chunk := response.read(1024 * 1024):
                        output.write(chunk)
                        digest.update(chunk)
                        sha1.update(chunk)
                        md5.update(chunk)
                        size += len(chunk)
                    output.flush()
                    os.fsync(output.fileno())
            checksum = digest.hexdigest()
            sha1_checksum = sha1.hexdigest()
            md5_checksum = md5.hexdigest()
            if size != request.expected_size_bytes:
                raise DownloadValidationError(
                    f"downloaded size {size} does not match {request.expected_size_bytes}"
                )
            if request.expected_sha256 is not None and checksum != request.expected_sha256:
                raise DownloadValidationError("downloaded SHA-256 mismatch")
            if request.expected_sha1 is not None and sha1_checksum != request.expected_sha1:
                raise DownloadValidationError("downloaded SHA-1 mismatch")
            if request.expected_md5 is not None and md5_checksum != request.expected_md5:
                raise DownloadValidationError("downloaded MD5 mismatch")
            os.replace(partial, destination)
            directory_fd = os.open(destination.parent, os.O_RDONLY)
            try:
                os.fsync(directory_fd)
            finally:
                os.close(directory_fd)
            return DownloadResult(
                destination=destination,
                sha256=checksum,
                size_bytes=size,
                source_generation=generation,
                attempts_used=attempt,
                reused=False,
                md5=md5_checksum,
                sha1=sha1_checksum,
            )
        except (DownloadValidationError, OSError, TimeoutError, URLError) as error:
            last_error = error
            if attempt < request.attempts:
                sleeper(min(2 ** (attempt - 1), 30))
        finally:
            try:
                partial.unlink(missing_ok=True)
            except OSError:
                pass
    raise DownloadError(
        f"download failed validation after {request.attempts} attempts: {last_error}"
    ) from last_error
