#!/usr/bin/env python3
"""Build and verify a lossless unprefixed public WGS VCF derivative."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import uuid


REFERENCES = {
    "challenge_alignment": (
        b"##reference=GCA_000001405.15_GRCh38_no_alt_analysis_set_"
        b"plus_hs38d1_maskedGRC_exclusions_v2_no_chr.fasta\n"
    ),
    "public_annotation_benchmark": (
        b"##reference=GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta\n"
    ),
}
SHA256 = re.compile(r"^[0-9a-f]{64}$")


class DerivativeError(RuntimeError):
    pass


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _canonical(value: object) -> bytes:
    return (
        json.dumps(
            value, sort_keys=True, separators=(",", ":"),
            ensure_ascii=True, allow_nan=False,
        )
        + "\n"
    ).encode("ascii")


def _transform(source, reference_line: bytes):
    inserted_reference = False
    header_complete = False
    records = 0
    source_records = hashlib.sha256()
    derived_records = hashlib.sha256()
    contigs = set()
    for line in source:
        if not line.endswith(b"\n") or b"\r" in line or b"\0" in line:
            raise DerivativeError("source VCF line encoding is invalid")
        if not header_complete:
            if line.startswith(b"##reference="):
                raise DerivativeError("source VCF already declares a reference")
            if line.startswith(b"##contig=<ID="):
                if not inserted_reference:
                    yield reference_line
                    inserted_reference = True
                prefix = b"##contig=<ID=chr"
                if not line.startswith(prefix):
                    raise DerivativeError("source header contig is not chr-prefixed")
                line = b"##contig=<ID=" + line[len(prefix):]
            elif line.startswith(b"#CHROM\t"):
                if not inserted_reference:
                    raise DerivativeError("source VCF has no contig declarations")
                header_complete = True
            elif not line.startswith(b"##"):
                raise DerivativeError("source VCF header is malformed")
            yield line
            continue
        fields = line.split(b"\t", 1)
        if len(fields) != 2 or not fields[0].startswith(b"chr"):
            raise DerivativeError("source VCF record contig is not chr-prefixed")
        mapped = fields[0][3:]
        if not mapped:
            raise DerivativeError("source VCF record contig is empty")
        derived = mapped + b"\t" + fields[1]
        records += 1
        contigs.add(mapped.decode("ascii"))
        source_records.update(line)
        derived_records.update(derived)
        yield derived
    if not header_complete or records == 0:
        raise DerivativeError("source VCF is incomplete")
    return {
        "record_count": records,
        "record_contigs": sorted(contigs),
        "source_record_stream_sha256": source_records.hexdigest(),
        "derived_record_stream_sha256": derived_records.hexdigest(),
    }


def build(
    source_path: Path,
    output_path: Path,
    manifest_path: Path,
    expected_source_sha256: str,
    expected_record_count: int,
    reference_role: str,
) -> None:
    if (
        not source_path.is_absolute()
        or not output_path.is_absolute()
        or not manifest_path.is_absolute()
        or output_path == manifest_path
        or output_path.exists()
        or manifest_path.exists()
        or output_path.is_symlink()
        or manifest_path.is_symlink()
        or SHA256.fullmatch(expected_source_sha256) is None
        or expected_record_count <= 0
        or reference_role not in REFERENCES
        or _sha256(source_path) != expected_source_sha256
    ):
        raise DerivativeError("derivative build binding is invalid")
    output_path.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
    temporary = output_path.parent / (".tmp-" + uuid.uuid4().hex + ".vcf.gz")
    descriptor = os.open(
        temporary,
        os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
        0o600,
    )
    process = None
    try:
        executable = next(
            (path for path in (Path("/usr/local/bin/bgzip"), Path("/usr/bin/bgzip"))
             if path.is_file()),
            None,
        )
        if executable is None:
            raise DerivativeError("pinned bgzip executable is unavailable")
        with os.fdopen(descriptor, "wb") as output:
            process = subprocess.Popen(
                (str(executable), "--threads", "1", "--stdout"),
                stdin=subprocess.PIPE,
                stdout=output,
                stderr=subprocess.PIPE,
                close_fds=True,
            )
            try:
                with gzip.open(source_path, "rb") as source:
                    transformed = _transform(source, REFERENCES[reference_role])
                    while True:
                        try:
                            process.stdin.write(next(transformed))
                        except StopIteration as stopped:
                            summary = stopped.value
                            break
                process.stdin.close()
                stderr = process.stderr.read(1024 * 1024 + 1)
                status = process.wait()
            except Exception:
                process.kill()
                process.wait()
                raise
            if status != 0 or stderr or summary["record_count"] != expected_record_count:
                raise DerivativeError("pinned bgzip derivative build failed")
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary, output_path)
        directory = os.open(output_path.parent, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
        output_sha256 = _sha256(output_path)
        implementation_sha256 = _sha256(Path(__file__))
        bgzip_version = subprocess.run(
            (str(executable), "--version"),
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        ).stdout.splitlines()[0]
        manifest = {
            "schema": "genosapiens.public-wgs-unprefixed-derivative",
            "version": 1,
            "classification": "public",
            "source": {
                "path": str(source_path),
                "sha256": expected_source_sha256,
                "size_bytes": source_path.stat().st_size,
                "record_count": expected_record_count,
            },
            "transformation": {
                "contig_rule": "remove-one-leading-chr-from-every-header-and-record-contig",
                "other_record_columns": "byte-identical",
                "reference_role": reference_role,
                "reference_header": REFERENCES[reference_role].decode("ascii").rstrip("\n"),
                "implementation_sha256": implementation_sha256,
                "bgzip_version": bgzip_version,
                "bgzip_threads": 1,
            },
            "closure": summary,
            "output": {
                "path": str(output_path),
                "sha256": output_sha256,
                "size_bytes": output_path.stat().st_size,
            },
        }
        content = _canonical(manifest)
        manifest_path.write_bytes(content)
        os.chmod(output_path, 0o444)
        os.chmod(manifest_path, 0o444)
    finally:
        if process is not None and process.poll() is None:
            process.kill()
            process.wait()
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def verify(source_path: Path, output_path: Path, manifest_path: Path) -> None:
    try:
        content = manifest_path.read_bytes()
        manifest = json.loads(content)
        if _canonical(manifest) != content:
            raise ValueError
        source = manifest["source"]
        output = manifest["output"]
        if (
            source["path"] != str(source_path)
            or output["path"] != str(output_path)
            or _sha256(source_path) != source["sha256"]
            or source_path.stat().st_size != source["size_bytes"]
            or _sha256(output_path) != output["sha256"]
            or output_path.stat().st_size != output["size_bytes"]
            or manifest["transformation"]["implementation_sha256"]
            != _sha256(Path(__file__))
        ):
            raise ValueError
        count = 0
        source_digest = hashlib.sha256()
        derived_digest = hashlib.sha256()
        with gzip.open(source_path, "rb") as original, gzip.open(output_path, "rb") as derived:
            role = manifest["transformation"]["reference_role"]
            transformed = _transform(original, REFERENCES[role])
            while True:
                try:
                    expected = next(transformed)
                except StopIteration as stopped:
                    summary = stopped.value
                    break
                actual = derived.readline()
                if actual != expected:
                    raise ValueError
                if not actual.startswith(b"#"):
                    count += 1
                    source_fields = expected.split(b"\t", 1)
                    source_line = b"chr" + source_fields[0] + b"\t" + source_fields[1]
                    source_digest.update(source_line)
                    derived_digest.update(actual)
            if derived.read(1):
                raise ValueError
        if (
            count != source["record_count"]
            or summary != manifest["closure"]
            or source_digest.hexdigest()
            != summary["source_record_stream_sha256"]
            or derived_digest.hexdigest()
            != summary["derived_record_stream_sha256"]
        ):
            raise ValueError
    except Exception:
        raise DerivativeError("unprefixed derivative verification failed") from None


def main() -> int:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    builder = subparsers.add_parser("build")
    for command in (builder, subparsers.add_parser("verify")):
        command.add_argument("source", type=Path)
        command.add_argument("output", type=Path)
        command.add_argument("manifest", type=Path)
    builder.add_argument("--source-sha256", required=True)
    builder.add_argument("--record-count", required=True, type=int)
    builder.add_argument(
        "--reference-role", required=True, choices=tuple(REFERENCES),
    )
    args = parser.parse_args()
    if args.command == "build":
        build(
            args.source, args.output, args.manifest,
            args.source_sha256, args.record_count, args.reference_role,
        )
    else:
        verify(args.source, args.output, args.manifest)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except DerivativeError as error:
        print(str(error), file=sys.stderr)
        raise SystemExit(1)
