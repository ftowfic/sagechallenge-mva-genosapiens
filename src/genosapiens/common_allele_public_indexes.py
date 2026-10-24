"""Canonical public clinical/exception projection and bounded SQLite lookup index."""

from __future__ import annotations

from contextlib import AbstractContextManager, contextmanager
from dataclasses import dataclass
import errno
import gzip
import hashlib
import json
import mmap
import os
from pathlib import Path
import re
import shutil
import sqlite3
import stat
import tempfile
from typing import Any, Iterable, Iterator, Literal, Mapping, Sequence
import uuid

from genosapiens import admission_links
from genosapiens.common_allele_exception_bundle import verify_exception_bundle
from genosapiens.bounded_wgs_config import PublicQueryLimits
from genosapiens.wgs_query_limits import QueryBudget


RECORD_SCHEMA = "genosapiens.common-allele-public-protection-record"
IDENTITY_SCHEMA = "genosapiens.common-allele-public-index-identity"
MANIFEST_SCHEMA = "genosapiens.common-allele-public-index-manifest"
SHA256 = re.compile(r"^[0-9a-f]{64}$")
DNA = re.compile(r"^[ACGT]+$")
CHROMOSOMES = {str(value) for value in range(1, 23)} | {"X", "Y", "MT"}
CACHE_KIB = 65536
MAX_ROW_BYTES = 8 * 1024 * 1024
CURSOR_BATCH_ROWS = 1024
QUERY_CTE_ROWS = 64
CACHE_DROP_WINDOW_BYTES = 64 * 1024 * 1024
RETENTION_SCOPE = 'shared-retention-evidence-v1'
RETENTION_EXACT_PROJECTION_VERSION = 'retention-exact-records/v1'


def retention_disease_ids(row):
    """Index explicit disease identifiers without interpreting disease names.

    This is a retrieval key, not a pathogenicity decision. The shared selector
    evaluates classifications and functional assertions after retrieval.
    """
    if row.get('record_type') != 'exact':
        return ()
    evidence = row.get('evidence', {})
    if row['source'] == 'clinvar':
        values = (evidence.get('condition_ids'),)
    elif row['source'] == 'open_targets':
        values = (evidence.get('disease_from_source_id'),)
    elif row['source'] == 'curated_functional':
        values = evidence['functional']['disease_ids']
    else:
        return ()
    result = set()
    for value in values:
        if value is None:
            continue
        if not isinstance(value, str):
            _fail()
        for namespace, number in re.findall(
                r'(?<![A-Za-z0-9_])(MONDO|OMIM|ORPHA|Orphanet)[:_]([0-9]+)(?![0-9A-Za-z])', value):
            result.add(('ORPHA' if namespace == 'Orphanet' else namespace) + ':' + number)
    return tuple(sorted(result))


def retention_coverage_record(*, source_bindings, curated_counts):
    """Public preparation closes all configured sources, including empty ones."""
    from genosapiens.curated_alleles import EXTRACTION_POLICY
    if (set(source_bindings) != {'clinvar', 'open_targets', 'curated', 'reference', 'accessions'}
            or any(not SHA256.fullmatch(value) for value in source_bindings.values())
            or set(curated_counts) != {'source_rows', 'assertions', 'gaps'}
            or any(type(v) is not int or v < 0 for v in curated_counts.values())):
        _fail()
    return dict(schema=RECORD_SCHEMA, version=1, record_type='scope', source='retention_coverage',
        scope=RETENTION_SCOPE, clinical_fields_version=2, disease_lookup_version=1,
        extraction_policy=EXTRACTION_POLICY,
        source_bindings=source_bindings, curated_counts=curated_counts)


class CommonAlleleIndexError(RuntimeError):
    """Static failure for a public protection projection/index."""


def _fail() -> None:
    raise CommonAlleleIndexError("common-allele public index operation failed")


def _canonical(value: object) -> bytes:
    try:
        content = (json.dumps(value, sort_keys=True, separators=(",", ":"),
                              ensure_ascii=True, allow_nan=False) + "\n").encode("ascii")
    except (TypeError, ValueError, UnicodeError, RecursionError, OverflowError):
        _fail()
    if len(content) > MAX_ROW_BYTES:
        _fail()
    return content


def _pairs(items: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in items:
        if key in result:
            _fail()
        result[key] = value
    return result


def _json(content: bytes, *, canonical: bool = False) -> dict[str, Any]:
    try:
        value = json.loads(content.decode("utf-8"), object_pairs_hook=_pairs,
                           parse_constant=lambda _value: _fail())
    except CommonAlleleIndexError:
        raise
    except (ValueError, TypeError, UnicodeError, RecursionError):
        _fail()
    if type(value) is not dict or canonical and _canonical(value) != content:
        _fail()
    return value


def _sha(value: object) -> str:
    if type(value) is not str or SHA256.fullmatch(value) is None:
        _fail()
    return value


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    descriptor = None
    offset = 0
    released = 0
    try:
        descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
        if not stat.S_ISREG(os.fstat(descriptor).st_mode):
            _fail()
        os.posix_fadvise(descriptor, 0, 0, os.POSIX_FADV_SEQUENTIAL)
        while chunk := os.read(descriptor, 1024 * 1024):
            digest.update(chunk)
            offset += len(chunk)
            if offset - released >= CACHE_DROP_WINDOW_BYTES:
                os.posix_fadvise(
                    descriptor,
                    released,
                    offset - released,
                    os.POSIX_FADV_DONTNEED,
                )
                released = offset
        os.posix_fadvise(descriptor, released, 0, os.POSIX_FADV_DONTNEED)
        released = offset
    except CommonAlleleIndexError:
        raise
    except (AttributeError, OSError):
        _fail()
    finally:
        if descriptor is not None:
            try:
                if released < offset:
                    os.posix_fadvise(
                        descriptor, released, 0, os.POSIX_FADV_DONTNEED
                    )
            except (AttributeError, OSError):
                pass
            os.close(descriptor)
    return digest.hexdigest()


def _discard_file_cache(descriptor: int, *, required: bool = True) -> None:
    """Release clean pages after one complete authenticated SQLite pass."""

    try:
        os.posix_fadvise(descriptor, 0, 0, os.POSIX_FADV_DONTNEED)
    except (AttributeError, OSError):
        if required:
            _fail()


def _regular(root: Path, relative: str, sha256: str, size: int | None = None,
             *, trust_published: bool = False) -> Path:
    path = root / relative
    try:
        metadata = path.lstat()
        if (type(trust_published) is not bool
                or trust_published and size is None
                or path.resolve(strict=True) != path
                or not path.is_relative_to(root)
                or not stat.S_ISREG(metadata.st_mode)
                or stat.S_IMODE(metadata.st_mode) & 0o022
                or size is not None and metadata.st_size != size
                or not trust_published and _file_sha256(path) != sha256):
            _fail()
    except (OSError, ValueError):
        _fail()
    return path


def _manifest_outputs(root: Path, relative: str, sha256: str) -> tuple[dict[str, Any], dict[str, dict[str, Any]]]:
    value = _json(_regular(root, relative, sha256).read_bytes())
    outputs = value.get("outputs")
    if value.get("classification") != "public" or type(outputs) is not list:
        _fail()
    by_path = {}
    for row in outputs:
        if (type(row) is not dict or set(row) != {"path", "retention", "sha256", "size_bytes"}
                or type(row["path"]) is not str or row["path"] in by_path):
            _fail()
        _sha(row["sha256"])
        if type(row["size_bytes"]) is not int or row["size_bytes"] <= 0:
            _fail()
        by_path[row["path"]] = row
    return value, by_path


@dataclass(frozen=True, slots=True)
class SequenceAlias:
    chromosome: str
    position: int
    reference: str
    alternate: str


@dataclass(frozen=True, slots=True)
class PublicIndexQuery:
    chromosome: str
    position: int
    reference: str
    alternate: str
    normalized_alias: SequenceAlias | None = None
    source_span: tuple[int, int] | None = None
    normalized_span: tuple[int, int] | None = None
    native_scope: str | None = None


@dataclass(frozen=True, slots=True)
class BuiltPublicIndexes:
    index_sha256: str
    projection_path: str
    projection_sha256: str
    sqlite_path: str
    sqlite_sha256: str
    manifest_path: str
    manifest_sha256: str
    record_count: int


@dataclass(frozen=True, slots=True)
class PublicIndexQueryResult:
    status: Literal["matches", "no_match", "retain_unknown", "unavailable"]
    records: tuple[Mapping[str, Any], ...]
    reasons: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class VerifiedPublicIndexes:
    index_sha256: str
    manifest_sha256: str
    closure_sha256: str
    projection_path: Path
    sqlite_path: Path
    record_count: int


class _MmapFasta(AbstractContextManager["_MmapFasta"]):
    def __init__(self, fasta: Path, fai: Path) -> None:
        self._stream = fasta.open("rb")
        self._map = mmap.mmap(self._stream.fileno(), 0, access=mmap.ACCESS_READ)
        self.index: dict[str, tuple[int, int, int, int]] = {}
        for line in fai.read_text(encoding="ascii").splitlines():
            fields = line.split("\t")
            if len(fields) < 5:
                _fail()
            name = fields[0]
            values = tuple(int(value) for value in fields[1:5])
            if name in self.index or min(values) <= 0:
                _fail()
            self.index[name] = values

    def fetch(self, chromosome: str, position: int, length: int) -> str:
        contig = "chrM" if chromosome == "MT" else "chr" + chromosome
        try:
            contig_length, offset, line_bases, line_width = self.index[contig]
        except KeyError:
            raise KeyError from None
        if position <= 0 or length <= 0 or position + length - 1 > contig_length:
            raise IndexError
        zero = position - 1
        chunks = []
        remaining = length
        while remaining:
            line, within = divmod(zero, line_bases)
            count = min(remaining, line_bases - within)
            start = offset + line * line_width + within
            chunks.append(self._map[start:start + count])
            zero += count
            remaining -= count
        return b"".join(chunks).decode("ascii")

    def __exit__(self, *_args: object) -> None:
        self._map.close()
        self._stream.close()


def normalize_public_alias(
    chromosome: str, position: int, reference: str, alternate: str,
    fetch_reference: Any,
) -> SequenceAlias:
    """Reference-check and left-normalize one publisher-native sequence allele."""

    if (chromosome not in CHROMOSOMES or type(position) is not int or position <= 0
            or type(reference) is not str or type(alternate) is not str
            or DNA.fullmatch(reference) is None or DNA.fullmatch(alternate) is None
            or reference == alternate):
        _fail()
    try:
        if fetch_reference(chromosome, position, len(reference)) != reference:
            _fail()
        while reference[-1] == alternate[-1]:
            if len(reference) > 1 and len(alternate) > 1:
                reference, alternate = reference[:-1], alternate[:-1]
            elif position > 1:
                preceding = fetch_reference(chromosome, position - 1, 1)
                reference, alternate = preceding + reference, preceding + alternate
                position -= 1
            else:
                break
        while len(reference) > 1 and len(alternate) > 1 and reference[0] == alternate[0]:
            reference, alternate, position = reference[1:], alternate[1:], position + 1
        if fetch_reference(chromosome, position, len(reference)) != reference:
            _fail()
        return SequenceAlias(chromosome, position, reference, alternate)
    except CommonAlleleIndexError:
        raise
    except (KeyError, IndexError, OSError, UnicodeError, ValueError, TypeError):
        _fail()


_DDL = """
PRAGMA foreign_keys=ON;
CREATE TABLE binding (closure_sha256 TEXT PRIMARY KEY NOT NULL);
CREATE TABLE records (
  sequence INTEGER PRIMARY KEY NOT NULL,
  record_type TEXT NOT NULL CHECK(record_type IN ('exact','interval','unknown','scope')),
  canonical_sha256 TEXT UNIQUE NOT NULL,
  canonical_row BLOB NOT NULL
);
CREATE TABLE exact_records (
  sequence INTEGER PRIMARY KEY REFERENCES records(sequence), source TEXT NOT NULL,
  chromosome TEXT NOT NULL, position INTEGER NOT NULL, reference TEXT NOT NULL, alternate TEXT NOT NULL,
  alias_chromosome TEXT NOT NULL, alias_position INTEGER NOT NULL,
  alias_reference TEXT NOT NULL, alias_alternate TEXT NOT NULL
);
CREATE INDEX exact_source_key ON exact_records(source,chromosome,position,reference,alternate);
CREATE INDEX exact_alias_key ON exact_records(source,alias_chromosome,alias_position,alias_reference,alias_alternate);
CREATE TABLE retention_diseases (
  disease_id TEXT NOT NULL, sequence INTEGER NOT NULL REFERENCES records(sequence),
  PRIMARY KEY(disease_id,sequence)
) WITHOUT ROWID;
CREATE INDEX retention_disease_sequence ON retention_diseases(sequence,disease_id);
CREATE TABLE intervals (
  sequence INTEGER PRIMARY KEY REFERENCES records(sequence), source TEXT NOT NULL,
  chromosome TEXT NOT NULL, start INTEGER NOT NULL, end INTEGER NOT NULL
);
CREATE INDEX interval_key ON intervals(source,chromosome,start,end);
CREATE TABLE unknown_loci (
  sequence INTEGER PRIMARY KEY REFERENCES records(sequence), source TEXT NOT NULL,
  chromosome TEXT, start INTEGER, end INTEGER
);
CREATE INDEX unknown_key ON unknown_loci(source,chromosome,start,end);
CREATE TABLE scopes (
  sequence INTEGER PRIMARY KEY REFERENCES records(sequence), source TEXT NOT NULL, scope TEXT NOT NULL
);
"""


class _ProjectionWriter(AbstractContextManager["_ProjectionWriter"]):
    def __init__(self, projection: Path, database: Path, closure_sha256: str) -> None:
        self.projection_path = projection
        self.database_path = database
        self._raw = projection.open("xb")
        self._gzip = gzip.GzipFile(filename="", fileobj=self._raw, mode="wb", mtime=0, compresslevel=9)
        self.connection = sqlite3.connect(database)
        self.connection.execute(f"PRAGMA cache_size=-{CACHE_KIB}")
        self.connection.execute("PRAGMA temp_store=FILE")
        self.connection.execute("PRAGMA journal_mode=DELETE")
        self.connection.executescript(_DDL)
        self.connection.execute("INSERT INTO binding VALUES (?)", (closure_sha256,))
        self.uncompressed_sha256 = hashlib.sha256()
        self.uncompressed_size = 0
        self.counts = {"exact": 0, "interval": 0, "unknown": 0, "scope": 0}
        self.source_counts: dict[str, int] = {}
        self.sequence = 0

    def emit(self, row: Mapping[str, Any]) -> None:
        content = _canonical(row)
        if row.get("schema") != RECORD_SCHEMA or row.get("version") != 1:
            _fail()
        record_type, source = row.get("record_type"), row.get("source")
        if record_type not in self.counts or type(source) is not str:
            _fail()
        self.sequence += 1
        self._gzip.write(content)
        self.uncompressed_sha256.update(content)
        self.uncompressed_size += len(content)
        row_sha = hashlib.sha256(content).hexdigest()
        self.connection.execute("INSERT INTO records VALUES (?,?,?,?)",
                                (self.sequence, record_type, row_sha, content))
        if record_type == "exact":
            allele, alias = row["allele"], row["normalized_alias"]
            self.connection.execute("INSERT INTO exact_records VALUES (?,?,?,?,?,?,?,?,?,?)", (
                self.sequence, source, allele["chromosome"], allele["position"],
                allele["reference"], allele["alternate"], alias["chromosome"],
                alias["position"], alias["reference"], alias["alternate"],
            ))
            self.connection.executemany('INSERT INTO retention_diseases VALUES (?,?)',
                ((disease, self.sequence) for disease in retention_disease_ids(row)))
        elif record_type == "interval":
            self.connection.execute("INSERT INTO intervals VALUES (?,?,?,?,?)", (
                self.sequence, source, row["chromosome"], row["start"], row["end"],
            ))
        elif record_type == "unknown":
            locus = row["affected_locus"]
            self.connection.execute("INSERT INTO unknown_loci VALUES (?,?,?,?,?)", (
                self.sequence, source, locus.get("chromosome"), locus.get("start"), locus.get("end"),
            ))
        else:
            self.connection.execute("INSERT INTO scopes VALUES (?,?,?)",
                                    (self.sequence, source, row["scope"]))
        self.counts[record_type] += 1
        self.source_counts[source] = self.source_counts.get(source, 0) + 1
        if self.sequence % 100_000 == 0:
            self.connection.commit()

    def __exit__(self, exc_type: object, _exc: object, _tb: object) -> None:
        try:
            if exc_type is None:
                self.connection.commit()
                self.connection.execute("ANALYZE")
                self.connection.execute("PRAGMA optimize")
                self.connection.commit()
        finally:
            self.connection.close()
            self._gzip.close()
            self._raw.close()


def _unknown(
    source: str, reason: str, locator: Mapping[str, Any],
    chromosome: str | None = None, start: int | None = None, end: int | None = None,
) -> dict[str, Any]:
    return {"schema": RECORD_SCHEMA, "version": 1, "record_type": "unknown",
            "source": source, "reason": reason, "source_locator": locator,
            "affected_locus": {"chromosome": chromosome, "start": start, "end": end}}


def _exact(
    source: str, chromosome: str, position: int, reference: str, alternate: str,
    alias: SequenceAlias, locator: Mapping[str, Any], evidence: Mapping[str, Any],
    *, source_records: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    return {"schema": RECORD_SCHEMA, "version": 1, "record_type": "exact",
            "source": source,
            "allele": {"chromosome": chromosome, "position": position,
                       "reference": reference, "alternate": alternate},
            "normalized_alias": {"chromosome": alias.chromosome, "position": alias.position,
                                 "reference": alias.reference, "alternate": alias.alternate},
            "source_locator": locator, "evidence": evidence,
            **({} if source_records is None else {"source_records": source_records})}


def _variant_id_locus(value: object) -> tuple[str | None, int | None, int | None]:
    if type(value) is not str:
        return None, None, None
    fields = value.split("_", 3)
    if len(fields) != 4 or not fields[1].isdecimal() or not fields[2]:
        return None, None, None
    chromosome, position, reference = fields[0], int(fields[1]), fields[2]
    if chromosome not in CHROMOSOMES or position <= 0:
        return None, None, None
    return chromosome, position, position + len(reference) - 1


def _publish_file(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if source.is_symlink() or destination.is_symlink() or destination.parent.is_symlink():
        _fail()
    temporary = None
    try:
        try:
            os.link(source, destination)
        except OSError as error:
            if error.errno != errno.EXDEV:
                raise
            # Separate Docker bind mounts can reject a hard link even when
            # their host paths share a filesystem. Copy in bounded chunks to
            # the destination mount, then commit with a local no-clobber link.
            with source.open("rb") as incoming, tempfile.NamedTemporaryFile(
                dir=destination.parent, prefix=".publish-", delete=False,
            ) as outgoing:
                temporary = Path(outgoing.name)
                shutil.copyfileobj(incoming, outgoing, 1024 * 1024)
                outgoing.flush()
                os.fsync(outgoing.fileno())
            if (temporary.stat().st_size != source.stat().st_size
                    or _file_sha256(temporary) != _file_sha256(source)):
                _fail()
            temporary.chmod(0o444)
            os.link(temporary, destination)
    except FileExistsError:
        if (destination.is_symlink() or not destination.is_file()
                or destination.stat().st_size != source.stat().st_size
                or _file_sha256(destination) != _file_sha256(source)):
            _fail()
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    os.chmod(destination, 0o444)
    with destination.open("rb") as published:
        os.fsync(published.fileno())
    directory = os.open(destination.parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(directory)
    finally:
        os.close(directory)


def _reference_files(
    references: Path, manifest_path: str, manifest_sha256: str,
) -> tuple[Path, Path, dict[str, str]]:
    _manifest, outputs = _manifest_outputs(references, manifest_path, manifest_sha256)
    fasta_rows = [row for path, row in outputs.items()
                  if path.endswith("/public-native.no-alt-plus-hs38d1.chr.fna")]
    fai_rows = [row for path, row in outputs.items()
                if path.endswith("/public-native.no-alt-plus-hs38d1.chr.fna.fai")]
    if len(fasta_rows) != 1 or len(fai_rows) != 1:
        _fail()
    fasta_row, fai_row = fasta_rows[0], fai_rows[0]
    fasta = _regular(references, fasta_row["path"], fasta_row["sha256"], fasta_row["size_bytes"])
    fai = _regular(references, fai_row["path"], fai_row["sha256"], fai_row["size_bytes"])
    return fasta, fai, {"fasta": fasta_row["sha256"], "fai": fai_row["sha256"]}


def _contract_artifacts(
    repository: Path, contract_path: str, expected_sha256: str,
    roles: set[str],
) -> tuple[dict[str, Any], dict[str, dict[str, Any]]]:
    content = _contract_bytes(repository, contract_path, expected_sha256)
    contract = _json(content)
    rows = contract.get("artifacts")
    if type(rows) is not list:
        _fail()
    artifacts = {row.get("role"): row for row in rows if type(row) is dict}
    if set(artifacts) != roles:
        _fail()
    return contract, artifacts


def _contract_bytes(repository: Path, relative: str, expected_sha256: str) -> bytes:
    """Bind repository configuration bytes without imposing resource file modes."""
    from genosapiens.bounded_json import read_bounded_bytes, BoundedJSONError
    try:
        content = read_bounded_bytes(repository / relative, reference_root=repository)
    except BoundedJSONError as error:
        raise CommonAlleleIndexError(f'Cannot read public configuration: {relative}') from error
    if hashlib.sha256(content).hexdigest() != _sha(expected_sha256):
        raise CommonAlleleIndexError(f'Public configuration checksum differs: {relative}')
    return content


@contextmanager
def _preparation_stage(name, progress):
    if progress is not None:
        progress(name)
    try:
        yield
    except (CommonAlleleIndexError, OSError, ValueError, RuntimeError) as error:
        raise CommonAlleleIndexError(f'Public index preparation failed during {name}: {error}') from error


def _projection_clinvar(
    writer: _ProjectionWriter, vcf: Path, artifact: Mapping[str, Any],
    reference: _MmapFasta,
) -> None:
    record_number = 0
    try:
        with gzip.open(vcf, "rb") as stream:
            for raw in stream:
                if len(raw) > MAX_ROW_BYTES:
                    _fail()
                if raw.startswith(b"#"):
                    continue
                record_number += 1
                locator = {"path": artifact["path"], "record_number": record_number,
                           "sha256": artifact["sha256"]}
                try:
                    text = raw.decode("utf-8")
                    parsed = admission_links.parse_clinvar_vcf_row(text)
                    columns, info = parsed["columns"], parsed["info"]
                    chrom, position, ref, alt = columns[0], int(columns[1]), columns[3], columns[4]
                    if "," in alt:
                        raise ValueError
                    alias = normalize_public_alias(chrom, position, ref, alt, reference.fetch)
                    writer.emit(_exact(
                        "clinvar", chrom, position, ref, alt, alias, locator,
                        {"variation_id": columns[2], "allele_id": info.get("ALLELEID"),
                         "clinical_significance": info.get("CLNSIG"),
                         "conflicting_significance": info.get("CLNSIGCONF"),
                         "review_status": info.get("CLNREVSTAT"),
                         "condition_ids": info.get("CLNDISDB"),
                         "condition_names": info.get("CLNDN"),
                         "gene_info": info.get("GENEINFO"),
                         "clinical_sources": info.get("CLNVI")},
                        source_records={"clinvar": parsed},
                    ))
                except Exception:
                    try:
                        fields = raw.decode("utf-8").rstrip("\r\n").split("\t")
                        chrom = fields[0] if fields[0] in CHROMOSOMES else None
                        start = int(fields[1]) if chrom is not None and fields[1].isdecimal() else None
                        end = start + len(fields[3]) - 1 if start is not None and len(fields) > 3 else None
                    except Exception:
                        chrom = start = end = None
                    writer.emit(_unknown("clinvar", "unassessable_source_record", locator,
                                         chrom, start, end))
    except (OSError, gzip.BadGzipFile, UnicodeError):
        _fail()


def _projection_open_targets(
    writer: _ProjectionWriter, variant_paths: Sequence[Path], evidence_paths: Sequence[Path],
    output_by_absolute: Mapping[str, Mapping[str, Any]], reference: _MmapFasta,
    scratch: Path,
) -> None:
    try:
        import duckdb
        if duckdb.__version__ != "1.4.2":
            _fail()
        temp_directory = scratch / "duckdb"
        temp_directory.mkdir()
        connection = duckdb.connect(database=":memory:", config={
            "threads": "2", "memory_limit": "4GB", "preserve_insertion_order": "false",
            "temp_directory": str(temp_directory), "autoinstall_known_extensions": "false",
            "autoload_known_extensions": "false",
        })
        connection.execute("SET enable_progress_bar=false")
        connection.execute("""
          CREATE TEMP TABLE variants AS
          SELECT variantId, chromosome, position, referenceAllele, alternateAllele,
                 filename AS variant_filename, file_row_number AS variant_row
          FROM read_parquet(?, filename=true, file_row_number=true)
        """, [[str(path) for path in variant_paths]])
        counts = connection.execute(
            "SELECT count(*),count(DISTINCT variantId),count(*)-count(variantId) FROM variants"
        ).fetchone()
        if counts != (7_538_243, 7_538_243, 0):
            _fail()
        connection.execute("CREATE UNIQUE INDEX variant_id_index ON variants(variantId)")
        eva_fields = tuple(admission_links._EVA_TEXT) + tuple(admission_links._EVA_LISTS) + ("score",)
        select_eva = ",".join(f'e."{name}"' for name in eva_fields)
        for evidence_path in evidence_paths:
            cursor = connection.execute(
                f"SELECT {select_eva}, e.file_row_number AS evidence_row, "
                "v.chromosome, v.position, v.referenceAllele, v.alternateAllele, "
                "v.variant_filename, v.variant_row "
                "FROM read_parquet(?, file_row_number=true) e "
                "LEFT JOIN variants v USING(variantId) ORDER BY e.file_row_number",
                [str(evidence_path)],
            )
            names = tuple(item[0] for item in cursor.description)
            evidence_output = output_by_absolute[str(evidence_path)]
            while batch := cursor.fetchmany(CURSOR_BATCH_ROWS):
                for values in batch:
                    joined = dict(zip(names, values, strict=True))
                    evidence = {name: joined[name] for name in eva_fields}
                    locator = {"path": evidence_output["path"],
                               "row_number": joined["evidence_row"],
                               "sha256": evidence_output["sha256"],
                               "assertion_id": evidence.get("id")}
                    try:
                        admission_links.parse_open_targets_eva_row(evidence)
                        variant = {name: joined[name] for name in (
                            "variantId", "chromosome", "position", "referenceAllele", "alternateAllele"
                        )}
                        chrom, position, ref, alt = admission_links.parse_open_targets_variant_row(variant)
                        alias = normalize_public_alias(chrom, position, ref, alt, reference.fetch)
                        variant_output = output_by_absolute[joined["variant_filename"]]
                        locator["variant"] = {"path": variant_output["path"],
                                              "row_number": joined["variant_row"],
                                              "sha256": variant_output["sha256"],
                                              "variant_id": evidence["variantId"]}
                        writer.emit(_exact(
                            "open_targets", chrom, position, ref, alt, alias, locator,
                            {"assertion_id": evidence["id"], "variant_id": evidence["variantId"],
                             "clinical_significances": evidence["clinicalSignificances"],
                             "confidence": evidence["confidence"], "target_id": evidence["targetId"],
                             "disease_from_source_id": evidence["diseaseFromSourceId"]},
                            source_records={"evidence": evidence, "variant": variant},
                        ))
                    except Exception:
                        chrom, start, end = _variant_id_locus(evidence.get("variantId"))
                        writer.emit(_unknown(
                            "open_targets", "unassessable_or_unjoined_assertion",
                            locator, chrom, start, end,
                        ))
        connection.close()
    except CommonAlleleIndexError:
        raise
    except Exception as error:
        raise CommonAlleleIndexError("common-allele Open Targets projection failed") from error


def _projection_exceptions(
    writer: _ProjectionWriter, records: Sequence[Mapping[str, Any]], reference: _MmapFasta,
    manifest_path: str, manifest_sha256: str,
) -> None:
    for source in records:
        locator = {"path": manifest_path, "sha256": manifest_sha256,
                   "record_id": source["record_id"]}
        if source["record_type"] == "exact_allele":
            alias = normalize_public_alias(
                source["chromosome"], source["position"], source["reference"],
                source["alternate"], reference.fetch,
            )
            writer.emit(_exact(
                "common_allele_exceptions", source["chromosome"], source["position"],
                source["reference"], source["alternate"], alias, locator,
                {"exception_record_id": source["record_id"],
                 "protection_reason": source["protection_reason"],
                 "curated_vus_exception": source["curated_vus_exception"]},
            ))
        elif source["record_type"] == "gene_interval":
            writer.emit({"schema": RECORD_SCHEMA, "version": 1, "record_type": "interval",
                         "source": "common_allele_exceptions", "chromosome": source["chromosome"],
                         "start": source["start"], "end": source["end"],
                         "source_locator": locator, "evidence": {
                             "exception_record_id": source["record_id"],
                             "gene_symbol": source["gene_symbol"], "target_id": source["target_id"],
                             "protection_reason": source["protection_reason"]}})
        elif source["record_type"] == "native_scope":
            writer.emit({"schema": RECORD_SCHEMA, "version": 1, "record_type": "scope",
                         "source": "common_allele_exceptions", "scope": source["scope"],
                         "source_locator": locator, "evidence": {
                             "exception_record_id": source["record_id"],
                             "protection_reason": source["protection_reason"]}})
        else:
            _fail()


def reconcile_projection_index(
    projection_path: str | Path, sqlite_path: str | Path, closure_sha256: str,
) -> tuple[int, str]:
    """Compare canonical rows and every query-table semantic field."""

    _sha(closure_sha256)
    projection = Path(projection_path)
    database = Path(sqlite_path)
    connection = None
    count = 0
    digest = hashlib.sha256()
    try:
        descriptor = os.open(database, os.O_RDONLY | os.O_NOFOLLOW)
        connection = sqlite3.connect(
            f"file:/proc/self/fd/{descriptor}?mode=ro&immutable=1", uri=True
        )
        connection.execute(f"PRAGMA cache_size=-{CACHE_KIB}")
        connection.execute("PRAGMA temp_store=FILE")
        connection.execute("PRAGMA trusted_schema=OFF")
        if connection.execute("PRAGMA integrity_check").fetchone() != ("ok",):
            _fail()
        _discard_file_cache(descriptor)
        if connection.execute("PRAGMA foreign_key_check").fetchone() is not None:
            _fail()
        if connection.execute("SELECT * FROM binding").fetchall() != [(closure_sha256,)]:
            _fail()
        indexed_cursor = connection.execute(
            "SELECT sequence,record_type,canonical_sha256,canonical_row FROM records ORDER BY sequence"
        )
        indexed = iter(indexed_cursor)
        projection_descriptor = os.open(projection, os.O_RDONLY | os.O_NOFOLLOW)
        try:
            os.posix_fadvise(
                projection_descriptor, 0, 0, os.POSIX_FADV_SEQUENTIAL
            )
            with os.fdopen(os.dup(projection_descriptor), "rb") as compressed:
                with gzip.GzipFile(fileobj=compressed, mode="rb") as stream:
                    for line in stream:
                        if len(line) > MAX_ROW_BYTES or not line.endswith(b"\n"):
                            _fail()
                        row = _json(line, canonical=True)
                        count += 1
                        digest.update(line)
                        expected = (
                            count,
                            row.get("record_type"),
                            hashlib.sha256(line).hexdigest(),
                            line,
                        )
                        if next(indexed, None) != expected:
                            _fail()
                        actual_diseases = tuple(value[0] for value in connection.execute(
                            'SELECT disease_id FROM retention_diseases WHERE sequence=? ORDER BY disease_id',
                            (count,)))
                        if actual_diseases != retention_disease_ids(row):
                            _fail()
            if next(indexed, None) is not None:
                _fail()
            _discard_file_cache(descriptor)
            _discard_file_cache(projection_descriptor)
        finally:
            try:
                indexed_cursor.close()
            finally:
                _discard_file_cache(descriptor, required=False)
                _discard_file_cache(projection_descriptor, required=False)
                os.close(projection_descriptor)
        typed_count = sum(connection.execute(f"SELECT count(*) FROM {table}").fetchone()[0]
                          for table in ("exact_records", "intervals", "unknown_loci", "scopes"))
        _discard_file_cache(descriptor)
        if typed_count != count:
            _fail()
        for table in ("exact_records", "intervals", "unknown_loci", "scopes"):
            if connection.execute(
                f"SELECT 1 FROM {table} WHERE source NOT IN "
                "('clinvar','open_targets','common_allele_exceptions','curated_functional',"
                "'curated_curation_gap','retention_coverage') LIMIT 1"
            ).fetchone() is not None:
                _fail()
        semantic_queries = (
            """SELECT count(*) FROM records r LEFT JOIN exact_records e USING(sequence)
               WHERE r.record_type='exact' AND (e.sequence IS NULL OR
                 e.source != json_extract(CAST(r.canonical_row AS TEXT),'$.source') OR
                 e.chromosome != json_extract(CAST(r.canonical_row AS TEXT),'$.allele.chromosome') OR
                 e.position != json_extract(CAST(r.canonical_row AS TEXT),'$.allele.position') OR
                 e.reference != json_extract(CAST(r.canonical_row AS TEXT),'$.allele.reference') OR
                 e.alternate != json_extract(CAST(r.canonical_row AS TEXT),'$.allele.alternate') OR
                 e.alias_chromosome != json_extract(CAST(r.canonical_row AS TEXT),'$.normalized_alias.chromosome') OR
                 e.alias_position != json_extract(CAST(r.canonical_row AS TEXT),'$.normalized_alias.position') OR
                 e.alias_reference != json_extract(CAST(r.canonical_row AS TEXT),'$.normalized_alias.reference') OR
                 e.alias_alternate != json_extract(CAST(r.canonical_row AS TEXT),'$.normalized_alias.alternate'))""",
            """SELECT count(*) FROM records r LEFT JOIN intervals i USING(sequence)
               WHERE r.record_type='interval' AND (i.sequence IS NULL OR
                 i.source != json_extract(CAST(r.canonical_row AS TEXT),'$.source') OR
                 i.chromosome != json_extract(CAST(r.canonical_row AS TEXT),'$.chromosome') OR
                 i.start != json_extract(CAST(r.canonical_row AS TEXT),'$.start') OR
                 i.end != json_extract(CAST(r.canonical_row AS TEXT),'$.end'))""",
            """SELECT count(*) FROM records r LEFT JOIN unknown_loci u USING(sequence)
               WHERE r.record_type='unknown' AND (u.sequence IS NULL OR
                 u.source != json_extract(CAST(r.canonical_row AS TEXT),'$.source') OR
                 u.chromosome IS NOT json_extract(CAST(r.canonical_row AS TEXT),'$.affected_locus.chromosome') OR
                 u.start IS NOT json_extract(CAST(r.canonical_row AS TEXT),'$.affected_locus.start') OR
                 u.end IS NOT json_extract(CAST(r.canonical_row AS TEXT),'$.affected_locus.end'))""",
            """SELECT count(*) FROM records r LEFT JOIN scopes s USING(sequence)
               WHERE r.record_type='scope' AND (s.sequence IS NULL OR
                 s.source != json_extract(CAST(r.canonical_row AS TEXT),'$.source') OR
                 s.scope != json_extract(CAST(r.canonical_row AS TEXT),'$.scope'))""",
            """SELECT count(*) FROM records r LEFT JOIN exact_records e USING(sequence)
                 LEFT JOIN intervals i USING(sequence) LEFT JOIN unknown_loci u USING(sequence)
                 LEFT JOIN scopes s USING(sequence)
               WHERE (e.sequence IS NOT NULL)+(i.sequence IS NOT NULL)+
                     (u.sequence IS NOT NULL)+(s.sequence IS NOT NULL) != 1""",
        )
        for query in semantic_queries:
            if connection.execute(query).fetchone() != (0,):
                _fail()
            _discard_file_cache(descriptor)
        return count, digest.hexdigest()
    except CommonAlleleIndexError:
        raise
    except (OSError, sqlite3.Error, gzip.BadGzipFile, EOFError):
        _fail()
    finally:
        if connection is not None:
            connection.close()
        try:
            os.close(descriptor)
        except (UnboundLocalError, OSError):
            pass


def build_common_allele_public_indexes(
    *, repository_root: str | Path, reference_root: str | Path,
    clinvar_contract_path: str, clinvar_contract_sha256: str,
    clinvar_manifest_path: str, clinvar_manifest_sha256: str,
    open_targets_contract_path: str, open_targets_contract_sha256: str,
    open_targets_manifest_path: str, open_targets_manifest_sha256: str,
    reference_contract_path: str, reference_contract_sha256: str,
    reference_manifest_path: str, reference_manifest_sha256: str,
    exception_manifest_path: str, exception_manifest_sha256: str,
    curated_contract_path: str, curated_bundle_sha256: str,
    accession_report_path: str, accession_report_sha256: str,
    progress=None,
) -> BuiltPublicIndexes:
    """Scan each complete pinned source and publish projection plus accelerator."""

    repository = Path(repository_root).resolve(strict=True)
    references = Path(reference_root).resolve(strict=True)
    with _preparation_stage('contract and manifest metadata', progress):
        clinvar_contract, clinvar_artifacts = _contract_artifacts(
            repository, clinvar_contract_path, clinvar_contract_sha256, {"vcf", "tbi"}
        )
        _contract_bytes(repository, open_targets_contract_path, open_targets_contract_sha256)
        _contract_bytes(repository, reference_contract_path, reference_contract_sha256)
        clinvar_manifest, clinvar_outputs = _manifest_outputs(
            references, clinvar_manifest_path, clinvar_manifest_sha256
        )
        open_targets_manifest, open_targets_outputs = _manifest_outputs(
            references, open_targets_manifest_path, open_targets_manifest_sha256
        )
        if (clinvar_manifest.get("transformation", {}).get("parameters_sha256")
                != clinvar_contract_sha256
                or open_targets_manifest.get("transformation", {}).get("parameters_sha256")
                != open_targets_contract_sha256):
            _fail()
    selected_artifacts: list[dict[str, Any]] = []
    with _preparation_stage('clinical source integrity', progress):
        clinvar_paths = {}
        for role, artifact in clinvar_artifacts.items():
            output = clinvar_outputs.get(artifact["path"])
            if output is None or (output["sha256"], output["size_bytes"]) != (
                    artifact["sha256"], artifact["size_bytes"]):
                _fail()
            clinvar_paths[role] = _regular(
                references, output["path"], output["sha256"], output["size_bytes"]
            )
            selected_artifacts.append({"role": f"clinvar_{role}", **{
                key: output[key] for key in ("path", "sha256", "size_bytes")
            }})
        variant_outputs = sorted(
            (row for path, row in open_targets_outputs.items()
             if "/output/variant/" in path and path.endswith(".parquet")),
            key=lambda row: row["path"],
        )
        evidence_outputs = sorted(
            (row for path, row in open_targets_outputs.items()
             if "/output/evidence_eva/" in path and path.endswith(".parquet")),
            key=lambda row: row["path"],
        )
        if len(variant_outputs) != 25 or len(evidence_outputs) != 8:
            _fail()
        variant_paths, evidence_paths = [], []
        output_by_absolute = {}
        for role, rows, paths in (("open_targets_variant", variant_outputs, variant_paths),
                                  ("open_targets_eva", evidence_outputs, evidence_paths)):
            for row in rows:
                path = _regular(references, row["path"], row["sha256"], row["size_bytes"])
                paths.append(path)
                output_by_absolute[str(path)] = row
                selected_artifacts.append({"role": role, **{
                    key: row[key] for key in ("path", "sha256", "size_bytes")
                }})
    with _preparation_stage('reference source integrity', progress):
        fasta, fai, reference_hashes = _reference_files(
            references, reference_manifest_path, reference_manifest_sha256
        )
        for role, path, sha256 in (("reference_fasta", fasta, reference_hashes["fasta"]),
                                   ("reference_fai", fai, reference_hashes["fai"])):
            selected_artifacts.append({"role": role, "path": path.relative_to(references).as_posix(),
                                       "sha256": sha256, "size_bytes": path.stat().st_size})
    with _preparation_stage('curated evidence preparation', progress):
        exception = verify_exception_bundle(
            reference_root=references, manifest_path=exception_manifest_path,
            expected_manifest_sha256=exception_manifest_sha256,
        )
        from genosapiens.disease_gene_evidence import read_pinned_disease_gene_evidence, _plain
        from genosapiens.curated_alleles import accession_references_from_report, project_curated_index
        from genosapiens.reference_concordance import load_reference_identity
        genes = read_pinned_disease_gene_evidence(contract_path=repository / curated_contract_path,
            reference_root=references, expected_bundle_sha256=curated_bundle_sha256)
        gene_output = _plain(genes.source_binding)['outputs']['gene_assertion']
        selected_artifacts.append({'role': 'curated_gene_assertions', **{
            key: gene_output[key] for key in ('path', 'sha256', 'size_bytes')}})
        report = _regular(references, accession_report_path, accession_report_sha256)
        if report.stat().st_size > MAX_ROW_BYTES:
            _fail()
        selected_artifacts.append(dict(role='curated_accessions', path=accession_report_path,
            sha256=accession_report_sha256, size_bytes=report.stat().st_size))
        public_identity = load_reference_identity(repository / reference_contract_path).reference_role(
            'public_annotation_benchmark').identity_id
        accessions = accession_references_from_report(report.read_text(encoding='ascii'),
            reference_identity=public_identity, reference_sha256=reference_manifest_sha256)
    selected_artifacts.sort(key=lambda row: (row["role"], row["path"]))
    parser_hashes = {
        "admission_links": _file_sha256(Path(admission_links.__file__)),
        "common_allele_public_indexes": _file_sha256(Path(__file__)),
    }
    closure = {
        "schema": "genosapiens.common-allele-public-source-closure", "version": 1,
        "contracts": {"clinvar": clinvar_contract_sha256,
                      "open_targets": open_targets_contract_sha256,
                      "reference": reference_contract_sha256},
        "manifests": {"clinvar": clinvar_manifest_sha256,
                      "open_targets": open_targets_manifest_sha256,
                      "reference": reference_manifest_sha256,
                      "exceptions": exception_manifest_sha256},
        "source_artifacts": selected_artifacts, "parser_sha256s": parser_hashes,
        "curated_source_binding": _plain(genes.source_binding),
        "limits": {"cache_kib": CACHE_KIB, "cursor_batch_rows": CURSOR_BATCH_ROWS,
                   "maximum_row_bytes": MAX_ROW_BYTES, "duckdb_memory": "4GB",
                   "duckdb_threads": 2},
    }
    closure_sha = hashlib.sha256(_canonical(closure)).hexdigest()
    staging_parent = references / "derived/common-allele-public-indexes"
    staging_parent.mkdir(parents=True, exist_ok=True)
    staging = staging_parent / (".staging-" + uuid.uuid4().hex)
    staging.mkdir(mode=0o700)
    projection_stage = staging / "protection.jsonl.gz"
    sqlite_stage = staging / "protection.sqlite"
    try:
        with _MmapFasta(fasta, fai) as public_reference:
            with _ProjectionWriter(projection_stage, sqlite_stage, closure_sha) as writer:
                with _preparation_stage('exception projection', progress):
                    _projection_exceptions(
                        writer, exception.records, public_reference,
                        exception_manifest_path, exception_manifest_sha256,
                    )
                with _preparation_stage('ClinVar projection', progress):
                    _projection_clinvar(
                        writer, clinvar_paths["vcf"], clinvar_artifacts["vcf"], public_reference
                    )
                with _preparation_stage('Open Targets projection', progress):
                    _projection_open_targets(
                        writer, variant_paths, evidence_paths, output_by_absolute,
                        public_reference, staging,
                    )
                if writer.source_counts != {
                    "common_allele_exceptions": 144, "clinvar": 4_467_990,
                    "open_targets": 3_998_459,
                } or writer.sequence != 8_466_593:
                    _fail()
                with _preparation_stage('curated allele projection', progress):
                    curated_counts = project_curated_index(writer,
                        (_plain(record.source_fields) for record in genes.records.values()
                         if record.table == 'gene_assertion'), source_sha256=gene_output['sha256'],
                        accession_references=accessions, fetch_reference=public_reference.fetch)
                    writer.emit(retention_coverage_record(source_bindings=dict(
                        clinvar=clinvar_manifest_sha256, open_targets=open_targets_manifest_sha256,
                        curated=gene_output['sha256'], reference=reference_manifest_sha256,
                        accessions=accession_report_sha256), curated_counts=curated_counts))
        with _preparation_stage('projection completion', progress):
            record_count, record_sha = reconcile_projection_index(
                projection_stage, sqlite_stage, closure_sha
            )
            if record_count != writer.sequence or record_sha != writer.uncompressed_sha256.hexdigest():
                _fail()
        projection_sha = _file_sha256(projection_stage)
        sqlite_sha = _file_sha256(sqlite_stage)
        outputs_identity = {
            "projection": {"sha256": projection_sha, "size_bytes": projection_stage.stat().st_size,
                           "uncompressed_sha256": record_sha,
                           "uncompressed_size_bytes": writer.uncompressed_size},
            "sqlite": {"sha256": sqlite_sha, "size_bytes": sqlite_stage.stat().st_size},
        }
        identity = {"schema": IDENTITY_SCHEMA, "version": 1, "classification": "public",
                    "closure": closure, "closure_sha256": closure_sha,
                    "record_count": record_count, "record_type_counts": writer.counts,
                    "source_record_counts": writer.source_counts, "outputs": outputs_identity}
        index_sha = hashlib.sha256(_canonical(identity)).hexdigest()
        prefix = f"derived/common-allele-public-indexes/{index_sha}"
        projection_relative = f"{prefix}/protection.jsonl.gz"
        sqlite_relative = f"{prefix}/protection.sqlite"
        manifest_relative = f"manifests/common-allele-public-indexes/{index_sha}.json"
        _publish_file(projection_stage, references / projection_relative)
        _publish_file(sqlite_stage, references / sqlite_relative)
        manifest = {
            "schema": MANIFEST_SCHEMA, "version": 1, "classification": "public",
            "qualification_status": "qualified", "controlled_use_authorized": False,
            "challenge_or_truth_inputs_used": False, "index_sha256": index_sha,
            "identity": identity,
            "outputs": {
                "projection": {"path": projection_relative, **outputs_identity["projection"]},
                "sqlite": {"path": sqlite_relative, **outputs_identity["sqlite"]},
            },
        }
        manifest_content = _canonical(manifest)
        manifest_stage = staging / "manifest.json"
        manifest_stage.write_bytes(manifest_content)
        _publish_file(manifest_stage, references / manifest_relative)
        # Publication links the authenticated staging files or verifies its copy.
        # Do not repeat whole-output hashing after that completed publication.
        if (references / manifest_relative).read_bytes() != manifest_content:
            _fail()
        return BuiltPublicIndexes(
            index_sha, projection_relative, projection_sha, sqlite_relative, sqlite_sha,
            manifest_relative, hashlib.sha256(manifest_content).hexdigest(), record_count,
        )
    finally:
        shutil.rmtree(staging, ignore_errors=True)


def verify_common_allele_public_indexes(
    *, reference_root: str | Path, manifest_path: str,
    expected_manifest_sha256: str, verify_sources: bool = True,
    reuse_qualified_reconciliation: bool = False,
) -> VerifiedPublicIndexes:
    if (type(verify_sources) is not bool
            or type(reuse_qualified_reconciliation) is not bool
            or verify_sources and reuse_qualified_reconciliation):
        _fail()
    references = Path(reference_root).resolve(strict=True)
    manifest_file = _regular(references, manifest_path, _sha(expected_manifest_sha256))
    manifest = _json(manifest_file.read_bytes(), canonical=True)
    if set(manifest) != {
        "schema", "version", "classification", "qualification_status",
        "controlled_use_authorized", "challenge_or_truth_inputs_used", "index_sha256",
        "identity", "outputs",
    } or manifest["schema"] != MANIFEST_SCHEMA or manifest["version"] != 1 \
            or manifest["classification"] != "public" \
            or manifest["qualification_status"] != "qualified" \
            or manifest["controlled_use_authorized"] is not False \
            or manifest["challenge_or_truth_inputs_used"] is not False:
        _fail()
    index_sha = _sha(manifest["index_sha256"])
    identity = manifest["identity"]
    if (type(identity) is not dict or identity.get("schema") != IDENTITY_SCHEMA
            or hashlib.sha256(_canonical(identity)).hexdigest() != index_sha):
        _fail()
    closure, closure_sha = identity.get("closure"), _sha(identity.get("closure_sha256"))
    if type(closure) is not dict or hashlib.sha256(_canonical(closure)).hexdigest() != closure_sha:
        _fail()
    parsers = closure.get("parser_sha256s")
    if closure.get("limits") != {
        "cache_kib": CACHE_KIB, "cursor_batch_rows": CURSOR_BATCH_ROWS,
        "maximum_row_bytes": MAX_ROW_BYTES, "duckdb_memory": "4GB", "duckdb_threads": 2,
    } or type(parsers) is not dict or set(parsers) != {
        "admission_links", "common_allele_public_indexes",
    } or any(type(value) is not str or SHA256.fullmatch(value) is None
             for value in parsers.values()):
        _fail()
    if not reuse_qualified_reconciliation and parsers != {
        "admission_links": _file_sha256(Path(admission_links.__file__)),
        "common_allele_public_indexes": _file_sha256(Path(__file__)),
    }:
        _fail()
    outputs = manifest["outputs"]
    if type(outputs) is not dict or set(outputs) != {"projection", "sqlite"}:
        _fail()
    projection = outputs["projection"]
    database = outputs["sqlite"]
    expected_prefix = f"derived/common-allele-public-indexes/{index_sha}/"
    if (type(projection) is not dict or type(database) is not dict
            or projection.get("path") != expected_prefix + "protection.jsonl.gz"
            or database.get("path") != expected_prefix + "protection.sqlite"
            or identity.get("outputs", {}).get("projection")
            != {key: projection[key] for key in projection if key != "path"}
            or identity.get("outputs", {}).get("sqlite")
            != {key: database[key] for key in database if key != "path"}):
        _fail()
    projection_path = _regular(
        references, projection["path"], _sha(projection["sha256"]), projection["size_bytes"],
        trust_published=reuse_qualified_reconciliation,
    )
    database_path = _regular(
        references, database["path"], _sha(database["sha256"]), database["size_bytes"],
        trust_published=reuse_qualified_reconciliation,
    )
    source_artifacts = closure.get("source_artifacts")
    if type(source_artifacts) is not list or not source_artifacts:
        _fail()
    for artifact in source_artifacts:
        if type(artifact) is not dict or set(artifact) != {"role", "path", "sha256", "size_bytes"}:
            _fail()
        _sha(artifact["sha256"])
        if type(artifact["size_bytes"]) is not int or artifact["size_bytes"] <= 0:
            _fail()
    if verify_sources:
        for artifact in source_artifacts:
            _regular(references, artifact["path"], _sha(artifact["sha256"]), artifact["size_bytes"])
    if reuse_qualified_reconciliation:
        # The qualified manifest records the completed projection/index join.
        # Reopening a process checks bindings and the small SQLite binding row;
        # it intentionally does not certify unchanged database pages.
        count = identity.get("record_count")
        if type(count) is not int or count <= 0:
            _fail()
        _sha(projection.get("uncompressed_sha256"))
        descriptor = None
        connection = None
        try:
            descriptor = os.open(database_path, os.O_RDONLY | os.O_NOFOLLOW)
            connection = sqlite3.connect(
                f"file:/proc/self/fd/{descriptor}?mode=ro&immutable=1", uri=True,
            )
            connection.execute("PRAGMA trusted_schema=OFF")
            if connection.execute("SELECT closure_sha256 FROM binding").fetchall() != [(closure_sha,)]:
                _fail()
            for table in ("records", "exact_records", "intervals", "unknown_loci", "scopes"):
                connection.execute(f"SELECT 1 FROM {table} LIMIT 1").fetchone()
        except (OSError, sqlite3.Error):
            _fail()
        finally:
            if connection is not None:
                connection.close()
            if descriptor is not None:
                os.close(descriptor)
    else:
        count, record_sha = reconcile_projection_index(projection_path, database_path, closure_sha)
        if (count != identity.get("record_count")
                or record_sha != projection.get("uncompressed_sha256")):
            _fail()
    return VerifiedPublicIndexes(
        index_sha, expected_manifest_sha256, closure_sha, projection_path, database_path, count
    )


class CommonAllelePublicIndexSession(AbstractContextManager["CommonAllelePublicIndexSession"]):
    """One authenticated read-only stage session with fixed SQLite limits."""

    def __init__(self, verified: VerifiedPublicIndexes, *, limits: PublicQueryLimits) -> None:
        if type(verified) is not VerifiedPublicIndexes or type(limits) is not PublicQueryLimits:
            _fail()
        self.limits = limits
        self.verified = verified
        self.disease_lookup = None
        self._descriptor = os.open(verified.sqlite_path, os.O_RDONLY | os.O_NOFOLLOW)
        self._identity = os.fstat(self._descriptor)
        self.connection = sqlite3.connect(
            f"file:/proc/self/fd/{self._descriptor}?mode=ro&immutable=1", uri=True
        )
        self.connection.execute(f"PRAGMA cache_size=-{CACHE_KIB}")
        self.connection.execute("PRAGMA temp_store=FILE")
        self.connection.execute("PRAGMA trusted_schema=OFF")
        self.connection.execute("PRAGMA query_only=ON")
        self.connection.execute("PRAGMA mmap_size=0")

    def assert_stable(self):
        """Check immutable ownership before returning cached or queried facts."""
        if self.disease_lookup is not None:
            self.disease_lookup.assert_stable()
        def signature(info):
            return (info.st_dev,info.st_ino,info.st_size,info.st_mtime_ns,info.st_ctime_ns)
        try:
            if (self.verified.sqlite_path.is_symlink()
                    or signature(os.fstat(self._descriptor)) != signature(self._identity)
                    or signature(self.verified.sqlite_path.stat()) != signature(self._identity)):
                _fail()
        except OSError:
            _fail()

    @staticmethod
    def _bound_query(query: PublicIndexQuery) -> tuple[object, ...]:
        if (
            type(query) is not PublicIndexQuery
            or query.chromosome not in CHROMOSOMES
            or type(query.position) is not int
            or query.position <= 0
            or DNA.fullmatch(query.reference) is None
            or DNA.fullmatch(query.alternate) is None
            or query.reference == query.alternate
        ):
            _fail()
        alias = query.normalized_alias or SequenceAlias(
            query.chromosome, query.position, query.reference, query.alternate,
        )
        if (
            type(alias) is not SequenceAlias
            or alias.chromosome not in CHROMOSOMES
            or type(alias.position) is not int
            or alias.position <= 0
            or DNA.fullmatch(alias.reference) is None
            or DNA.fullmatch(alias.alternate) is None
            or alias.reference == alias.alternate
        ):
            _fail()
        source_start, source_end = query.source_span or (
            query.position, query.position + len(query.reference) - 1,
        )
        alias_start, alias_end = query.normalized_span or (
            alias.position, alias.position + len(alias.reference) - 1,
        )
        if (
            any(type(value) is not int for value in (
                source_start, source_end, alias_start, alias_end,
            ))
            or min(source_start, source_end, alias_start, alias_end) <= 0
            or source_end < source_start
            or alias_end < alias_start
            or query.native_scope is not None
            and type(query.native_scope) is not str
        ):
            _fail()
        return (
            query.chromosome, query.position, query.reference, query.alternate,
            alias.chromosome, alias.position, alias.reference, alias.alternate,
            source_start, source_end, alias_start, alias_end, query.native_scope,
        )

    def query_many(
        self,
        queries: tuple[PublicIndexQuery, ...],
        *,
        budget: QueryBudget | None = None,
        _include_retention: bool = False,
    ) -> tuple[PublicIndexQueryResult, ...]:
        """Resolve a logical batch through fixed-size set-based SQLite queries."""

        self.assert_stable()
        if type(queries) is not tuple or type(_include_retention) is not bool:
            _fail()
        budget = budget or QueryBudget(
            self.limits, max_row_bytes=MAX_ROW_BYTES, origin="clinical_index",
        )
        if budget.limits != self.limits:
            _fail()
        budget.batch(len(queries))
        bound = tuple(self._bound_query(query) for query in queries)
        records: list[list[Mapping[str, Any]]] = [[] for _query in queries]
        sources = "'clinvar','open_targets','common_allele_exceptions'"
        if _include_retention:
            sources += ",'curated_functional'"
        with budget.at_origin("clinical_index"), budget.sqlite(self.connection):
            for offset in range(0, len(bound), QUERY_CTE_ROWS):
                rows = bound[offset:offset + QUERY_CTE_ROWS]
                values = ",".join("(?,?,?,?,?,?,?,?,?,?,?,?,?,?)" for _row in rows)
                parameters: list[object] = []
                for index, row in enumerate(rows, offset):
                    parameters.extend((index, *row))
                cursor = self.connection.execute(f"""
                  WITH input(
                    query_ordinal,chromosome,position,reference,alternate,
                    alias_chromosome,alias_position,alias_reference,alias_alternate,
                    source_start,source_end,alias_start,alias_end,native_scope
                  ) AS (VALUES {values}),
                  matched(query_ordinal,sequence) AS (
                    SELECT q.query_ordinal,e.sequence FROM input q JOIN exact_records e ON
                      e.source IN ({sources})
                      AND e.chromosome=q.chromosome AND e.position=q.position
                      AND e.reference=q.reference AND e.alternate=q.alternate
                    UNION
                    SELECT q.query_ordinal,e.sequence FROM input q JOIN exact_records e ON
                      e.source IN ({sources})
                      AND e.alias_chromosome=q.alias_chromosome
                      AND e.alias_position=q.alias_position
                      AND e.alias_reference=q.alias_reference
                      AND e.alias_alternate=q.alias_alternate
                    UNION
                    SELECT q.query_ordinal,i.sequence FROM input q JOIN intervals i ON
                      i.source IN ('clinvar','open_targets','common_allele_exceptions')
                      AND i.chromosome=q.chromosome
                      AND i.start<=q.source_end AND i.end>=q.source_start
                    UNION
                    SELECT q.query_ordinal,i.sequence FROM input q JOIN intervals i ON
                      i.source IN ('clinvar','open_targets','common_allele_exceptions')
                      AND i.chromosome=q.alias_chromosome
                      AND i.start<=q.alias_end AND i.end>=q.alias_start
                    UNION
                    SELECT q.query_ordinal,s.sequence FROM input q JOIN scopes s ON
                      s.scope=q.native_scope
                    UNION
                    SELECT q.query_ordinal,u.sequence FROM input q JOIN unknown_loci u ON
                      u.source IN ('clinvar','open_targets','common_allele_exceptions')
                      AND u.chromosome=q.chromosome
                      AND u.start<=q.source_end AND u.end>=q.source_start
                    UNION
                    SELECT q.query_ordinal,u.sequence FROM input q JOIN unknown_loci u ON
                      u.source IN ('clinvar','open_targets','common_allele_exceptions')
                      AND u.chromosome=q.alias_chromosome
                      AND u.start<=q.alias_end AND u.end>=q.alias_start
                  )
                  SELECT m.query_ordinal,m.sequence,length(r.canonical_row),r.canonical_row
                  FROM matched m CROSS JOIN records r ON r.sequence=m.sequence
                  ORDER BY m.query_ordinal,m.sequence
                """, parameters)
                for ordinal, _sequence, size, content in cursor:
                    # Charge SQLite's reported BLOB length before JSON decoding.
                    budget.charge(size, "projection")
                    records[ordinal].append(_json(content, canonical=True))
        results = []
        for matched in records:
            values = tuple(matched)
            if any(row["record_type"] == "unknown" for row in values):
                results.append(PublicIndexQueryResult(
                    "retain_unknown", values,
                    ("applicable_public_source_record_unassessable",),
                ))
            elif values:
                results.append(PublicIndexQueryResult(
                    "matches", values, ("authenticated_public_protection_match",),
                ))
            else:
                results.append(PublicIndexQueryResult(
                    "no_match", (), ("no_protection_found_in_complete_pinned_sources",),
                ))
        return tuple(results)

    def retention_coverage(self):
        """Reject an obsolete or incomplete selection index before empty lookup."""
        self.assert_stable()
        if not hasattr(self, '_retention_coverage'):
            rows = self.connection.execute('''SELECT r.sequence, length(r.canonical_row) FROM scopes s
                JOIN records r USING(sequence) WHERE s.source='retention_coverage'
                AND s.scope=? LIMIT 2''', (RETENTION_SCOPE,)).fetchall()
            if len(rows) != 1 or not 0 < rows[0][1] <= MAX_ROW_BYTES:
                raise CommonAlleleIndexError('selection index needs current clinical and curated source closure')
            content, = self.connection.execute('SELECT canonical_row FROM records WHERE sequence=?',
                (rows[0][0],)).fetchone()
            row = _json(content, canonical=True)
            expected = retention_coverage_record(source_bindings=row.get('source_bindings', {}),
                curated_counts=row.get('curated_counts', {}))
            # Allele queries can use the published format that predates the
            # optional disease lookup. Preserve its actual capability metadata.
            if 'disease_lookup_version' not in row:
                expected.pop('disease_lookup_version')
            if row != expected:
                _fail()
            self._retention_coverage = row
        # Do not expose mutable cached metadata to consumers.
        return _json(_canonical(self._retention_coverage), canonical=True)

    def retention_many(self, queries, *, budget=None):
        self.retention_coverage()
        return self.query_many(queries, budget=budget, _include_retention=True)

    def iter_retention_exact_records(self):
        """Stream exact public assertions once for public scope production.

        The projection version is independent of runtime query implementation.
        Fetch at most 64 MiB and 1,024 rows before yielding each bounded batch.
        Oversized blobs never cross the SQLite boundary into Python.
        """
        self.retention_coverage()
        batch_rows = min(CURSOR_BATCH_ROWS, (64 * 1024**2) // MAX_ROW_BYTES)
        try:
            cursor = self.connection.execute('''SELECT e.source,
                e.chromosome,e.position,e.reference,e.alternate,
                e.alias_chromosome,e.alias_position,e.alias_reference,e.alias_alternate,
                r.canonical_sha256,length(r.canonical_row),
                CASE WHEN length(r.canonical_row)<=? THEN r.canonical_row END
                FROM exact_records e JOIN records r USING(sequence)
                WHERE e.source IN ('clinvar','open_targets','curated_functional')
                ORDER BY e.sequence''', (MAX_ROW_BYTES,))
            while batch := cursor.fetchmany(batch_rows):
                self.assert_stable()
                for source, chrom, pos, ref, alt, ac, ap, ar, aa, digest, size, content in batch:
                    if not 0 < size <= MAX_ROW_BYTES or content is None \
                            or hashlib.sha256(content).hexdigest() != digest:
                        _fail()
                    row = _json(content, canonical=True)
                    if row.get('record_type') != 'exact' or row.get('source') != source \
                            or row.get('allele') != dict(chromosome=chrom, position=pos,
                                reference=ref, alternate=alt) \
                            or row.get('normalized_alias') != dict(chromosome=ac, position=ap,
                                reference=ar, alternate=aa):
                        _fail()
                    yield row
            self.assert_stable()
        except sqlite3.Error:
            raise CommonAlleleIndexError('Exact retention projection is incomplete.') from None

    def retention_for_diseases(self, disease_ids, *, budget=None):
        """Read exact assertions by frozen disease identity with a bounded query.

        No gene metadata or genomic window is required. Unknown disease names
        are not promoted to identifier equivalence.
        """
        embedded = self.retention_coverage().get('disease_lookup_version') == 1
        if not embedded and self.disease_lookup is None:
            raise CommonAlleleIndexError('Published index has no disease lookup. Prepare that lookup before disease queries.')
        table = 'retention_diseases' if embedded else 'disease_lookup.retention_diseases'
        if (type(disease_ids) is not tuple or not disease_ids
                or any(type(value) is not str or not re.fullmatch(
                    r'(MONDO|OMIM|ORPHA):[0-9]+', value) for value in disease_ids)):
            _fail()
        budget = budget or QueryBudget(self.limits, max_row_bytes=MAX_ROW_BYTES, origin='clinical_index')
        if budget.limits != self.limits:
            _fail()
        budget.batch(len(disease_ids))
        found = {}
        with budget.at_origin('clinical_index'), budget.sqlite(self.connection):
            for disease in sorted(set(disease_ids)):
                cursor = self.connection.execute(f'''SELECT d.sequence,length(r.canonical_row),r.canonical_row
                    FROM {table} d JOIN records r ON r.sequence=d.sequence
                    WHERE d.disease_id=? ORDER BY d.sequence''', (disease,))
                for sequence, size, content in cursor:
                    if sequence in found:
                        continue
                    budget.charge(size, 'projection')
                    row = _json(content, canonical=True)
                    if disease not in retention_disease_ids(row):
                        _fail()
                    found[sequence] = row
        self.assert_stable()
        return tuple(found[key] for key in sorted(found))

    def query(
        self, chromosome: str, position: int, reference: str, alternate: str, *,
        normalized_alias: SequenceAlias | None = None,
        source_span: tuple[int, int] | None = None,
        normalized_span: tuple[int, int] | None = None,
        native_scope: str | None = None,
        budget: QueryBudget | None = None,
    ) -> PublicIndexQueryResult:
        return self.query_many((PublicIndexQuery(
            chromosome, position, reference, alternate,
            normalized_alias, source_span, normalized_span, native_scope,
        ),), budget=budget)[0]

    def __exit__(self, *_args: object) -> None:
        self.connection.close()
        final = os.fstat(self._descriptor)
        os.close(self._descriptor)
        if (self._identity.st_size, self._identity.st_mtime_ns, self._identity.st_ctime_ns) != (
                final.st_size, final.st_mtime_ns, final.st_ctime_ns):
            _fail()


def unavailable_public_index_result() -> PublicIndexQueryResult:
    return PublicIndexQueryResult("unavailable", (), ("required_public_index_unavailable",))
