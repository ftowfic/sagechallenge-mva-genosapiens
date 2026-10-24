"""Strict JSON and typed Parquet codecs for task 9.1, always controlled/no-LLM.

JSON operates only on supplied bytes/values. Parquet accesses exactly the local
file the caller selects, uses an in-memory DuckDB, and never resolves references.
Inventory, retention, storage authorization, and source authentication belong to
the caller. No artifact is made LLM-readable by serialization or validation.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import fields, is_dataclass
from enum import Enum
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import stat
from types import UnionType
from typing import Literal, TypeVar, Union, get_args, get_origin

from genosapiens.allele_identity import AlleleIdentity, identity_to_dict
from genosapiens.disease_gene_evidence import NATIVE_LOCUS_TABLES, PUBLIC_TABLES, SCHEMA as GENE_SCHEMA
from genosapiens.disease_profiles import SCHEMA as PROFILE_SCHEMA
from genosapiens.evidence_model import (
    ENTITY_TYPES, SCHEMA, EvidenceModelError, _Controlled, _hints, _typed,
    _validation_session,
)
from genosapiens.rare_disease import TABLE_SCHEMAS


DUCKDB_VERSION = "1.4.2"
_T = TypeVar("_T", bound=_Controlled)
_IDENTITIES = {cls.__dataclass_fields__["kind"].default.value: cls for cls in get_args(AlleleIdentity)}


def _plain(value):
    if value is None or type(value) in (str, int, float, bool):
        return value
    if isinstance(value, Enum):
        return value.value
    if type(value) in _IDENTITIES.values():
        return identity_to_dict(value)
    if is_dataclass(value):
        return {item.name: _plain(getattr(value, item.name)) for item in fields(value)}
    if isinstance(value, Mapping):
        return {key: _plain(child) for key, child in value.items()}
    if type(value) in (tuple, list):
        return [_plain(child) for child in value]
    return value


def _canonical(value):
    return _encode_plain(_plain(value))


def _encode_plain(value):
    """Encode a plain tree already produced by the validating caller."""
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False)


def _owned_json_default(value):
    """Expose one immutable node; the JSON encoder traverses scalar leaves."""
    if isinstance(value,Enum):return value.value
    if type(value) in _IDENTITIES.values():return identity_to_dict(value)
    if is_dataclass(value):return {item.name:getattr(value,item.name) for item in fields(value)}
    if isinstance(value,Mapping):return dict(value)
    raise TypeError('unsupported owned JSON node')


def _owned_canonical_bytes(value):
    """Canonical bytes for a private already-validated tree; no authority cache.

    Shallow dataclass/proxy views let the C JSON encoder traverse primitive
    leaves without first building another complete plain tree in Python.
    Ordinary evidence codecs retain their independent strict validation path.
    """
    return (json.dumps(value,sort_keys=True,separators=(',',':'),ensure_ascii=True,
        allow_nan=False,default=_owned_json_default)+'\n').encode('ascii')


def _pairs(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise EvidenceModelError("duplicate JSON object key")
        result[key] = value
    return result


def _nonfinite(_value):
    raise EvidenceModelError("nonfinite JSON number")


def _loads(content):
    return json.loads(content, object_pairs_hook=_pairs, parse_constant=_nonfinite)


def _entity_type(cls):
    if cls not in ENTITY_TYPES:
        raise EvidenceModelError("unsupported evidence entity type")


def _type_contract(hint):
    """Pin fields AND nullability/enum vocabulary, not only DuckDB storage types."""
    origin, args = get_origin(hint), get_args(hint)
    if is_dataclass(hint):
        return {
            "record": hint.__module__ + "." + hint.__name__,
            "fields": [
                [item.name, _type_contract(_hints(hint)[item.name]), item.init]
                for item in fields(hint)
            ],
        }
    if isinstance(hint, type) and issubclass(hint, Enum):
        return {"enum": hint.__module__ + "." + hint.__name__, "values": [item.value for item in hint]}
    if origin is Literal:
        return {"literal": list(args)}
    if origin:
        return {"container": str(origin), "types": ["..." if arg is Ellipsis else _type_contract(arg) for arg in args]}
    return str(hint)


def schema_sha256(entity_type: type[_T]) -> str:
    """Runtime schema pin; callers may additionally freeze this in a run manifest."""
    _entity_type(entity_type)
    return _schema_sha256(entity_type)


@lru_cache(maxsize=len(ENTITY_TYPES))
def _schema_sha256(entity_type):
    # Only schema_sha256 admits keys. Never retain evidence instances or arbitrary
    # runtime types; the supported frozen entity vocabulary bounds this cache.
    _entity_type(entity_type)
    return hashlib.sha256(_canonical({
        "schema": SCHEMA, "contract": _type_contract(entity_type), "duckdb": DUCKDB_VERSION,
        "parquet": "typed-columns; tagged-identity-slots; integer-or-real-struct; raw-json-varchar; row-ordinal",
        "public_schemas": [PROFILE_SCHEMA, GENE_SCHEMA],
        "native_source_tables": NATIVE_LOCUS_TABLES,
        "source_tables": {name: TABLE_SCHEMAS[name] for name in ("phenotype_assertion", *PUBLIC_TABLES)},
    }).encode()).hexdigest()


def _metadata(cls):
    return {
        "schema": SCHEMA, "entity_type": cls.__name__, "schema_sha256": schema_sha256(cls),
        "classification": "controlled", "llm_policy": "no-llm", "duckdb_version": DUCKDB_VERSION,
    }


def _owned_json_mapping(value):
    """Plain envelope for a private, strictly validated snapshot.

    Callers must own the validation operation and its copied value. External
    to_json/from_json always establish a fresh isolated validation boundary.
    This helper carries no validation token and never caches caller objects.
    """
    return {**_metadata(type(value)), 'data': _plain(value)}


def _owned_json_envelope(value):
    """Private envelope retaining the validated object for the JSON encoder."""
    return {**_metadata(type(value)),'data':value}


def to_json(entity: _T) -> bytes:
    """Canonical UTF-8 JSON envelope; preserves every null, record, and array order."""
    try:
        _entity_type(type(entity))
        # An external encoding is a new strict operation, even inside a caller's
        # larger validation scope. The copied children are discarded on exit.
        with _validation_session(isolated=True):
            value = _typed(entity, type(entity))
            return (_encode_plain(_owned_json_mapping(value)) + "\n").encode("utf-8")
    except (ValueError, TypeError, RuntimeError, AttributeError, OverflowError, RecursionError):
        raise EvidenceModelError("cannot encode controlled evidence JSON") from None


def from_json(content: bytes | str, entity_type: type[_T]) -> _T:
    """Reject missing/unknown fields, duplicate keys, nonfinite values and wrong pins."""
    try:
        _entity_type(entity_type)
        if type(content) not in (bytes, str):
            raise EvidenceModelError("JSON bytes or text required")
        if type(content) is bytes:
            content = content.decode("utf-8", errors="strict")
        data = _loads(content)
        return _from_mapping(data, entity_type)
    except (ValueError, TypeError, RuntimeError, AttributeError, OverflowError, RecursionError):
        raise EvidenceModelError("invalid controlled evidence JSON") from None


def _from_mapping(data, entity_type):
    """Strict typed envelope already parsed by a canonical outer decoder."""
    _entity_type(entity_type)
    metadata = _metadata(entity_type)
    if type(data) is not dict or set(data) != {*metadata, 'data'} or any(
            data[key] != value for key, value in metadata.items()):
        raise EvidenceModelError('evidence envelope differs from its pin')
    return _typed(data['data'], entity_type, decode=True)


def _optional(hint):
    args = get_args(hint)
    if get_origin(hint) in (UnionType, Union) and type(None) in args and len(args) == 2:
        return next(item for item in args if item is not type(None))
    return hint


def _raw_json(hint):
    return hint is object or get_origin(hint) is Mapping and get_args(hint) == (str, object)


def _quote(name):
    return '"' + name.replace('"', '""') + '"'


def _sql_type(hint):
    hint = _optional(hint)
    origin, args = get_origin(hint), get_args(hint)
    if hint == AlleleIdentity:
        return "STRUCT(" + ", ".join([
            '"kind" VARCHAR',
            *(_quote(kind) + " " + _sql_type(cls) for kind, cls in _IDENTITIES.items()),
        ]) + ")"
    if _raw_json(hint):
        return "VARCHAR"
    if is_dataclass(hint):
        return "STRUCT(" + ", ".join(
            _quote(item.name) + " " + _sql_type(_hints(hint)[item.name]) for item in fields(hint)
        ) + ")"
    if origin is tuple:
        if len(args) == 2 and args[1] is Ellipsis:
            return _sql_type(args[0]) + "[]"
        # Fixed tuples (e.g. confidence bounds) use structs, preserving arity/types.
        return "STRUCT(" + ", ".join(_quote("f" + str(i)) + " " + _sql_type(arg) for i, arg in enumerate(args)) + ")"
    if origin is Mapping:
        return "MAP(" + _sql_type(args[0]) + ", " + _sql_type(args[1]) + ")"
    if origin is Literal or isinstance(hint, type) and issubclass(hint, Enum):
        return "VARCHAR"
    if hint is float:
        # Upstream numeric records accept both 4 and 4.0 and retain raw identity.
        return 'STRUCT("integer_value" BIGINT, "real_value" DOUBLE)'
    return {str: "VARCHAR", int: "BIGINT", bool: "BOOLEAN"}[hint]


def _sql_value(value, hint, *, decode=False):
    """Translate typed values to/from Parquet structs, never a whole-entity blob."""
    if value is None:
        return None
    hint = _optional(hint)
    origin, args = get_origin(hint), get_args(hint)
    if hint == AlleleIdentity:
        if decode:
            if not isinstance(value, Mapping) or set(value) != {"kind", *_IDENTITIES}:
                raise EvidenceModelError("invalid native identity struct")
            kind = value["kind"]
            if kind not in _IDENTITIES or value[kind] is None or any(
                value[other] is not None for other in _IDENTITIES if other != kind
            ):
                raise EvidenceModelError("invalid native identity slots")
            result = _sql_value(value[kind], _IDENTITIES[kind], decode=True)
            if result["kind"] != kind:
                raise EvidenceModelError("native identity tags disagree")
            return result
        kind = value.kind.value
        return {"kind": kind, **{
            slot: _sql_value(value, cls) if slot == kind else None for slot, cls in _IDENTITIES.items()
        }}
    if _raw_json(hint):
        if not decode:
            return _canonical(value)
        result = _loads(value)
        if _canonical(result) != value:
            raise EvidenceModelError("noncanonical raw source JSON")
        return result
    if is_dataclass(hint):
        return {
            item.name: _sql_value(value[item.name] if decode else getattr(value, item.name),
                                  _hints(hint)[item.name], decode=decode)
            for item in fields(hint)
        }
    if origin is tuple:
        if len(args) == 2 and args[1] is Ellipsis:
            return [_sql_value(item, args[0], decode=decode) for item in value]
        if decode:
            return [_sql_value(value["f" + str(i)], arg, decode=True) for i, arg in enumerate(args)]
        return {"f" + str(i): _sql_value(item, arg) for i, (item, arg) in enumerate(zip(value, args, strict=True))}
    if origin is Mapping:
        return {key: _sql_value(child, args[1], decode=decode) for key, child in value.items()}
    if hint is float:
        if decode:
            if (value["integer_value"] is None) == (value["real_value"] is None):
                raise EvidenceModelError("invalid numeric value slots")
            return value["integer_value"] if value["integer_value"] is not None else value["real_value"]
        return {"integer_value": value if type(value) is int else None,
                "real_value": value if type(value) is float else None}
    return value.value if isinstance(value, Enum) else value


def _connection():
    import duckdb

    if duckdb.__version__ != DUCKDB_VERSION:
        raise EvidenceModelError("Parquet codec requires its pinned DuckDB version")
    return duckdb.connect(":memory:", config={
        "threads": "1", "temp_directory": "",
        "autoinstall_known_extensions": "false", "autoload_known_extensions": "false",
    })


def _path(path, *, writing):
    if type(path) is not str and not isinstance(path, Path):
        raise EvidenceModelError("explicit local file path required")
    path = Path(path)
    if not path.is_absolute() or any(token in str(path) for token in (":", "*", "?", "[", "]", "{", "}")):
        raise EvidenceModelError("explicit local file path required")
    if path.is_symlink() or (writing and (path.exists() or not path.parent.is_dir())):
        raise EvidenceModelError("invalid local output file")
    if not writing and not stat.S_ISREG(path.stat().st_mode):
        raise EvidenceModelError("local input must be a regular file")
    return str(path)


def _create_table(connection, entity_type):
    definitions = ['"_ordinal" BIGINT'] + [
        _quote(item.name) + " " + _sql_type(_hints(entity_type)[item.name]) for item in fields(entity_type)
    ]
    connection.execute("CREATE TABLE evidence (" + ", ".join(definitions) + ")")
    return tuple((row[0], row[1]) for row in connection.execute("DESCRIBE evidence").fetchall())


def _decode_row(row, entity_type):
    return _typed({
        item.name: _sql_value(value, _hints(entity_type)[item.name], decode=True)
        for item, value in zip(fields(entity_type), row, strict=True)
    }, entity_type, decode=True)


def write_parquet(path: str | Path, records: tuple[_T, ...], entity_type: type[_T]) -> None:
    """Write a new, caller-selected local file; no overwrite, mkdir, or sidecars.

    Typed columns and nested structs are checked for exact semantic round trips
    before publication. Arbitrary raw source fields alone use canonical JSON.
    """
    try:
        from duckdb import Value

        _entity_type(entity_type)
        path = _path(path, writing=True)
        records = _typed(records, tuple[entity_type, ...])
        with _connection() as connection:
            column_types = {name: connection.sqltype(dtype) for name, dtype in _create_table(connection, entity_type)}
            if records:
                placeholders = ", ".join("?" for _ in range(len(fields(entity_type)) + 1))
                # Explicit types prevent NULL inference from nested empty maps/lists.
                connection.executemany("INSERT INTO evidence VALUES (" + placeholders + ")", [
                    [index, *(Value(_sql_value(getattr(record, item.name), _hints(entity_type)[item.name]), column_types[item.name])
                              for item in fields(entity_type))]
                    for index, record in enumerate(records)
                ])
            rows = connection.execute('SELECT * FROM evidence ORDER BY "_ordinal"').fetchall()
            for record, row in zip(records, rows, strict=True):
                if to_json(_decode_row(row[1:], entity_type)) != to_json(record):
                    raise EvidenceModelError("lossy Parquet value coercion")
            metadata = ", ".join(
                "'" + key + "': '" + value.replace("'", "''") + "'" for key, value in _metadata(entity_type).items()
            )
            connection.execute(
                "COPY evidence TO ? (FORMAT PARQUET, COMPRESSION ZSTD, KV_METADATA {" + metadata + "})", [path],
            )
    except Exception:
        # DuckDB/OS exceptions can embed a path or controlled SQL value.
        raise EvidenceModelError("cannot write controlled evidence Parquet") from None


def read_parquet(path: str | Path, entity_type: type[_T]) -> tuple[_T, ...]:
    """Validate schema/pins/types before decoding. No coercing SELECTs or repairs."""
    try:
        _entity_type(entity_type)
        path = _path(path, writing=False)
        with _connection() as connection:
            expected = _create_table(connection, entity_type)
            pairs = connection.execute("SELECT key, value FROM parquet_kv_metadata(?)", [path]).fetchall()
            metadata = _pairs((bytes(key).decode("utf-8"), bytes(value).decode("utf-8")) for key, value in pairs)
            if metadata != _metadata(entity_type):
                raise EvidenceModelError("Parquet metadata differs from pinned schema")
            query = "SELECT * FROM read_parquet(?, hive_partitioning = false, union_by_name = false)"
            actual = tuple((row[0], row[1]) for row in connection.execute("DESCRIBE " + query, [path]).fetchall())
            if actual != expected:
                raise EvidenceModelError("Parquet column types differ from pinned schema")
            rows = connection.execute(query + ' ORDER BY "_ordinal"', [path]).fetchall()
            if any(type(row[0]) is not int or row[0] != index for index, row in enumerate(rows)):
                raise EvidenceModelError("invalid Parquet row ordering")
            return tuple(_decode_row(row[1:], entity_type) for row in rows)
    except Exception:
        raise EvidenceModelError("invalid controlled evidence Parquet") from None
