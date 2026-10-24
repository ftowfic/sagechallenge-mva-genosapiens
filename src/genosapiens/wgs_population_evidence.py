"""Complete population policy in the normalized evidence owner.

The ordered context and exception arrays retain the original logical policy
digest. A disposable disk index supports bounded context lookup and preserves
global exception ordinals. Member evaluation uses the ordinary policy kernel.
"""
import hashlib
import json
import sqlite3

from genosapiens.evidence_model import _sha, _typed
from genosapiens.evidence_serialization import _canonical
from genosapiens.population_policy import (
    ContextKey, ContextPolicy, CuratedException, PreparedPopulationPolicy,
    evaluate_population_member_inputs,
)
from genosapiens.wgs_evidence_store import EvidenceStoreError
from genosapiens.wgs_population_policy import _LogicalDigest

MAX_POLICY_ROWS = 1_000_000
MAX_POLICY_WINDOW_BYTES = 64 * 1024**2


def _content(value):
    content = _canonical(value)
    if len(content) > MAX_POLICY_WINDOW_BYTES:
        raise EvidenceStoreError("population policy window exceeds capacity", capacity=True)
    return content


class _PopulationPolicyOwner:
    def _initialize(self, connection, name):
        self.connection, self.name, self.closed = connection, name, False
        self.prefix = "population_" + hashlib.sha256(name.encode()).hexdigest()
        self.tables = []
        for suffix, definition in (
            ("contexts", "ordinal INTEGER UNIQUE,key TEXT PRIMARY KEY,content TEXT"),
            ("exceptions", "ordinal INTEGER PRIMARY KEY,key TEXT,content TEXT UNIQUE"),
        ):
            table = self.prefix + "_" + suffix
            self.connection.execute(f"CREATE TABLE {table} ({definition})")
            self.tables.append(table)
        self.connection.execute(f"CREATE INDEX {self.prefix}_exception_key ON {self.prefix}_exceptions(key,ordinal)")
        self.context_count = self.exception_count = 0

    def _contexts(self, rows, digest, *, decode=False):
        for row in rows:
            context = _typed(row, ContextPolicy, decode=decode)
            self.context_count += 1
            if self.context_count > MAX_POLICY_ROWS:
                raise EvidenceStoreError("population context count exceeds capacity", capacity=True)
            self.connection.execute(f"INSERT INTO {self.prefix}_contexts VALUES (?,?,?)",
                (self.context_count, _canonical(context.context), _content(context)))
            digest.add_context(context)
            yield context

    def _exceptions(self, rows, digest, *, decode=False):
        for row in rows:
            exception = _typed(row, CuratedException, decode=decode)
            key = _canonical(exception.context)
            if not self.connection.execute(f"SELECT 1 FROM {self.prefix}_contexts WHERE key=?", (key,)).fetchone():
                raise EvidenceStoreError("population exception context is missing")
            ordinal = self.exception_count
            self.exception_count += 1
            if self.exception_count > MAX_POLICY_ROWS:
                raise EvidenceStoreError("population exception count exceeds capacity", capacity=True)
            self.connection.execute(f"INSERT INTO {self.prefix}_exceptions VALUES (?,?,?)",
                (ordinal, key, _content(exception)))
            digest.add_exception(exception)
            yield exception

    def _require_open(self):
        if self.closed:
            raise EvidenceStoreError("population evidence owner is closed")

    @property
    def reference(self):
        self._require_open()
        return {"name": self.name, "kind": "population-policy-evidence-v1", "sha256": self.sha256,
                "policy_id": self.policy_id, "policy_sha256": self.policy_sha256}

    def inputs(self, key):
        self._require_open()
        key = _typed(key, ContextKey | None)
        row = self.connection.execute(f"SELECT content FROM {self.prefix}_contexts WHERE key=?",
                                      (_canonical(key),)).fetchone()
        if row is None:
            return None, ()
        settings = _typed(json.loads(row[0]), ContextPolicy, decode=True)
        size, exceptions = len(row[0]), []
        for ordinal, content in self.connection.execute(
                f"SELECT ordinal,content FROM {self.prefix}_exceptions WHERE key=? ORDER BY ordinal", (_canonical(key),)):
            size += len(content)
            if size > MAX_POLICY_WINDOW_BYTES:
                raise EvidenceStoreError("population policy window exceeds capacity", capacity=True)
            exceptions.append((ordinal, _typed(json.loads(content), CuratedException, decode=True)))
        return settings, tuple(exceptions)

    def evaluate(self, key, alleles):
        settings, exceptions = self.inputs(key)
        return evaluate_population_member_inputs(key, alleles, settings=settings, indexed_exceptions=exceptions)

    def iter_contexts(self):
        self._require_open()
        for (content,) in self.connection.execute(f"SELECT content FROM {self.prefix}_contexts ORDER BY ordinal"):
            yield _typed(json.loads(content), ContextPolicy, decode=True)

    def iter_exceptions(self):
        self._require_open()
        for (content,) in self.connection.execute(f"SELECT content FROM {self.prefix}_exceptions ORDER BY ordinal"):
            yield _typed(json.loads(content), CuratedException, decode=True)

    def close(self):
        for table in reversed(self.tables):
            self.connection.execute("DROP TABLE " + table)
        self.tables.clear()
        self.closed = True

    def __enter__(self):
        self._require_open()
        return self

    def __exit__(self, *_):
        self.close()


class CompletePopulationPolicy(_PopulationPolicyOwner):
    def __init__(self, writer, prepared_policy, *, name):
        from genosapiens.wgs_population_input import GeneratedPopulationPolicy
        generated = type(prepared_policy) is GeneratedPopulationPolicy
        if not generated and type(prepared_policy) is not PreparedPopulationPolicy:
            raise EvidenceStoreError("population policy owner is invalid")
        if generated:
            prepared_policy._require()
            if not prepared_policy.sealed:
                raise EvidenceStoreError("population policy input has not globally closed")
        self.tables = []
        try:
            self._initialize(writer.connection, name)
            self.policy_id = prepared_policy.policy_id if generated else prepared_policy.policy.policy_id
            self.policy_sha256 = prepared_policy.policy_sha256
            digest = _LogicalDigest(self.policy_id)
            contexts = writer.add_sequence(name + ".contexts", "population-policy-contexts-v1",
                self._contexts(prepared_policy.iter_contexts() if generated else prepared_policy.policy.contexts, digest))
            exceptions = writer.add_sequence(name + ".exceptions", "population-policy-exceptions-v1",
                self._exceptions(prepared_policy.iter_exceptions() if generated else prepared_policy.policy.exceptions, digest))
            if digest.finish() != self.policy_sha256:
                raise EvidenceStoreError("logical population policy digest changed")
            self.sha256 = writer.add(name, "population-policy-evidence-v1", {
                "schema": "genosapiens.wgs-population-policy-evidence", "version": 1,
                "policy_id": self.policy_id, "policy_sha256": self.policy_sha256,
                "context_count": self.context_count, "exception_count": self.exception_count,
                "contexts": {"name": name + ".contexts", "sha256": contexts},
                "exceptions": {"name": name + ".exceptions", "sha256": exceptions},
            })
        except BaseException as error:
            writer.failed = True
            self.close()
            if isinstance(error, (ValueError, TypeError, KeyError, sqlite3.Error)):
                raise EvidenceStoreError("invalid complete population policy") from None
            raise


class VerifiedPopulationPolicy(_PopulationPolicyOwner):
    def __init__(self, reader, reference, *, expected_policy_sha256):
        self.tables = []
        try:
            _sha(expected_policy_sha256)
            if reference.kind != "population-policy-evidence-v1":
                raise EvidenceStoreError("foreign population evidence root")
            header = reader.reconstruct(reference)
            if (set(header) != {"schema", "version", "policy_id", "policy_sha256", "context_count",
                    "exception_count", "contexts", "exceptions"}
                    or header["schema"] != "genosapiens.wgs-population-policy-evidence"
                    or type(header["version"]) is not int or header["version"] != 1
                    or type(header["policy_id"]) is not str or not header["policy_id"].strip()
                    or header["policy_sha256"] != expected_policy_sha256
                    or any(type(header[key]) is not int or not 0 <= header[key] <= MAX_POLICY_ROWS
                           for key in ("context_count", "exception_count"))):
                raise EvidenceStoreError("invalid population evidence header")
            self._initialize(reader.connection, reference.name)
            self.policy_id, self.policy_sha256 = header["policy_id"], expected_policy_sha256
            digest = _LogicalDigest(self.policy_id)
            for key, function in (("contexts", self._contexts), ("exceptions", self._exceptions)):
                child = header[key]
                if set(child) != {"name", "sha256"} or child["name"] != reference.name + "." + key:
                    raise EvidenceStoreError("population policy child binding changed")
                child_reference = reader.reference(child["name"], kind="population-policy-" + key + "-v1")
                if child_reference.sha256 != child["sha256"]:
                    raise EvidenceStoreError("population policy child digest changed")
                for _ in function(reader.iter_array(child_reference, max_item_bytes=MAX_POLICY_WINDOW_BYTES), digest, decode=True):
                    pass
            if ((self.context_count, self.exception_count) != (header["context_count"], header["exception_count"])
                    or digest.finish() != self.policy_sha256):
                raise EvidenceStoreError("complete logical population policy changed")
            self.sha256 = reference.sha256
        except BaseException as error:
            self.close()
            if isinstance(error, (ValueError, TypeError, KeyError, sqlite3.Error)):
                raise EvidenceStoreError("invalid complete population policy") from None
            raise
