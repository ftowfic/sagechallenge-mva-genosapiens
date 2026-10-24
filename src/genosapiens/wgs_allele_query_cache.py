"""Case-bound resource facts below context report and policy construction.

Only private new SQLite state is writable. Cached facts do not decide admission,
population eligibility, clinical significance or hypothesis membership. Resource
owners remain responsible for their strict source/adaptation contracts.
"""
from __future__ import annotations

from collections import OrderedDict
from collections.abc import Mapping
from dataclasses import dataclass, fields, is_dataclass
import hashlib
import json
import os
from pathlib import Path
import sqlite3
import stat
import sys
from types import MappingProxyType
import zlib

from genosapiens.evidence_serialization import _plain
from genosapiens.wgs_ledger import canonical


class CaseAlleleQueryCacheError(RuntimeError):
    """Static private-cache failure, with no scientific payload in diagnostics."""


_BINDINGS = {'schema','version','consumer_executable_sha256','case_input_sha256',
    'reference_role','reference_identity_id','reference_manifest_sha256',
    'population_sources','clinical_sources','query_semantics_sha256','adapter_sha256'}
_PROFILE_BINDINGS={'profile','profile_policy_sha256','profile_catalog_sha256','regulatory_index_sha256'}
_STATUSES = frozenset({'positive','no_match','uncertain','unsupported'})
_ROLES = frozenset({'population_joint','population_all_sites','population_adapted',
    'clinical_query','clinical_association','clinical_adapted'})
_MEMORY_BYTES = 64*1024*1024
_MEMORY_ENTRIES = 1024


def _sha(value):
    return type(value) is str and len(value)==64 and all(c in '0123456789abcdef' for c in value)


def _freeze(value):
    if type(value) is dict:
        return MappingProxyType({k:_freeze(v) for k,v in value.items()})
    if type(value) is list:
        return tuple(_freeze(v) for v in value)
    return value


def _paths(value):
    if isinstance(value,Path):
        return str(value)
    if isinstance(value,dict):
        return {k:_paths(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):
        return [_paths(v) for v in value]
    return value


def adapter_identity():
    root=Path(__file__).parent
    return hashlib.sha256(canonical({name:hashlib.sha256((root/name).read_bytes()).hexdigest()
        for name in ('wgs_allele_query_cache.py','population_annotation.py','wgs_clinical_session.py',
                     'population_clinical_session.py','wgs_query_limits.py',
                     'clinical_annotation.py','admission_links.py')})).hexdigest()


def population_source_bindings(*,joint,all_sites,reference,limits):
    return _paths({key:_plain(value) for key,value in
        dict(joint=joint,all_sites=all_sites,reference=reference,limits=limits).items()})


def clinical_source_bindings(session):
    return _paths({key:_plain(value) for key,value in dict(indexes=session.indexes,
        open_targets=session.open_targets,clinvar=session.clinvar,reference=session.reference,
        limits=session.execution_config.public_queries,
        raw_record_owner=str(session.scratch_path/'raw-record-cache.sqlite')).items()})


def _binding(content):
    try:
        data=json.loads(content)
        if type(data) is not dict:
            raise ValueError
        expected=_BINDINGS | _PROFILE_BINDINGS if data.get('version')==2 else _BINDINGS
        if (type(content) is not bytes or not 0<len(content)<=4*1024*1024
                or canonical(data)!=content or set(data)!=expected
                or data['schema']!='genosapiens.wgs-case-allele-query-cache' or data['version'] not in (1,2)
                or type(data['version']) is not int
                or data['version']==2 and data['profile'] not in {'focused','exhaustive'}
                or any(not _sha(v) for k,v in data.items() if k.endswith('_sha256'))
                or any(type(data[k]) is not str or not data[k] for k in ('reference_role','reference_identity_id'))
                or any(type(data[k]) is not dict or not data[k] for k in ('population_sources','clinical_sources'))):
            raise ValueError
        return data
    except (TypeError,ValueError,RecursionError):
        raise CaseAlleleQueryCacheError('case query cache binding is invalid') from None


@dataclass(frozen=True,slots=True,repr=False)
class QueryFact:
    role: str
    key_sha256: str
    status: str
    value: object
    size_bytes: int
    sha256: str


@dataclass(frozen=True,slots=True,repr=False)
class _MemoryFact:
    fact: QueryFact
    key_bytes: bytes
    size_bytes: int
    population: tuple | None = None


def _resident_bytes(*values):
    """Conservative retained object charge, including proxy backing tables.

    Charge each reachable object once per entry, its allocator allowance, and
    mapping table slots. Shared objects in separate entries are charged again.
    Whole-container measurements still cover transient decoder/SQLite/IPC work.
    """
    pending=list(values);seen=set();size=128
    while pending:
        value=pending.pop()
        if id(value) in seen:continue
        seen.add(id(value))
        size+=sys.getsizeof(value)+64
        if value is None or type(value) in (str,bytes,int,float,bool):continue
        if isinstance(value,Mapping):
            size+=256+72*len(value)
            pending.extend(value.keys());pending.extend(value.values())
        elif is_dataclass(value) and not isinstance(value,type):
            size+=sys.getsizeof(getattr(value,'__dict__',{}))
            pending.extend(getattr(value,item.name) for item in fields(value))
        elif type(value) in (tuple,list):
            pending.extend(value)
    return size


class CaseAlleleQueryCache:
    """One disk owner and a 64 MiB/1,024-entry immutable decoded-fact LRU.

    The memory bound charges canonical payloads plus retained plain/typed object
    overhead. SQLite/transient overhead also needs whole-process measurements.
    Cold writes persist without entering the LRU; a subsequent read admits the
    fact. Larger valid facts bypass memory retention, never validation.
    """
    def __init__(self,path,*,binding,execution_config):
        self.parents=_binding(binding)
        self.binding,self.binding_sha256=binding,hashlib.sha256(binding).hexdigest()
        self.current_adapter_identity=adapter_identity()
        self.config=execution_config
        self.path=Path(path)
        self.connection=None
        self.memory=OrderedDict()
        self.memory_bytes=0
        self.pending_writes=0
        self.max_fact_bytes=execution_config.public_queries.max_response_bytes
        self.metrics=dict(hits=0,misses=0,writes=0,memory_hits=0,max_memory_bytes=0,max_memory_entries=0,
            population_typed_decodes=0,population_typed_hits=0,population_typed_writes=0,population_adaptations=0)
        try:
            path=self.path
            if not path.is_absolute() or any(p.is_symlink() for p in (path,*path.parents)):
                raise ValueError
            parent=path.parent.stat()
            if parent.st_uid!=os.getuid() or stat.S_IMODE(parent.st_mode)!=0o700 or not stat.S_ISDIR(parent.st_mode):
                raise ValueError
            existing=path.exists()
            if existing:
                metadata=path.lstat()
                if (not stat.S_ISREG(metadata.st_mode) or metadata.st_uid!=os.getuid()
                        or stat.S_IMODE(metadata.st_mode)!=0o600 or metadata.st_nlink!=1):
                    raise ValueError
            else:
                fd=os.open(path,os.O_WRONLY|os.O_CREAT|os.O_EXCL|os.O_NOFOLLOW,0o600)
                os.close(fd)
            db=sqlite3.connect(path)
            self.connection=db
            db.enable_load_extension(False)
            db.execute('PRAGMA trusted_schema=OFF')
            db.execute('PRAGMA mmap_size=0')
            db.execute('PRAGMA temp_store=FILE')
            db.execute(f'PRAGMA cache_size=-{max(1,execution_config.sqlite.cache_bytes//1024)}')
            db.execute('PRAGMA journal_mode=DELETE')
            db.execute('PRAGMA synchronous=FULL')
            if not existing:
                db.executescript('''
                    PRAGMA user_version=1;
                    CREATE TABLE binding(singleton INTEGER PRIMARY KEY CHECK(singleton=1),content BLOB NOT NULL);
                    CREATE TABLE facts(role TEXT NOT NULL,key_sha256 TEXT NOT NULL,key_bytes BLOB NOT NULL,
                        status TEXT NOT NULL,content BLOB NOT NULL,size_bytes INTEGER NOT NULL,sha256 TEXT NOT NULL,
                        PRIMARY KEY(role,key_sha256));
                ''')
                db.execute('INSERT INTO binding VALUES (1,?)',(binding,))
                db.commit()
            if (db.execute('PRAGMA user_version').fetchone()!=(1,)
                    or db.execute('SELECT content FROM binding WHERE singleton=1').fetchone()!=(binding,)):
                raise ValueError
            self.data_version=db.execute('PRAGMA data_version').fetchone()[0]
        except (OSError,sqlite3.Error,ValueError):
            self.close()
            raise CaseAlleleQueryCacheError('case query cache is unavailable or foreign') from None

    def _active(self):
        if self.connection is None:
            raise CaseAlleleQueryCacheError('case query cache is closed')
        if self.connection.execute('PRAGMA data_version').fetchone()[0]!=self.data_version:
            self.memory.clear()
            self.memory_bytes=0
            raise CaseAlleleQueryCacheError('case query cache changed outside its owner')

    def _key(self,role,key):
        if role not in _ROLES:
            raise CaseAlleleQueryCacheError('case query cache role is invalid')
        content=canonical(_plain(key))
        if not 0<len(content)<=self.max_fact_bytes:
            raise CaseAlleleQueryCacheError('case query identity exceeds capacity')
        return content,hashlib.sha256(content).hexdigest()

    def _remember(self,fact,key_bytes,*,population=None):
        identity=(fact.role,fact.key_sha256)
        size=fact.size_bytes+len(key_bytes)+_resident_bytes(fact,key_bytes,population)
        previous=self.memory.pop(identity,None)
        if previous is not None:self.memory_bytes-=previous.size_bytes
        while self.memory and (self.memory_bytes>_MEMORY_BYTES or
                size<=_MEMORY_BYTES and (len(self.memory)>=_MEMORY_ENTRIES or self.memory_bytes+size>_MEMORY_BYTES)):
            _,old=self.memory.popitem(last=False)
            self.memory_bytes-=old.size_bytes
        if size>_MEMORY_BYTES:return
        self.memory[identity]=_MemoryFact(fact,key_bytes,size,population)
        self.memory_bytes+=size
        self.metrics['max_memory_bytes']=max(self.metrics['max_memory_bytes'],self.memory_bytes)
        self.metrics['max_memory_entries']=max(self.metrics['max_memory_entries'],len(self.memory))

    def get(self,role,key):
        self._active()
        key_bytes,key_sha=self._key(role,key)
        identity=(role,key_sha)
        found=self.memory.get(identity)
        if found is not None:
            fact=found.fact
            if found.key_bytes!=key_bytes:
                raise CaseAlleleQueryCacheError('case query cache key conflicts')
            self.memory.move_to_end(identity)
            self.metrics['hits']+=1
            self.metrics['memory_hits']+=1
            return fact
        try:
            row=self.connection.execute('SELECT key_bytes,status,content,size_bytes,sha256 FROM facts '
                'WHERE role=? AND key_sha256=?',identity).fetchone()
            if row is None:
                self.metrics['misses']+=1
                return None
            saved_key,status,blob,size,digest=row
            if (saved_key!=key_bytes or status not in _STATUSES or type(size) is not int
                    or not 0<size<=self.max_fact_bytes or not _sha(digest)
                    or type(blob) is not bytes or len(blob)>self.max_fact_bytes+65536):
                raise ValueError
            decoder=zlib.decompressobj()
            content=decoder.decompress(blob,self.max_fact_bytes+1)
            if (len(content)!=size or not decoder.eof or decoder.unused_data or decoder.unconsumed_tail
                    or hashlib.sha256(content).hexdigest()!=digest):
                raise ValueError
            value=json.loads(content)
            if canonical(value)!=content:
                raise ValueError
            fact=QueryFact(role,key_sha,status,_freeze(value),size,digest)
            self._remember(fact,key_bytes)
            self.metrics['hits']+=1
            return fact
        except (sqlite3.Error,ValueError,TypeError,zlib.error,RecursionError):
            raise CaseAlleleQueryCacheError('case query cache fact is corrupt') from None

    def population_adaptation(self,key):
        """Return strict immutable population copies under this exact bound key.

        There is one LRU entry for the plain fact and optional typed view. Disk
        reopens and evictions always perform a fresh strict typed decode.
        """
        from genosapiens.evidence_model import PopulationEvidence,_typed
        fact=self.get('population_adapted',key)
        if fact is None:return None
        identity=(fact.role,fact.key_sha256)
        entry=self.memory.get(identity)
        if entry is not None and entry.population is not None:
            self.metrics['population_typed_hits']+=1
            return fact,entry.population
        try:
            rows=_typed(_plain(fact.value),tuple[PopulationEvidence,...],decode=True)
        except (ValueError,TypeError,RuntimeError,RecursionError):
            raise CaseAlleleQueryCacheError('case population fact is corrupt') from None
        self.metrics['population_typed_decodes']+=1
        key_bytes,_=self._key('population_adapted',key)
        self._remember(fact,key_bytes,population=rows)
        return fact,rows

    def _put_population_adaptation(self,key,*,status,rows):
        """Return the adapter's strict copied view when creating a cold fact.

        Generic disk/memory insertion still authenticates canonical bytes and
        conflicts. Reopened/evicted facts still use the ordinary strict decoder.
        Cold copies belong to the synchronous adapter call, not the LRU. The
        first reuse strictly decodes the persisted fact before retaining it.
        """
        from genosapiens.evidence_model import PopulationEvidence,_typed,_validation_session
        with _validation_session():
            try:
                copied=_typed(rows,tuple[PopulationEvidence,...])
            except (ValueError,TypeError,RuntimeError,RecursionError):
                raise CaseAlleleQueryCacheError('case population fact is corrupt') from None
            fact=self._put_fact('population_adapted',key,status=status,value=copied,population=copied)
            self.metrics['population_typed_writes']+=1
            return fact,copied

    def put(self,role,key,*,status,value):
        return self._put_fact(role,key,status=status,value=value)

    def _put_fact(self,role,key,*,status,value,population=None):
        # Only _put_population_adaptation supplies its strictly copied view.
        # Persist cold facts without charging/retaining a memory entry that may
        # never be reused. Existing facts have already crossed the read boundary.
        self._active()
        key_bytes,key_sha=self._key(role,key)
        content=canonical(_plain(value))
        if status not in _STATUSES or not 0<len(content)<=self.max_fact_bytes:
            raise CaseAlleleQueryCacheError('case query fact exceeds capacity or has invalid status')
        digest=hashlib.sha256(content).hexdigest()
        previous=self.get(role,key)
        if previous is not None:
            if (previous.status,previous.size_bytes,previous.sha256)!=(status,len(content),digest):
                raise CaseAlleleQueryCacheError('case query facts conflict')
            if population is not None:self._remember(previous,key_bytes,population=population)
            return previous
        try:
            self.connection.execute('INSERT INTO facts VALUES (?,?,?,?,?,?,?)',
                (role,key_sha,key_bytes,status,zlib.compress(content,level=1),len(content),digest))
            self.pending_writes+=1
            if self.pending_writes>=self.config.sqlite.cursor_batch_rows:
                self.flush()
            fact=QueryFact(role,key_sha,status,_freeze(json.loads(content)),len(content),digest)
            self.metrics['writes']+=1
            return fact
        except (sqlite3.Error,ValueError,TypeError):
            raise CaseAlleleQueryCacheError('case query fact write failed') from None

    def population_results(self,requests,*,session):
        """The existing owner rebuilds the current ordered report from facts."""
        self._active()
        if (self.parents['adapter_sha256']!=self.current_adapter_identity
                or self.parents['reference_manifest_sha256']!=session.reference.manifest_sha256
                or self.parents['population_sources']!=population_source_bindings(joint=session.joint,
                    all_sites=session.all_sites,reference=session.reference,limits=session.limits)):
            raise CaseAlleleQueryCacheError('population cache session binding changed')
        result=session.query(requests,_case_cache=self)
        self.flush()
        return result

    def clinical_results(self,requests,*,session,_budget=None):
        """Share strict clinical facts; role/batch reports stay in the session."""
        self._active()
        if (self.parents['adapter_sha256']!=self.current_adapter_identity
                or self.parents['reference_manifest_sha256']!=session.reference.manifest_sha256
                or self.parents['clinical_sources']!=clinical_source_bindings(session)):
            raise CaseAlleleQueryCacheError('clinical cache session binding changed')
        result=session._case_query_results(requests,cache=self,budget=_budget)
        self.flush()
        return result

    def flush(self):
        if self.pending_writes:
            try:
                self.connection.commit()
                self.pending_writes=0
            except sqlite3.Error:
                raise CaseAlleleQueryCacheError('case query transaction failed') from None

    def close(self):
        self.memory.clear()
        self.memory_bytes=0
        if self.connection is not None:
            try:
                self.flush()
            finally:
                self.connection.close()
                self.connection=None

    def __enter__(self):return self
    def __exit__(self,kind,*_):
        if kind is not None and self.connection is not None:
            self.connection.rollback()
            self.pending_writes=0
        self.close()
