"""Read scientific JSON through a private disk index without a full object graph.

The parser retains a 64 KiB input buffer, one scalar, and container nesting.
Views retain node identities. Scalars load on demand. SQLite uses a 2 MiB cache.
Scientific owners still validate complete schemas, scope, ordering, and meaning.
"""
from collections.abc import Mapping
import copy
import json
import os
from pathlib import Path
import re
import sqlite3
import stat
import tempfile

from genosapiens.track2_scientific_storage import ScientificStorageError

BUFFER_BYTES = 64 * 1024
SPECIAL = re.compile(br'["\\]')
DELIMITER = re.compile(br'[\s,\[\]{}:]')


class _Document:
    def __init__(self, path):
        path=Path(path).absolute()
        if any(p.is_symlink() for p in (path,*path.parents)):
            raise ScientificStorageError('A scientific input contains a symbolic link.')
        descriptor=os.open(path,os.O_RDONLY|os.O_NOFOLLOW|os.O_NONBLOCK)
        if not stat.S_ISREG(os.fstat(descriptor).st_mode):
            os.close(descriptor)
            raise ScientificStorageError('A scientific input is not a regular file.')
        self.stream=os.fdopen(descriptor,'rb')
        self.signature=self._signature()
        self.temporary=tempfile.TemporaryDirectory(prefix='track2-json-')
        self.database=sqlite3.connect(Path(self.temporary.name)/'records.sqlite')
        (Path(self.temporary.name)/'records.sqlite').chmod(0o600)
        self.database.execute('PRAGMA cache_size=-2048')
        self.database.execute('PRAGMA temp_store=FILE')
        self.database.execute('CREATE TABLE nodes (id INTEGER PRIMARY KEY, parent INTEGER, ordinal INTEGER, name TEXT, kind TEXT, start INTEGER, length INTEGER)')
        self.database.execute('CREATE UNIQUE INDEX positions ON nodes(parent,ordinal)')
        self.database.execute('CREATE UNIQUE INDEX names ON nodes(parent,name) WHERE name IS NOT NULL')
        self.buffer=b'';self.position=0;self.base=0
        self.next_id=0

    def _signature(self):
        info=os.fstat(self.stream.fileno())
        return info.st_dev,info.st_ino,info.st_size,info.st_mtime_ns,info.st_ctime_ns

    def close(self):
        for name in ('database','stream','temporary','source_owner'):
            owner=getattr(self,name,None)
            if owner is not None:
                owner.cleanup() if name in ('temporary','source_owner') else owner.close()
                setattr(self,name,None)

    def __del__(self):
        try:self.close()
        except OSError:pass

    def _fill(self):
        if self.position==len(self.buffer):
            self.base+=self.position;self.position=0
            self.buffer=self.stream.read(BUFFER_BYTES)
        return bool(self.buffer)

    def _peek(self):
        return self.buffer[self.position] if self._fill() else None

    def _space(self):
        while self._peek() in (9,10,13,32):self.position+=1

    def _expect(self, value):
        self._space()
        if self._peek()!=value:raise ScientificStorageError('Scientific JSON has invalid structure.')
        self.position+=1

    def _raw(self, start, length):
        if self._signature()!=self.signature:
            raise ScientificStorageError('A scientific input changed during reading.')
        data=os.pread(self.stream.fileno(),length,start)
        if len(data)!=length:raise ScientificStorageError('A scientific input is incomplete.')
        return data

    def _string(self):
        start=self.base+self.position
        self.position+=1
        while self._fill():
            match=SPECIAL.search(self.buffer,self.position)
            if match is None:
                self.position=len(self.buffer);continue
            self.position=match.end()
            if match.group()==b'"':
                length=self.base+self.position-start
                # Validate escapes and UTF-8 for every scalar, even unused fields.
                return start,length,json.loads(self._raw(start,length))
            if not self._fill():break
            self.position+=1  # An escaped quote cannot close the string.
        raise ScientificStorageError('Scientific JSON has an incomplete string.')

    def _parse(self, parent=None, ordinal=0, name=None):
        self._space();start=self.base+self.position;token=self._peek()
        if token is None:raise ScientificStorageError('Scientific JSON is incomplete.')
        self.next_id+=1;identity=self.next_id
        kind='object' if token==123 else 'array' if token==91 else 'scalar'
        try:
            self.database.execute('INSERT INTO nodes VALUES (?,?,?,?,?,?,?)',
                (identity,parent,ordinal,name,kind,start,0))
        except sqlite3.IntegrityError:
            raise ScientificStorageError('Scientific JSON has duplicate object keys.') from None
        if kind!='scalar':
            self.position+=1;self._space();index=0;closing=125 if kind=='object' else 93
            if self._peek()!=closing:
                while True:
                    key=None
                    if kind=='object':
                        self._space()
                        if self._peek()!=34:raise ScientificStorageError('Scientific JSON requires object keys.')
                        _,_,key=self._string();self._expect(58)
                    self._parse(identity,index,key);index+=1;self._space()
                    if self._peek()==closing:break
                    self._expect(44)
            self._expect(closing)
        elif token==34:
            self._string()
        else:
            while self._fill():
                match=DELIMITER.search(self.buffer,self.position)
                self.position=match.start() if match else len(self.buffer)
                if match:break
            length=self.base+self.position-start
            raw=self._raw(start,length)
            value=json.loads(raw)
            if isinstance(value,(dict,list,str)) or raw in (b'NaN',b'Infinity',b'-Infinity'):
                raise ScientificStorageError('Scientific JSON has an invalid scalar.')
        self.database.execute('UPDATE nodes SET length=? WHERE id=?',(self.base+self.position-start,identity))
        return identity,kind

    def read(self):
        identity,kind=self._parse();self._space()
        if self._peek() is not None:raise ScientificStorageError('Scientific JSON has trailing content.')
        self.database.commit()
        if self._signature()!=self.signature:raise ScientificStorageError('A scientific input changed during parsing.')
        self.buffer=b''
        return self.value(identity,kind)

    def value(self, identity, kind, start=None, length=None):
        if kind=='object':return _Object(self,identity)
        if kind=='array':return _Array(self,identity)
        if start is None:start,length=self.database.execute('SELECT start,length FROM nodes WHERE id=?',(identity,)).fetchone()
        return json.loads(self._raw(start,length))

    def children(self, identity):
        cursor=self.database.execute('SELECT id,name,kind,start,length FROM nodes WHERE parent=? ORDER BY ordinal',(identity,))
        for child,name,kind,start,length in cursor:
            yield name,self.value(child,kind,start,length)


class _Object(dict):
    """Read-only dict-compatible view. Only object keys remain in memory."""
    def __new__(cls,document,identity=None):
        # dataclasses.asdict reconstructs containers through their type.
        # Such copies are ordinary records, never live disk views.
        return dict(document) if identity is None else super().__new__(cls)
    def __init__(self, document, identity):
        self.document=document;self.identity=identity
        super().__init__((row[0],None) for row in document.database.execute(
            'SELECT name FROM nodes WHERE parent=? ORDER BY ordinal',(identity,)))

    def __iter__(self):return super().__iter__()
    def __getitem__(self,key):
        row=self.document.database.execute('SELECT id,kind,start,length FROM nodes WHERE parent=? AND name=?',
            (self.identity,key)).fetchone()
        if row is None:raise KeyError(key)
        return self.document.value(*row)
    def get(self,key,default=None):return self[key] if key in self else default
    def items(self):return self.document.children(self.identity)
    def values(self):return (value for _,value in self.items())
    def copy(self):return dict(self.items())
    def __or__(self,other):return self.copy() | other
    def __ror__(self,other):return other | self.copy()
    def __eq__(self,other):
        return isinstance(other,Mapping) and self.keys()==other.keys() and all(value==other[key] for key,value in self.items())
    def __ne__(self,other):return not self==other
    def __deepcopy__(self,memo):return {key:copy.deepcopy(value,memo) for key,value in self.items()}
    def __repr__(self):return '<scientific JSON object>'
    def _readonly(self,*args,**kwargs):raise TypeError('Scientific JSON views are read-only.')
    __setitem__=__delitem__=clear=pop=popitem=setdefault=update=__ior__=_readonly


class _Array(list):
    """Read-only list-compatible view with one member decoded at a time."""
    def __new__(cls,document,identity=None):
        return list(document) if identity is None else super().__new__(cls)
    def __init__(self,document,identity):
        self.document=document;self.identity=identity
        self.length=document.database.execute('SELECT count(*) FROM nodes WHERE parent=?',(identity,)).fetchone()[0]
    def __len__(self):return self.length
    def __iter__(self):return (value for _,value in self.document.children(self.identity))
    def __getitem__(self,index):
        if isinstance(index,slice):return [self[i] for i in range(*index.indices(self.length))]
        if index<0:index+=self.length
        row=self.document.database.execute('SELECT id,kind,start,length FROM nodes WHERE parent=? AND ordinal=?',
            (self.identity,index)).fetchone()
        if row is None:raise IndexError(index)
        return self.document.value(*row)
    def __eq__(self,other):
        return isinstance(other,(list,tuple)) and len(self)==len(other) and all(a==b for a,b in zip(self,other))
    def __ne__(self,other):return not self==other
    def __contains__(self,value):return any(row==value for row in self)
    def count(self,value):return sum(row==value for row in self)
    def index(self,value,start=0,stop=None):
        for index in range(*slice(start,stop).indices(self.length)):
            if self[index]==value:return index
        raise ValueError('The value is absent from the scientific array.')
    def copy(self):return list(self)
    def __add__(self,other):return list(self)+other
    def __radd__(self,other):return other+list(self)
    def __reversed__(self):return (self[i] for i in range(self.length-1,-1,-1))
    def __deepcopy__(self,memo):return [copy.deepcopy(value,memo) for value in self]
    def __repr__(self):return '<scientific JSON array>'
    def _readonly(self,*args,**kwargs):raise TypeError('Scientific JSON views are read-only.')
    __setitem__=__delitem__=append=extend=insert=pop=remove=clear=reverse=sort=__iadd__=__imul__=_readonly


def read_scientific(path):
    document=_Document(path)
    try:return document.read()
    except BaseException:
        document.close()
        raise
