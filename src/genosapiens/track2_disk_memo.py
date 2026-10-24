"""Private invocation memoization with one decoded value and a bounded cache."""
import json
from pathlib import Path
import sqlite3
import tempfile
from genosapiens.track2_scientific_storage import json_chunks


class DiskMemo:
    def __init__(self):
        self.temp=tempfile.TemporaryDirectory(prefix='track2-evidence-memo-')
        path=Path(self.temp.name)/'memo.sqlite';self.db=sqlite3.connect(path);path.chmod(0o600)
        self.db.execute('PRAGMA cache_size=-2048');self.db.execute('PRAGMA temp_store=FILE')
        self.db.execute('CREATE TABLE memo (key TEXT PRIMARY KEY,value BLOB)')
    def __enter__(self):return self
    def __exit__(self,*args):self.close()
    def close(self):self.db.close();self.temp.cleanup()
    def _key(self,key):return json.dumps(key,sort_keys=True,separators=(',',':'))
    def __contains__(self,key):return self.db.execute('SELECT 1 FROM memo WHERE key=?',(self._key(key),)).fetchone() is not None
    def __getitem__(self,key):
        row=self.db.execute('SELECT value FROM memo WHERE key=?',(self._key(key),)).fetchone()
        if row is None:raise KeyError(key)
        return json.loads(row[0])
    def __setitem__(self,key,value):self.db.execute('INSERT OR REPLACE INTO memo VALUES (?,?)',(self._key(key),b''.join(json_chunks(value))))
    def get(self,key,default=None):return self[key] if key in self else default
