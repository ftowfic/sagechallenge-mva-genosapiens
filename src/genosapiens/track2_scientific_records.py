"""Private disk-backed collections for scientific grouping and ordered replay."""
import json
import os
from pathlib import Path
import sqlite3
import tempfile

from genosapiens.track2_scientific_storage import Records,json_chunks,publish_json
from genosapiens.track2_scientific_reader import read_scientific


def stored_json(value):
    """Retain reusable JSON views over a complete temporary scientific product."""
    temporary=tempfile.TemporaryDirectory(prefix='track2-records-')
    try:
        path=Path(temporary.name)/'records.json'
        publish_json(path,value)
        result=read_scientific(path)
        result.document.source_owner=temporary
        return result
    except BaseException:
        temporary.cleanup()
        raise


def ordered_records(rows,key):
    """Stable external ordering with one record and a 2 MiB SQLite cache.

    Each key position contains comparable SQLite scalars. Scientific owners
    supply their existing ordering keys. No record is dropped or pooled.
    """
    with tempfile.TemporaryDirectory(prefix='track2-order-') as temporary:
        connection=sqlite3.connect(Path(temporary)/'order.sqlite')
        (Path(temporary)/'order.sqlite').chmod(0o600)
        connection.execute('PRAGMA cache_size=-2048');connection.execute('PRAGMA temp_store=FILE')
        try:
            with (Path(temporary)/'data.jsonl').open('w+b') as stream:
                os.fchmod(stream.fileno(),0o600)
                columns=None
                for index,row in enumerate(rows):
                    values=tuple(key(row))
                    if columns is None:
                        columns=','.join('key'+str(i) for i in range(len(values)))
                        if not columns:raise ValueError('Scientific ordering requires a key.')
                        connection.execute('CREATE TABLE records (ordinal INTEGER,start INTEGER,length INTEGER,'+columns+')')
                    start=stream.tell()
                    for chunk in json_chunks(row):stream.write(chunk)
                    parameters=(index,start,stream.tell()-start,*values)
                    connection.execute('INSERT INTO records VALUES ('+','.join('?' for _ in parameters)+')',parameters)
                stream.flush();connection.commit()
                if columns is not None:
                    for start,length in connection.execute('SELECT start,length FROM records ORDER BY '+columns+',ordinal'):
                        stream.seek(start)
                        yield json.loads(stream.read(length))
        finally:connection.close()
