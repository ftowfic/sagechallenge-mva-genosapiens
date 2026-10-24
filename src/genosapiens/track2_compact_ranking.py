"""Store evidence once. Sort exact scalar keys and stable record references."""
from fractions import Fraction
import json
from pathlib import Path
import sqlite3
import tempfile

from genosapiens.track2_scientific_storage import Records,json_chunks
from genosapiens.track2_scientific_records import stored_json


def encoded(values):
    return json.dumps([{'fraction':[v.numerator,v.denominator]} if isinstance(v,Fraction) else v for v in values],separators=(',',':'))


def compare(left,right):
    def decode(value):
        return tuple(Fraction(*v['fraction']) if isinstance(v,dict) else v for v in json.loads(value))
    a,b=decode(left),decode(right)
    return (a>b)-(a<b)


class CompactRankStore:
    """One evidence record per insertion, with a 2 MiB disk-backed sort cache."""
    def __init__(self):
        self.temp=tempfile.TemporaryDirectory(prefix='track2-compact-ranks-')
        path=Path(self.temp.name)/'ranks.sqlite';self.db=sqlite3.connect(path);path.chmod(0o600)
        self.db.execute('PRAGMA cache_size=-2048');self.db.execute('PRAGMA temp_store=FILE')
        self.db.create_collation('exact',compare)
        self.db.execute('CREATE TABLE evidence (id INTEGER PRIMARY KEY,data BLOB)')
        self.db.execute('CREATE TABLE keys (group_id TEXT,entity TEXT,record_id INTEGER,score TEXT COLLATE exact,display TEXT COLLATE exact)')
        self.db.execute('CREATE TABLE unique_rows (identity TEXT PRIMARY KEY)')
        self.serialized_records=0

    def __enter__(self):return self
    def __exit__(self,*args):self.db.close();self.temp.cleanup()

    def add(self,row, *, unique=None):
        if unique is not None:
            try:self.db.execute('INSERT INTO unique_rows VALUES (?)',(json.dumps(unique),))
            except sqlite3.IntegrityError:raise ValueError('A derived objective score occurs more than once.') from None
        cursor=self.db.execute('INSERT INTO evidence(data) VALUES (?)',(b''.join(json_chunks(row)),))
        self.serialized_records+=1
        return cursor.lastrowid

    def key(self,group,entity,record,score,display):
        self.db.execute('INSERT INTO keys VALUES (?,?,?,?,?)',(group,entity,record,encoded(score),encoded(display)))

    def finish(self, *, representatives):
        self.db.commit()
        selection=('SELECT * FROM (SELECT *,ROW_NUMBER() OVER (PARTITION BY group_id,entity ORDER BY display,record_id) AS chosen FROM keys) WHERE chosen=1'
            if representatives else 'SELECT * FROM keys')
        self.db.execute('CREATE TABLE ranked AS SELECT group_id,record_id,display,'
            'RANK() OVER (PARTITION BY group_id ORDER BY score) AS rank,'
            'COUNT(*) OVER (PARTITION BY group_id,score) AS ties FROM ('+selection+')')
        self.db.execute('CREATE INDEX ranked_record ON ranked(record_id)')

    def rows(self,group, *, derived=False):
        for data,rank,ties in self.db.execute('SELECT e.data,r.rank,r.ties FROM ranked r JOIN evidence e ON e.id=r.record_id WHERE r.group_id=? ORDER BY r.display COLLATE exact,r.record_id',(group,)):
            row=json.loads(data)
            if derived:row.update(method_rank=rank,tie_count=ties)
            else:row.update(rank=rank)
            yield row

    def complete(self):
        for data,rank,ties in self.db.execute('SELECT e.data,r.rank,r.ties FROM evidence e LEFT JOIN ranked r ON e.id=r.record_id ORDER BY e.id'):
            row=json.loads(data)
            if rank is not None:row.update(method_rank=rank,tie_count=ties)
            yield row


def select_compact_tables(scores, *, tissue):
    from genosapiens.track2_report_model import presentation_policy
    from genosapiens.track2_tissue_ranking import exact_score
    policy=presentation_policy(tissue=tissue);result={}
    with CompactRankStore() as store:
        for row in scores:
            record=store.add(row)
            for table,metric in policy['metrics'].items():
                admitted=(row['lincs_nominated'] or row['mechanism_nominated']) if table=='combined' else row[table+'_nominated']
                if not admitted:continue
                value=exact_score(row,metric) if tissue else row['half_points'][metric]
                score=(value is None,-value if value is not None else Fraction(0))
                display=(*score,row['name'].casefold(),row['canonical_active_ingredient_id'],row['objective_id'])
                store.key(table,row['canonical_active_ingredient_id'],record,score,display)
        store.finish(representatives=True)
        for table,metric in policy['metrics'].items():
            ranking=stored_json(Records(dict(row,metric=metric) for row in store.rows(table)))
            rows=ranking[:policy['main_table_limit']]
            cutoff=(exact_score(rows[-1],metric) if tissue else rows[-1]['half_points'][metric]) if rows else None
            ties=stored_json(Records(dict(canonical_active_ingredient_id=row['canonical_active_ingredient_id'],objective_id=row['objective_id'],rank=row['rank'])
                for i,row in enumerate(ranking) if i>=len(rows) and (exact_score(row,metric) if tissue else row['half_points'][metric])==cutoff))
            result[table]=dict(rows=rows,ranking=ranking,cutoff_ties=ties,metric=metric)
            result[table].update(cutoff_exact=rows[-1]['exact_scores'][metric] if rows else None) if tissue else result[table].update(cutoff_half_points=cutoff)
    return result
