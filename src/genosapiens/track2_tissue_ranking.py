"""Exact rational ranking with bounded external sorting and unchanged mechanism points."""
from fractions import Fraction
from itertools import groupby
import json
from pathlib import Path
import sqlite3
import tempfile

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_scientific_storage import Records, json_chunks
from genosapiens.track2_scientific_records import stored_json
from genosapiens.track2_tissue_records import validate_tissue_record


def exact_score(row,metric):
    value=row['exact_scores'][metric]
    return Fraction(value['numerator'],value['denominator']) if value['denominator'] else None


def tissue_report_score(rows,objective,coverage):
    from genosapiens.track2_report_ranking import score_objective
    rows=list(rows); tissue=rows[0]['tissue_score']
    validate_tissue_record(tissue,'score')
    require(all(r['tissue_score']==tissue for r in rows), 'A report objective has inconsistent tissue scores.')
    require(all(tissue[k]==objective[k] for k in ('gene_id','scenario_id','objective_id')),
        'Tissue report evidence belongs to another objective.')
    base=score_objective([{k:v for k,v in r.items() if k!='tissue_score'} | {'cells':[]} for r in rows],objective,coverage)
    mechanism=Fraction(base['half_points']['S_M'],2)
    lincs=Fraction(tissue['numerator'],tissue['denominator']) if tissue['denominator'] else None
    total=(lincs or Fraction(0))+mechanism
    base.pop('half_points');base.pop('counted_cells')
    base.update(schema='genosapiens.track2-report-objective-score/v2',tissue_score=tissue,
        components=dict(I_L=int(base['lincs_nominated']),I_M=base['components']['I_M'],M=base['components']['M'],
            P=tissue['P'],C=tissue['C'],N=tissue['N']),
        E=int(base['lincs_nominated'])+int(base['mechanism_nominated']),
        scores=dict(S_L=float(lincs) if lincs is not None else None,S_M=float(mechanism),S=float(total)),
        exact_scores={key:dict(numerator=value.numerator,denominator=value.denominator) if value is not None
            else dict(numerator=0,denominator=0) for key,value in [('S_L',lincs),('S_M',mechanism),('S',total)]},
        provisional=tissue['provisional'] or any(r['state']=='failed' for r in base['coverage']))
    return base


def ordered_exact(rows,key):
    """SQLite compares complete rational keys without rounding or retaining all rows."""
    def encode(values):
        return json.dumps([{'fraction':[v.numerator,v.denominator]} if isinstance(v,Fraction) else v for v in values])
    def decode(value):
        return tuple(Fraction(*v['fraction']) if isinstance(v,dict) else v for v in json.loads(value))
    def compare(a,b):
        left,right=decode(a),decode(b)
        return (left>right)-(left<right)
    with tempfile.TemporaryDirectory(prefix='track2-tissue-ranks-') as temporary:
        connection=sqlite3.connect(Path(temporary)/'sort.sqlite')
        connection.execute('PRAGMA cache_size=-2048');connection.execute('PRAGMA temp_store=FILE')
        connection.create_collation('rational',compare)
        connection.execute('CREATE TABLE rows (ordinal INTEGER, key TEXT COLLATE rational, start INTEGER, size INTEGER)')
        try:
            with (Path(temporary)/'data.jsonl').open('w+b') as stream:
                for index,row in enumerate(rows):
                    start=stream.tell()
                    for chunk in json_chunks(row): stream.write(chunk)
                    connection.execute('INSERT INTO rows VALUES (?,?,?,?)',(index,encode(key(row)),start,stream.tell()-start))
                stream.flush();connection.commit()
                for start,size in connection.execute('SELECT start,size FROM rows ORDER BY key,ordinal'):
                    stream.seek(start);yield json.loads(stream.read(size))
        finally: connection.close()


def select_tissue_tables(scores):
    from genosapiens.track2_compact_ranking import select_compact_tables
    return select_compact_tables(scores,tissue=True)
