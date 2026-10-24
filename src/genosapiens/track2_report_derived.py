"""Derive report strengths from retained conditions without matrix access."""
from contextlib import contextmanager
import json
from pathlib import Path
import sqlite3
import tempfile

from genosapiens.track2_derived_ranking import rank_methods
from genosapiens.track2_derived_strength import METHODS, aggregate_strength
from genosapiens.track2_frozen_calibration import canonical, identity, query_content
from genosapiens.track2_operator_inputs import require


@contextmanager
def condition_index(rows):
    """Store each unique retained observation once in a bounded SQLite cache."""
    with tempfile.TemporaryDirectory(prefix='track2-derived-conditions-') as temporary:
        connection=sqlite3.connect(Path(temporary)/'conditions.sqlite')
        connection.execute('PRAGMA cache_size=-2048')
        connection.execute('PRAGMA temp_store=FILE')
        connection.execute('CREATE TABLE conditions (query TEXT, signature TEXT, ingredient TEXT, gene TEXT, method TEXT, tissue TEXT, record BLOB, PRIMARY KEY (query,signature))')
        try:
            for row in rows:
                observation={k:v for k,v in row.items() if k not in {'objective_id','scenario_id'}}
                data=canonical(observation)
                key=row['query_sha256'],row['sig_id']
                prior=connection.execute('SELECT record FROM conditions WHERE query=? AND signature=?',key).fetchone()
                if prior:
                    require(prior[0]==data, 'Repeated retained observations disagree.')
                else:
                    connection.execute('INSERT INTO conditions VALUES (?,?,?,?,?,?,?)',
                        (*key,row['medication_entity_key'],row['gene_id'],row['method'],row['tissue_id'],data))
            connection.execute('CREATE INDEX ingredient_scope ON conditions (ingredient,gene,method)')
            connection.commit()
            yield connection
        finally:connection.close()


def derive_rankings(*, conditions, evaluated_scores, objectives, historical_lincs, queries=()):
    """Retain all evaluated objective/method scores, including nonnominees.

    This operation computes no random backgrounds. Calibrated enrichment is an
    explicit subsequent operation with independently admitted resource identities.
    """
    objectives={row['objective_id']:row for row in objectives}
    query_keys={}
    for query in queries:
        key=identity(query_content(query))
        require(query_keys.setdefault(query['query_sha256'],key)==key, 'A retained query identity has conflicting content.')
    representatives={}
    for row in historical_lincs:
        tissue=row.get('tissue_score')
        require(tissue is not None, 'Derived tissue scoring cannot use legacy cell representatives.')
        for method in METHODS:
            representatives[tissue['medication_entity_key'],method]=row['objective_id']
    from genosapiens.track2_disk_memo import DiskMemo
    from genosapiens.track2_scientific_records import stored_json
    from genosapiens.track2_scientific_storage import Records
    nominated=set()
    with condition_index(conditions) as connection, DiskMemo() as aggregates, DiskMemo() as fingerprints:
        def scores():
            for source in evaluated_scores:
                oid=source['objective_id'];objective=objectives.get(oid)
                require(objective is not None and all(source[k]==objective[k] for k in ('gene_id','scenario_id')),
                        'A derived source score has another objective identity.')
                for method in METHODS:
                    coverage=[r for r in source['coverage'] if r['method']==method]
                    if not coverage:continue
                    tissues=sorted({r['tissue_id'] for r in coverage})
                    selected={r['query_sha256'] for r in coverage if r.get('query_sha256')}
                    def observations():
                        for (raw,) in connection.execute('SELECT record FROM conditions WHERE ingredient=? AND gene=? AND method=? ORDER BY query,signature',
                                (source['medication_entity_key'],source['gene_id'],method)):
                            row=json.loads(raw)
                            if row['query_sha256'] in selected and row['tissue_id'] in tissues:yield row
                    incomplete=any(r.get('incomplete_jobs') for r in coverage)
                    cache_key=(source['medication_entity_key'],source['gene_id'],method,objective['desired_modulation'],
                        tuple(tissues),tuple(sorted(selected)),incomplete)
                    if cache_key not in aggregates:
                        aggregates[cache_key]=aggregate_strength(observations(),method=method,direction=objective['desired_modulation'],
                            declared_tissues=tissues,incomplete=incomplete)
                    score=aggregates[cache_key]
                    if query_keys and cache_key not in fingerprints:
                        from genosapiens.track2_calibration_bank import CONDITION_FIELDS
                        scientific={}
                        for observation in observations():
                            require(observation['query_sha256'] in query_keys, 'A retained condition has no scientific query.')
                            value={field:observation[field] for field in CONDITION_FIELDS}
                            value['query_key']=query_keys[observation['query_sha256']]
                            key=value['query_key'],value['sig_id']
                            require(scientific.setdefault(key,value)==value, 'Scientific query aliases disagree.')
                        fingerprints[cache_key]=identity([scientific[key] for key in sorted(scientific)])
                    condition_scope_sha256=fingerprints.get(cache_key)
                    for field in ('n_conditions','n_valid','n_undefined'):
                        if all(field in r for r in coverage):
                            # Unsupported directions have no valid derived coefficients.
                            if field != 'n_conditions' and objective['desired_modulation'] != 'increase':continue
                            require(score[field] == sum(r[field] for r in coverage),
                                    'Derived coverage differs from the retained condition inventory.')
                    row=dict(score,medication_entity_key=source['medication_entity_key'],
                        canonical_active_ingredient_id=source['canonical_active_ingredient_id'],name=source['name'],
                        objective_id=oid,scenario_id=source['scenario_id'],gene_id=source['gene_id'],
                        condition_scope_sha256=condition_scope_sha256,
                        numerator=None,denominator=None,percentile=None,entry_key=None,distribution_sha256=None)
                    if any(r['support'] or r['weaker_nomination'] for r in coverage):
                        nominated.add((row['medication_entity_key'],oid,method))
                    yield row
        rows=stored_json(Records(scores()))
    return rank_methods(rows,source_selected=representatives,nominated=nominated)


def report_conditions(root):
    """Read the supported native report's retained job streams only."""
    from genosapiens.track2_figure_stream import condition_rows
    root=Path(root)
    for path in sorted((root/'lincs').glob('*/*/science/conditions.jsonl')):
        require(path.resolve()==path.absolute(), 'A retained condition stream contains a symbolic link.')
        yield from condition_rows(path)


def calibrate_rankings(derived, *, bank, resource_contract, expected_resource_sha256,
                       expected_fda_snapshot_sha256, expected_source_release):
    """Use exact admitted scopes, leaving missing entries explicitly uncalibrated."""
    import hashlib
    from genosapiens.track2_calibration_bank import _digest, _path, admit_entry, percentile
    from genosapiens.track2_derived_strength import policy
    from genosapiens.track2_report_calibration import _resource
    _digest(expected_fda_snapshot_sha256)
    _, bindings = _resource(resource_contract, expected_resource_sha256)
    module_root = Path(__file__).resolve().parent
    code = {name:hashlib.sha256((module_root/name).read_bytes()).hexdigest() for name in
            ('track2_derived_strength.py','track2_frozen_calibration.py','track2_lincs_scoring.py')}
    from genosapiens.track2_disk_memo import DiskMemo
    from genosapiens.track2_calibration_bank import AdmittedBackground
    from dataclasses import asdict
    with DiskMemo() as scopes:
        for key,binding in bindings.items():
            if not _path(bank,'entries/'+key).exists():continue
            background=admit_entry(bank,expected_entry_key=key,expected_distribution_sha256=binding['distribution_sha256'])
            scope=json.loads(_path(Path(bank)/'entries'/key,'scope.json').read_bytes())
            if (scope['fda_snapshot_sha256'] != expected_fda_snapshot_sha256
                    or scope['source_release'] != expected_source_release or scope['policy'] != policy()
                    or scope['scoring_code'] != code):
                continue
            matching=(scope['medication_entity_key'],scope['method'],scope['direction'],identity(scope['conditions']))
            require(matching not in scopes, 'Several resource entries claim the same report scope.')
            scopes[matching]=asdict(background)
        from genosapiens.track2_scientific_records import stored_json
        from genosapiens.track2_scientific_storage import Records
        representatives={};nominated=set()
        def calibrated_rows():
            for original in derived['complete_export']:
                row=dict(original)
                for field in ('numerator','denominator','percentile','entry_key','distribution_sha256'):
                    row[field]=None
                for field in ('lower_count','equal_count','degenerate_background','tie_count'):
                    row.pop(field,None)
                key=row['medication_entity_key'],row['objective_id'],row['method']
                if row['source_nominated']:nominated.add(key)
                if row['source_representative']:representatives[key[0],key[2]]=key[1]
                scope_key=row['medication_entity_key'],row['method'],row['direction'],row.get('condition_scope_sha256')
                if row['status']=='complete' and scope_key in scopes:
                    row.update(percentile(AdmittedBackground(**scopes[scope_key]),row['strength']))
                yield row
        rows=stored_json(Records(calibrated_rows()))
    return rank_methods(rows,source_selected=representatives,nominated=nominated)
