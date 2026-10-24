"""SQLite ingestion of release-bound public discovery metadata."""
from collections import defaultdict
import csv
import gzip
import json
import math
from pathlib import Path
import sqlite3

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_published_genotype import canonical, sha256

DDL = '''
PRAGMA user_version=1;
CREATE TABLE sources(source_id TEXT PRIMARY KEY, release TEXT, manifest_sha256 TEXT, schema TEXT, lineage_id TEXT);
CREATE TABLE genes(gene_id TEXT PRIMARY KEY, approved_symbol TEXT, aliases_json TEXT, source_id TEXT);
CREATE TABLE gene_aliases(alias TEXT, gene_id TEXT, resolution TEXT, PRIMARY KEY(alias,gene_id));
CREATE TABLE medications(ingredient_key TEXT PRIMARY KEY, data TEXT NOT NULL);
CREATE TABLE products(application TEXT, product_number TEXT, ingredient_key TEXT, data TEXT, PRIMARY KEY(application,product_number,ingredient_key));
CREATE TABLE compound_identity(provider TEXT, compound_id TEXT, ingredient_key TEXT, resolution TEXT, reason TEXT, source_json TEXT, data TEXT, PRIMARY KEY(provider,compound_id));
CREATE INDEX compound_ingredient ON compound_identity(ingredient_key,provider);
CREATE TABLE mechanisms(mechanism_id TEXT PRIMARY KEY, action TEXT, description TEXT, target_name TEXT, target_type TEXT, targets_json TEXT, references_json TEXT, source_json TEXT);
CREATE TABLE mechanism_origins(mechanism_id TEXT, partition TEXT, row_number INTEGER, PRIMARY KEY(partition,row_number));
CREATE TABLE mechanism_genes(mechanism_id TEXT, gene_id TEXT, group_warning TEXT, PRIMARY KEY(mechanism_id,gene_id));
CREATE INDEX mechanism_gene ON mechanism_genes(gene_id,mechanism_id);
CREATE TABLE mechanism_compounds(mechanism_id TEXT, chembl_id TEXT, PRIMARY KEY(mechanism_id,chembl_id));
CREATE INDEX mechanism_compound ON mechanism_compounds(chembl_id,mechanism_id);
CREATE TABLE lincs_signatures(sig_id TEXT PRIMARY KEY, pert_id TEXT, pert_type TEXT, pert_iname TEXT, gene_id TEXT, cell_id TEXT, pert_dose REAL, pert_dose_unit TEXT, pert_time REAL, pert_time_unit TEXT, distil_nsample REAL, distil_cc_q75 REAL, tas REAL, mapping_reason TEXT, source_json TEXT);
CREATE INDEX signature_gene ON lincs_signatures(gene_id,cell_id,pert_type);
CREATE INDEX signature_cell ON lincs_signatures(cell_id,pert_type,pert_id);
CREATE TABLE lincs_features(pr_gene_id TEXT PRIMARY KEY, symbol TEXT, measured_landmark INTEGER);
CREATE TABLE lincs_cells(cell_id TEXT PRIMARY KEY, annotations_json TEXT);
CREATE TABLE coverage(source TEXT, reason TEXT, count INTEGER, examples_json TEXT, PRIMARY KEY(source,reason));
'''


def encoded(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def parquet_rows(paths, required):
    import duckdb
    with duckdb.connect() as connection:
        for path in sorted(paths):
            cursor = connection.execute('SELECT * FROM read_parquet(?)', [str(path)])
            fields = [c[0] for c in cursor.description]
            require(set(required) <= set(fields), 'A required Parquet schema differs: '+Path(path).name)
            number = 0
            while rows := cursor.fetchmany(1024):
                for row in rows:
                    number += 1
                    yield Path(path).name, number, dict(zip(fields, row))


def tsv_rows(path, required):
    with gzip.open(path, 'rt') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        require(set(required) <= set(reader.fieldnames or []), 'A required LINCS metadata schema differs: '+Path(path).name)
        yield from reader


def finite(value):
    try:
        parsed = float(value)
        return parsed if math.isfinite(parsed) and parsed != -666 else None
    except (TypeError, ValueError):
        return None


def ingest_catalog(database, *, genes, molecules, mechanisms, lincs, fda, policy, sources, audit_path):
    """Stream each supplied source into a new catalog. Never publish here."""
    from genosapiens.track2_medication_catalog import resolve_compound
    from genosapiens.track2_lincs_selection import resolve_perturbation_gene
    require(not Path(database).exists(), 'Catalog construction needs a new staging file.')
    connection = sqlite3.connect(database)
    counts = defaultdict(lambda: [0, []])
    def exclusion(source, reason, identifier):
        value = counts[(source, reason)]
        value[0] += 1
        if len(value[1]) < 10:
            value[1].append(identifier)
    try:
        connection.executescript(DDL)
        for row in sources:
            connection.execute('INSERT INTO sources VALUES(?,?,?,?,?)',
                tuple(row[k] for k in ['source_id','release','manifest_sha256','schema','lineage_id']))
        approved, aliases = defaultdict(set), defaultdict(set)
        for _, _, row in genes:
            gene, symbol = row['id'], row['approvedSymbol']
            require(isinstance(gene, str) and gene and isinstance(symbol, str) and symbol,
                'A target row has no gene identity.')
            names = sorted({item['label'] for field in ['symbolSynonyms','obsoleteSymbols'] for item in row.get(field) or []})
            connection.execute('INSERT INTO genes VALUES(?,?,?,?)', (gene, symbol, encoded(names), 'open_targets'))
            approved[symbol.upper()].add(gene)
            for alias in names:
                aliases[alias.upper()].add(gene)
        gene_index = {'approved': approved, 'aliases': aliases}
        for alias in sorted(set(approved) | set(aliases)):
            matches = approved.get(alias) or aliases[alias]
            for gene in sorted(matches):
                connection.execute('INSERT INTO gene_aliases VALUES(?,?,?)',
                    (alias, gene, 'resolved' if len(matches)==1 else 'ambiguous'))
        gene_ids = {g for values in approved.values() for g in values}
        for row in fda['medications'].values():
            connection.execute('INSERT INTO medications VALUES(?,?)', (row['ingredient_key'], encoded(row)))
        for row in fda['products']:
            connection.execute('INSERT INTO products VALUES(?,?,?,?)',
                (row['application'], row['product_number'], row['ingredient_key'], encoded(row)))
        for row in fda['exclusions']:
            exclusion(row['source'], row['reason'], row['identifier'])
        with Path(audit_path).open('w') as audit:
            def compound(provider, row):
                resolution = resolve_compound(provider, row, fda, policy['identity'])
                audit.write(encoded(resolution)+'\n')
                connection.execute('INSERT INTO compound_identity VALUES(?,?,?,?,?,?,?)',
                    (provider, resolution['compound_id'], resolution['ingredient_key'], resolution['resolution'],
                     resolution['reason'], encoded(row), encoded(resolution)))
                if resolution['resolution'] != 'included':
                    exclusion(provider, resolution['reason'], resolution['compound_id'])
            for _, _, row in molecules:
                compound('open_targets', row)
            for row in tsv_rows(lincs['pert_info'], ['pert_id','pert_type','pert_iname','canonical_smiles','inchi_key']):
                if row['pert_type'] == 'trt_cp':
                    compound('lincs', row)
        for partition, number, row in mechanisms:
            require(isinstance(row['chemblIds'], list) and isinstance(row['targets'], list)
                and isinstance(row['references'], list), 'A mechanism has invalid source arrays.')
            mechanism_id = 'mechanism:'+sha256(canonical({'release': sources[0]['release'], 'record': row}))
            connection.execute('INSERT OR IGNORE INTO mechanisms VALUES(?,?,?,?,?,?,?,?)',
                (mechanism_id, row['actionType'], row['mechanismOfAction'], row['targetName'], row['targetType'],
                    encoded(row['targets']), encoded(row['references']), encoded(row)))
            connection.execute('INSERT INTO mechanism_origins VALUES(?,?,?)', (mechanism_id, partition, number))
            group = (row['targetType'] or '').upper() in {'PROTEIN FAMILY','PROTEIN COMPLEX','PROTEIN COMPLEX GROUP'} or len(row['targets']) > 1
            warning = 'group_target_applicability_unconfirmed' if group else None
            if not row['targets']:
                exclusion('open_targets', 'unmapped_target', mechanism_id)
            for gene in sorted(set(row['targets'])):
                if gene not in gene_ids:
                    exclusion('open_targets', 'target_missing_from_release', gene)
                    continue
                connection.execute('INSERT OR IGNORE INTO mechanism_genes VALUES(?,?,?)', (mechanism_id, gene, warning))
            for chembl in sorted(set(row['chemblIds'])):
                connection.execute('INSERT OR IGNORE INTO mechanism_compounds VALUES(?,?)', (mechanism_id, chembl))
        for row in tsv_rows(lincs['cell_info'], ['cell_id']):
            connection.execute('INSERT INTO lincs_cells VALUES(?,?)', (row['cell_id'], encoded(row)))
        for row in tsv_rows(lincs['gene_info'], ['pr_gene_id','pr_gene_symbol','pr_is_lm']):
            require(row['pr_is_lm'] in {'0','1'}, 'The measured feature flag is invalid.')
            connection.execute('INSERT INTO lincs_features VALUES(?,?,?)',
                (row['pr_gene_id'], row['pr_gene_symbol'], int(row['pr_is_lm'])))
        connection.execute('CREATE TABLE _metrics(sig_id TEXT PRIMARY KEY, data TEXT)')
        for row in tsv_rows(lincs['sig_metrics'], ['sig_id','pert_id','pert_type','pert_iname','distil_nsample','distil_cc_q75','tas']):
            connection.execute('INSERT INTO _metrics VALUES(?,?)', (row['sig_id'], encoded(row)))
        for row in tsv_rows(lincs['sig_info'], ['sig_id','pert_id','pert_type','pert_iname','cell_id','pert_dose','pert_dose_unit','pert_time','pert_time_unit']):
            found = connection.execute('SELECT data FROM _metrics WHERE sig_id=?', (row['sig_id'],)).fetchone()
            metrics = json.loads(found[0]) if found else {}
            conflict=bool(metrics) and any(metrics[k]!=row[k] for k in ['pert_id','pert_type','pert_iname'])
            mapping = resolve_perturbation_gene(row, gene_index)
            if conflict:
                # Keep both original assertions. A contradictory join cannot
                # supply genetic identity or quality-qualified drug conditions.
                mapping=dict(gene_id=None,reason='signature_metadata_identity_conflict')
                exclusion('lincs','signature_metadata_identity_conflict',row['sig_id'])
            if not conflict and row['pert_type'] in {'trt_oe','trt_sh.cgs'} and mapping['gene_id'] is None:
                exclusion('lincs', mapping['reason'], row['sig_id'])
            connection.execute('INSERT INTO lincs_signatures VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)', (
                row['sig_id'], row['pert_id'], row['pert_type'], row['pert_iname'], mapping['gene_id'],
                row['cell_id'], finite(row['pert_dose']), row['pert_dose_unit'], finite(row['pert_time']), row['pert_time_unit'],
                None if conflict else finite(metrics.get('distil_nsample')),
                None if conflict else finite(metrics.get('distil_cc_q75')),
                None if conflict else finite(metrics.get('tas')),
                mapping['reason'], encoded(dict(metadata=row, metrics=metrics))))
        connection.execute('DROP TABLE _metrics')
        coverage = [dict(source=s, reason=r, count=n, examples=examples) for (s,r),(n,examples) in sorted(counts.items())]
        for row in coverage:
            connection.execute('INSERT INTO coverage VALUES(?,?,?,?)',
                (row['source'], row['reason'], row['count'], encoded(row['examples'])))
        connection.commit()
        return coverage
    finally:
        connection.close()
