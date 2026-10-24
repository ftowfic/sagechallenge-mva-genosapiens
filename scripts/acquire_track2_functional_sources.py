#!/usr/bin/env python3
"""Freeze reviewed human UniProt variant annotations without case inputs."""
import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import re
import sqlite3
import sys
import time
import urllib.parse
import urllib.request

QUERY = 'organism_id:9606 AND reviewed:true AND ft_variant:*'
FIELDS = 'accession,id,gene_names,organism_id,xref_ensembl,ft_variant,cc_disease,cc_alternative_products,sequence,version,sequence_version'
SCHEMA = 'genosapiens.track2-uniprot-functional-catalog/v1'


def canonical(value):
    return (json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n').encode()


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def acquire(output, *, maximum_pages=100, maximum_bytes=512*1024*1024):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'src'))
    from genosapiens.track2_protein_loss import enrich_entry, POLICY as ALTERATION_POLICY
    root = Path(output).absolute()
    if root.exists() or root.resolve() != root:
        raise ValueError('Choose a new real public resource directory.')
    root.mkdir(parents=True)
    (root/'pages').mkdir()
    database = sqlite3.connect(root/'catalog.sqlite')
    database.executescript('CREATE TABLE entries(accession TEXT PRIMARY KEY, content BLOB NOT NULL);'
        'CREATE TABLE membership(gene TEXT NOT NULL, accession TEXT NOT NULL, PRIMARY KEY(gene, accession));'
        'CREATE TABLE genes(gene TEXT PRIMARY KEY, content BLOB NOT NULL);')
    url = 'https://rest.uniprot.org/uniprotkb/search?' + urllib.parse.urlencode(
        dict(query=QUERY, format='json', size=500, fields=FIELDS))
    pages, visited, release, total, count, received = [], set(), None, None, 0, 0
    started = datetime.now(timezone.utc).isoformat()
    (root/'progress.json').write_bytes(canonical(dict(status='incomplete', release=None,
        records=0, expected_records=None, pages=0, bytes_received=0)))
    try:
        while url:
            parsed = urllib.parse.urlsplit(url)
            if (parsed.scheme, parsed.netloc, parsed.path) != (
                    'https', 'rest.uniprot.org', '/uniprotkb/search') or url in visited:
                raise ValueError('Unexpected or repeated public pagination URL.')
            if len(pages) >= maximum_pages:
                raise ValueError('Public acquisition reached its page bound.')
            visited.add(url)
            for attempt in range(3):
                try:
                    request = urllib.request.Request(url, headers={'Accept':'application/json',
                        'User-Agent':'Genosapiens-public-functional-acquisition/1.0'})
                    with urllib.request.urlopen(request, timeout=45) as response:
                        raw = response.read(32*1024*1024+1)
                        headers = dict(response.headers.items())
                        current_release = response.headers.get('X-UniProt-Release')
                        current_total = int(response.headers['X-Total-Results'])
                        link = response.headers.get('Link', '')
                    break
                except (OSError, TimeoutError):
                    if attempt == 2:
                        raise
                    time.sleep(2*(attempt+1))
            received += len(raw)
            if len(raw) > 32*1024*1024 or received > maximum_bytes:
                raise ValueError('Public acquisition reached its byte bound.')
            if not current_release or current_total < 1:
                raise ValueError('The public source lacks a release or total count.')
            if release is None:
                release, total = current_release, current_total
            if (release, total) != (current_release, current_total):
                raise ValueError('The public source changed during acquisition.')
            rows = json.loads(raw)['results']
            if not rows:
                raise ValueError('Unexpected empty public source page.')
            name = f'pages/{len(pages)+1:04d}.json.gz'
            (root/name).write_bytes(gzip.compress(raw, mtime=0))
            pages.append(dict(path=name, content_sha256=digest(raw), size_bytes=len(raw),
                records=len(rows), url=url, response_headers=headers))
            for row in rows:
                if (row['entryType'] != 'UniProtKB reviewed (Swiss-Prot)'
                        or row['organism']['taxonId'] != 9606):
                    raise ValueError('The public page contains an out-of-scope entry.')
                accession = row['primaryAccession']
                if any(f['type'] == 'Natural variant' and not f.get('alternativeSequence')
                       for f in row.get('features', [])):
                    if not re.fullmatch(r'[A-Z0-9]+', accession):
                        raise ValueError('Unexpected protein source accession.')
                    request = urllib.request.Request('https://rest.uniprot.org/uniprotkb/'+accession+'.txt',
                        headers={'Accept':'text/plain', 'User-Agent':'Genosapiens-public-functional-acquisition/2.0'})
                    with urllib.request.urlopen(request, timeout=45) as response:
                        flat = response.read(2*1024*1024+1)
                        flat_release = response.headers.get('X-UniProt-Release')
                    received += len(flat)
                    if len(flat) > 2*1024*1024 or received > maximum_bytes:
                        raise ValueError('Supplementary protein evidence reached its byte bound.')
                    row = enrich_entry(row, flat.decode('utf-8'), json_release=release, flat_release=flat_release)
                database.execute('INSERT INTO entries VALUES (?,?)', (accession, canonical(row)))
                genes = set()
                for crossref in row.get('uniProtKBCrossReferences', []):
                    if crossref['database'] != 'Ensembl':
                        continue
                    for prop in crossref.get('properties', []):
                        if prop['key'] == 'GeneId' and re.fullmatch(r'ENSG[0-9]{11}(?:\.[0-9]+)?', prop['value']):
                            genes.add(prop['value'].split('.')[0])
                database.executemany('INSERT INTO membership VALUES (?,?)',
                    [(gene, accession) for gene in sorted(genes)])
            count += len(rows)
            database.commit()
            (root/'progress.json').write_bytes(canonical(dict(status='incomplete', release=release,
                records=count, expected_records=total, pages=len(pages), bytes_received=received)))
            print(f'Frozen public functional entries: {count}/{total}.', flush=True)
            match = re.search(r'<([^>]+)>;\s*rel="next"', link)
            url = match[1] if match else None
        if count != total:
            raise ValueError('Public pagination ended before the declared total.')
        genes = {}
        for (gene,) in database.execute('SELECT DISTINCT gene FROM membership ORDER BY gene').fetchall():
            entries = [json.loads(row[0]) for row in database.execute(
                'SELECT e.content FROM membership m JOIN entries e USING(accession) WHERE gene=? ORDER BY accession',
                (gene,))]
            content = canonical(dict(gene_id=gene, entries=entries))
            genes[gene] = dict(sha256=digest(content), entries=len(entries))
            database.execute('INSERT INTO genes VALUES (?,?)', (gene, gzip.compress(content, mtime=0)))
        database.commit()
        unmapped = database.execute('SELECT count(*) FROM entries WHERE accession NOT IN '
            '(SELECT accession FROM membership)').fetchone()[0]
        database.close()
        manifest = dict(schema=SCHEMA, classification='public', status='complete', source='UniProtKB',
            release=release, query=QUERY, fields=FIELDS, alteration_policy=ALTERATION_POLICY, started_at=started,
            completed_at=datetime.now(timezone.utc).isoformat(), records=count, unmapped_entries=unmapped,
            genes=genes, pages=pages, database='catalog.sqlite',
            limitations=['Coverage includes only reviewed human entries with natural variant annotations.',
                'Entries without an Ensembl gene mapping remain outside gene lookup.',
                'Source coverage does not establish an absence of contrary evidence in all literature.',
                'Protein annotations do not establish genotype effects or treatment response.'])
        raw = canonical(manifest)
        (root/'manifest.json').write_bytes(raw)
        (root/'progress.json').write_bytes(canonical(dict(status='complete', manifest_sha256=digest(raw))))
        print(json.dumps(dict(manifest_sha256=digest(raw), release=release, records=count,
            genes=len(genes), unmapped_entries=unmapped)), flush=True)
        return manifest
    finally:
        database.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    acquire(args.output)
