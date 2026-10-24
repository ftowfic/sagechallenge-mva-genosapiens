"""Revise report gene metadata while preserving scientific analysis products."""
import argparse
import copy
import json
from pathlib import Path
import shutil

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_report_gene_metadata import FILENAME, SCHEMA
from genosapiens.track2_scientific_storage import file_identity, publish_json


def prepare_revision(source, metadata_package, destination):
    """Copy a retained report and change only its gene display context."""
    source, metadata_package, destination = (Path(p).absolute() for p in (source, metadata_package, destination))
    require(source.resolve() == source and destination.resolve() == destination
            and not destination.exists() and not destination.is_relative_to(source),
            'Choose a new report revision directory outside the source report.')
    model = json.loads((source/'presentation.json').read_bytes())
    before = copy.deepcopy(model)
    # This presentation-only operation uses the existing compact projection.
    # It does not open the source Parquet partitions, discovery catalog, or matrix.
    metadata_path = metadata_package/FILENAME
    data = json.loads(metadata_path.read_bytes())
    require(data['schema'] == SCHEMA and data['classification'] == 'public',
            'The revision requires public gene metadata.')
    mappings = {}
    for hypothesis in model['context']['hypotheses']:
        for gene in hypothesis['genes']:
            gid = gene['gene_id']
            rows = data['records'].get(gid, [])
            require(all(row['gene_id'] == gid and row['source_record_id'] == gid for row in rows),
                    'A revision gene record has a different stable identifier.')
            if len(rows) == 1:
                record = dict(rows[0], mapping_state='exact_stable_identifier')
            elif rows:
                record = dict(gene_id=gid, mapping_state='ambiguous', alternatives=rows,
                    source_manifest_sha256=data['source_manifest_sha256'])
            else:
                record = None
            mappings[gid] = record
            gene['reference_metadata'] = record
            gene['metadata_state'] = record['mapping_state'] if record else 'unavailable'
            if record and record.get('symbol'):
                gene['symbol'] = record['symbol']
    def ignore(path, names):
        if Path(path) == source:
            return set(names) & {'manifest.json','artifact-index.jsonl','artifact-index-parts'}
        return []
    shutil.copytree(source, destination, ignore=ignore)
    destination.chmod(0o700)
    for path in destination.rglob('*'):
        if path.is_dir():
            path.chmod(0o700)
    for path in [destination/'presentation.json', destination/'supplements/context.json']:
        path.unlink()
    publish_json(destination/'presentation.json', model)
    publish_json(destination/'supplements/context.json', model['context'])
    publish_json(destination/'supplements/gene-metadata.json', dict(schema=SCHEMA, classification='public',
        source_manifest_sha256=data['source_manifest_sha256'], records={k:[v] if v else [] for k,v in mappings.items()},
        source_companion=file_identity(metadata_path)))
    proof = dict(schema='genosapiens.track2-report-metadata-revision/v1', state='prepared',
        source_report=str(source), source_manifest=file_identity(source/'manifest.json'),
        source_presentation=file_identity(source/'presentation.json'), metadata_companion=file_identity(metadata_path),
        scientific_analysis='reused_without_reanalysis', human_output_acceptance=False,
        unchanged_model_fields=[key for key in before if key not in {'context','figures','supplements'}])
    require(all(before[key] == model[key] for key in proof['unchanged_model_fields']),
            'A gene metadata revision changed scientific report data.')
    publish_json(destination/'revision.json', proof)
    return destination


def finish_revision(destination):
    """Refresh overview provenance and render the updated presentation."""
    from genosapiens.track2_report_model import check_model
    from genosapiens.track2_report_overviews import render_overviews
    from genosapiens.track2_report_paper import paper_markdown
    root = Path(destination)
    proof = json.loads((root/'revision.json').read_bytes())
    require(proof['state'] in {'prepared','complete'}, 'The metadata revision is not prepared.')
    model = json.loads((root/'presentation.json').read_bytes())
    for path in (root/'overviews').glob('*'):
        if path.is_file():
            path.chmod(0o600)
    overview = render_overviews(root/'overviews', model)
    for figure in model['figures']:
        if figure['main']:
            figure['provenance'] = overview['binding']
    for title, path in [('Gene metadata and source binding','supplements/gene-metadata.json'),
                        ('Report revision provenance','revision.json')]:
        supplement = dict(title=title, path=path)
        if supplement not in model['supplements']:
            model['supplements'].append(supplement)
    check_model(model)
    (root/'presentation.json').unlink()
    publish_json(root/'presentation.json', model)
    title = (root/'report.md').read_text().splitlines()[0]
    text = ''.join(paper_markdown(model, title, json.loads((root/'readiness.json').read_bytes())))
    (root/'report.md').write_text(text)
    proof['state'] = 'complete'
    from genosapiens import track2_report_paper
    proof['renderer'] = file_identity(track2_report_paper.__file__)
    proof['updated_files'] = {str(p.relative_to(root)):file_identity(p) for p in
        [root/'report.md', root/'presentation.json', root/'supplements/context.json',
         root/'supplements/gene-metadata.json', root/'overviews/manifest.json']}
    publish_json(root/'revision.json', proof)
    publish_json(root/'manifest.json', dict(schema=proof['schema'], state='complete',
        source_manifest=proof['source_manifest'], revision=file_identity(root/'revision.json'),
        scientific_analysis=proof['scientific_analysis']))
    return root/'report.md'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path)
    parser.add_argument('--metadata-package', type=Path)
    parser.add_argument('--destination', type=Path, required=True)
    parser.add_argument('--finish', action='store_true')
    args = parser.parse_args()
    if args.finish:
        print(finish_revision(args.destination))
    else:
        if args.source is None or args.metadata_package is None:
            parser.error('Supply --source and --metadata-package.')
        print(prepare_revision(args.source, args.metadata_package, args.destination))


if __name__ == '__main__':
    main()
