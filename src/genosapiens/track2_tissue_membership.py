"""Release-bound primary-site membership from retained catalog annotations."""
import json

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_tissue_records import NORMALIZATION_POLICY, schema_name, validate_tissue_record

UNRESOLVED = frozenset(('', '-666', 'unknown', 'na', 'n/a', 'none', 'null', 'not available'))


def membership_record(cell_id, annotations, *, source_release, metadata_sha256):
    raw = annotations.get('primary_site')
    require(raw is None or isinstance(raw, str), 'LINCS primary_site must be source text or absent.')
    normalized = ' '.join((raw or '').split()).casefold()
    unresolved = normalized in UNRESOLVED
    record = dict(schema=schema_name('membership'), source_release=source_release,
        metadata_sha256=metadata_sha256, normalization_policy=NORMALIZATION_POLICY, cell_id=cell_id,
        raw_primary_site=raw, tissue_id=None if unresolved else normalized,
        unresolved_reason='missing_primary_site' if raw is None else 'unresolved_primary_site' if unresolved else None,
        annotations=dict(annotations))
    return validate_tissue_record(record, 'membership')


def catalog_membership(connection, *, source_release, metadata_sha256):
    # Include cells referenced by signatures but absent from the cell metadata.
    rows = connection.execute('''SELECT c.cell_id, c.annotations_json FROM lincs_cells c
        UNION ALL SELECT DISTINCT s.cell_id, '{}' FROM lincs_signatures s
        WHERE NOT EXISTS (SELECT 1 FROM lincs_cells c WHERE c.cell_id=s.cell_id)
        ORDER BY cell_id''')
    return [membership_record(cell, json.loads(annotations), source_release=source_release,
        metadata_sha256=metadata_sha256) for cell, annotations in rows]
