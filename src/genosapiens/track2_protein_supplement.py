"""Bind supplementary representations to unchanged frozen source projections."""
import copy
import json
from pathlib import Path

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_published_genotype import sha256
from genosapiens.track2_protein_loss import enrich_entry

PATH = 'config/track2/protein-representation-supplement.json'
PIN_SHA256 = '4945457cba6feebff7ec49ebffd26d01afd005a40993dd15272cc19a7f2016bb'
POLICY = 'track2-same-release-exact-protein-projection/1.0.0'


def installed_supplement():
    raw = (Path(__file__).resolve().parents[2]/PATH).read_bytes()
    require(sha256(raw) == PIN_SHA256, 'The protein representation supplement changed.')
    data = json.loads(raw)
    require(data['schema'] == 'genosapiens.track2-uniprot-representation-supplement/v1',
            'Unsupported protein representation supplement.')
    return data


def exact_projection(retained, current):
    """Every retained field must agree. Extra current fields supply provenance."""
    if isinstance(retained, dict):
        return isinstance(current, dict) and all(k in current and exact_projection(v, current[k]) for k,v in retained.items())
    if isinstance(retained, list):
        if not isinstance(current, list):
            return False
        remaining = list(current)
        for value in retained:
            matches = [i for i, candidate in enumerate(remaining) if exact_projection(value, candidate)]
            if len(matches) != 1:
                return False
            remaining.pop(matches[0])
        return True
    return type(retained) is type(current) and retained == current


def supplement_entry(entry, release, supplement):
    record = supplement['records'].get(entry['primaryAccession'])
    if not record or 'track2_terminal_loss' in entry:
        return entry
    require(sha256(record['json_text'].encode()) == record['json_sha256']
            and sha256(record['flat_text'].encode()) == record['flat_sha256'],
            'Protein supplement source bytes changed.')
    current = json.loads(record['json_text'])
    accession = entry['primaryAccession']
    require(record['json_url'] == 'https://rest.uniprot.org/uniprotkb/'+accession+'.json'
            and record['flat_url'] == 'https://rest.uniprot.org/uniprotkb/'+accession+'.txt'
            and current['primaryAccession'] == accession, 'Protein supplement source identity differs.')
    compatible = (release == record['json_release'] == record['flat_release']
                  and exact_projection(entry, current))
    proof = dict(policy=POLICY, state='compatible_source_projection' if compatible else 'unresolved',
        reason=None if compatible else 'protein_supplement_release_or_retained_projection_differs',
        frozen_release=release, supplemental_release=record['json_release'],
        frozen_entry_version=entry.get('entryAudit', {}).get('entryVersion'),
        supplemental_entry_version=current.get('entryAudit', {}).get('entryVersion'),
        json_sha256=record['json_sha256'], flat_sha256=record['flat_sha256'],
        json_url=record['json_url'], flat_url=record['flat_url'])
    result = copy.deepcopy(entry)
    if compatible:
        # Keep the original selected feature/comment domain. Audit fields alone
        # come from the authenticated additional representation.
        result['entryAudit'] = current['entryAudit']
        result = enrich_entry(result, record['flat_text'], json_release=record['json_release'],
                              flat_release=record['flat_release'])
    result['track2_representation_proof'] = proof
    return result
