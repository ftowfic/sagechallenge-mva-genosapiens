"""Explicit UniProt terminal-loss evidence, separate from nucleotide identity."""
import copy
import re

from genosapiens.track2_published_genotype import canonical, sha256

SCHEMA = 'genosapiens.track2-protein-alteration/v1'
POLICY = 'track2-explicit-terminal-protein-loss/1.0.0'


def _flat_features(text):
    blocks = []
    current = None
    for line in text.splitlines():
        header = re.match(r'^FT   (\w+)\s+(\S+)\s*$', line)
        if header:
            current = [header[1], header[2], []]
            blocks.append(current)
        elif line.startswith('FT ') and current:
            current[2].append(line[5:].strip())
    result = {}
    for kind, position, lines in blocks:
        if kind != 'VARIANT':
            continue
        content = ' '.join(lines)
        identifier = re.search(r'/id="([^"]+)"', content)
        note = re.search(r'/note="([^"]+)"', content)
        if identifier and note:
            result.setdefault(identifier[1], []).append(dict(position=position, note=note[1], raw=content))
    return result


def normalize_terminal_loss(entry, text, *, json_release, flat_release):
    """Retain ambiguous records as unresolved. Never infer deletion from emptiness."""
    audit = entry.get('entryAudit', {})
    entry_version = re.search(r'^DT .* entry version (\d+)\.', text, re.M)
    sequence_version = re.search(r'^DT .* sequence version (\d+)\.', text, re.M)
    accession = re.search(r'^AC\s+(\w+);', text, re.M)
    sequence_text = text.split('\nSQ ', 1)[-1].split('\n', 1)[-1].split('//', 1)[0] if '\nSQ ' in text else ''
    sequence = re.sub(r'\s|\d', '', sequence_text)
    versions = bool(json_release and json_release == flat_release and entry_version and sequence_version
        and audit.get('entryVersion') == int(entry_version[1])
        and audit.get('sequenceVersion') == int(sequence_version[1])
        and accession and accession[1] == entry['primaryAccession'])
    features = _flat_features(text)
    records = {}
    for feature in entry.get('features', []):
        if feature['type'] != 'Natural variant' or feature.get('alternativeSequence'):
            continue
        key = feature.get('featureId')
        candidates = features.get(key, [])
        location = feature['location']
        left, right = location['start'].get('value'), location['end'].get('value')
        reason = None
        if not versions:
            reason = 'protein_source_versions_disagree_or_unavailable'
        elif not sequence or sequence != entry.get('sequence', {}).get('value'):
            reason = 'protein_source_sequences_disagree_or_unavailable'
        elif len(candidates) != 1 or not re.match(r'^Missing(?: \(|$)', candidates[0]['note']):
            reason = 'explicit_protein_loss_meaning_unavailable'
        elif (not isinstance(left, int) or not isinstance(right, int) or not 1 <= left <= right <= len(sequence)
              or location['start'].get('modifier') != 'EXACT' or location['end'].get('modifier') != 'EXACT'
              or candidates[0]['position'] != (str(left) if left == right else f'{left}..{right}')):
            reason = 'protein_source_boundaries_disagree_or_unavailable'
        elif right != len(sequence):
            reason = 'internal_protein_deletion_requires_separate_evidence'
        records[key] = dict(schema=SCHEMA, policy=POLICY, state='unresolved' if reason else 'explicit_terminal_loss',
            reason=reason, feature_id=key, start=left, end=right, accession=entry['primaryAccession'],
            entry_version=audit.get('entryVersion'), sequence_version=audit.get('sequenceVersion'),
            sequence_sha256=sha256(sequence.encode()), missing_sequence=sequence[left-1:right] if reason is None else None,
            source_release=json_release, flat_release=flat_release, flat_sha256=sha256(text.encode()),
            source_feature=candidates[0] if len(candidates) == 1 else None)
    return records


def enrich_entry(entry, text, *, json_release, flat_release):
    result = copy.deepcopy(entry)
    result['track2_terminal_loss'] = dict(schema=SCHEMA, policy=POLICY, json_release=json_release,
        flat_release=flat_release, flat_text=text, flat_sha256=sha256(text.encode()),
        records=normalize_terminal_loss(entry, text, json_release=json_release, flat_release=flat_release))
    return result


def terminal_loss_record(entry, feature):
    retained = entry.get('track2_terminal_loss')
    if not retained:
        return None
    from genosapiens.track2_operator_inputs import require
    require(retained['schema'] == SCHEMA and retained['policy'] == POLICY
        and sha256(retained['flat_text'].encode()) == retained['flat_sha256'], 'Protein alteration source changed.')
    replay = normalize_terminal_loss(entry, retained['flat_text'], json_release=retained['json_release'],
                                     flat_release=retained['flat_release'])
    require(replay == retained['records'], 'Protein alteration interpretation changed.')
    return replay.get(feature.get('featureId'))


def compatible_stop(record, *, start, end, reference, alternate, terms, sequence):
    """Call only after exact transcript, protein, isoform, and residue checks."""
    if not record or record['state'] != 'explicit_terminal_loss':
        return None
    if (set(terms or ()) != {'stop_gained'} or alternate != '*' or start != end
            or start != record['start'] or record['end'] != len(sequence)
            or len(reference) != 1 or sequence[start-1:start] != reference
            or record['sequence_sha256'] != sha256(sequence.encode())
            or record['missing_sequence'] != sequence[start-1:]):
        return None
    return dict(kind='protein_effect_compatibility', alteration=record,
        reference_residue=reference, stop_position=start, nucleotide_identity='not_established',
        patient_assay_transfer='not_established')
