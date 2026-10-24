"""Scoped source observations, separate from sequence genotype interpretation."""
from collections import Counter
from dataclasses import replace
import re

from genosapiens.allele_identity import ReferenceRole
from genosapiens.evidence_serialization import _plain
from genosapiens.supplied_structural import POLICY, parse_structural_record, structural_envelope, structural_definitions
from genosapiens.track2_published_genotype import canonical, sha256
from genosapiens.variant_normalization import source_allele_id
from genosapiens.vcf_qc import _parse_info


SCHEMA = 'genosapiens.track2-scoped-structural-evidence/v1'
GEOMETRY = 'nominal-or-supplied-endpoint-envelope;unbounded-is-unresolved/v1'
OUTCOMES = ('scoped', 'outside_scope', 'non_carrier', 'unsupported', 'unresolved')
LIMITATIONS = [
    'Structural observations come only from supplied calls. The builder does not call new variants.',
    'Nominal intervals do not establish exact affected bases. Missing endpoint bounds remain unknown.',
    'Structural genotype construction, mechanism qualification, and medication matching remain unassessed.',
    'Overlap does not establish causality, a complete genotype, therapeutic direction, or readiness.',
]


def source_context(header):
    return (structural_definitions(header.raw_lines),
        next((line.strip().split('=',1)[1] for line in header.raw_lines
            if line.startswith('##fileformat=')), None))


def observe(record, context, regions_by_contig, bound):
    """Keep one PASS observation per original ALT and selected sample."""
    if record.filter_value != 'PASS':
        return []
    if all(re.fullmatch(r'[ACGTN]+', alt) or alt == '.' for alt in record.alternates):
        return []
    sample = next(s for s in record.sample_calls if s.sample_name == bound['sample'])
    definitions, version = context
    mapping = bound['reference_mapping']
    chromosome = mapping['contig_map'][record.contig]
    regions = regions_by_contig.get(chromosome,())
    result = []
    for number, alternate in enumerate(record.alternates, 1):
        if re.fullmatch(r'[ACGTN]+', alternate) or alternate == '.':
            continue
        tokens = re.split(r'[/|]', sample.genotype or '.')
        carriage = ('carried' if str(number) in tokens else
            'unknown' if '.' in tokens else 'non_carrier')
        identity = None
        reasons = []
        try:
            # Another sample's malformed optional fields cannot erase this call.
            calls = parse_structural_record(replace(record,sample_calls=(sample,)), definitions,
                mapping['source_lengths'][record.contig], role=ReferenceRole(bound['reference_role']),
                reference_id=bound['reference_identity_id'], vcf_version=version)
            identity = calls[0][1]
        except (ValueError, TypeError, KeyError, OverflowError):
            reasons.append('unsupported-structural-source-format')
        matches = []
        geometry_state = 'unresolved'
        envelope = None
        if identity is None:
            outcome = 'unsupported'
        else:
            interval = identity.interval
            envelope = structural_envelope(identity)
            geometry_state = 'bounded_possible_interval' if envelope else 'nominal_only_unbounded_uncertainty'
            for region in regions:
                nominal = interval.start - 1 < region['end'] and region['start'] < interval.end
                possible = envelope is not None and envelope.start - 1 < region['end'] and region['start'] < envelope.end
                if nominal or possible:
                    matches.append({'interval':region,
                        'overlap':'nominal_overlap' if nominal else 'possible_overlap'})
            matches = sorted({canonical(m):m for m in matches}.values(),key=canonical)
            outcome = ('scoped' if matches else 'outside_scope' if not regions or envelope is not None else 'unresolved')
            if outcome == 'unresolved':
                reasons.append('endpoint-envelope-unbounded')
        if carriage == 'non_carrier':
            outcome = 'non_carrier'
        result.append({'observation_id':source_allele_id(record.source_record_id,number),
            'source_record_id':record.source_record_id,'record_number':record.record_number,
            'allele_number':number,'alternate':alternate,'sample':sample.sample_name,
            'raw_line_sha256':record.raw_line_sha256,'source_line':record.raw_line,
            'filter':record.filter_value,'call':_plain(sample),
            'info':_parse_info(record.raw_line.split('\t')[7]),
            'carriage':carriage,'identity':_plain(identity),'geometry_state':geometry_state,
            'possible_interval':_plain(envelope),'matches':matches,'outcome':outcome,'reasons':reasons})
    return result


def component(observations, header, bound):
    observations = sorted(observations,key=lambda r:(r['record_number'],r['allele_number']))
    counts = Counter(r['outcome'] for r in observations)
    from dataclasses import asdict
    from genosapiens.variant_retention import (
        SourceFact, StructuralFact, LookupCoverage, RetentionFacts, select_variant,
    )
    retained=[]
    for row in observations:
        if row['outcome'] != 'scoped' or row['identity'] is None:
            continue
        source=SourceFact('supplied-structural',bound['source_sha256'],row['source_record_id'])
        facts=RetentionFacts(row['observation_id'],None,'assess',(row['observation_id'],),(),(),(),
            (StructuralFact(source,row['observation_id'],True,('structural_interpretation_not_assessed',)),),
            (LookupCoverage(source.dataset,source.release_sha256,True,'complete'),))
        retained.append(asdict(select_variant(facts,established_genes=frozenset(),catalog_sha256=bound['scope_sha256'])))
    return {'schema':SCHEMA,'source_sha256':bound['source_sha256'],'sample':bound['sample'],
        'scope_sha256':bound['scope_sha256'],'reference_mapping':bound['reference_mapping'],
        'parser_policy':POLICY,'geometry_policy':GEOMETRY,'header':_plain(header),
        'retention':_plain(retained),
        'observations':observations,'coverage':{'state':'complete_scoped_screening',
            'counting_unit':'source_ALT_and_selected_sample','supplied':len(observations),
            **{name:counts[name] for name in OUTCOMES},
            'unresolved_geometry':sum(r['geometry_state'] in ('unresolved','nominal_only_unbounded_uncertainty') for r in observations),
            'unknown_genotype':sum(r['carriage']=='unknown' for r in observations)},
        'interpretation':{'genotype':'not_assessed','mechanism':'not_assessed','medication':'not_assessed'},
        'limitations':LIMITATIONS}


def validate_component(value, selection, scope):
    """Replay retained observations and coverage without access to raw case files."""
    from genosapiens.track2_codec import decode
    from genosapiens.track2_vcf_selection import require, prepare_intervals
    from genosapiens.vcf_context import prepare_header, FormatLayoutCache
    from genosapiens.vcf_qc import VcfQcStreamHeader, _parse_stream_record
    require(value['schema']==SCHEMA and value['source_sha256']==selection['source_sha256']
        and value['sample']==selection['sample'] and value['scope_sha256']==sha256(canonical(scope))
        and value['reference_mapping']['reference_sha256']==selection['reference_sha256'],
        'structural evidence bindings differ')
    header=decode(VcfQcStreamHeader,value['header'])
    require(header.input_sha256==selection['source_sha256'])
    bound={'sample':value['sample'],'source_sha256':value['source_sha256'],
        'scope_sha256':value['scope_sha256'],'reference_mapping':value['reference_mapping'],
        'reference_role':value['reference_mapping']['source_role'],
        'reference_identity_id':value['reference_mapping']['source_identity_id']}
    records={}
    for row in value['observations']:
        key=row['record_number']
        require(key not in records or records[key]==row['source_line'])
        records[key]=row['source_line']
    replay=[]
    context=source_context(header)
    regions=prepare_intervals(scope).structural_by_contig
    prepared=prepare_header(header);format_cache=FormatLayoutCache()
    for key,line in sorted(records.items()):
        record=_parse_stream_record(line.rstrip('\r\n'),record_number=key,header=header,
            raw_line_sha256=sha256(line.encode('ascii')),raw_line=line,
            prepared_header=prepared,format_cache=format_cache)
        replay.extend(observe(record,context,regions,bound))
    require(value==component(replay,header,bound),'structural evidence or coverage changed')
    require(not ({r['observation_id'] for r in replay}&{r['allele_id'] for r in selection['selected_alleles']}),
        'structural identity entered sequence selection')


def render_structural(value):
    """Both reports use one deterministic view of the bound observations."""
    def cell(text):
        return str(text).replace('&','&amp;').replace('<','&lt;').replace('>','&gt;').replace('|','&#124;').replace('\n',' ')
    coverage=value['coverage']
    labels={'scoped':'Scope overlap','outside_scope':'Outside scope','non_carrier':'Explicit non-carrier',
        'unsupported':'Unsupported format','unresolved':'Unresolved scope'}
    lines=['<!-- structural-evidence:start -->','## Supplied structural observations','',
        'This section screens supplied structural calls against the frozen diagnosis scope.',
        'It counts each original alternate allele for the selected sample once. Repeated interval hits do not increase counts.',
        'The tables include outside-scope observations and unsupported formats for coverage accounting.',
        'An empty sequence result does not establish the absence of structural observations.','',
        'Component SHA-256: `'+sha256(canonical(value))+'`.','',
        '| Coverage outcome | Observations |','|---|---|',
        f"| Supplied structural or unsupported alternate alleles | {coverage['supplied']} |"]
    lines += [f'| {labels[key]} | {coverage[key]} |' for key in OUTCOMES]
    lines += [f"| Geometry lacks complete endpoint bounds | {coverage['unresolved_geometry']} |",
        f"| Genotype is unknown or partly missing | {coverage['unknown_genotype']} |",'',
        'The last two counts describe uncertainty. They can overlap other outcomes.','',
        '| Source observation | Alternate | Nominal source interval | Source genotype | Copy number | Geometry | Outcome | Scope hits |',
        '|---|---|---|---|---|---|---|---|']
    for row in value['observations']:
        identity=row['identity']
        source=row['source_line'].split('\t')
        interval=(f"{source[0]}:{identity['interval']['start']}-{identity['interval']['end']}"
            if identity else f'{source[0]}:{source[1]} (interval not assessed)')
        dosage=identity.get('dosage') if identity else None
        copy_number=dosage['copy_number'] if dosage else None
        geometry=('Supplied endpoint bounds' if row['geometry_state']=='bounded_possible_interval'
            else 'Endpoint bounds unknown' if identity else 'Unresolved')
        hits=', '.join(sorted({m['overlap'].replace('_',' ') for m in row['matches']})) or 'None established'
        lines.append('| '+' | '.join(cell(v) for v in (row['observation_id'],row['alternate'],interval,
            row['call']['genotype'] or 'Unknown',copy_number if copy_number is not None else 'Unknown',
            geometry,labels[row['outcome']],hits))+' |')
    if not value['observations']:
        lines += ['','The supplied source contains no structural or unsupported alternate alleles.']
    lines += ['',*value['limitations'],'',
        'The retained component contains source records, interval matches, endpoint offsets, call fields, and unsupported-format reasons.','',
        '### Structural terms','',
        '| Term | Meaning |','|---|---|',
        '| Structural observation | A supplied non-sequence alternate allele and its selected sample call. |',
        '| DEL, DUP, CNV | Source symbols for deletion, duplication, and copy-number variation. |',
        '| Nominal interval | The source POS through END coordinates, with both endpoints included. Exact affected bases remain uncertain. |',
        '| Endpoint bounds | Supplied CIPOS and CIEND offsets around the two reported endpoints. Missing offsets remain unknown. |',
        '| Possible overlap | Only the interval allowed by supplied endpoint bounds overlaps the diagnosis scope. |',
        '| Copy number | The supplied CN value. It does not establish a normal baseline or a disease effect. |',
        '| Source genotype | The supplied GT allele indexes. A dot means that the source omits that allele call. |',
        '| Explicit non-carrier | The complete supplied genotype excludes this alternate allele. |',
        '| Scope | The fixed gene and regulatory intervals used for screening. |',
        '| SHA-256 | A hash that identifies exact component contents. |','',
        'Source identifiers, field names, symbols, and genotype values retain their exact technical forms.',
        '<!-- structural-evidence:end -->','']
    return '\n'.join(lines)


def validate_report_section(text, value):
    if (text.count('<!-- structural-evidence:start -->') != 1
            or text.count('<!-- structural-evidence:end -->') != 1
            or render_structural(value) not in text):
        raise ValueError('report structural section disagrees with bound evidence')
