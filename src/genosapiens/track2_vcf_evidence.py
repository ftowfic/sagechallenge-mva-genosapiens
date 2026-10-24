"""Bound interpretation limits for the current text VCF adapter."""
from genosapiens.vcf_contracts import ReferenceDeclaration, VcfReadLimits


NUMERIC = {'adapter':'genosapiens.vcf_qc/text-vcf-v1', 'encoding':'ascii',
    'source_identity':'original-bytes', 'qual_conversion':'decimal-to-python-binary64',
    'raw_tokens':'QUAL-and-INFO-preserved', 'library_rendering':'not-used'}


def retrieval(strategy):
    if strategy not in {'streaming','indexed'}:
        raise ValueError('VCF retrieval strategy is invalid.')
    return {'method':'streaming' if strategy=='streaming' else 'custom-selection-index',
        'index_build':'not-used' if strategy=='streaming' else 'complete-source-pass',
        'source_authentication':'complete-source-digest',
        'selection_domain':'source-event-and-normalized-spans',
        'exclusion_accounting':'complete-source-occurrences', 'native_tbi':'deferred'}


def evidence(bound, strategy):
    return {'schema':'genosapiens.track2-vcf-evidence/v1', 'numeric':dict(NUMERIC),
        'reference':{key:bound[key] for key in ('reference_sha256','reference_identity_id',
            'reference_role','concordance_sha256','reference_mapping','reference_declarations')},
        'read_limits':bound['read_limits'], 'pid_contract_sha256':bound['pid_contract_sha256'],
        'retrieval':retrieval(strategy)}


def validate(value, selection, structural):
    if (set(value)!={'schema','numeric','reference','read_limits','pid_contract_sha256','retrieval'}
            or value['schema']!='genosapiens.track2-vcf-evidence/v1'
            or value['numeric']!=NUMERIC or value['retrieval']!=retrieval(selection['strategy'])):
        raise ValueError('VCF evidence contract is invalid.')
    VcfReadLimits(**value['read_limits'])
    reference=value['reference']
    if (set(reference)!={'reference_sha256','reference_identity_id','reference_role',
            'concordance_sha256','reference_mapping','reference_declarations'}
            or reference['reference_sha256']!=selection['reference_sha256']
            or reference['reference_mapping']!=structural['reference_mapping']):
        raise ValueError('VCF reference evidence differs.')
    mapping=reference['reference_mapping']
    if (reference['reference_identity_id']!=mapping['source_identity_id']
            or reference['reference_role']!=mapping['source_role']):
        raise ValueError('VCF operational reference identity differs.')
    for declaration in reference['reference_declarations']:
        ReferenceDeclaration(**declaration)
    from genosapiens.track2_diagnosis_package import _digest
    if not _digest(reference['concordance_sha256']) or (value['pid_contract_sha256'] is not None
            and not _digest(value['pid_contract_sha256'])):
        raise ValueError('VCF evidence digest is invalid.')


def source_numbers(record):
    """Keep source text independently from Python's numeric conversion."""
    return {'raw_qual':record.raw_quality, 'typed_qual':record.quality,
        'raw_info':record.raw_line.split('\t',8)[7], **NUMERIC}


def validate_calls(parts):
    """Reject stale source projections even when their package was rehashed."""
    import json
    from genosapiens.evidence_model import AlleleEvidence, _typed
    from genosapiens.track2_vcf_selection import replay_sequence_quality
    from genosapiens.vcf_qc import _parse_quality
    for row in parts['evidence']:
        if not row['evidence_id'].startswith('scoped-allele:'):continue
        allele=_typed(json.loads(row['assertion']),AlleleEvidence,decode=True)
        for call in allele.calls or ():
            fields=call.source.source_fields
            n=fields.get('numeric_evidence',{})
            if (set(n)!={'raw_qual','typed_qual','raw_info',*NUMERIC}
                    or any(n[key]!=value for key,value in NUMERIC.items())
                    or n['typed_qual']!=_parse_quality(n['raw_qual'])
                    or n['typed_qual']!=call.quality or not isinstance(n['raw_info'],str)):
                raise ValueError('VCF numeric source evidence differs.')
            replay_sequence_quality(call)


def render(parts):
    """Render all selected calls, including calls without a genotype scenario."""
    def cell(value):
        return str(value if value is not None else 'unresolved').replace('|','\\|').replace('\n',' ')
    item=parts['selection']['vcf_evidence'];ref=item['reference'];method=item['retrieval']
    lines=['## VCF evidence contracts','',
        'The text adapter preserves source QUAL and INFO tokens. Python converts typed QUAL values to binary64.',
        'Binary64 conversion can lose decimal precision. This adapter does not use library rendering.','',
        'The reader requires ASCII input. This requirement does not establish conformance for other files.',
        '',
        'Retrieval method: `'+method['method']+'`. Index construction: `'+method['index_build']+'`.',
        'Selection checks source-event spans and normalized spans. Exclusion counts cover the complete source.',
        'Native TBI retrieval remains deferred. A region probe cannot establish complete selection.','',
        'The installed reference and accepted concordance establish operational reference identity.',
        'Unresolved header declarations remain metadata. They do not supply contig aliases.',
        'Reference identity: `'+cell(ref['reference_identity_id'])+'`. Role: `'+cell(ref['reference_role'])+'`.',
        'Concordance SHA-256: `'+ref['concordance_sha256']+'`.','',
        'PID contract SHA-256: `'+cell(item['pid_contract_sha256'])+'`.','',
        '| Reference declaration | Resolution |','|---|---|',
        *('| '+cell(r['raw_value'])+' | '+cell(r['resolution'])+' |' for r in ref['reference_declarations']),
        '', 'AD and DP remain separate measurements. The balance policy divides ALT depth by FORMAT/DP.',
        'An unexplained AD total above DP prevents automatic quality acceptance.','',
        'PGT/PID grouping requires an applicable producer contract. Standard GT/PS evidence remains independent.',
        'The tables preserve exact source tokens, states, and reason identifiers.','',
        '| Source allele | GT | AD | DP | Balance | Quality | AD/DP state | Reasons |',
        '|---|---|---|---|---|---|---|---|']
    measurements=[];phases=[];formats=[];numbers=[]
    for row in parts['evidence']:
        if not row['evidence_id'].startswith('scoped-allele:'):continue
        import json
        allele=json.loads(row['assertion'])
        for call in allele['calls'] or ():
            fields=call['source']['source_fields'];q=fields['automatic_quality'];key=fields['source_allele_id']
            measurements.append((key,q['genotype'],q['allele_depths'],q['depth'],q['allele_balance'],
                q['state'],q['ad_dp_state'],q['reasons']))
            phases.append((key,q['phase_state'],q['phase_set'],q['pgt'],
                q['pid'],q['pid_contract_sha256'],q['phase_reasons']))
            formats.extend((key,f['key'],f['state'],f['raw_value']) for f in q['format_fields'])
            n=fields['numeric_evidence'];numbers.append((key,n['raw_qual'],n['typed_qual'],n['raw_info']))
    lines+=['| '+' | '.join(cell(v) for v in row)+' |' for row in sorted(measurements,key=lambda r:r[0])]
    for headings,rows in ((('Source allele','Phase state','PS','PGT','PID','Contract SHA-256','Reasons'),phases),
            (('Source allele','FORMAT key','State','Raw value'),formats),
            (('Source allele','Raw QUAL','Typed QUAL','Raw INFO'),numbers)):
        lines+=['','| '+' | '.join(headings)+' |','|'+'---|'*len(headings),
            *('| '+' | '.join(cell(v) for v in row)+' |' for row in sorted(rows,key=lambda r:r[:2]))]
    lines+=['','### VCF terms','',
        '| Term | Meaning |','|---|---|',
        '| AD | Read depths for the reference allele and each alternate allele. |',
        '| DP | The sample read depth reported in FORMAT. |',
        '| GT | The reported genotype. |',
        '| PS | A standard phase-set identifier. |',
        '| PGT/PID | Producer-specific phase evidence and its block identifier. |',
        '| QUAL | The source record quality value. |',
        '| Binary64 | The floating-point number format used for typed QUAL. |',
        '| TBI | A native genomic region index. |',
        '| SHA-256 | A digest that binds exact bytes. |']
    return lines
