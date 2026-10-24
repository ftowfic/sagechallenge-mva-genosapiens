"""Shared record parser for supplied symbolic structural calls.

This extracts the established Track 1 format rules without importing its
execution inventory. It produces nominal source geometry, never affected bases.
"""
import math
import re

from genosapiens.allele_identity import (
    Breakpoint, CopyNumberEvidence, GenomicInterval, StructuralAllele, StructuralVariantType,
)


POLICY = 'single-symbolic-DEL-DUP-CNV-nominal-pos-end/v1'
INTERVAL_SEMANTICS = 'nominal-vcf-pos-end-context'


def parse_structural_record(record, definitions, reference_length, *, role, reference_id,
                            vcf_version):
    """Return sample/identity pairs or reject unsupported source semantics.

    The caller retains rejected source observations. Missing optional values
    remain unknown. A non-carrier still has a source call, not a carried event.
    """
    from genosapiens.vcf_qc import _parse_info
    info = _parse_info(record.raw_line.rstrip('\r\n').split('\t')[7])
    end = info.get('END')
    symbol = record.alternates[0] if len(record.alternates) == 1 else None
    svtype = info.get('SVTYPE')
    if (vcf_version not in ('VCFv4.2', 'VCFv4.3')
            or symbol not in ('<DEL>', '<DUP>', '<CNV>') or svtype != symbol[1:-1]
            or definitions.get(('INFO', 'SVTYPE')) != ('1', 'String')
            or definitions.get(('INFO', 'END')) != ('1', 'Integer')
            or type(end) is not str or not end.isdigit()
            or not record.position <= int(end) <= reference_length):
        raise ValueError('unsupported-structural-source-format')
    interval = GenomicInterval(record.contig, record.position, int(end))
    points = []
    for name, position in (('CIPOS', record.position), ('CIEND', int(end))):
        if name not in info:
            continue
        if definitions.get(('INFO', name)) != ('2', 'Integer'):
            raise ValueError('unsupported-structural-source-format')
        if info[name] == '.':
            continue
        if type(info[name]) is not str or re.fullmatch(r'-?[0-9]+,-?[0-9]+', info[name]) is None:
            raise ValueError('unsupported-structural-source-format')
        bounds = tuple(int(v) for v in info[name].split(','))
        point = Breakpoint(record.contig, position, confidence_interval=bounds)
        if position + bounds[1] > reference_length:
            raise ValueError('unsupported-structural-source-format')
        points.append(point)
    if 'IMPRECISE' in info and (definitions.get(('INFO', 'IMPRECISE')) != ('0', 'Flag')
                              or info['IMPRECISE'] is not None):
        raise ValueError('unsupported-structural-source-format')
    calls = []
    for sample in record.sample_calls:
        fields = dict(sample.class_format_fields)
        value = fields.get('CN')
        count = None
        if 'CN' in fields:
            contract = definitions.get(('FORMAT', 'CN'))
            if contract not in {('1', 'Integer'), ('1', 'Float')}:
                raise ValueError('unsupported-structural-source-format')
            if value not in (None, '.'):
                if contract[1] == 'Integer' and not value.isdigit():
                    raise ValueError('unsupported-structural-source-format')
                count = int(value) if contract[1] == 'Integer' else float(value)
                if not math.isfinite(count) or count < 0:
                    raise ValueError('unsupported-structural-source-format')
        identity = StructuralAllele(role, reference_id, StructuralVariantType(svtype), interval,
            tuple(points), CopyNumberEvidence(count, None, sample.genotype, sample.ploidy))
        calls.append((sample, identity))
    return tuple(calls)


def structural_definitions(lines):
    """Read declared types through the shared VCF metadata parser."""
    from genosapiens.vcf_qc import _structured_attributes
    definitions = {}
    for raw in lines:
        line = raw.strip()
        if not line.startswith(('##INFO=<', '##FORMAT=<')) or not line.endswith('>'):
            continue
        attributes = _structured_attributes(line.split('=<',1)[1][:-1])
        if attributes is None or 'ID' not in attributes:
            continue
        key = (line[2:].split('=',1)[0],attributes['ID'])
        definitions[key] = ((None, None) if key in definitions else
            (attributes.get('Number'),attributes.get('Type')))
    return definitions


def structural_envelope(identity):
    """Use both supplied endpoint bounds, or preserve unbounded uncertainty."""
    if identity.interval is None or len(identity.breakpoints) != 2:
        return None
    left, right = identity.breakpoints
    if left.confidence_interval is None or right.confidence_interval is None:
        return None
    return GenomicInterval(identity.interval.contig,
        identity.interval.start + left.confidence_interval[0],
        identity.interval.end + right.confidence_interval[1])
