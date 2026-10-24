"""Source-specific projections for shared selection from authenticated indexes.

The index owner establishes source release and reference compatibility. These
adapters validate requested records and preserve the original classification.
No source field or record establishes that compatibility by itself.
"""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
import json
import re

from genosapiens.variant_retention import (
    ClinicalFact, ExactAllele, FunctionalFact, LookupCoverage, RetentionError, SourceFact,
)


@dataclass(frozen=True, slots=True)
class RetentionEvidence:
    allele: ExactAllele | None
    clinical: tuple[ClinicalFact, ...]
    functional: tuple[FunctionalFact, ...]
    coverage: tuple[LookupCoverage, ...]
    mapping_reason: str | None = None


def _values(value):
    if value in (None, '', '.'):
        return ()
    if type(value) is not str:
        raise RetentionError('clinical source classification is invalid')
    return (value,)


def clinical_from_index(record, *, reference_identity, reference_sha256, source_release_sha256):
    """Read a complete normalized exact index row without any gene-field gate."""
    if (not isinstance(record, Mapping) or record.get('record_type') != 'exact'
            or record.get('source') not in {'clinvar', 'open_targets'}):
        raise RetentionError('clinical selection needs an exact source index record')
    try:
        key, raw, locator = record['normalized_alias'], record['evidence'], record['source_locator']
        allele = ExactAllele(reference_identity, reference_sha256, key['chromosome'],
            key['position'], key['reference'], key['alternate'])
        record_id = json.dumps(dict(locator), sort_keys=True, separators=(',', ':'))
        if record['source'] == 'clinvar':
            # Missing fields indicate an obsolete index, not a successful lookup.
            if not {'condition_ids', 'condition_names', 'gene_info'} <= raw.keys():
                raise RetentionError('clinical selection index lacks current condition fields')
            source = SourceFact('ClinVar VCF', source_release_sha256, record_id)
            classifications = list(_values(raw.get('clinical_significance')))
            conflict = _values(raw.get('conflicting_significance'))
            # CLNSIGCONF exposes positive and contrary submitted classifications
            # with counts. Preserve the aggregate and the raw conflict string.
            for value in conflict:
                for item in value.split('|'):
                    match = re.fullmatch(r'([A-Za-z_/ ]+)\(([1-9][0-9]*)\)', item)
                    if match:
                        classifications.append(match[1])
            conditions = _values(raw['condition_ids']) or _values(raw['condition_names'])
            return ClinicalFact(source, allele, 'exact', tuple(classifications), conditions,
                review_status=raw.get('review_status'), conflicts=conflict)
        source = SourceFact('Open Targets evidence_eva', source_release_sha256, record_id)
        labels = raw['clinical_significances']
        if labels is not None and (not isinstance(labels, (list, tuple))
                or any(type(value) is not str for value in labels)):
            raise RetentionError('clinical source classification is invalid')
        return ClinicalFact(source, allele, 'exact', tuple(labels or ()),
            _values(raw['disease_from_source_id']), gene_id=raw['target_id'],
            review_status=raw['confidence'])
    except RetentionError:
        raise
    except (KeyError, TypeError, ValueError):
        raise RetentionError('clinical selection source record is malformed') from None


def clinical_from_source(record, *, exact_allele):
    """Project a source row after the clinical session verifies its exact match."""
    if type(exact_allele) is not ExactAllele:
        raise RetentionError('clinical source lacks an authenticated sequence match')
    raw = record.source_fields['record']
    key = dict(chromosome=exact_allele.contig, position=exact_allele.position,
               reference=exact_allele.reference, alternate=exact_allele.alternate)
    if record.source_dataset == 'ClinVar VCF':
        info = raw['info']
        evidence = dict(clinical_significance=info.get('CLNSIG'),
            conflicting_significance=info.get('CLNSIGCONF'), review_status=info.get('CLNREVSTAT'),
            condition_ids=info.get('CLNDISDB'), condition_names=info.get('CLNDN'), gene_info=info.get('GENEINFO'))
        source = 'clinvar'
    elif record.source_dataset == 'Open Targets evidence_eva':
        evidence = dict(clinical_significances=raw['clinicalSignificances'], confidence=raw['confidence'],
            target_id=raw['targetId'], disease_from_source_id=raw['diseaseFromSourceId'])
        source = 'open_targets'
    else:
        raise RetentionError('clinical source does not contain an allele assertion')
    return clinical_from_index(dict(record_type='exact', source=source, normalized_alias=key,
        evidence=evidence, source_locator={'source_record_id':record.source_record_id,
            'source_sha256':record.source_sha256}), reference_identity=exact_allele.reference_identity,
        reference_sha256=exact_allele.reference_sha256, source_release_sha256=record.public_manifest_sha256)
