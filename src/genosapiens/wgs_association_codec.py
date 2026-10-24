"""One bounded compact codec for discovery rows; ordinary evidence is unchanged.

Tables preserve first occurrence and all proof multiplicity. Live owners retain
the independently validated typed row and its single canonical encoding; neither
a supplied digest nor a serialized trust flag can manufacture that ownership.
"""
from dataclasses import dataclass, fields
import base64
import hashlib
import zlib

from genosapiens.allele_identity import RepeatAllele
from genosapiens.evidence_model import IdentityMapping, SourcedRecord, _typed_value, _validation_session
from genosapiens.evidence_serialization import _loads, _owned_canonical_bytes
from genosapiens.wgs_shards import WgsScientificError

MAX_LOGICAL_BYTES = 32 * 1024**2
MAX_PART_BYTES = 4 * 1024**2
MAX_FRAGMENT_BYTES = 3 * 1024**2 - 8 * 1024
MAX_PART_HEADER_BYTES = 4096
MAX_PARTS = (MAX_LOGICAL_BYTES + 65536 + MAX_FRAGMENT_BYTES - 1) // MAX_FRAGMENT_BYTES
ROW_SCHEMA = 'genosapiens.wgs-association-row'
CODEC = 'canonical-row-v1/zlib-1/base64'
_OWNED = object()
_IDENTITIES = {'IdentityMapping': IdentityMapping, 'RepeatAllele': RepeatAllele}


class AssociationCodecError(WgsScientificError):
    """Static, payload-free failure at the compact discovery boundary."""


@dataclass(frozen=True, slots=True, repr=False)
class OwnedAssociationRow:
    row: object
    canonical_bytes: bytes
    sha256: str
    _proof: tuple

    def __post_init__(self):
        self.require(self.row)

    @property
    def size_bytes(self):
        return len(self.canonical_bytes)

    def require(self, row):
        proof = self._proof
        if (type(proof) is not tuple or len(proof) != 4 or proof[0] is not _OWNED
                or row is not self.row or proof[1] is not row
                or proof[2] is not self.canonical_bytes or proof[3] != self.sha256):
            raise AssociationCodecError('foreign compact association owner')
        return self


def _own(row, content):
    if type(content) is not bytes or not 0 < len(content) <= MAX_LOGICAL_BYTES:
        raise AssociationCodecError('compact association exceeds logical capacity', category='response_capacity')
    digest = hashlib.sha256(content).hexdigest()
    return OwnedAssociationRow(row, content, digest, (_OWNED, row, content, digest))


def _document(row):
    tables = {'records': [], 'identities': []}
    keys = {name: {} for name in tables}
    seen = {name: {} for name in tables}

    def intern(name, value):
        cached = seen[name].get(id(value))
        if cached is not None and cached[0] is value:
            return cached[1]
        if name == 'identities':
            tag = type(value).__name__
            if _IDENTITIES.get(tag) is not type(value):
                raise AssociationCodecError('invalid association identity type')
            entry = {'tag': tag, 'value': value}
        else:
            if type(value) is not SourcedRecord:
                raise AssociationCodecError('invalid association source record')
            entry = value
        key = _owned_canonical_bytes(entry)
        index = keys[name].get(key)
        if index is None:
            index = len(tables[name])
            keys[name][key] = index
            tables[name].append(entry)
        seen[name][id(value)] = (value, index)
        return index

    def values(value, excluded):
        return {field.name: getattr(value, field.name) for field in fields(value) if field.name not in excluded}

    associations = []
    for proof in row.associations:
        associations.append({**values(proof, {'source_refs'}),
            'source_refs': [intern('records', source) for source in proof.source_refs]})
    links = []
    for proof in row.link_proofs:
        links.append({**values(proof, {'native', 'source_refs', 'provenance'}),
            'native': intern('identities', proof.native),
            'source_refs': [intern('records', source) for source in proof.source_refs],
            'provenance': intern('records', proof.provenance)})
    repeats = []
    for proof in row.repeat_proofs:
        repeats.append({**values(proof, {'native_identity', 'input_source'}),
            'native_identity': intern('identities', proof.native_identity),
            'input_source': intern('records', proof.input_source)})
    return dict(schema=ROW_SCHEMA, version=1, classification='controlled', llm_policy='no-llm',
        **tables, row={**values(row, {'associations', 'link_proofs', 'repeat_proofs'}),
            'associations': associations, 'link_proofs': links, 'repeat_proofs': repeats})


@_validation_session(isolated=True)
def encode_association_row(row):
    """Validate a typed row and encode its canonical compact document once."""
    from genosapiens.wgs_association_projection import AssociationRow
    try:
        validated = _typed_value(row, AssociationRow)
        return _own(row, _owned_canonical_bytes(_document(validated)))
    except AssociationCodecError:
        raise
    except (ValueError, TypeError, KeyError, RuntimeError, UnicodeError, AttributeError, RecursionError):
        raise AssociationCodecError('invalid compact association row') from None


@_validation_session(isolated=True)
def decode_association_row(content):
    """Independently validate the sole shipping codec and retain its byte owner."""
    from genosapiens import candidate_admission as admission
    from genosapiens.wgs_association_projection import AssociationRow
    try:
        if type(content) is not bytes or not 0 < len(content) <= MAX_LOGICAL_BYTES:
            raise ValueError
        value = _loads(content)
        if (type(value) is not dict or set(value) !=
                {'schema', 'version', 'classification', 'llm_policy', 'row', 'records', 'identities'}
                or value['schema'] != ROW_SCHEMA or type(value['version']) is not int or value['version'] != 1
                or value['classification'] != 'controlled' or value['llm_policy'] != 'no-llm'
                or type(value['row']) is not dict
                or any(type(value[name]) is not list for name in ('records', 'identities'))):
            raise ValueError
        records = value['records']
        for record in records:
            _typed_value(record, SourcedRecord, decode=True)
        identities = []
        for entry in value['identities']:
            if (type(entry) is not dict or set(entry) != {'tag', 'value'}
                    or type(entry['tag']) is not str or entry['tag'] not in _IDENTITIES):
                raise ValueError
            kind = _IDENTITIES[entry['tag']]
            _typed_value(entry['value'], kind, decode=True)
            identities.append((kind, entry['value']))

        def reference(table, index, kind=None):
            if type(index) is not int or not 0 <= index < len(table):
                raise ValueError
            entry = table[index]
            if kind is not None:
                if entry[0] is not kind:
                    raise ValueError
                return entry[1]
            return entry

        def rows(name, kind):
            group = value['row'][name]
            if type(group) is not list:
                raise ValueError
            names = {field.name for field in fields(kind)}
            for item in group:
                if type(item) is not dict or set(item) != names:
                    raise ValueError
                yield item

        associations = []
        for item in rows('associations', admission._Association):
            if type(item['source_refs']) is not list:
                raise ValueError
            associations.append({**item, 'source_refs': [reference(records, index) for index in item['source_refs']]})
        links = []
        for item in rows('link_proofs', admission.links.AssociationProof):
            if type(item['source_refs']) is not list:
                raise ValueError
            links.append({**item, 'native': reference(identities, item['native'], IdentityMapping),
                'source_refs': [reference(records, index) for index in item['source_refs']],
                'provenance': reference(records, item['provenance'])})
        repeats = []
        for item in rows('repeat_proofs', admission.repeat.RepeatAssociationProof):
            repeats.append({**item, 'native_identity': reference(identities, item['native_identity'], RepeatAllele),
                'input_source': reference(records, item['input_source'])})
        row = _typed_value({**value['row'], 'associations': associations,
            'link_proofs': links, 'repeat_proofs': repeats}, AssociationRow, decode=True)
        # This enforces first occurrence, table uniqueness/use, all exact fields,
        # canonical JSON spelling and order. It is independent readback work.
        if _owned_canonical_bytes(_document(row)) != content:
            raise ValueError
        return _own(row, content)
    except (ValueError, TypeError, KeyError, RuntimeError, UnicodeError, AttributeError, RecursionError):
        raise AssociationCodecError('invalid serialized compact association row') from None


@dataclass(frozen=True, slots=True, repr=False)
class AssociationPart:
    retained_ordinal: int
    part_index: int
    part_count: int
    logical_bytes: int
    logical_sha256: str
    payload: str
    schema: str = 'genosapiens.wgs-association-part'
    version: int = 1
    codec: str = CODEC

    def __post_init__(self):
        from genosapiens.evidence_model import _sha
        try:
            _sha(self.logical_sha256)
            if (self.schema != 'genosapiens.wgs-association-part' or type(self.version) is not int or self.version != 1
                    or self.codec != CODEC or type(self.retained_ordinal) is not int or not 1 <= self.retained_ordinal < 2**63
                    or type(self.part_count) is not int or not 1 <= self.part_count <= MAX_PARTS
                    or type(self.part_index) is not int or not 0 <= self.part_index < self.part_count
                    or type(self.logical_bytes) is not int or not 0 < self.logical_bytes <= MAX_LOGICAL_BYTES
                    or type(self.payload) is not str or not 0 < len(self.payload) <= 4 * ((MAX_FRAGMENT_BYTES + 2) // 3)
                    or not self.payload.isascii()):
                raise ValueError
        except (ValueError, TypeError, RuntimeError):
            raise AssociationCodecError('invalid physical association part') from None

    @property
    def commitment(self):
        return self.retained_ordinal, self.part_count, self.logical_bytes, self.logical_sha256


def iter_association_parts(owner):
    """Compress once, then own each bounded canonical physical line once."""
    from genosapiens.wgs_shards import encode_row
    if type(owner) is not OwnedAssociationRow:
        raise AssociationCodecError('compact association owner required')
    owner.require(owner.row)
    compressed = zlib.compress(owner.canonical_bytes, level=1)
    count = (len(compressed) + MAX_FRAGMENT_BYTES - 1) // MAX_FRAGMENT_BYTES
    for index in range(count):
        fragment = compressed[index * MAX_FRAGMENT_BYTES:(index + 1) * MAX_FRAGMENT_BYTES]
        part = AssociationPart(owner.row.retained_ordinal, index, count, owner.size_bytes,
            owner.sha256, base64.b64encode(fragment).decode('ascii'))
        encoded = encode_row(part)
        if (len(encoded.content) > MAX_PART_BYTES
                or len(encoded.content) - len(part.payload) > MAX_PART_HEADER_BYTES):
            raise AssociationCodecError('association part exceeds physical capacity', category='response_capacity')
        yield encoded


def read_association_rows(parts):
    """Reassemble ordered groups with bounded decompression and strict EOF."""
    source = iter(parts)
    pending = decoder = output = None
    previous = expected = 0
    try:
        for part in source:
            if type(part) is not AssociationPart:
                raise AssociationCodecError('physical association part required')
            part.__post_init__()
            if pending is None:
                if part.part_index != 0 or part.retained_ordinal <= previous:
                    raise AssociationCodecError('association part group order differs')
                pending, expected = part.commitment, 0
                decoder, output = zlib.decompressobj(), bytearray()
            if part.commitment != pending or part.part_index != expected:
                raise AssociationCodecError('incomplete or foreign association part group')
            try:
                fragment = base64.b64decode(part.payload, validate=True)
                if (not 0 < len(fragment) <= MAX_FRAGMENT_BYTES
                        or base64.b64encode(fragment).decode('ascii') != part.payload
                        or part.part_index + 1 < part.part_count and len(fragment) != MAX_FRAGMENT_BYTES):
                    raise ValueError
                remaining = part.logical_bytes - len(output)
                chunk = decoder.decompress(fragment, remaining + 1)
                if len(chunk) > remaining or decoder.unconsumed_tail or decoder.unused_data:
                    raise ValueError
                output.extend(chunk)
                del chunk, fragment
                expected += 1
                if expected < part.part_count:
                    if decoder.eof:
                        raise ValueError
                    continue
                if not decoder.eof or len(output) != part.logical_bytes:
                    raise ValueError
                content = bytes(output)
                output = None
                if hashlib.sha256(content).hexdigest() != part.logical_sha256:
                    raise ValueError
                owned = decode_association_row(content)
                if owned.row.retained_ordinal != part.retained_ordinal:
                    raise ValueError
            except (ValueError, TypeError, RuntimeError, zlib.error):
                raise AssociationCodecError('invalid compressed association group') from None
            previous = part.retained_ordinal
            pending = decoder = None
            yield owned
            del owned, content
        if pending is not None:
            raise AssociationCodecError('incomplete association part stream')
    finally:
        if hasattr(source, 'close'):
            source.close()
