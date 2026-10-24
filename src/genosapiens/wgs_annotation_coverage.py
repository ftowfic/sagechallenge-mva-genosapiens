"""Bounded, authenticated reference-name coverage for annotation outcomes."""
from dataclasses import dataclass, field
import hashlib
from pathlib import Path
import re

from genosapiens.bounded_json import read_bounded_bytes

REASON = 'contig-outside-pinned-vep-reference'
POLICY = 'missing-response-outside-pinned-vep-reference-v1'
ALIAS_REASON = 'nonprimary-contig-alias-unresolved'
ALIAS_POLICY = 'missing-response-unresolved-nonprimary-alias-v1'
REASON_POLICIES = {REASON: POLICY, ALIAS_REASON: ALIAS_POLICY}
MAX_INDEX_BYTES = 1024 * 1024


def _name_tokens(name):
    bare = name[3:] if name.startswith('chr') else name
    return {name, bare, bare.split('.', 1)[0], *re.findall(r'(?:GL|KI|JH|KN|KQ|KV|KZ|ML|MU)[0-9]+', bare)}


def _primary_name(contig):
    bare = contig[3:] if contig.startswith('chr') else contig
    accession = re.fullmatch(r'NC_0*([0-9]+)(?:\.[0-9]+)?', bare)
    return bool(bare in {'X', 'Y', 'M', 'MT'}
        or re.fullmatch(r'(?:[1-9]|1[0-9]|2[0-2])', bare)
        or accession and (1 <= int(accession[1]) <= 24 or int(accession[1]) == 12920))


@dataclass(frozen=True, slots=True, repr=False)
class ReferenceCoverage:
    content: bytes
    sha256: str
    names: frozenset[str] = field(init=False)
    alias_tokens: frozenset[str] = field(init=False)

    def __post_init__(self):
        if (type(self.content) is not bytes or not 0 < len(self.content) <= MAX_INDEX_BYTES
                or hashlib.sha256(self.content).hexdigest() != self.sha256):
            raise ValueError('annotation reference index binding is invalid')
        names = set()
        try:
            for line in self.content.decode('ascii').splitlines():
                name, length, offset, bases, width = line.split('\t')
                if (not name or name in names or any(c.isspace() for c in name)
                        or int(length) <= 0 or int(offset) < 0
                        or int(bases) <= 0 or int(width) < int(bases)):
                    raise ValueError
                names.add(name)
        except (ValueError, UnicodeError):
            raise ValueError('annotation reference index is invalid') from None
        if not names:
            raise ValueError('annotation reference index is empty')
        object.__setattr__(self, 'names', frozenset(names))
        object.__setattr__(self, 'alias_tokens', frozenset(token for name in names for token in _name_tokens(name)))

    def absent(self, contig):
        # Conservative aliases prevent uncertain naming from becoming absence.
        if _primary_name(contig):
            return False
        # A possible version/accession alias is uncertainty, never proof of absence.
        return _name_tokens(contig).isdisjoint(self.alias_tokens)

    def unassessed_reason(self, contig):
        """Keep absence separate from an unverified nonprimary alias match."""
        if self.absent(contig):
            return REASON
        bare = contig[3:] if contig.startswith('chr') else contig
        if not _primary_name(contig) and contig not in self.names and bare not in self.names:
            return ALIAS_REASON
        return None


def load_reference_coverage(fasta_path, *, expected_sha256, expected_bytes):
    if type(expected_bytes) is not int or not 0 < expected_bytes <= MAX_INDEX_BYTES:
        raise ValueError('annotation reference index size is invalid')
    content = read_bounded_bytes(Path(str(fasta_path) + '.fai'), maximum_bytes=MAX_INDEX_BYTES)
    if len(content) != expected_bytes:
        raise ValueError('annotation reference index size differs')
    return ReferenceCoverage(content, expected_sha256)
