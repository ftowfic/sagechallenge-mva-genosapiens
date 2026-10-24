"""Operation-owned immutable header maps and finite FORMAT layout caches."""

from collections import OrderedDict
from dataclasses import dataclass
import sys
from types import MappingProxyType

from genosapiens.vcf_contracts import VcfContractError


@dataclass(frozen=True, slots=True)
class PreparedVcfHeader:
    header: object
    contig_lengths: object
    format_definitions: frozenset[str]
    field_declarations: object
    pgt_pid_declared: bool


def prepare_header(header):
    from genosapiens.vcf_qc import VcfQcStreamHeader, _structured_attributes
    if type(header) is not VcfQcStreamHeader:
        raise VcfContractError("invalid prepared VCF header")
    definitions = {}
    for raw in header.raw_lines:
        line = raw.rstrip("\r\n")
        if line.startswith(("##FORMAT=<", "##INFO=<")):
            attributes = _structured_attributes(line.split("=<", 1)[1][:-1])
            if attributes is not None and "ID" in attributes:
                namespace = line[2:].split("=", 1)[0]
                definitions[namespace, attributes["ID"]] = (attributes.get("Number"), attributes.get("Type"))
    return PreparedVcfHeader(header, MappingProxyType({c.contig: c.length for c in header.declared_contigs}),
        frozenset(header.format_definitions), MappingProxyType(definitions),
        {"PGT", "PID"}.issubset(header.format_definitions))


class FormatLayoutCache:
    """One read owner retains at most 128 layouts and 256 KiB of charged storage."""

    def __init__(self, *, max_entries=128, max_bytes=256 * 1024):
        if any(type(value) is not int or value < 0 for value in (max_entries, max_bytes)):
            raise VcfContractError("invalid FORMAT cache limits")
        self.max_entries, self.max_bytes = max_entries, max_bytes
        self._entries = OrderedDict()
        self.retained_bytes = 0
        self.builds = 0

    def positions(self, keys, *, allowed=None):
        if (type(keys) is not tuple or len(set(keys)) != len(keys)
                or any(type(key) is not str or not key or any(c.isspace() for c in key) for key in keys)
                or "GT" in keys and keys[0] != "GT"
                or allowed is not None and any(key not in allowed for key in keys)):
            raise VcfContractError("invalid FORMAT layout")
        if keys in self._entries:
            self._entries.move_to_end(keys)
            return self._entries[keys][0]
        positions = {key: index for index, key in enumerate(keys)}
        self.builds += 1
        # Charge retained keys, values, mapping, and a conservative node allowance.
        size = (512 + sys.getsizeof(keys) + sys.getsizeof(positions)
                + sum(sys.getsizeof(key) + sys.getsizeof(index) + 64 for key, index in positions.items()))
        positions = MappingProxyType(positions)
        if self.max_entries and size <= self.max_bytes:
            while self._entries and (len(self._entries) >= self.max_entries
                                     or self.retained_bytes + size > self.max_bytes):
                _, (_, old_size) = self._entries.popitem(last=False)
                self.retained_bytes -= old_size
            self._entries[keys] = positions, size
            self.retained_bytes += size
        return positions

    def __len__(self):
        return len(self._entries)
