"""Portable submission naming inputs and upstream-attested reference queries.

The producer validates the installed reference. Saved base observations are
bound to that reference and the producing case; they are not FASTA membership
proofs. Replay authenticates the original metadata bytes and uses the original
contig renderer, then reruns every representation decision from these inputs.
"""
from dataclasses import replace
import hashlib
from pathlib import Path, PurePosixPath
import tempfile

from genosapiens.admission_links import PinnedPublicSource
from genosapiens.bounded_json import read_bounded_bytes, MAX_GOVERNING_JSON_BYTES
from genosapiens.evidence_serialization import _loads
from genosapiens.manifest import PublicResourceManifest
from genosapiens.reference_concordance import (
    ReferenceIdentityContract, ReferenceRoleBinding, verify_reference_manifest_contract,
    verify_reference_construction,
)
from genosapiens.submission_rendering import load_submission_contig_renderer
from genosapiens.track1_submission import _OpenSubmissionReference, SubmissionReferenceBinding
from genosapiens.wgs_evidence_store import EvidenceStoreError

METADATA_KIND = "submission-reference-metadata-v1"
MAX_METADATA_BYTES = 16 * 1024**2
MAX_QUERY_BYTES = 64 * 1024**2


def _safe_relative(value):
    if (type(value) is not str or not value or "\\" in value
            or PurePosixPath(value).is_absolute() or ".." in PurePosixPath(value).parts
            or str(PurePosixPath(value)) != value):
        raise EvidenceStoreError("invalid retained reference metadata path")
    return value


def _metadata(contract_bytes, manifest_bytes, binding):
    if (type(binding) is not SubmissionReferenceBinding
            or hashlib.sha256(contract_bytes).hexdigest() != binding.contract_sha256
            or hashlib.sha256(manifest_bytes).hexdigest() != binding.manifest_sha256):
        raise EvidenceStoreError("retained reference metadata identity changed")
    contract = ReferenceIdentityContract.from_mapping(_loads(contract_bytes))
    manifest = PublicResourceManifest.from_mapping(_loads(manifest_bytes))
    verify_reference_manifest_contract(contract,manifest,contract_sha256=binding.contract_sha256)
    outputs = {item.path:item for item in manifest.outputs}
    public = contract.reference_role("public_annotation_benchmark")
    def output(role):
        return outputs[contract.expected_output(role).path]
    expected = SubmissionReferenceBinding(contract.contract_id,binding.contract_sha256,binding.manifest_sha256,
        public.role,public.identity_id,output(public.fasta_output_role).sha256,
        output(public.fai_output_role).sha256,output("contig_aliases").sha256)
    if binding != expected:
        raise EvidenceStoreError("retained reference role binding changed")
    roles = {"contig_aliases","construction_report",*(r.fai_output_role for r in contract.reference_roles)}
    files = tuple(sorted((output(role) for role in roles),key=lambda item:item.path))
    if len({item.path for item in files}) != len(files):
        raise EvidenceStoreError("duplicate retained reference metadata path")
    for item in files:
        _safe_relative(item.path)
    return contract,files


class SubmissionReferenceMetadata:
    def __init__(self, writer, pin, opened, *, name):
        self.name,self.binding = name,opened.binding
        try:
            if (type(pin) is not PinnedPublicSource or pin.contract_sha256 != self.binding.contract_sha256
                    or pin.manifest_sha256 != self.binding.manifest_sha256):
                raise EvidenceStoreError("foreign submission reference metadata")
            contract_bytes = read_bounded_bytes(pin.contract_path)
            manifest_bytes = read_bounded_bytes(pin.manifest_path,reference_root=pin.root)
            _,files = _metadata(contract_bytes,manifest_bytes,self.binding)
            size = len(contract_bytes)+len(manifest_bytes)
            def rows():
                nonlocal size
                for item in files:
                    content = read_bounded_bytes(pin.root/item.path,reference_root=pin.root,
                        expected_size_bytes=item.size_bytes)
                    size += len(content)
                    if size > MAX_METADATA_BYTES:
                        raise EvidenceStoreError("reference metadata exceeds capacity",capacity=True)
                    if hashlib.sha256(content).hexdigest() != item.sha256:
                        raise EvidenceStoreError("reference metadata checksum changed")
                    yield {"path":item.path,"bytes_hex":content.hex()}
            digest = writer.add_sequence(name+".files","submission-reference-files-v1",rows())
            self.sha256 = writer.add(name,METADATA_KIND,{
                "schema":"genosapiens.wgs-submission-reference-metadata","version":1,
                "binding":self.binding,"contract_hex":contract_bytes.hex(),"manifest_hex":manifest_bytes.hex(),
                "byte_count":size,"files":{"name":name+".files","sha256":digest}})
        except BaseException:
            writer.failed = True
            raise


def _decode_hex(value):
    if type(value) is not str or len(value) > 2*MAX_GOVERNING_JSON_BYTES:
        raise EvidenceStoreError("reference metadata file exceeds capacity",capacity=True)
    try:
        content = bytes.fromhex(value)
    except ValueError:
        raise EvidenceStoreError("invalid retained reference metadata bytes") from None
    if content.hex() != value:
        raise EvidenceStoreError("noncanonical retained reference metadata bytes")
    return content


def _no_queries(*_):
    raise EvidenceStoreError("submission reference query has no retained context owner")


class VerifiedSubmissionReferenceMetadata:
    def __init__(self, reader, root, *, expected_binding):
        self.temporary = None
        try:
            if root.kind != METADATA_KIND or type(expected_binding) is not SubmissionReferenceBinding:
                raise EvidenceStoreError("foreign submission reference metadata root")
            value = reader.reconstruct(root,max_bytes=2*MAX_METADATA_BYTES)
            from genosapiens.evidence_serialization import _canonical
            if (type(value) is not dict or set(value) != {"schema","version","binding","contract_hex",
                    "manifest_hex","byte_count","files"}
                    or value["schema"] != "genosapiens.wgs-submission-reference-metadata"
                    or type(value["version"]) is not int or value["version"] != 1
                    or _canonical(value["binding"]) != _canonical(expected_binding)):
                raise EvidenceStoreError("invalid retained submission reference metadata")
            contract_bytes,manifest_bytes = _decode_hex(value["contract_hex"]),_decode_hex(value["manifest_hex"])
            contract,files = _metadata(contract_bytes,manifest_bytes,expected_binding)
            self.temporary = tempfile.TemporaryDirectory(prefix="reference-metadata-",dir=reader.temporary.name)
            directory = Path(self.temporary.name)
            contract_path = directory/"identity.json"
            manifest_path = directory/"manifests/grch38"/(contract.resource_id+".json")
            def save(relative, content):
                path = directory/_safe_relative(relative)
                path.parent.mkdir(parents=True,exist_ok=True)
                with path.open("xb") as handle:
                    handle.write(content)
            save("identity.json",contract_bytes)
            save(manifest_path.relative_to(directory).as_posix(),manifest_bytes)
            source = value["files"]
            if type(source) is not dict or set(source) != {"name","sha256"}:
                raise EvidenceStoreError("invalid retained reference file binding")
            file_root = reader.reference(source["name"],kind="submission-reference-files-v1")
            if source["sha256"] != file_root.sha256:
                raise EvidenceStoreError("retained reference file binding changed")
            count,size = 0,len(contract_bytes)+len(manifest_bytes)
            for row in reader.iter_array(file_root):
                if (count == len(files) or type(row) is not dict or set(row) != {"path","bytes_hex"}
                        or row["path"] != files[count].path):
                    raise EvidenceStoreError("retained reference metadata coverage changed")
                content = _decode_hex(row["bytes_hex"])
                item = files[count]
                if len(content) != item.size_bytes or hashlib.sha256(content).hexdigest() != item.sha256:
                    raise EvidenceStoreError("retained reference metadata checksum changed")
                size += len(content)
                if size > MAX_METADATA_BYTES:
                    raise EvidenceStoreError("reference metadata exceeds capacity",capacity=True)
                save(item.path,content)
                count += 1
            if count != len(files) or type(value["byte_count"]) is not int or value["byte_count"] != size:
                raise EvidenceStoreError("retained reference metadata coverage changed")
            verify_reference_construction(contract,PublicResourceManifest.from_mapping(_loads(manifest_bytes)),
                directory,contract_sha256=expected_binding.contract_sha256)
            self.pin = PinnedPublicSource(contract_path,expected_binding.contract_sha256,
                directory,manifest_path,expected_binding.manifest_sha256)
            self.renderers = {}
            for role in contract.reference_roles:
                self.renderers[(role.role,role.identity_id)] = load_submission_contig_renderer(
                    identity_path=contract_path,reference_manifest_path=manifest_path,reference_root=directory,
                    reference_binding=ReferenceRoleBinding(role.role,role.identity_id))
            public_renderer = self.renderers[(expected_binding.role,expected_binding.identity_id)]
            names = {}
            for native,(canonical,_) in public_renderer._names.items():
                if canonical in names:
                    raise EvidenceStoreError("ambiguous retained public reference name")
                names[canonical] = native
            self.reference = _OpenSubmissionReference(contract,expected_binding,names,_no_queries)
        except BaseException:
            self.close()
            raise

    def close(self):
        if self.temporary is not None:
            self.temporary.cleanup()
            self.temporary = None

    def __enter__(self):
        return self

    def __exit__(self,*_):
        self.close()


class ReferenceQueries:
    """One context's exact query occurrences, with conflicting duplicates refused."""
    def __init__(self, writer, reference, *, name, maximum_count):
        self.writer,self.connection,self.name = writer,writer.connection,name
        self.reference,self.maximum_count = reference,maximum_count
        self.table = "reference_queries_"+hashlib.sha256(name.encode()).hexdigest()
        self.count,self.closed = 0,False
        self.connection.execute(f"CREATE TABLE {self.table} (ordinal INTEGER PRIMARY KEY,contig TEXT,position INTEGER,length INTEGER,bases TEXT)")
        self.connection.execute(f"CREATE INDEX {self.table}_key ON {self.table}(contig,position,length)")
        self.opened = replace(reference,fetch_public=self.fetch)

    def fetch(self, contig, position, length):
        if self.closed:
            raise EvidenceStoreError("retained reference queries are closed")
        if (type(contig) is not str or not contig or len(contig) > 16384
                or type(position) is not int or not 0 < position < 2**63
                or type(length) is not int or not 0 < length <= MAX_QUERY_BYTES
                or self.count == self.maximum_count):
            raise EvidenceStoreError("reference query exceeds capacity",capacity=True)
        try:
            bases = self.reference.fetch_public(contig,position,length)
        except (ValueError,KeyError,IndexError,OSError):
            bases = None
        if bases is not None and (type(bases) is not str or len(bases) != length):
            raise EvidenceStoreError("invalid reference query result")
        key = contig,position,length
        prior = self.connection.execute(f"SELECT bases FROM {self.table} WHERE contig=? AND position=? AND length=? LIMIT 1",key).fetchone()
        if prior is not None and prior[0] != bases:
            raise EvidenceStoreError("conflicting retained reference query")
        self.connection.execute(f"INSERT INTO {self.table} VALUES (?,?,?,?,?)",(self.count,*key,bases))
        self.count += 1
        if bases is None:
            raise ValueError("retained reference query unavailable")
        return bases

    def rows(self):
        if self.closed:
            raise EvidenceStoreError("retained reference queries are closed")
        for contig,position,length,bases in self.connection.execute(
                f"SELECT contig,position,length,bases FROM {self.table} ORDER BY ordinal"):
            yield {"contig":contig,"position":position,"length":length,"bases":bases}

    def publish(self):
        self.sha256 = self.writer.add_sequence(self.name,"submission-reference-queries-v1",self.rows())

    def close(self):
        if not self.closed:
            self.connection.execute("DROP TABLE "+self.table)
            self.closed = True


def replay_query_source(reader, root, reference):
    """A sequential source that refuses missing, extra or differently requested bases."""
    rows = iter(reader.iter_array(root))
    def fetch(contig, position, length):
        row = next(rows,None)
        if (type(row) is not dict or set(row) != {"contig","position","length","bases"}
                or type(row["position"]) is not int or type(row["length"]) is not int
                or (row["contig"],row["position"],row["length"]) != (contig,position,length)):
            raise EvidenceStoreError("retained reference query coverage changed")
        bases = row["bases"]
        if bases is None:
            raise ValueError("retained reference query unavailable")
        if type(bases) is not str or len(bases) != length:
            raise EvidenceStoreError("invalid retained reference query result")
        return bases
    def finish():
        if next(rows,None) is not None:
            raise EvidenceStoreError("unused retained reference queries")
    return replace(reference,fetch_public=fetch),finish,rows
