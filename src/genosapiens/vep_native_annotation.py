"""Pinned VEP structural effects and consequence-free native geometry context.

Internal backend: the main runner verifies resources, reference bindings, network,
and publication storage before calling this module. No acquisition or allele
conversion occurs here. All records and subprocess diagnostics are controlled.
"""

from collections import Counter
from collections.abc import Mapping
from dataclasses import asdict, dataclass, replace
from decimal import Decimal
import json
import math
import os
from pathlib import Path
import re
import subprocess
from types import MappingProxyType
from urllib.parse import quote

from genosapiens.allele_identity import ReferenceRole, RepeatAllele, StructuralAllele, StructuralVariantType, identity_to_dict
from genosapiens.evidence_model import AlleleEvidence, SourcedRecord
from genosapiens.evidence_serialization import _plain
from genosapiens.vep_annotation import (
    VepAnnotationError, _Controlled, _digest, _file_identity, _json,
    _preflight_resources, _transcript,
)


SCHEMA = "genosapiens.vep-native-annotation"
INTERVAL_RECIPE_SCHEMA = "genosapiens.native-interval-recipe"
INTERVAL_RECIPE_DATASET = "genosapiens.native-interval-recipe"
INTERVAL_RECIPE_CONTRACT = {
    "schema": INTERVAL_RECIPE_SCHEMA,
    "fields": ["schema", "record_id", "reference", "event_type", "event_subtype", "interval_semantics", "contig", "start", "end"],
    "reference_fields": ["role", "identity_id", "manifest_sha256"],
    "endpoint_fields": ["min", "max"],
    "interval_semantics": ["affected_reference_bases", "breakpoint_envelope", "unknown"],
    "coordinates": "1-based-closed", "supported_exact_subtypes": ["DEL", "DUP", "INV", "CNV"],
    "source_fields": ["recipe_utf8"], "contributors": "all-identity-mappings-and-all-call-sources",
}
BIN_SIZE = 1000000
SV_CLASSES = {"DEL": "deletion", "DUP": "duplication", "INV": "inversion", "CNV": "copy_number_variation"}
NATIVE_FLAGS = (
    "--cache", "--offline", "--cache_version", "116", "--species", "homo_sapiens",
    "--assembly", "GRCh38", "--format", "ensembl", "--json", "--max_sv_size", "-1",
    "--no_stats", "--no_progress", "--terms", "SO", "--hgvs", "--hgvsg",
    "--protein", "--protein_version", "--gene_version", "--transcript_version",
    "--symbol", "--numbers", "--domains", "--canonical", "--mane", "--biotype",
    "--uniprot", "--tsl", "--appris", "--variant_class", "--sift", "b",
    "--polyphen", "b", "--allele_number", "--check_ref", "--use_given_ref",
)
SOURCE_REVISIONS = {
    "ensembl-vep": "2cb0bbe216bb31c75de8f8000e2da7ff4fb7b451",
    "ensembl": "0d852313c420a1d5128b6774635125ebbac03350",
    "ensembl-variation": "db7394a02523a6e519cb0d89f84b2889bac12737",
    "ensembl-io": "6afb5dc27a5ae6881d75959153fbe6e9a4a7e788",
    "ensembl-compara": "aa58db6179591cfa35cdd6d5a14767355d63433a",
    "ensembl-funcgen": "90049ea7ee4d8ae3a6d298dca46d6c6ab20538c4",
}
SOURCE_HASHES = {
    "modules/Bio/EnsEMBL/VEP/Config.pm": "523ef2348fc0d4e16258f4a053b1cb35ce2dd038bc349430e84c9fafa8bfa6b9",
    "modules/Bio/EnsEMBL/VEP/Parser.pm": "67da534caa71976d7152c135124e631338a4dfd7111d49f25a979f5818ebf5a2",
    "modules/Bio/EnsEMBL/VEP/Parser/VEP_input.pm": "1fe20af54abe77f904993a9cd1ecb886af30391ee12f75762c7a3b5931be1dbf",
    "modules/Bio/EnsEMBL/VEP/OutputFactory.pm": "f4e645db89724b13273c3ce65b9a790c67b36fe829016aa50e595da6638fd90f",
    "modules/Bio/EnsEMBL/VEP/OutputFactory/JSON.pm": "a80e981fcc902f8d00ded16508994965144d325c536d1cdde5fe75b4d7d165af",
    "modules/Bio/EnsEMBL/VEP/Utils.pm": "80e75270ec014a8b9d3ef43c7ba8872dadf12a62b0d3e88cbf2a290ff1dbc050",
    "modules/Bio/EnsEMBL/VEP/CacheDir.pm": "34d31b5bf840329d6533934fc8d04a4c30a5109b389495cc42a346dd53a40b47",
    "modules/Bio/EnsEMBL/VEP/AnnotationSource/Cache/Transcript.pm": "88a5090edf2284f92b3426e21e165a9b06262709cb9215c6f1ccb5e003cbeaff",
    "modules/Bio/EnsEMBL/VEP/AnnotationSource/Cache/BaseSerialized.pm": "1bda2b7488454db4159a2ba51a44b7f26f406800b4823e0af144c4501082c3c3",
    "Bio/EnsEMBL/TranscriptMapper.pm": "a9769ac71d1877964cc0e5daebb94393f76737d1a57de53f83da5ca8d64b4218",
    "Bio/EnsEMBL/Variation/BaseTranscriptVariation.pm": "36dcf3f23107d317012932cb7f1565b91a23a8b62ab24fe10ea1d8f661f9470f",
    "Bio/EnsEMBL/Variation/Utils/VariationEffect.pm": "6dcef8278c2e24b37eb3389a1e9d00f593336549f2d9ce72ca11271476c16faf",
}
PROTOCOL = {
    "schema": SCHEMA, "coordinates": "1-based-closed", "cache_bin_size": BIN_SIZE,
    "query_fields": ["query_id", "contig", "start", "end", "ordinal"],
    "roles": ["affected_interval", "structural_interval", "breakpoint_anchor", "breakpoint_uncertainty", "repeat_locus"],
    "scopes": ["structural_effect", "transcript_context"],
    "context_fields": ["schema", "query", "cache_bin_size", "cache_versions", "cache_files", "transcripts"],
    "context_transcript_fields": ["transcript_id", "gene_id", "protein_id", "gene_symbol", "gene_symbol_source",
        "hgnc_id", "biotype", "canonical", "raw_fields", "attributes", "transcript_geometry", "exons", "introns",
        "cds", "overlapping_exons", "overlapping_introns", "overlapping_cds", "mappings", "peptide_length",
        "protein_features", "domain_overlaps"],
    "native_identity_fields": ["id", "input", "seq_region_name", "start", "end", "strand", "assembly_name",
        "variant_class", "allele_string"],
    "joined_claims": {"identifiers": "raw-stable-cores-and-explicit-versions",
        "domains": "pinned-word-token-capture-carry-and-numberify-projection; positive-claim-multiplicity",
        "protein_positions": "mapper-first-last-coordinates; gap-null; peptide-length-plus-mapped-stop"},
}

# Only cached transcript objects are inspected. In particular, no VariationFeature
# is created for a context query and no variant consequence predicate is called.
CACHE_HELPER = r'''use strict;
use warnings;
use JSON::PP;
use Digest::SHA;
use Bio::EnsEMBL::VEP::Config;
use Bio::EnsEMBL::VEP::CacheDir;
use Bio::EnsEMBL::TranscriptMapper;

my ($root, $fasta, $input, $output) = @ARGV;
die "invalid helper arguments\n" unless @ARGV == 4;
my $json = JSON::PP->new->canonical->ascii;
my $config = Bio::EnsEMBL::VEP::Config->new({offline => 1, cache => 1,
  cache_version => 116, species => 'homo_sapiens', assembly => 'GRCh38',
  dir => '.', dir_cache => $root, fasta => $fasta, no_stats => 1, no_progress => 1,
  warning_file => 'context.warnings'});
my $cache = Bio::EnsEMBL::VEP::CacheDir->new({config => $config, root_dir => $root});
my ($source) = grep {$_->isa('Bio::EnsEMBL::VEP::AnnotationSource::Cache::Transcript')}
  @{$cache->get_all_AnnotationSources};
die "unsupported cache layout\n" unless $source && $source->{cache_region_size} == 1000000
  && $source->serializer_type eq 'storable';
my $versions = $cache->version_data;
die "unsupported cache version\n" unless $versions->{assembly} eq 'GRCh38.p14'
  && $versions->{genebuild} eq 'GENCODE50' && $versions->{gencode} eq 'GENCODE 50';

sub versioned {
  my ($id, $version) = @_;
  return undef unless defined $id && $id ne '-';
  return $id =~ /\.[0-9]+$/ || !defined($version) ? $id : "$id.$version";
}
sub blocks {
  my ($features) = @_;
  my $rank = 0;
  return [map {my $f = $_; {start => 0 + $f->start, end => 0 + $f->end,
    strand => 0 + $f->strand, rank => ++$rank, total => scalar(@$features),
    stable_id => $f->can('stable_id') ? $f->stable_id : undef,
    version => $f->can('version') ? $f->version : undef}} @$features];
}
sub overlap {
  my ($f, $q) = @_;
  return $f->{start} <= $q->{end} && $f->{end} >= $q->{start};
}

open my $in, '<', $input or die "input unavailable\n";
open my $out, '>', $output or die "output unavailable\n";
my %ids;
while (my $line = <$in>) {
  my $q = $json->decode($line);
  die "invalid query\n" unless ref($q) eq 'HASH'
    && join(',', sort keys %$q) eq 'contig,end,ordinal,query_id,start'
    && !ref($q->{query_id}) && $q->{query_id} =~ /^native-query-[0-9a-f]{64}$/
    && !$ids{$q->{query_id}}++ && !ref($q->{contig}) && $q->{contig} =~ /^[A-Za-z0-9_.-]+$/;
  for my $key (qw(start end ordinal)) {
    die "invalid query coordinate\n" if !defined($q->{$key}) || ref($q->{$key}) || $q->{$key} !~ /^[0-9]+$/;
  }
  die "invalid query span\n" unless $q->{start} >= 1 && $q->{end} >= $q->{start};
  my (%records, %fingerprints, @files);
  for my $bin (int(($q->{start} - 1) / 1000000) .. int(($q->{end} - 1) / 1000000)) {
    my $relative = $q->{contig}.'/'.($bin * 1000000 + 1).'-'.(($bin + 1) * 1000000).'.gz';
    my $file = $source->get_dump_file_name($q->{contig}, $bin * 1000000 + 1, ($bin + 1) * 1000000);
    my $sha;
    if (-f $file) {
      open my $fh, '<', $file or die "cache unavailable\n";
      binmode $fh;
      $sha = Digest::SHA->new(256)->addfile($fh)->hexdigest;
      close $fh;
    }
    push @files, {bin => $bin, path => $relative, sha256 => $sha};
    my $features = $source->get_features_by_regions_uncached([[$q->{contig}, $bin]]);
    for my $tr (@$features) {
      next unless $tr->slice->seq_region_name eq $q->{contig}
        && $tr->start <= $q->{end} && $tr->end >= $q->{start};
      my $tl = $tr->translation;
      my $vc = $tr->{_variation_effect_feature_cache} || {};
      my $length = defined($vc->{peptide}) ? length($vc->{peptide}) : undef;
      my $mapper = Bio::EnsEMBL::TranscriptMapper->new($tr);
      my %maps;
      for my $method (qw(genomic2cdna genomic2cds genomic2pep)) {
        $maps{$method} = [map {{class => ref($_), start => 0 + $_->start, end => 0 + $_->end}}
          $mapper->$method($q->{start}, $q->{end}, $tr->strand)];
      }
      my $domains = exists($vc->{protein_features}) ? [map {
        my $f = $_;
        my $a = $f->analysis;
        {start => 0 + $f->start, end => 0 + $f->end, hseqname => $f->hseqname,
          analysis => {map {$_ => $a->$_()} qw(display_label db db_version program program_version)}}
      } @{$vc->{protein_features}}] : undef;
      my @pep = grep {$_->{class} eq 'Bio::EnsEMBL::Mapper::Coordinate'
        && defined($length) && $_->{start} >= 1 && $_->{start} <= $length}
        @{$maps{genomic2pep}};
      my $overlaps = defined($domains) && defined($length) ? [grep {
        my $d = $_;
        scalar grep {$d->{start} <= ($_->{end} < $length ? $_->{end} : $length)
          && $d->{end} >= $_->{start}} @pep
      } @$domains] : undef;
      my $record = {
        transcript_id => versioned($tr->stable_id, $tr->version),
        gene_id => versioned($tr->{_gene_stable_id}, $tr->{_gene_version}),
        protein_id => $tl ? versioned($tl->stable_id, $tl->version) : undef,
        gene_symbol => $tr->{_gene_symbol}, gene_symbol_source => $tr->{_gene_symbol_source},
        hgnc_id => $tr->{_gene_hgnc_id}, biotype => $tr->biotype,
        canonical => defined($tr->is_canonical) ? 0 + $tr->is_canonical : undef,
        raw_fields => {map {$_ => $tr->{$_}} qw(_gene_stable_id _gene_version _gene_symbol
          _gene_symbol_source _gene_hgnc_id _protein _swissprot _trembl _uniparc _uniprot_isoform _refseq)},
        attributes => [map {{code => $_->code, value => $_->value}} @{$tr->get_all_Attributes}],
        transcript_geometry => {contig => $tr->slice->seq_region_name, start => 0 + $tr->start,
          end => 0 + $tr->end, strand => 0 + $tr->strand},
        exons => blocks($tr->get_all_Exons), introns => blocks($tr->get_all_Introns),
        cds => $tl ? blocks($tr->get_all_CDS) : [], mappings => \%maps,
        peptide_length => $length, protein_features => $domains, domain_overlaps => $overlaps,
      };
      for my $kind (qw(exons introns cds)) {
        $record->{'overlapping_'.$kind} = [grep {overlap($_, $q)} @{$record->{$kind}}];
      }
      # Compare the entire retained source projection, not just the stable ID.
      my $key = $tr->stable_id;
      my $fingerprint = $json->encode($record);
      die "conflicting cached transcript copies\n" if exists($fingerprints{$key}) && $fingerprints{$key} ne $fingerprint;
      $fingerprints{$key} = $fingerprint;
      $records{$key} = $record;
    }
    $source->clean_cache;
  }
  print $out $json->encode({schema => 'genosapiens.vep-native-annotation', query => $q,
    cache_bin_size => 1000000, cache_versions => $versions, cache_files => \@files,
    transcripts => [map {$records{$_}} sort keys %records]})."\n";
}
close $in;
close $out or die "output incomplete\n";
'''


def configuration():
    return {"schema": SCHEMA, "schema_sha256": _digest(_json(PROTOCOL)),
            "flags": list(NATIVE_FLAGS), "helper_sha256": _digest(CACHE_HELPER),
            "cache_bin_size": BIN_SIZE, "source_revisions": dict(SOURCE_REVISIONS),
            "source_sha256": dict(SOURCE_HASHES),
            "interval_recipe_contract": INTERVAL_RECIPE_CONTRACT,
            "interval_recipe_contract_sha256": _digest(_json(INTERVAL_RECIPE_CONTRACT)),
            "effect_context_join": "exact-versioned-transcript; source-separated-domains; no-extra-transcript-row"}


def _interval_recipe(text):
    data = _strict_json(text)
    if type(data) is not dict or set(data) != set(INTERVAL_RECIPE_CONTRACT["fields"]):
        raise ValueError
    if data["schema"] != INTERVAL_RECIPE_SCHEMA:
        raise ValueError
    for key in ("record_id", "contig", "event_type", "event_subtype"):
        _text(data[key])
    if data["event_type"] not in {v.value for v in StructuralVariantType} or data["interval_semantics"] not in INTERVAL_RECIPE_CONTRACT["interval_semantics"]:
        raise ValueError
    reference = data["reference"]
    if type(reference) is not dict or set(reference) != set(INTERVAL_RECIPE_CONTRACT["reference_fields"]):
        raise ValueError
    if (reference["role"] not in {v.value for v in ReferenceRole}
            or re.fullmatch(r"[a-z0-9][a-z0-9._-]*", _text(reference["identity_id"])) is None
            or re.fullmatch(r"[0-9a-f]{64}", _text(reference["manifest_sha256"])) is None):
        raise ValueError
    for key in ("start", "end"):
        endpoint = data[key]
        if type(endpoint) is not dict or set(endpoint) != {"min", "max"}:
            raise ValueError
        for value in endpoint.values():
            if value is not None:
                _integer(value)
        if endpoint["min"] is not None and endpoint["max"] is not None and endpoint["min"] > endpoint["max"]:
            raise ValueError
    if data["start"]["min"] is not None and data["end"]["max"] is not None and data["start"]["min"] > data["end"]["max"]:
        raise ValueError
    return data


def native_interval_recipe_source(recipe_utf8: bytes, *, public_manifest_sha256: str | None = None) -> SourcedRecord:
    """Bind a closed source recipe's actual bytes, not a caller precision boolean.

    This is a source-format adapter, not a constructor from a StructuralAllele.
    It neither rewrites existing sources nor excuses other contributing records.
    Caller-specific VCF semantics are not defined by this recipe format.
    """
    try:
        if type(recipe_utf8) is not bytes:
            raise ValueError
        text = recipe_utf8.decode("utf-8")
        data = _interval_recipe(text)
        return SourcedRecord(INTERVAL_RECIPE_DATASET, "current", data["record_id"], _digest(text),
                             public_manifest_sha256, {"recipe_utf8": text})
    except (ValueError, TypeError, KeyError, AttributeError, OverflowError, RecursionError):
        raise VepAnnotationError("invalid native interval recipe") from None


def native_interval_precision(allele) -> SourcedRecord:
    """Intersect every bound source claim; unknown formats cannot grant precision."""
    if type(allele) is not AlleleEvidence:
        raise VepAnnotationError("native structural interval identity required")
    native = next(i for i in allele.identities if i.identity_id == allele.native_identity_id)
    identity = native.identity
    if type(identity) is not StructuralAllele:
        raise VepAnnotationError("native structural interval identity required")
    contributors = [("identity", i.identity_id, i.stage.value, i, i.source) for i in sorted(allele.identities, key=lambda i: i.identity_id)]
    contributors += [("call", n, None, None, c.source) for n, c in enumerate(allele.calls or ())]
    assessments = []
    for kind, origin, stage, mapping, source in contributors:
        recipe, reason = None, None
        try:
            if (source.source_dataset != INTERVAL_RECIPE_DATASET or source.source_release != "current"
                    or source.source_fields is None or set(source.source_fields) != {"recipe_utf8"}):
                raise ValueError
            text = source.source_fields["recipe_utf8"]
            if type(text) is not str or _digest(text) != source.source_sha256:
                raise ValueError
            recipe = _interval_recipe(text)
            if recipe["record_id"] != source.source_record_id:
                raise ValueError
        except (ValueError, TypeError, KeyError, AttributeError, OverflowError, RecursionError):
            reason = "source-interval-contract-unrecognized-or-invalid"
        if reason is None:
            expected_reference = {"role": identity.reference_role.value, "identity_id": identity.reference_identity_id,
                                  "manifest_sha256": native.reference_manifest_sha256}
            if type(identity) is not StructuralAllele or identity.interval is None:
                reason = "native-affected-interval-not-established"
            elif recipe["interval_semantics"] != "affected_reference_bases":
                reason = "source-affected-interval-not-established"
            elif recipe["event_type"] != identity.sv_type.value or recipe["event_subtype"] != identity.sv_type.value or recipe["event_subtype"] not in SV_CLASSES:
                reason = "source-exact-subtype-not-supported-or-conflicting"
            elif recipe["reference"] != expected_reference or recipe["contig"] != identity.interval.contig:
                reason = "source-interval-reference-conflict"
            elif any(recipe[k]["min"] is None or recipe[k]["max"] is None for k in ("start", "end")):
                reason = "source-endpoint-precision-unknown"
            elif any(recipe[k]["min"] != recipe[k]["max"] for k in ("start", "end")):
                reason = "source-endpoints-imprecise"
            elif (recipe["start"]["min"], recipe["end"]["min"]) != (identity.interval.start, identity.interval.end):
                reason = "source-affected-interval-conflict"
        if mapping is not None:
            original = mapping.identity
            if (type(original) is not StructuralAllele
                    or any(b.confidence_interval not in (None, (0, 0)) for b in original.breakpoints)):
                reason = reason or "contributing-identity-precision-not-established"
            elif (type(identity) is not StructuralAllele or original.interval != identity.interval
                    or original.sv_type != identity.sv_type or original.reference_role != identity.reference_role
                    or original.reference_identity_id != identity.reference_identity_id
                    or mapping.reference_manifest_sha256 != native.reference_manifest_sha256):
                reason = reason or "contributing-identity-interval-conflict"
        assessments.append({"origin_kind": kind, "origin": origin, "stage": stage,
            "source": _plain(source), "identity_mapping": _plain(mapping) if mapping is not None else None,
            "recipe": recipe, "status": "established" if reason is None else "context_only", "reason": reason})
    fields = {"schema": INTERVAL_RECIPE_SCHEMA, "contract_sha256": _digest(_json(INTERVAL_RECIPE_CONTRACT)),
        "allele_id": allele.allele_id, "native_identity_id": allele.native_identity_id,
        "status": "established" if assessments and all(a["reason"] is None for a in assessments) else "context_only",
        "contributors": assessments}
    digest = _digest(_json(fields))
    return SourcedRecord("genosapiens.native-interval-precision", "current", "interval-precision-" + digest,
                        digest, None, fields)


def _verify_sources():
    root = Path("/opt/vep/src/ensembl-vep")
    paths = [root / path for path in SOURCE_HASHES]
    paths += [root / ".version" / name for name in SOURCE_REVISIONS if name != "ensembl-vep"]
    _preflight_resources(files=tuple(paths))
    for name, digest in SOURCE_HASHES.items():
        if _file_identity(root / name)[0] != digest:
            raise VepAnnotationError("native annotation source pin mismatch")
    for name, revision in SOURCE_REVISIONS.items():
        if name != "ensembl-vep" and (root / ".version" / name).read_text().split() != ["release", "116", "sub", revision]:
            raise VepAnnotationError("native annotation source version mismatch")


@dataclass(frozen=True, slots=True, repr=False)
class NativeQuery(_Controlled):
    query_id: str
    role: str
    annotation_scope: str
    contig: str
    start: int
    end: int
    ordinal: int
    source_geometry: Mapping
    sv_type: str | None = None
    reason: str | None = None
    interval_evidence: SourcedRecord | None = None

    def __post_init__(self):
        try:
            if (type(self.query_id) is not str or re.fullmatch(r"native-query-[0-9a-f]{64}", self.query_id) is None
                    or self.role not in PROTOCOL["roles"] or self.annotation_scope not in PROTOCOL["scopes"]
                    or not isinstance(self.source_geometry, Mapping)
                    or self.interval_evidence is not None and type(self.interval_evidence) is not SourcedRecord
                    or (self.sv_type in SV_CLASSES) != (self.annotation_scope == "structural_effect")
                    or self.reason is None and re.fullmatch(r"[A-Za-z0-9_.-]+", self.contig) is None
                    or _integer(self.start) > _integer(self.end)):
                raise ValueError
            _integer(self.ordinal, 0)
            _text(self.contig)
            if self.reason is not None:
                _text(self.reason)
            object.__setattr__(self, "source_geometry", MappingProxyType(dict(self.source_geometry)))
        except (ValueError, TypeError):
            raise VepAnnotationError("invalid native annotation query") from None

    def transport(self):
        return {"query_id": self.query_id, "contig": self.contig, "start": self.start,
                "end": self.end, "ordinal": self.ordinal}

    def line(self):
        return f"{self.contig}\t{self.start}\t{self.end}\t{self.sv_type}\t+\t{self.query_id}\n"


def native_queries(request, coordinate_maps, reference, cache_contigs):
    """Preserve each source geometry independently; never derive an event span."""
    native = next(i for i in request.allele.identities if i.identity_id == request.identity_id)
    identity = native.identity
    if type(identity) not in (StructuralAllele, RepeatAllele):
        raise VepAnnotationError("native geometry identity required")
    coordinate_map = coordinate_maps.get((identity.reference_role.value, identity.reference_identity_id))
    common_reason = None
    if coordinate_map is None:
        common_reason = "source-reference-not-coordinate-bound"
    elif native.reference_manifest_sha256 != coordinate_map[0]:
        common_reason = "source-reference-manifest-mismatch"
    geometries = []
    precision = native_interval_precision(request.allele) if type(identity) is StructuralAllele else None
    sv_type = identity.sv_type.value if type(identity) is StructuralAllele else None
    if type(identity) is RepeatAllele:
        geometries.append(("repeat_locus", "transcript_context", 0, asdict(identity.interval), asdict(identity.interval)))
    else:
        if identity.interval is not None:
            precise = precision.source_fields["status"] == "established"
            geometries.append(("affected_interval" if precise else "structural_interval",
                "structural_effect" if precise else "transcript_context", 0,
                asdict(identity.interval), asdict(identity.interval)))
        for ordinal, point in enumerate(identity.breakpoints):
            anchor = {"contig": point.contig, "start": point.position, "end": point.position}
            geometries.append(("breakpoint_anchor", "transcript_context", ordinal, anchor, asdict(point)))
            if point.confidence_interval is not None:
                lo, hi = point.confidence_interval
                window = dict(anchor, start=point.position + lo, end=point.position + hi)
                geometries.append(("breakpoint_uncertainty", "transcript_context", ordinal, window, asdict(point)))
    queries = []
    for role, scope, ordinal, geometry, source_geometry in geometries:
        contig, start, end = geometry["contig"], geometry["start"], geometry["end"]
        reason = common_reason
        if reason is None:
            contig = coordinate_map[1].get(contig)
            if contig is None or contig not in reference.index or re.fullmatch(r"[A-Za-z0-9_.-]+", contig) is None:
                reason = "contig-outside-resource-alias-coverage"
                contig = geometry["contig"]
            elif not 1 <= start <= end <= reference.index[contig].length:
                reason = "geometry-span-out-of-bounds"
            elif contig not in cache_contigs:
                reason = "contig-outside-cache-coverage"
        key = "native-query-" + _digest(_json([SCHEMA, request.allele.allele_id, request.identity_id,
            native.reference_manifest_sha256, identity_to_dict(identity), role, scope, ordinal, geometry,
            precision.source_sha256 if precision else None]))
        queries.append(NativeQuery(key, role, scope, contig, start, end, ordinal,
                                   source_geometry, sv_type if scope == "structural_effect" else None, reason, precision))
    return tuple(queries)


def _strict_json(text):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError
            result[key] = value
        return result
    def invalid(_):
        raise ValueError
    def validate(value):
        if type(value) is float and not math.isfinite(value):
            raise ValueError
        if type(value) is str:
            value.encode("utf-8", errors="strict")
        if type(value) is dict:
            for key, child in value.items():
                validate(key)
                validate(child)
        elif type(value) is list:
            for child in value:
                validate(child)
    value = json.loads(text, object_pairs_hook=pairs, parse_constant=invalid)
    validate(value)
    return value


def _integer(value, minimum=1):
    if type(value) is not int or not minimum <= value < 2**63:
        raise ValueError
    return value


def _text(value):
    if type(value) is not str or not value or any(ord(c) < 32 or ord(c) == 127 for c in value):
        raise ValueError
    return value


def _terms(value):
    if type(value) is not list or not value or any(re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_]*", _text(v)) is None for v in value):
        raise ValueError
    return value


def _interval(value, strand=False):
    if type(value) is not dict or _integer(value["start"]) > _integer(value["end"]):
        raise ValueError
    if strand and (type(value["strand"]) is not int or value["strand"] not in (-1, 1)):
        raise ValueError


def _bound_source(source, query, raw, kind):
    identifier = "native-" + kind + "-" + _digest(_json([
        SCHEMA, source.source_dataset, source.source_release, source.public_manifest_sha256,
        query.query_id, kind, raw,
    ]))
    return replace(source, source_record_id=identifier, source_fields={
        "annotation_binding": source.source_fields, "query_id": query.query_id,
        "role": query.role, "annotation_scope": query.annotation_scope,
        "query_geometry": {k: v for k, v in query.transport().items() if k != "query_id"},
        "source_geometry": _plain(query.source_geometry), "raw_record": raw,
        "interval_evidence": _plain(query.interval_evidence),
    })


def _format_transcript(raw, source, query, context=False):
    """Reuse typed identifier formatting, but never label JSON/cache evidence CSQ."""
    projection = {}
    for key, target in {
        "gene_id": "Gene", "transcript_id": "Feature", "protein_id": "ENSP", "gene_symbol": "SYMBOL",
        "gene_symbol_source": "SYMBOL_SOURCE", "hgnc_id": "HGNC_ID", "biotype": "BIOTYPE",
        "mane_select": "MANE_SELECT", "mane_plus_clinical": "MANE_PLUS_CLINICAL",
        "appris": "APPRIS", "exon": "EXON", "intron": "INTRON", "hgvsc": "HGVSc",
        "hgvsp": "HGVSp", "hgvsg": "HGVSg", "amino_acids": "Amino_acids",
    }.items():
        value = raw.get(key)
        if value not in (None, "", ".", "-"):
            projection[target] = _text(value)
    for key, target in {"swissprot": "SWISSPROT", "trembl": "TREMBL", "uniparc": "UNIPARC",
                        "uniprot_isoform": "UNIPROT_ISOFORM", "mane": "MANE", "flags": "FLAGS"}.items():
        value = raw.get(key)
        if value is not None:
            if type(value) is not list:
                raise ValueError
            projection[target] = "&".join(_text(v) for v in value)
    if raw.get("canonical") is not None:
        if type(raw["canonical"]) is not int or raw["canonical"] not in (0, 1):
            raise ValueError
        projection["CANONICAL"] = "YES" if raw["canonical"] else ""
    if raw.get("tsl") is not None:
        projection["TSL"] = str(_integer(raw["tsl"]))
    if context:
        attributes = raw["attributes"]
        if type(attributes) is not list or type(raw["raw_fields"]) is not dict:
            raise ValueError
        statuses, mane, mane_values = [], [], {}
        for attribute in attributes:
            code, value = _text(attribute["code"]), attribute["value"]
            if value not in (None, ""):
                _text(value)
            if code in ("MANE_Select", "MANE_Plus_Clinical"):
                mane.append(code)
                mane_values.setdefault(code, set()).add(value)
            elif code in ("TSL", "appris") or code.startswith("cds_"):
                statuses.append(code + ":" + (value or ""))
        projection["MANE"] = "&".join(mane)
        for code, values in mane_values.items():
            if len(values) == 1:
                projection["MANE_SELECT" if code == "MANE_Select" else "MANE_PLUS_CLINICAL"] = next(iter(values)) or ""
        for key in ("swissprot", "trembl", "uniparc", "uniprot_isoform"):
            value = raw["raw_fields"].get("_" + key)
            if value not in (None, "-", ""):
                projection[key.upper()] = _text(value).replace(",", "&")
        # Mapping ranges live only in raw context, not mutation-position fields.
        for key in ("EXON", "INTRON", "HGVSc", "HGVSp", "HGVSg", "Amino_acids"):
            projection.pop(key, None)
    else:
        projection["Consequence"] = "&".join(_terms(raw["consequence_terms"]))
        projection["VARIANT_CLASS"] = SV_CLASSES[query.sv_type]
        for prefix in ("cdna", "cds", "protein"):
            for suffix in ("start", "end"):
                if raw.get(prefix + "_" + suffix) is not None:
                    _integer(raw[prefix + "_" + suffix])
            lo, hi = raw.get(prefix + "_start"), raw.get(prefix + "_end")
            if lo is not None and hi is not None and lo > hi:
                raise ValueError
            if prefix == "protein" and (lo is not None or hi is not None):
                projection["Protein_position"] = str(lo) if lo == hi else f"{lo or '?'}-{hi or '?'}"
        for key in ("bp_overlap", "distance"):
            if raw.get(key) is not None:
                _integer(raw[key], 0)
        if raw.get("percentage_overlap") is not None and (type(raw["percentage_overlap"]) not in (int, float)
                                                          or not 0 <= raw["percentage_overlap"] <= 100):
            raise ValueError
    if raw.get("strand") is not None and (type(raw["strand"]) is not int or raw["strand"] not in (-1, 1)):
        raise ValueError
    bound = _bound_source(source, query, raw, "context" if context else "effect")
    typed = _transcript({k: quote(v, safe="") for k, v in projection.items()}, bound, bound.source_record_id)
    if context and statuses:
        typed = replace(typed, transcript_status=(typed.transcript_status or ()) + tuple(statuses))
    if context and any(len(values) > 1 for values in mane_values.values()):
        typed = replace(typed, mane_refseq=None)
    domains = raw.get("domain_overlaps" if context else "domains")
    if domains is not None:
        if type(domains) is not list or any(type(d) is not dict for d in domains):
            raise ValueError
        if not context:
            for domain in domains:
                # Pinned VEP JSON can coerce a domain name to a number. Keep the
                # original JSON value; do not reconstruct a missing label suffix.
                for key in ("db", "name"):
                    if type(domain[key]) is str:
                        _text(domain[key])
                    elif type(domain[key]) not in (int, float):
                        raise ValueError
        typed = replace(typed, domain_overlaps=tuple(_bound_source(source, query,
            {"transcript_id": raw["transcript_id"], "ordinal": i, "domain": d},
            "reference-domain" if context else "domain") for i, d in enumerate(domains)))
    if not context:
        if raw["impact"] not in ("HIGH", "MODERATE", "LOW", "MODIFIER"):
            raise ValueError
        predictions = [_bound_source(source, query, {"transcript_id": raw["transcript_id"], "impact": raw["impact"]}, "prediction")]
        for tool in ("sift", "polyphen"):
            score, prediction = raw.get(tool + "_score"), raw.get(tool + "_prediction")
            if score is not None and (type(score) not in (int, float) or not 0 <= score <= 1):
                raise ValueError
            if prediction is not None:
                _text(prediction)
            if score is not None or prediction is not None:
                predictions.append(_bound_source(source, query, {"transcript_id": raw["transcript_id"],
                    tool + "_score": score, tool + "_prediction": prediction}, "prediction"))
        typed = replace(typed, predictions=tuple(predictions))
    return typed


def _validate_context(raw, query):
    if type(raw) is not dict or not set(PROTOCOL["context_transcript_fields"]) <= set(raw):
        raise ValueError
    geometry = raw["transcript_geometry"]
    _interval(geometry, strand=True)
    if geometry["contig"] != query.contig or geometry["start"] > query.end or geometry["end"] < query.start:
        raise ValueError
    for kind in ("exons", "introns", "cds"):
        blocks = raw[kind]
        if type(blocks) is not list or type(raw["overlapping_" + kind]) is not list:
            raise ValueError
        for rank, block in enumerate(blocks, 1):
            _interval(block, strand=True)
            if (_integer(block["rank"]) != rank or _integer(block["total"]) != len(blocks)
                    or not geometry["start"] <= block["start"] <= block["end"] <= geometry["end"]
                    or block["strand"] != geometry["strand"]):
                raise ValueError
            if block.get("stable_id") is not None:
                _text(block["stable_id"])
            if block.get("version") not in (None, ""):
                _integer(block["version"])
        overlaps = [b for b in blocks if b["start"] <= query.end and b["end"] >= query.start]
        if _json(raw["overlapping_" + kind]) != _json(overlaps):
            raise ValueError
    maps = raw["mappings"]
    if type(maps) is not dict:
        raise ValueError
    for method in ("genomic2cdna", "genomic2cds", "genomic2pep"):
        if type(maps[method]) is not list:
            raise ValueError
        for mapping in maps[method]:
            _interval(mapping)
            if mapping["class"] not in ("Bio::EnsEMBL::Mapper::Gap", "Bio::EnsEMBL::Mapper::Coordinate"):
                raise ValueError
    length = raw["peptide_length"]
    if length is not None:
        _integer(length, 0)
    domains = raw["protein_features"]
    if domains is not None:
        if type(domains) is not list:
            raise ValueError
        for domain in domains:
            _interval(domain)
            _text(domain["hseqname"])
            if type(domain["analysis"]) is not dict:
                raise ValueError
            for value in domain["analysis"].values():
                if value is not None:
                    _text(value)
    expected = None
    if domains is not None and length is not None:
        pep = [p for p in maps["genomic2pep"] if p["class"] == "Bio::EnsEMBL::Mapper::Coordinate" and 1 <= p["start"] <= length]
        expected = [d for d in domains if any(d["start"] <= min(p["end"], length) and d["end"] >= p["start"] for p in pep)]
    if _json(raw["domain_overlaps"]) != _json(expected):
        raise ValueError


def parse_native_output(content, *, expected, source, context=False, cache_entries=None):
    """Strict JSONL identity accounting. Unknown source keys remain raw."""
    try:
        if source.source_sha256 != _digest(content):
            raise ValueError
        queries = {q.query_id: q for q in expected}
        if len(queries) != len(expected) or any(q.reason is not None for q in expected):
            raise ValueError
        results = {}
        for line in content.splitlines():
            record = _strict_json(line)
            if type(record) is not dict:
                raise ValueError
            key = record["query"]["query_id"] if context else record["id"]
            if type(key) is not str or key not in queries or key in results:
                raise ValueError
            query = queries[key]
            if context:
                if (query.annotation_scope != "transcript_context" or record["schema"] != SCHEMA
                        or _integer(record["cache_bin_size"]) != BIN_SIZE or record["query"] != query.transport()
                        or any(type(record["query"][k]) is not int for k in ("start", "end", "ordinal"))):
                    raise ValueError
                versions = record["cache_versions"]
                if type(versions) is not dict or any(versions.get(k) != v for k, v in {
                        "assembly": "GRCh38.p14", "genebuild": "GENCODE50", "gencode": "GENCODE 50"}.items()):
                    raise ValueError
                files = []
                for b in range((query.start - 1) // BIN_SIZE, (query.end - 1) // BIN_SIZE + 1):
                    path = f"{query.contig}/{b * BIN_SIZE + 1}-{(b + 1) * BIN_SIZE}.gz"
                    digest = cache_entries.get(path)
                    if digest is not None and re.fullmatch(r"[0-9a-f]{64}", _text(digest)) is None:
                        raise ValueError
                    files.append({"bin": b, "path": path, "sha256": digest})
                if record["cache_files"] != files or any(type(f["bin"]) is not int for f in record["cache_files"]):
                    raise ValueError
                rows, non_transcripts = record["transcripts"], []
                if rows and not any(f["sha256"] is not None for f in files):
                    raise ValueError
            else:
                if (query.annotation_scope != "structural_effect" or record["input"] != query.line().rstrip("\n")
                        or record["seq_region_name"] != query.contig or _integer(record["start"]) != query.start
                        or _integer(record["end"]) != query.end or type(record["strand"]) is not int or record["strand"] != 1
                        or record["assembly_name"] != "GRCh38" or record["variant_class"] != SV_CLASSES[query.sv_type]
                        or record["allele_string"] != SV_CLASSES[query.sv_type]):
                    raise ValueError
                rows = record.get("transcript_consequences", [])
                non_transcripts = record.get("intergenic_consequences", [])
                if type(non_transcripts) is not list or (not rows and not non_transcripts) or rows and non_transcripts:
                    raise ValueError
            if type(rows) is not list:
                raise ValueError
            transcripts, seen = [], set()
            for raw in rows:
                if type(raw) is not dict:
                    raise ValueError
                tid = _text(raw["transcript_id"])
                if re.fullmatch(r"ENST[0-9]{11}(?:\.[1-9][0-9]*)?", tid) is None:
                    raise ValueError
                stable = tid.split(".", 1)[0]
                if stable in seen:
                    raise ValueError
                seen.add(stable)
                if context:
                    _validate_context(raw, query)
                elif _integer(raw["allele_num"]) != 1 or raw["variant_allele"] != SV_CLASSES[query.sv_type]:
                    raise ValueError
                transcripts.append(_format_transcript(raw, source, query, context))
            for raw in non_transcripts:
                if (type(raw) is not dict or raw["variant_allele"] != SV_CLASSES[query.sv_type]
                        or _integer(raw["allele_num"]) != 1 or _terms(raw["consequence_terms"]) != ["intergenic_variant"]
                        or raw["impact"] != "MODIFIER"):
                    raise ValueError
            results[key] = (tuple(sorted(transcripts, key=lambda t: t.consequence_id)),
                            tuple(_bound_source(source, query, r, "intergenic") for r in non_transcripts), record)
        if set(results) != set(queries):
            raise ValueError
        return results
    except (ValueError, TypeError, KeyError, AttributeError, IndexError, OverflowError, RecursionError):
        raise VepAnnotationError("invalid or incomplete native annotation output") from None


def _validate_joined_claims(effect, context):
    """Reject incompatible claims before attaching reference context to an effect."""
    try:
        for key, prefix in (("gene_id", "ENSG"), ("protein_id", "ENSP")):
            identifiers = []
            for record in (effect, context):
                value = record.get(key)
                if value in (None, "", ".", "-"):
                    identifiers.append(None)
                    continue
                match = re.fullmatch(r"(" + prefix + r"[0-9]{11})(?:\.([1-9][0-9]*))?", _text(value))
                if match is None:
                    raise ValueError
                identifiers.append(match.groups())
            a, b = identifiers
            if a is not None and b is not None and (a[0] != b[0] or a[1] is not None and b[1] is not None and a[1] != b[1]):
                raise ValueError

        mappings = context["mappings"]["genomic2pep"]
        coordinate = "Bio::EnsEMBL::Mapper::Coordinate"
        lo = mappings[0]["start"] if mappings and mappings[0]["class"] == coordinate else None
        hi = mappings[-1]["end"] if mappings and mappings[-1]["class"] == coordinate else None
        if lo is not None and hi is not None:
            lo, hi = sorted((lo, hi))  # VEP format_coords also orders the bounds.
        length = context["peptide_length"]
        for key, mapped in (("protein_start", lo), ("protein_end", hi)):
            value = effect.get(key)
            if value is not None:
                _integer(value)
                if mapped != value or length is not None and value > length + 1:
                    raise ValueError

        def token(value):
            # The exporter numberifies word tokens, including 00123 and 1e3.
            # Decimal comparisons avoid conflating distinct large integer names.
            if type(value) in (int, float) or type(value) is str and re.fullmatch(r"[0-9]+(?:[eE][0-9]+)?", value):
                result = Decimal(str(value))
                if not result.is_finite():
                    raise ValueError
                return result
            return _text(value)

        features = context["protein_features"]
        if features is not None and effect.get("domains"):
            available = Counter()
            if lo is not None and hi is not None:
                captures = None
                for feature in features:
                    if feature["start"] > hi or feature["end"] < lo:
                        continue
                    # OutputFactory builds the label; JSON.pm takes its first
                    # unanchored word:word match (e.g. PDB-ENSP -> ENSP).
                    label = ":".join(v for v in (feature["analysis"].get("display_label"), feature["hseqname"]) if v)
                    label = re.sub(r"[\s;=]", "_", label)
                    match = re.search(r"(\w+):(\w+)", label, flags=re.ASCII)
                    # JSON.pm does not test the match result: a failed match
                    # reuses the previous Perl captures within this loop.
                    if match is not None:
                        captures = match.groups()
                    if captures is not None and all(v != "0" for v in captures):
                        available[tuple(token(v) for v in captures)] += 1
            reported = Counter((token(d["db"]), token(d["name"])) for d in effect["domains"])
            if reported - available:
                raise ValueError
    except (ValueError, TypeError, KeyError, AttributeError, ArithmeticError):
        raise VepAnnotationError("incompatible native and cache transcript claims") from None


def run_native_annotation(inputs, *, coordinate_maps, reference, release, work, env, executable,
                          timeout, buffer_size, binding, versions, execution_config=None):
    """Run within the main runner's verified private staging directory only."""
    try:
        from genosapiens.vep_annotation import _annotation_child, _annotation_content
        cache_dir = release / "homo_sapiens/116_GRCh38"
        cache_contigs = {p.name for p in cache_dir.iterdir() if p.is_dir()}
        inventory = _strict_json(_annotation_content(release / "cache.inventory.json", execution_config).decode())
        prefix = "homo_sapiens/116_GRCh38/"
        cache_entries = {e["path"][len(prefix):]: e["sha256"] for e in inventory["entries"] if e["path"].startswith(prefix)}
        by_index = {index: native_queries(req, coordinate_maps, reference, cache_contigs) for index, req in inputs}
        all_queries = tuple(q for queries in by_index.values() for q in queries)
        if len({q.query_id for q in all_queries}) != len(all_queries):
            raise ValueError
        valid = tuple(q for q in all_queries if q.reason is None)
        effect_context = {q.query_id: replace(q,
            query_id="native-query-" + _digest(_json([q.query_id, "reference-context"])),
            annotation_scope="transcript_context", sv_type=None)
            for q in valid if q.annotation_scope == "structural_effect"}
        execution_queries = valid + tuple(effect_context.values())
        artifacts, proof, results = [], {"configuration": configuration(), "versions": versions, "artifacts": {}}, {}
        if valid:
            _verify_sources()
        for context in (False, True):
            queries = tuple(sorted((q for q in execution_queries if (q.annotation_scope == "transcript_context") == context),
                                   key=lambda q: (q.contig, q.start, q.end, q.query_id)))
            if not queries:
                continue
            name = "context" if context else "native"
            input_path, output_path = work / (name + ".input.jsonl" if context else name + ".input.tsv"), work / (name + ".output.jsonl")
            text = "".join(_json(q.transport()) if context else q.line() for q in queries)
            input_path.write_text(text, encoding="ascii")
            os.chmod(input_path, 0o600)
            artifacts.append(input_path)
            input_hash = _file_identity(input_path)[0]
            helper = None
            if context:
                helper = work / "native-cache-context.pl"
                helper.write_text(CACHE_HELPER, encoding="ascii")
                os.chmod(helper, 0o600)
                artifacts.append(helper)
                command = ["/usr/bin/perl", str(helper), str(release), str(reference.fasta_path), str(input_path), str(output_path)]
            else:
                command = [str(executable), *NATIVE_FLAGS, "--buffer_size", str(buffer_size), "--dir", str(work),
                           "--dir_plugins", str(work), "--dir_cache", str(release), "--fasta", str(reference.fasta_path),
                           "--input_file", str(input_path), "--output_file", str(output_path), "--warning_file", str(work / "native.warnings")]
            completed = _annotation_child(command, work=work, env=env, timeout=timeout,
                name=name, execution_config=execution_config,
                output_files=(output_path, work / (name + ".warnings")))
            for suffix, data in (("stdout", completed.stdout), ("stderr", completed.stderr)):
                path = work / (name + "." + suffix)
                path.write_bytes(data)
                os.chmod(path, 0o600)
                artifacts.append(path)
            for path in (output_path, work / (name + ".warnings")):
                if path.exists():
                    os.chmod(path, 0o600)
                    artifacts.append(path)
            if (_file_identity(input_path)[0] != input_hash
                    or helper is not None and _file_identity(helper)[0] != _digest(CACHE_HELPER)):
                raise ValueError
            for path in artifacts:
                digest = _file_identity(path)[0]
                if proof["artifacts"].setdefault(path.name, digest) != digest:
                    raise ValueError
            content = _annotation_content(output_path, execution_config).decode("utf-8")
            if _digest(content) != proof["artifacts"][output_path.name]:
                raise ValueError
            source = SourcedRecord("Ensembl VEP cache context" if context else "Ensembl VEP native", "116",
                name + "-output", _digest(content), binding["bundle_manifest_sha256"], {
                    "schema": SCHEMA, "configuration_sha256": binding["configuration_sha256"],
                    "runtime_image": binding["runtime_image"], "upstream_image_digest": binding["upstream_image_digest"],
                    "executable_sha256": binding["executable_sha256"], "contract_sha256": binding["contract_sha256"],
                    "cache_inventory_sha256": binding["cache_inventory_sha256"],
                    "fasta_sha256": binding["fasta_sha256"], "fai_sha256": binding["fai_sha256"],
                    "resource_reference_identity_id": binding["resource_reference_identity_id"],
                    "input_evidence_sha256": binding["input_evidence_sha256"],
                    "tool_version": binding["tool_version"], "source_revisions": dict(SOURCE_REVISIONS),
                    "source_sha256": dict(SOURCE_HASHES), "versions": versions,
                    "helper_sha256": _digest(CACHE_HELPER), "schema_sha256": _digest(_json(PROTOCOL)),
                    "transport_sha256": input_hash, "artifacts": dict(proof["artifacts"]),
                })
            parsed = parse_native_output(content, expected=queries, source=source, context=context, cache_entries=cache_entries)
            results.update({key: (*value, source) for key, value in parsed.items()})
        if valid:
            _verify_sources()
        # Also detect a later subprocess modifying an earlier transport/output.
        if any(_file_identity(p)[0] != proof["artifacts"][p.name] for p in artifacts):
            raise ValueError
        outcomes = {}
        for index, request in inputs:
            records, transcripts, non_transcripts = [], [], []
            for query in by_index[index]:
                reference_context = None
                if query.reason:
                    source = SourcedRecord("Ensembl VEP native assessment", "116", query.query_id,
                        _digest(_json([query.transport(), _plain(query.source_geometry), query.reason])),
                        binding["bundle_manifest_sha256"], {"schema": SCHEMA, "versions": versions,
                            "configuration_sha256": binding["configuration_sha256"],
                            "cache_inventory_sha256": binding["cache_inventory_sha256"]})
                    raw, count = None, None
                else:
                    tx, other, raw, source = results[query.query_id]
                    if query.query_id in effect_context:
                        cq = effect_context[query.query_id]
                        ctx, _, context_raw, context_source = results[cq.query_id]
                        by_transcript = {t.ensembl_transcript_id: t for t in ctx if t.ensembl_transcript_id is not None}
                        enriched, matched = [], set()
                        for effect in tx:
                            cached = by_transcript.get(effect.ensembl_transcript_id)
                            if cached is None:
                                enriched.append(effect)
                                continue
                            _validate_joined_claims(effect.source.source_fields["raw_record"],
                                                   cached.source.source_fields["raw_record"])
                            matched.add(cached.ensembl_transcript_id)
                            identifier = "native-enriched-" + _digest(_json([effect.consequence_id, cached.consequence_id]))
                            merged_source = replace(effect.source, source_record_id=identifier, source_fields={
                                **effect.source.source_fields, "effect_source_record_id": effect.source.source_record_id,
                                "reference_context": _plain(cached.source),
                            })
                            domains = effect.domain_overlaps
                            if cached.domain_overlaps is not None:
                                # The original JSON domains remain in raw_record.
                                # Use the complete reference features, not two
                                # counts of the same cache-derived observation.
                                domains = cached.domain_overlaps
                            enriched.append(replace(effect, consequence_id=identifier, source=merged_source, domain_overlaps=domains))
                        tx = tuple(enriched)
                        context_record = _bound_source(context_source, cq, context_raw, "geometry-context")
                        reference_context = _plain(replace(context_record, source_fields={
                            **context_record.source_fields, "transcript_count": len(ctx),
                            "joined_transcripts": sorted(matched),
                            "unjoined_context_transcripts": [t.source.source_fields["raw_record"]["transcript_id"]
                                for t in ctx if t.ensembl_transcript_id not in matched],
                        }))
                    transcripts.extend(tx)
                    non_transcripts.extend(other)
                    count = len(tx)
                fields = {"role": query.role, "annotation_scope": query.annotation_scope,
                    "status": "unassessed" if query.reason else "completed", "reason": query.reason,
                    "query_id": query.query_id, "query_geometry": {k: v for k, v in query.transport().items() if k != "query_id"},
                    "source_geometry": _plain(query.source_geometry), "transcript_count": count,
                    "lookup_result": None if query.reason else "matched" if count else "completed-no-overlap",
                    "raw_result": raw, "annotation_binding": source.source_fields,
                    "interval_evidence": _plain(query.interval_evidence), "reference_context": reference_context,
                    "native_identity_id": request.identity_id, "allele_id": request.allele.allele_id}
                records.append(replace(source, source_record_id=query.query_id, source_fields=fields))
            completed_count = sum(q.reason is None for q in by_index[index])
            status = "unassessed" if not completed_count else "partially_annotated" if completed_count < len(records) else (
                "annotated" if any(q.annotation_scope == "structural_effect" for q in by_index[index]) else "annotated_context")
            reasons = {q.reason for q in by_index[index] if q.reason}
            reason = None if not reasons else next(iter(reasons)) if len(reasons) == 1 else "multiple-geometry-assessment-gaps"
            outcomes[index] = (status, reason, tuple(transcripts), tuple(non_transcripts), tuple(records))
        return outcomes, proof, tuple(artifacts)
    except VepAnnotationError:
        raise
    except Exception:
        raise VepAnnotationError("offline native VEP annotation failed") from None
