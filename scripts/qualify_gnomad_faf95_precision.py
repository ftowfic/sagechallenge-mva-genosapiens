#!/usr/bin/env python3
"""Qualify the focused gnomAD FAF95 float boundary through Hail and BCF."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

from genosapiens.gnomad_faf95_adapter import (
    authenticate_gnomad_faf95_precision_contract,
    load_gnomad_faf95_precision_contract,
    qualification_report,
    verify_gnomad_faf95_precision_qualification,
)
from genosapiens.gnomad_projection import (
    gnomad_vcf_metadata, load_gnomad_projection_contract,
    project_hail_projection_to_vcf,
)


DESCRIPTOR = "config/public-resources/gnomad-joint-faf95-precision-v1.json"
RUNTIME_IMAGE = (
    "genosapiens/reference-builder:disease-gene-v1@sha256:"
    "8813259c447025909c0982293b20eeb1493a0f9b60724d507c9081d904a796ac"
)


def _run(arguments: tuple[str, ...]) -> str:
    completed = subprocess.run(arguments, check=True, capture_output=True, text=True)
    if completed.stderr:
        raise RuntimeError("focused gnomAD precision command emitted diagnostics")
    return completed.stdout


def _dataset(hl, ancestries, maximum):
    by_ancestry = {}
    for index, ancestry in enumerate(ancestries):
        by_ancestry[ancestry] = hl.struct(
            ac=index + 1, an=(index + 1) * 1000,
            homozygote_count=0, hemizygote_count=0,
            faf95=maximum / 2, faf99=maximum / 3,
        )
    return hl.struct(
        adjusted=hl.struct(ac=10, an=6000, homozygote_count=0, hemizygote_count=0),
        ancestries=hl.struct(**by_ancestry),
        faf=hl.struct(faf95=maximum / 2, faf99=maximum / 3),
        grpmax=hl.struct(AC=500, AF=1 / 6, AN=3000,
                         homozygote_count=0, gen_anc="nfe"),
        fafmax=hl.struct(faf95_max=maximum, faf95_max_gen_anc="nfe",
                         faf99_max=maximum / 2, faf99_max_gen_anc="nfe"),
        filters=hl.empty_set(hl.tstr),
    )


def _build_boundary_rows(contract, projection, temporary: Path):
    import hail as hl

    hl.init(master="local[2]", tmp_dir=str(temporary), local_tmpdir=str(temporary),
            log=str(temporary / "hail.log"), quiet=True)
    try:
        table = hl.utils.range_table(len(contract.boundary_values))
        values = hl.literal([float(item.value) for item in contract.boundary_values])
        maximum = values[table.idx]
        datasets = {
            dataset.name: _dataset(hl, contract.ancestry_order, maximum)
            for dataset in projection.datasets
        }
        # The exact matching ancestry is nfe, index 4 in the authenticated order.
        nfe = datasets["joint"].ancestries.nfe
        datasets["joint"] = datasets["joint"].annotate(
            ancestries=datasets["joint"].ancestries.annotate(
                nfe=nfe.annotate(ac=contract.boundary_ac, an=contract.boundary_an)
            )
        )
        table = table.annotate(
            locus=hl.locus("chr1", table.idx + 1, reference_genome="GRCh38"),
            alleles=hl.array(["A", "C"]), **datasets,
            region_flags=hl.struct(**{flag: False for flag in projection.region_flags}),
            discordance=hl.struct(
                overall=hl.struct(p_value=1.0, odds_ratio=1.0),
                ancestries=hl.struct(**{
                    ancestry: hl.struct(p_value=1.0, odds_ratio=1.0)
                    for ancestry in contract.ancestry_order
                }),
                cochran_mantel_haenszel_test=hl.struct(p_value=1.0, chisq=0.0),
                stat_union=hl.struct(p_value=1.0, stat_test_name="synthetic-boundary",
                                     gen_ancs=hl.array(["nfe"])),
                is_discrepant=False,
            ),
        ).key_by("locus", "alleles")
        projected = project_hail_projection_to_vcf(table, hl, projection)
        vcf = temporary / "boundary.vcf.bgz"
        hl.export_vcf(projected, str(vcf), metadata=gnomad_vcf_metadata(projection), tabix=True)
    finally:
        hl.stop()

    bcf = temporary / "boundary.bcf"
    _run(("bcftools", "view", "--no-version", "-Ob", "-o", str(bcf), str(vcf)))
    _run(("bcftools", "index", "--csi", str(bcf)))
    header = _run(("bcftools", "view", "--header-only", "--no-version", str(bcf)))
    expected_header = (
        '##INFO=<ID=JT_FAFMAX95,Number=1,Type=Float,'
        'Description="gnomAD joint maximum ancestry 95% filtering allele frequency">'
    )
    if expected_header not in header:
        raise RuntimeError("focused gnomAD precision header mismatch")
    output = _run((
        "bcftools", "query", "--format",
        "%POS\\t%FILTER\\t%INFO/JT_FAFMAX95\\t%INFO/JT_FAFMAX95_ANC\\t"
        "%INFO/JT_AC_ANC\\t%INFO/JT_AN_ANC\\n", str(bcf),
    ))
    rows = []
    for expected, line in zip(contract.boundary_values, output.splitlines(), strict=True):
        position, filter_value, token, ancestry, acs, ans = line.split("\t")
        if int(position) != len(rows) + 1:
            raise RuntimeError("focused gnomAD precision row order mismatch")
        index = contract.ancestry_order.index(ancestry)
        ac_values, an_values = acs.split(","), ans.split(",")
        rows.append({
            "label": expected.label, "source_value": expected.value,
            "bcf_token": token, "filter": filter_value, "ancestry": ancestry,
            "ancestry_ac": ac_values[index], "ancestry_an": an_values[index],
        })
    return tuple(rows), hl.__version__, _run(("bcftools", "--version-only")).strip()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repository-root", default=".")
    parser.add_argument("--reference-root", default="/mnt/reference")
    arguments = parser.parse_args()
    repository = Path(arguments.repository_root).resolve(strict=True)
    references = Path(arguments.reference_root).resolve(strict=True)
    descriptor_path = repository / DESCRIPTOR
    descriptor_content = descriptor_path.read_bytes()
    descriptor_sha256 = hashlib.sha256(descriptor_content).hexdigest()
    contract = load_gnomad_faf95_precision_contract(
        descriptor_content, descriptor_sha256
    )
    authenticate_gnomad_faf95_precision_contract(
        contract, repository_root=repository, reference_root=references
    )
    projection = load_gnomad_projection_contract(
        repository / contract.projection_contract_path
    )
    with tempfile.TemporaryDirectory(prefix="gnomad-faf95-precision-") as directory:
        rows, hail_version, bcftools_version = _build_boundary_rows(
            contract, projection, Path(directory)
        )
    content = qualification_report(
        contract, boundary_rows=rows, runtime_image=RUNTIME_IMAGE,
        hail_version=hail_version, bcftools_version=bcftools_version,
    )
    report_sha256 = hashlib.sha256(content).hexdigest()
    verify_gnomad_faf95_precision_qualification(content, report_sha256, contract)
    relative = Path("manifests/gnomad/4.1/faf95-precision") / f"{report_sha256}.json"
    output = references / relative
    output.parent.mkdir(parents=True, exist_ok=True)
    try:
        with output.open("xb") as destination:
            destination.write(content)
    except FileExistsError:
        if output.read_bytes() != content:
            raise RuntimeError("focused gnomAD precision report path collision") from None
    print(json.dumps({
        "descriptor_sha256": descriptor_sha256,
        "report_path": relative.as_posix(), "report_sha256": report_sha256,
        "boundary_rows": rows,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
