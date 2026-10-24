"""Command-line entry points for containerized Genosapiens workflows."""

from __future__ import annotations

import argparse
from dataclasses import asdict
from collections.abc import Sequence
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import sys

from genosapiens import __version__
from genosapiens.challenge_fixture import (
    ChallengeFixtureError,
    bootstrap_challenge_fixture,
    challenge_fixture_result_as_json,
)
from genosapiens.clinvar import (
    ClinvarError,
    bootstrap_clinvar,
    clinvar_result_as_json,
)
from genosapiens.config import (
    ConfigError,
    ReferenceBuilderConfig,
    Track1PipelineConfig,
    config_as_json,
    load_config,
)
from genosapiens.download import DownloadError
from genosapiens.hpo_ic import HpoICError, publish_hpo_ic
from genosapiens.vcf_header import VcfHeaderError, sanitize_header_file
from genosapiens.deepvariant import (
    DeepVariantError,
    bootstrap_deepvariant_hg002,
    deepvariant_result_as_json,
)
from genosapiens.grch38_reference import (
    Grch38BootstrapError,
    bootstrap_grch38_reference,
    bootstrap_result_as_json,
)
from genosapiens.giab import (
    GiabError,
    bootstrap_giab_hg002,
    giab_result_as_json,
)
from genosapiens.gnomad_hail import (
    GnomadHailError,
    hail_report_as_json,
    run_bounded_public_scan,
    run_hail_smoke,
)
from genosapiens.gnomad_all_sites import (
    GnomadAllSitesError,
    bootstrap_gnomad_all_sites,
    gnomad_all_sites_result_as_json,
)
from genosapiens.gnomad_projection import GnomadProjectionError
from genosapiens.gnomad_pilot import (
    pilot_report_as_json,
    run_gnomad_pilot,
)
from genosapiens.gnomad_publication import (
    full_publication_result_as_json,
    publish_gnomad_joint_frequency,
)
from genosapiens.gnomad_validation import (
    validate_gnomad_pilot,
    validation_report_as_json,
)
from genosapiens.manifest import ManifestError
from genosapiens.ontology_bundle import (
    OntologyBundleError,
    bootstrap_ontology_bundle,
    ontology_result_as_json,
)
from genosapiens.open_targets import (
    OpenTargetsError,
    bootstrap_open_targets,
    open_targets_result_as_json,
)
from genosapiens.rare_disease import (
    RareDiseaseError,
    bootstrap_rare_disease_evidence,
    rare_disease_result_as_json,
)
from genosapiens.reference_concordance import (
    ReferenceConcordanceError,
    document_installed_contig_order_adjudication,
    evaluate_installed_reference_concordance,
    require_rescue_reference_concordance,
    rescue_authorization_as_json,
    write_concordance_report,
    write_contig_order_adjudication,
)
from genosapiens.storage import (
    CapacityPlan,
    InsufficientCapacity,
    StorageError,
    gb_to_bytes,
    initialize_layout,
    report_as_json,
    require_capacity,
)
from genosapiens.vep_bundle import (
    VepBundleError,
    bootstrap_vep_bundle,
    vep_result_as_json,
)
from genosapiens.variant_class_resources import (
    VariantClassResourceError,
    bootstrap_variant_class_resources,
    variant_class_result_as_json,
)


CHALLENGE_REFERENCE_IDENTITY_PATH = (
    Path(__file__).resolve().parents[2]
    / "config/reference-identities/challenge-grch38-masked-v2-hs38d1.json"
)
GIAB_HG002_CONTRACT_PATH = (
    Path(__file__).resolve().parents[2]
    / "config/public-resources/giab-hg002-v4.2.1.json"
)
DEEPVARIANT_HG002_CONTRACT_PATH = (
    Path(__file__).resolve().parents[2]
    / "config/public-resources/deepvariant-hg002-novaseq-35x-v1.0.0.json"
)
CHALLENGE_FIXTURE_CONTRACT_PATH = (
    Path(__file__).resolve().parents[2]
    / "config/public-resources/challenge-shaped-fixture.json"
)
ONTOLOGY_BUNDLE_CONTRACT_PATH = (
    Path(__file__).resolve().parents[2]
    / "config/public-resources/core-ontologies.json"
)
VEP_BUNDLE_CONTRACT_PATH = (
    Path(__file__).resolve().parents[2]
    / "config/public-resources/ensembl-vep-116-grch38.json"
)
CLINVAR_CONTRACT_PATH = (
    Path(__file__).resolve().parents[2]
    / "config/public-resources/clinvar-grch38-2026-08-22.json"
)
OPEN_TARGETS_CONTRACT_PATH = (
    Path(__file__).resolve().parents[2]
    / "config/public-resources/open-targets-platform-26.06.json"
)
RARE_DISEASE_CONTRACT_PATH = (
    Path(__file__).resolve().parents[2]
    / "config/public-resources/rare-disease-evidence-2026-08-30.json"
)
VARIANT_CLASS_RESOURCES_CONTRACT_PATH = (
    Path(__file__).resolve().parents[2]
    / "config/public-resources/variant-class-resources.json"
)
GNOMAD_JOINT_PROJECTION_CONTRACT_PATH = (
    Path(__file__).resolve().parents[2]
    / "config/public-resources/gnomad-joint-projection.json"
)
GNOMAD_PILOT_VALIDATION_CONTRACT_PATH = (
    Path(__file__).resolve().parents[2]
    / "config/public-resources/gnomad-chr21-chr22-validation.json"
)
GNOMAD_JOINT_PUBLICATION_CONTRACT_PATH = (
    Path(__file__).resolve().parents[2]
    / "config/public-resources/gnomad-joint-publication.json"
)
GNOMAD_ALL_SITES_CONTRACT_PATH = (
    Path(__file__).resolve().parents[2]
    / "config/public-resources/gnomad-all-sites-v4.1.json"
)


def _role_parser(prog: str, role: str) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=prog)
    parser.add_argument("--version", action="version", version=__version__)
    commands = parser.add_subparsers(dest="command", required=True)
    if role == 'track1-pipeline':
        commands.add_parser('run', help='run the complete focused pipeline with VCF, HPO and output paths')
    validate = commands.add_parser("validate-config", help=f"validate {role} JSON")
    validate.add_argument("path")
    show = commands.add_parser("show-config", help=f"validate and print {role} JSON")
    show.add_argument("path")
    if role == "reference-builder":
        initialize = commands.add_parser(
            "init-layout", help="initialize the reference filesystem layout"
        )
        initialize.add_argument("path", help="reference-builder configuration JSON")
        capacity = commands.add_parser(
            "preflight-capacity", help="measure capacity and enforce the reserve"
        )
        capacity.add_argument("path", help="reference-builder configuration JSON")
        capacity.add_argument("--operation", required=True)
        capacity.add_argument("--source-gb", default="0")
        capacity.add_argument("--output-gb", default="0")
        capacity.add_argument("--scratch-gb", default="0")
        bootstrap = commands.add_parser(
            "bootstrap-grch38",
            help="build and validate the pinned GRCh38 plus hs38d1 bundle",
        )
        bootstrap.add_argument("path", help="reference-builder configuration JSON")
        bootstrap.add_argument("--image-digest", required=True)
        bootstrap.add_argument("--retrieved-on")
        giab = commands.add_parser(
            "bootstrap-giab-hg002",
            help="install and validate GIAB HG002 NIST v4.2.1",
        )
        giab.add_argument("path", help="reference-builder configuration JSON")
        giab.add_argument("--image-digest", required=True)
        giab.add_argument("--retrieved-on")
        deepvariant = commands.add_parser(
            "bootstrap-deepvariant-hg002",
            help="install and validate the public HG002 NovaSeq 35x caller VCF",
        )
        deepvariant.add_argument("path", help="reference-builder configuration JSON")
        deepvariant.add_argument("--image-digest", required=True)
        deepvariant.add_argument("--retrieved-on")
        fixture = commands.add_parser(
            "bootstrap-challenge-fixture",
            help="build the challenge-shaped fixture from public HG002 alleles",
        )
        fixture.add_argument("path", help="reference-builder configuration JSON")
        fixture.add_argument("--image-digest", required=True)
        fixture.add_argument("--generated-on")
        ontologies = commands.add_parser(
            "bootstrap-ontologies",
            help="install and validate the pinned core ontology references",
        )
        ontologies.add_argument("path", help="reference-builder configuration JSON")
        ontologies.add_argument("--image-digest", required=True)
        ontologies.add_argument("--retrieved-on")
        hpo_ic = commands.add_parser(
            "build-hpo-ic", help="build immutable disease-count IC from the pinned public HPO corpus",
        )
        hpo_ic.add_argument("path", help="reference-builder configuration JSON")
        hpo_ic.add_argument("--image-digest", required=True)
        hpo_ic.add_argument("--retrieved-on")
        vep = commands.add_parser(
            "bootstrap-vep",
            help="install and validate the pinned offline Ensembl VEP bundle",
        )
        vep.add_argument("path", help="reference-builder configuration JSON")
        vep.add_argument("--image-digest", required=True)
        vep.add_argument("--retrieved-on")
        clinvar = commands.add_parser(
            "bootstrap-clinvar",
            help="install and validate the dated ClinVar GRCh38 VCF",
        )
        clinvar.add_argument("path", help="reference-builder configuration JSON")
        clinvar.add_argument("--image-digest", required=True)
        clinvar.add_argument("--retrieved-on")
        clinvar.add_argument("--manifest-path", help="fresh immutable manifest under manifests/clinvar")
        open_targets = commands.add_parser(
            "bootstrap-open-targets",
            help="install and validate the pinned Open Targets Platform subset",
        )
        open_targets.add_argument("path", help="reference-builder configuration JSON")
        open_targets.add_argument("--image-digest", required=True)
        open_targets.add_argument("--retrieved-on")
        rare_disease = commands.add_parser(
            "bootstrap-rare-disease-evidence",
            help="install direct Orphadata/GenCC and build cross-source projections",
        )
        rare_disease.add_argument(
            "path", help="reference-builder configuration JSON"
        )
        rare_disease.add_argument("--image-digest", required=True)
        rare_disease.add_argument("--retrieved-on")
        class_resources = commands.add_parser(
            "bootstrap-variant-class-resources",
            help="install population and annotation resources for non-small variants",
        )
        class_resources.add_argument(
            "path", help="reference-builder configuration JSON"
        )
        class_resources.add_argument("--image-digest", required=True)
        class_resources.add_argument("--retrieved-on")
        hail_smoke = commands.add_parser(
            "hail-smoke",
            help="run a bounded synthetic Hail job with confined local scratch",
        )
        hail_smoke.add_argument("path", help="reference-builder configuration JSON")
        hail_smoke.add_argument("--run-id", default="runtime-smoke")
        hail_smoke.add_argument("--workers", type=int, default=2)
        hail_smoke.add_argument("--rows", type=int, default=32)
        gnomad_scan = commands.add_parser(
            "scan-gnomad-bounded",
            help="stream a bounded interval from the pinned public gnomAD Hail table",
        )
        gnomad_scan.add_argument("path", help="reference-builder configuration JSON")
        gnomad_scan.add_argument("--run-id", default="bounded-public-scan")
        gnomad_scan.add_argument("--workers", type=int, default=2)
        gnomad_pilot = commands.add_parser(
            "run-gnomad-pilot",
            help="materialize and measure the pinned chr21/chr22 gnomAD pilot",
        )
        gnomad_pilot.add_argument(
            "path", help="reference-builder configuration JSON"
        )
        gnomad_pilot.add_argument("--image-digest", required=True)
        gnomad_pilot.add_argument("--run-id", default="chromosome-21-22-pilot")
        gnomad_pilot.add_argument("--workers", type=int, default=16)
        gnomad_validation = commands.add_parser(
            "validate-gnomad-pilot",
            help="compare the staged chr21/chr22 pilot with fixed official sentinels",
        )
        gnomad_validation.add_argument(
            "path", help="reference-builder configuration JSON"
        )
        gnomad_validation.add_argument("--image-digest", required=True)
        gnomad_validation.add_argument("--workers", type=int, default=2)
        gnomad_publication = commands.add_parser(
            "publish-gnomad-joint",
            help=(
                "publish the fixed canonical reduced gnomAD BCF resource after "
                "accepted chr21/chr22 validation"
            ),
        )
        gnomad_publication.add_argument(
            "path", help="reference-builder configuration JSON"
        )
        gnomad_publication.add_argument("--image-digest", required=True)
        gnomad_publication.add_argument("--retrieved-on")
        gnomad_all_sites = commands.add_parser(
            "bootstrap-gnomad-all-sites",
            help="download, index, and publish the fixed gnomAD v4.1 all-sites AN lookups",
        )
        gnomad_all_sites.add_argument(
            "path", help="reference-builder configuration JSON"
        )
        gnomad_all_sites.add_argument("--image-digest", required=True)
        gnomad_all_sites.add_argument("--retrieved-on")
    if role == "track1-pipeline":
        sanitize = commands.add_parser(
            "sanitize-header",
            help="sanitize an entire header-only input into a new controlled no-LLM artifact",
        )
        sanitize.add_argument("path", help="Track 1 pipeline configuration JSON")
        sanitize.add_argument("--input", required=True)
        sanitize.add_argument("--output", required=True)
        validate_reference = commands.add_parser(
            "validate-reference-concordance",
            help="stream and account for every VCF REF against the installed reference",
        )
        validate_reference.add_argument(
            "path", help="Track 1 pipeline configuration JSON"
        )
        validate_reference.add_argument("--input", required=True)
        validate_reference.add_argument("--revalidated-input")
        validate_reference.add_argument("--reference-manifest", required=True)
        validate_reference.add_argument("--output", required=True)
        validate_reference.add_argument("--order-adjudication")
        validate_reference.add_argument("--source-report")
        validate_reference.add_argument("--source-reference-manifest")
        validate_reference.add_argument("--source-reference-root")
        adjudicate_order = commands.add_parser(
            "adjudicate-contig-order",
            help="bind an order-only discrepancy to a non-sequence-changing correction",
        )
        adjudicate_order.add_argument(
            "path", help="Track 1 pipeline configuration JSON"
        )
        adjudicate_order.add_argument("--source-report", required=True)
        adjudicate_order.add_argument("--source-input", required=True)
        adjudicate_order.add_argument("--revalidated-input", required=True)
        adjudicate_order.add_argument("--reference-manifest", required=True)
        adjudicate_order.add_argument("--source-reference-manifest")
        adjudicate_order.add_argument("--source-reference-root")
        adjudicate_order.add_argument(
            "--resolution",
            choices=(
                "vcf-header-reordered-without-record-change",
                "reconstructed-reference-order-corrected",
            ),
            required=True,
        )
        adjudicate_order.add_argument("--output", required=True)
        authorize = commands.add_parser(
            "authorize-rescue",
            help="require a current accepted complete reference concordance report",
        )
        authorize.add_argument("path", help="Track 1 pipeline configuration JSON")
        authorize.add_argument("--concordance-report", required=True)
        authorize.add_argument("--input", required=True)
        authorize.add_argument("--revalidated-input")
        authorize.add_argument("--reference-manifest", required=True)
        authorize.add_argument("--order-adjudication")
        authorize.add_argument("--source-report")
        authorize.add_argument("--source-reference-manifest")
        authorize.add_argument("--source-reference-root")
        authorize.add_argument("--operator-override", action="store_true")
    return parser


def _require_distinct_output(output: Path, *inputs: str | Path | None) -> None:
    for value in inputs:
        if value is None:
            continue
        input_path = Path(value).resolve()
        same_artifact = output == input_path
        if not same_artifact:
            try:
                same_artifact = output.samefile(input_path)
            except OSError:
                pass
        if same_artifact:
            raise ReferenceConcordanceError(
                "controlled report output must not replace an input artifact"
            )




def _run_role(
    parser: argparse.ArgumentParser,
    config_type: (
        type[ReferenceBuilderConfig]
        | type[Track1PipelineConfig]
    ),
    argv: Sequence[str] | None,
) -> int:
    args = parser.parse_args(argv)
    try:
        config = load_config(args.path, config_type)
        if args.command == "init-layout":
            assert isinstance(config, ReferenceBuilderConfig)
            paths = initialize_layout(config.reference_root)
            print(f"initialized: {config.reference_root} ({len(paths)} directories)")
            return 0
        if args.command == "preflight-capacity":
            assert isinstance(config, ReferenceBuilderConfig)
            report = require_capacity(
                config.reference_root,
                CapacityPlan(
                    operation=args.operation,
                    source_bytes=gb_to_bytes(args.source_gb),
                    output_bytes=gb_to_bytes(args.output_gb),
                    scratch_bytes=gb_to_bytes(args.scratch_gb),
                    reserve_bytes=config.minimum_free_gb * 1_000_000_000,
                ),
            )
            print(report_as_json(report))
            return 0
        if args.command == "bootstrap-grch38":
            assert isinstance(config, ReferenceBuilderConfig)
            result = bootstrap_grch38_reference(
                identity_path=CHALLENGE_REFERENCE_IDENTITY_PATH,
                reference_root=config.reference_root,
                image_digest=args.image_digest,
                retrieved_on=args.retrieved_on,
            )
            print(bootstrap_result_as_json(result), end="")
            return 0
        if args.command == "bootstrap-giab-hg002":
            assert isinstance(config, ReferenceBuilderConfig)
            result = bootstrap_giab_hg002(
                contract_path=GIAB_HG002_CONTRACT_PATH,
                reference_identity_path=CHALLENGE_REFERENCE_IDENTITY_PATH,
                reference_root=config.reference_root,
                image_digest=args.image_digest,
                retrieved_on=args.retrieved_on,
            )
            print(giab_result_as_json(result), end="")
            return 0
        if args.command == "bootstrap-deepvariant-hg002":
            assert isinstance(config, ReferenceBuilderConfig)
            result = bootstrap_deepvariant_hg002(
                contract_path=DEEPVARIANT_HG002_CONTRACT_PATH,
                reference_identity_path=CHALLENGE_REFERENCE_IDENTITY_PATH,
                reference_root=config.reference_root,
                image_digest=args.image_digest,
                retrieved_on=args.retrieved_on,
            )
            print(deepvariant_result_as_json(result), end="")
            return 0
        if args.command == "bootstrap-challenge-fixture":
            assert isinstance(config, ReferenceBuilderConfig)
            result = bootstrap_challenge_fixture(
                contract_path=CHALLENGE_FIXTURE_CONTRACT_PATH,
                giab_contract_path=GIAB_HG002_CONTRACT_PATH,
                reference_identity_path=CHALLENGE_REFERENCE_IDENTITY_PATH,
                reference_root=config.reference_root,
                image_digest=args.image_digest,
                generated_on=args.generated_on,
            )
            print(challenge_fixture_result_as_json(result), end="")
            return 0
        if args.command == "bootstrap-ontologies":
            assert isinstance(config, ReferenceBuilderConfig)
            result = bootstrap_ontology_bundle(
                contract_path=ONTOLOGY_BUNDLE_CONTRACT_PATH,
                reference_root=config.reference_root,
                image_digest=args.image_digest,
                retrieved_on=args.retrieved_on,
                scratch_root=config.scratch_root,
                reserve_bytes=config.minimum_free_gb * 1_000_000_000,
            )
            print(ontology_result_as_json(result), end="")
            return 0
        if args.command == "bootstrap-vep":
            assert isinstance(config, ReferenceBuilderConfig)
            result = bootstrap_vep_bundle(
                contract_path=VEP_BUNDLE_CONTRACT_PATH,
                reference_root=config.reference_root,
                image_digest=args.image_digest,
                retrieved_on=args.retrieved_on,
                scratch_root=config.scratch_root,
                reserve_bytes=config.minimum_free_gb * 1_000_000_000,
            )
            print(vep_result_as_json(result), end="")
            return 0
        if args.command == "bootstrap-clinvar":
            assert isinstance(config, ReferenceBuilderConfig)
            result = bootstrap_clinvar(
                contract_path=CLINVAR_CONTRACT_PATH,
                reference_identity_path=CHALLENGE_REFERENCE_IDENTITY_PATH,
                reference_root=config.reference_root,
                image_digest=args.image_digest,
                manifest_path=args.manifest_path,
                retrieved_on=args.retrieved_on,
                scratch_root=config.scratch_root,
                reserve_bytes=config.minimum_free_gb * 1_000_000_000,
            )
            print(clinvar_result_as_json(result), end="")
            return 0
        if args.command == "build-hpo-ic":
            assert isinstance(config, ReferenceBuilderConfig)
            result = publish_hpo_ic(
                contract_path=ONTOLOGY_BUNDLE_CONTRACT_PATH,
                reference_root=config.reference_root, image_digest=args.image_digest,
                retrieved_on=args.retrieved_on,
                reserve_bytes=config.minimum_free_gb * 1_000_000_000,
            )
            print(json.dumps(asdict(result), indent=2, sort_keys=True))
            return 0
        if args.command == "bootstrap-open-targets":
            assert isinstance(config, ReferenceBuilderConfig)
            result = bootstrap_open_targets(
                contract_path=OPEN_TARGETS_CONTRACT_PATH,
                reference_root=config.reference_root,
                image_digest=args.image_digest,
                retrieved_on=args.retrieved_on,
                scratch_root=config.scratch_root,
                reserve_bytes=config.minimum_free_gb * 1_000_000_000,
            )
            print(open_targets_result_as_json(result), end="")
            return 0
        if args.command == "bootstrap-rare-disease-evidence":
            assert isinstance(config, ReferenceBuilderConfig)
            result = bootstrap_rare_disease_evidence(
                contract_path=RARE_DISEASE_CONTRACT_PATH,
                reference_root=config.reference_root,
                image_digest=args.image_digest,
                retrieved_on=args.retrieved_on,
                scratch_root=config.scratch_root,
                reserve_bytes=config.minimum_free_gb * 1_000_000_000,
            )
            print(rare_disease_result_as_json(result), end="")
            return 0
        if args.command == "bootstrap-variant-class-resources":
            assert isinstance(config, ReferenceBuilderConfig)
            result = bootstrap_variant_class_resources(
                contract_path=VARIANT_CLASS_RESOURCES_CONTRACT_PATH,
                reference_root=config.reference_root,
                image_digest=args.image_digest,
                retrieved_on=args.retrieved_on,
                scratch_root=config.scratch_root,
                reserve_bytes=config.minimum_free_gb * 1_000_000_000,
            )
            print(variant_class_result_as_json(result), end="")
            return 0
        if args.command == "hail-smoke":
            assert isinstance(config, ReferenceBuilderConfig)
            report = run_hail_smoke(
                reference_root=config.reference_root,
                scratch_root=config.scratch_root,
                run_id=args.run_id,
                workers=args.workers,
                rows=args.rows,
            )
            print(hail_report_as_json(report), end="")
            return 0
        if args.command == "scan-gnomad-bounded":
            assert isinstance(config, ReferenceBuilderConfig)
            report = run_bounded_public_scan(
                contract_path=GNOMAD_JOINT_PROJECTION_CONTRACT_PATH,
                reference_root=config.reference_root,
                scratch_root=config.scratch_root,
                run_id=args.run_id,
                workers=args.workers,
            )
            print(hail_report_as_json(report), end="")
            return 0
        if args.command == "run-gnomad-pilot":
            assert isinstance(config, ReferenceBuilderConfig)
            report = run_gnomad_pilot(
                contract_path=GNOMAD_JOINT_PROJECTION_CONTRACT_PATH,
                reference_root=config.reference_root,
                scratch_root=config.scratch_root,
                run_id=args.run_id,
                workers=args.workers,
                reserve_bytes=config.minimum_free_gb * 1_000_000_000,
                image_digest=args.image_digest,
            )
            print(pilot_report_as_json(report), end="")
            if not report.gate.allowed:
                parser.error(
                    "gnomAD continuation gate blocked: "
                    + ", ".join(report.gate.reasons)
                )
            return 0
        if args.command == "validate-gnomad-pilot":
            assert isinstance(config, ReferenceBuilderConfig)
            report = validate_gnomad_pilot(
                validation_contract_path=GNOMAD_PILOT_VALIDATION_CONTRACT_PATH,
                projection_contract_path=GNOMAD_JOINT_PROJECTION_CONTRACT_PATH,
                reference_root=config.reference_root,
                scratch_root=config.scratch_root,
                workers=args.workers,
                image_digest=args.image_digest,
            )
            print(validation_report_as_json(report), end="")
            return 0
        if args.command == "publish-gnomad-joint":
            assert isinstance(config, ReferenceBuilderConfig)
            result = publish_gnomad_joint_frequency(
                publication_contract_path=GNOMAD_JOINT_PUBLICATION_CONTRACT_PATH,
                projection_contract_path=GNOMAD_JOINT_PROJECTION_CONTRACT_PATH,
                validation_contract_path=GNOMAD_PILOT_VALIDATION_CONTRACT_PATH,
                reference_root=config.reference_root,
                scratch_root=config.scratch_root,
                reserve_bytes=config.minimum_free_gb * 1_000_000_000,
                image_digest=args.image_digest,
                retrieved_on=args.retrieved_on,
            )
            print(full_publication_result_as_json(result), end="")
            return 0
        if args.command == "bootstrap-gnomad-all-sites":
            assert isinstance(config, ReferenceBuilderConfig)
            result = bootstrap_gnomad_all_sites(
                contract_path=GNOMAD_ALL_SITES_CONTRACT_PATH,
                reference_root=config.reference_root,
                scratch_root=config.scratch_root,
                reserve_bytes=config.minimum_free_gb * 1_000_000_000,
                image_digest=args.image_digest,
                retrieved_on=args.retrieved_on,
            )
            print(gnomad_all_sites_result_as_json(result), end="")
            return 0
        if args.command == "sanitize-header":
            assert isinstance(config, Track1PipelineConfig)
            try:
                _require_distinct_output(Path(args.output).resolve(), args.input, args.path)
            except (OSError, RuntimeError):
                raise VcfHeaderError("header-output-invalid") from None
            sanitize_header_file(args.input, args.output, controlled_root=config.work_root)
            print('{"classification":"controlled","llm_policy":"no-llm","status":"sanitized"}')
            return 0
        if args.command == "validate-reference-concordance":
            assert isinstance(config, Track1PipelineConfig)
            output_path = Path(args.output).resolve()
            work_root = Path(config.work_root).resolve()
            if not output_path.is_relative_to(work_root):
                raise ReferenceConcordanceError(
                    "concordance report output must be under controlled work_root"
                )
            _require_distinct_output(
                output_path,
                args.path,
                args.input,
                args.revalidated_input,
                args.reference_manifest,
                args.order_adjudication,
                args.source_report,
                args.source_reference_manifest,
            )
            report = evaluate_installed_reference_concordance(
                identity_path=CHALLENGE_REFERENCE_IDENTITY_PATH,
                input_path=args.revalidated_input or args.input,
                source_input_path=(args.input if args.revalidated_input else None),
                source_report_path=args.source_report,
                reference_manifest_path=args.reference_manifest,
                reference_root=config.reference_root,
                source_reference_manifest_path=args.source_reference_manifest,
                source_reference_root=args.source_reference_root,
                order_adjudication_path=args.order_adjudication,
            )
            report_sha256 = write_concordance_report(
                output_path, report, allowed_root=work_root
            )
            print(
                json.dumps(
                    {
                        "accepted": report.accepted,
                        "report_sha256": report_sha256,
                        "schema": report.schema,
                    },
                    indent=2,
                    sort_keys=True,
                )
            )
            return 0
        if args.command == "adjudicate-contig-order":
            assert isinstance(config, Track1PipelineConfig)
            output_path = Path(args.output).resolve()
            work_root = Path(config.work_root).resolve()
            if not output_path.is_relative_to(work_root):
                raise ReferenceConcordanceError(
                    "contig-order adjudication output must be under controlled work_root"
                )
            _require_distinct_output(
                output_path,
                args.path,
                args.source_report,
                args.source_input,
                args.revalidated_input,
                args.reference_manifest,
                args.source_reference_manifest,
            )
            adjudication = document_installed_contig_order_adjudication(
                identity_path=CHALLENGE_REFERENCE_IDENTITY_PATH,
                source_report_path=args.source_report,
                source_input_path=args.source_input,
                revalidated_input_path=args.revalidated_input,
                reference_manifest_path=args.reference_manifest,
                reference_root=config.reference_root,
                source_reference_manifest_path=args.source_reference_manifest,
                source_reference_root=args.source_reference_root,
                resolution=args.resolution,
            )
            adjudication_sha256 = write_contig_order_adjudication(
                output_path, adjudication, allowed_root=work_root
            )
            print(
                json.dumps(
                    {
                        "accepted": adjudication.accepted,
                        "adjudication_sha256": adjudication_sha256,
                        "schema": adjudication.schema,
                    },
                    indent=2,
                    sort_keys=True,
                )
            )
            return 0
        if args.command == "authorize-rescue":
            assert isinstance(config, Track1PipelineConfig)
            authorization = require_rescue_reference_concordance(
                identity_path=CHALLENGE_REFERENCE_IDENTITY_PATH,
                report_path=args.concordance_report,
                input_path=args.input,
                validated_input_path=args.revalidated_input,
                reference_manifest_path=args.reference_manifest,
                reference_root=config.reference_root,
                order_adjudication_path=args.order_adjudication,
                source_report_path=args.source_report,
                source_reference_manifest_path=args.source_reference_manifest,
                source_reference_root=args.source_reference_root,
                operator_override=args.operator_override,
            )
            print(rescue_authorization_as_json(authorization), end="")
            return 0
    except (
        ConfigError,
        VcfHeaderError,
        StorageError,
        InsufficientCapacity,
        ReferenceConcordanceError,
        DownloadError,
        ManifestError,
        Grch38BootstrapError,
        GiabError,
        DeepVariantError,
        ChallengeFixtureError,
        OntologyBundleError,
        HpoICError,
        VepBundleError,
        ClinvarError,
        OpenTargetsError,
        RareDiseaseError,
        VariantClassResourceError,
        GnomadHailError,
        GnomadAllSitesError,
        GnomadProjectionError,
        OSError,
    ) as error:
        if args.command == "sanitize-header":
            parser.error("controlled header operation failed")
        if args.command in {
            "validate-reference-concordance",
            "adjudicate-contig-order",
            "authorize-rescue",
        }:
            reason = "controlled reference-concordance operation failed"
            if args.command == "authorize-rescue" and args.operator_override:
                reason += "; operator override cannot bypass reference concordance"
            parser.error(reason)
        parser.error(str(error))
    if args.command == "show-config":
        print(config_as_json(config))
    else:
        print(f"valid: {config.schema}")
    return 0


def reference_main(argv: Sequence[str] | None = None) -> int:
    return _run_role(
        _role_parser("genosapiens-reference", "reference-builder"),
        ReferenceBuilderConfig,
        argv,
    )


def track1_main(argv: Sequence[str] | None = None) -> int:
    arguments = list(sys.argv[1:] if argv is None else argv)
    if arguments and arguments[0] == 'run':
        from genosapiens.wgs_operator_invocation import run_operator_command
        return run_operator_command(arguments[1:])
    return _run_role(
        _role_parser("genosapiens-track1", "track1-pipeline"),
        Track1PipelineConfig,
        arguments,
    )




def _root_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="genosapiens")
    parser.add_argument("--version", action="version", version=__version__)
    commands = parser.add_subparsers(dest="role", required=True)
    commands.add_parser("reference", help="reference bootstrap commands")
    commands.add_parser("track1", help="Track 1 case analysis commands")
    commands.add_parser("track2", help="Track 2 public evidence and research reporting")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    arguments = list(sys.argv[1:] if argv is None else argv)
    parser = _root_parser()
    if not arguments:
        parser.print_help()
        return 0
    if arguments[0] == "reference":
        return reference_main(arguments[1:])
    if arguments[0] == "track1":
        return track1_main(arguments[1:])
    if arguments[0] == "track2":
        from genosapiens.track2_cli import main as track2_main
        return track2_main(arguments[1:])
    return parser.parse_args(arguments) or 0
