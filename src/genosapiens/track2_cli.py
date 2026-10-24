"""Current source-aware Track 2 evidence and rendering commands."""

import argparse
import hashlib
import json
import sys
from pathlib import Path

from genosapiens.llm_guard import LLMGuardError
from genosapiens.track2_codec import decode
from genosapiens.track2_manual_import import ManualEvidenceImport
from genosapiens.track2_policy import load_track2_policies
from genosapiens.track2_openrouter import ReportLanguageError
from genosapiens.track2_report_tables import (
    ReportTableInputs, build_report_tables, report_tables_as_json, report_tables_as_markdown,
)

ROOT = Path(__file__).resolve().parents[2]


def policies():
    return load_track2_policies(
        ROOT / "config/track2/source-domain-matrix-v1.json",
        ROOT / "config/track2/decision-policy-g0-g7-v1.json",
    )


def main(argv=None) -> int:
    arguments = list(sys.argv[1:] if argv is None else argv)
    if arguments and arguments[0] == 'run':
        from genosapiens.track2_operator import run_operator_command
        from genosapiens.track2_operator_inputs import OperatorError, failure_diagnostic
        from genosapiens.track2_storage import Track2StorageError
        try:
            return run_operator_command(arguments[1:])
        except (OperatorError, Track2StorageError) as error:
            print(str(error), file=sys.stderr)
            return 1
        except (OSError, ValueError, RuntimeError, KeyError, TypeError, AttributeError) as error:
            print(failure_diagnostic(error),file=sys.stderr)
            return 1
    parser = argparse.ArgumentParser(prog="genosapiens track2")
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("validate-policies", help="validate the current source and scientific gate policies")
    commands.add_parser('run', help='run the installed diagnosis-guided operator workflow')
    discovery=commands.add_parser('discover',help='execute frozen native discovery jobs offline')
    discovery.add_argument('plan',type=Path)
    discovery.add_argument('reference_root',type=Path)
    discovery.add_argument('output',type=Path)
    discovery.add_argument('--invocation-id')
    discovery.add_argument('--lincs-workers',type=int,choices=(1,2),default=2)
    admission = commands.add_parser("validate-input", help="validate the explicitly selected input authority (G0)")
    admission.add_argument("bundle", type=Path)
    admission.add_argument("--origin", required=True,
                           choices=("published_patient_genotype", "track1_case_handoff", "diagnosis_guided_vcf"))
    admission.add_argument("--draft", action="store_true")
    admission.add_argument("--allow-synthetic", action="store_true", help="independently synthetic inputs only")
    queries = commands.add_parser("lincs-queries", help="freeze independent local genetic queries before drug scoring")
    for name in ("bundle", "plan", "sources", "output"):
        queries.add_argument(name, type=Path)
    check = commands.add_parser('validate-research-report', help='verify report outputs and required LINCS figures')
    check.add_argument('output', type=Path)
    language = commands.add_parser('rewrite-report-language',
        help='create an OpenRouter language draft from a public report, using mandatory ZDR routing')
    language.add_argument('source', type=Path)
    language.add_argument('output', type=Path)
    language.add_argument('--openrouter-key-file', required=True, type=Path,
        help='runtime secret file containing the API key; never stored in report artifacts')
    language.add_argument('--model', default='nvidia/nemotron-3.5-super-vl')
    language.add_argument('--account-tier', required=True, help='actual API plan/tier for the methods disclosure')
    language.add_argument('--account-privacy-confirmed', required=True, action='store_true',
        help='attest that OpenRouter input/output logging and content use are off, with no conflicting credit terms')
    language.add_argument('--max-tokens', type=int, default=8192)
    language.add_argument('--timeout', type=int, default=60)
    manual = commands.add_parser("validate-manual-import", help="check a reviewed extraction against source bytes")
    manual.add_argument("record", type=Path)
    manual.add_argument("source", type=Path)
    dossier = commands.add_parser("input-dossier", help="validate and render the selected input package")
    dossier.add_argument("bundle", type=Path)
    dossier.add_argument("output", type=Path)
    dossier.add_argument("--draft", action="store_true", help="render for review without admitting the input")
    dossier.add_argument("--origin", choices=("published_patient_genotype", "diagnosis_guided_vcf"),
                         default="published_patient_genotype")
    dossier.add_argument("--allow-synthetic", action="store_true", help="independently synthetic inputs only")
    render = commands.add_parser("render-tables", help="render current structured scientific records")
    render.add_argument("input", type=Path)
    render.add_argument("output", type=Path)
    synth = commands.add_parser("synthesize", help="replay frozen scientific inputs into an offline report")
    synth.add_argument("bundle", type=Path)
    synth.add_argument("plan", type=Path)
    synth.add_argument("envelope", type=Path)
    synth.add_argument("output", type=Path)
    synth.add_argument("--draft", action="store_true")
    synth.add_argument("--input-snapshot", type=Path)
    synth.add_argument("--snapshot-sha256")
    synth.add_argument("--allow-synthetic", action="store_true", help="independently synthetic inputs only")
    private = commands.add_parser("freeze-input-snapshot", help="bind a diagnosis package to an independent public envelope")
    for name in ("bundle", "plan", "envelope", "runtime", "output"):
        private.add_argument(name, type=Path)
    private.add_argument("--draft", action="store_true")
    private.add_argument("--allow-synthetic", action="store_true", help="independently synthetic inputs only")
    stage = commands.add_parser('stage-diagnosis-input', help='copy the exact immutable package closure for synthesis')
    stage.add_argument('bundle', type=Path)
    stage.add_argument('output', type=Path)
    stage.add_argument('--input-sha256', required=True)
    stage.add_argument('--review-sha256', required=True)
    stage.add_argument('--allow-synthetic', action='store_true', help='independently synthetic inputs only')
    for command in (dossier, synth, private, stage):
        command.add_argument('--controlled-store', type=Path,
            help='private JSON with reference_root, inventory_path and case_id; required for controlled output')
    frozen = commands.add_parser("freeze", help="bind reviewed public evidence and runtime")
    frozen.add_argument("bundle", type=Path)
    frozen.add_argument("plan", type=Path)
    frozen.add_argument("scope", type=Path)
    frozen.add_argument("runtime", type=Path)
    frozen.add_argument("reference_root", type=Path)
    frozen.add_argument("--date", required=True)
    frozen.add_argument("--draft", action="store_true")
    args = parser.parse_args(argv)
    try:
        store = None
        if getattr(args, 'controlled_store', None) is not None:
            from genosapiens.track2_storage import ControlledStore
            store = ControlledStore(**json.loads(args.controlled_store.read_bytes()))
        policy = policies()
        if args.command == 'discover':
            from genosapiens.track2_discovery import discover_file
            result=discover_file(args.plan,args.reference_root,args.output,invocation_id=args.invocation_id,lincs_workers=args.lincs_workers)
            print('Native discovery status: '+result['completion_state']+'.')
            return 0 if result['completion_state']=='complete' else 2
        elif args.command == "validate-policies":
            print(json.dumps({"valid": True, "policy_version": policy.decision_policy["policy_version"]}))
        elif args.command == "validate-input":
            from genosapiens.track2_input import load_analysis_input
            bundle = load_analysis_input(args.bundle, input_origin=args.origin, require_review=not args.draft,
                                         allow_synthetic=args.allow_synthetic)
            # Never print case identifiers or source checksums to an operator log.
            print(json.dumps({"input_report_ready": bundle.input_report_ready,
                              "input_origin": args.origin, "scientific_readiness": "not_assessed"}))
        elif args.command == "lincs-queries":
            from genosapiens.track2_lincs import freeze_queries
            print(json.dumps(freeze_queries(args.bundle, args.plan, args.sources, args.output)))
        elif args.command == 'validate-research-report':
            from genosapiens.track2_synthesis import validate_research_report
            from genosapiens.track2_terminal import assessment_outcome
            result=validate_research_report(args.output)
            outcome=assessment_outcome(result)
            print(json.dumps({**result,'terminal_outcome':outcome.as_dict()}))
            return outcome.exit_code
        elif args.command == 'rewrite-report-language':
            from genosapiens.track2_openrouter import rewrite_report
            print(json.dumps(rewrite_report(args.source, args.output, key_file=args.openrouter_key_file,
                model=args.model, max_tokens=args.max_tokens, timeout=args.timeout,
                account_tier=args.account_tier, account_privacy_confirmed=args.account_privacy_confirmed)))
        elif args.command == "synthesize":
            from genosapiens.track2_synthesis import synthesize
            from genosapiens.track2_terminal import assessment_outcome
            result = synthesize(args.bundle, args.plan, args.envelope, args.output, draft=args.draft,
                allow_synthetic=args.allow_synthetic, input_snapshot_path=args.input_snapshot,
                expected_snapshot_sha256=args.snapshot_sha256, controlled_store=store)
            outcome=assessment_outcome(result,require_review=not args.draft)
            print(json.dumps({**{key: result[key] for key in (
                "input_report_ready", "research_report_ready", "candidate_objective_met", "case_matching", "submission_ready",
                'assessment_status')}, 'terminal_outcome':outcome.as_dict()}))
            return outcome.exit_code
        elif args.command == "freeze-input-snapshot":
            from genosapiens.track2_synthesis import freeze_input_snapshot
            freeze_input_snapshot(args.bundle, args.plan, args.envelope, args.output,
                runtime=json.loads(args.runtime.read_bytes()), draft=args.draft, allow_synthetic=args.allow_synthetic,
                controlled_store=store)
            print(json.dumps({"private_input_snapshot_written": True}))
        elif args.command == 'stage-diagnosis-input':
            from genosapiens.track2_storage import stage_diagnosis_input
            stage_diagnosis_input(args.bundle, args.output, expected_input_sha256=args.input_sha256,
                expected_review_sha256=args.review_sha256, allow_synthetic=args.allow_synthetic,
                controlled_store=store)
            print(json.dumps({'minimal_package_staged': True}))
        elif args.command == "freeze":
            from genosapiens.track2_synthesis import freeze
            print(json.dumps(freeze(args.bundle, args.plan, args.scope, json.loads(args.runtime.read_bytes()),
                                    args.reference_root, args.date, draft=args.draft)))
        elif args.command == "validate-manual-import":
            record = decode(ManualEvidenceImport, json.loads(args.record.read_text()))
            if hashlib.sha256(args.source.read_bytes()).hexdigest() != record.document_checksum:
                raise ValueError("manual source checksum mismatch")
            print(json.dumps({"valid": True, "import_id": record.import_id}))
        elif args.command == "render-tables":
            inputs = decode(ReportTableInputs, json.loads(args.input.read_text()))
            tables = build_report_tables(inputs, decision_policies={policy.decision_policy["policy_version"]: policy.decision_policy})
            args.output.mkdir(parents=True, exist_ok=False)
            (args.output / "tables.json").write_text(report_tables_as_json(tables))
            (args.output / "tables.md").write_text(report_tables_as_markdown(tables))
        elif args.command == "input-dossier":
            from genosapiens.track2_input import load_analysis_input
            from genosapiens.track2_published_genotype import render_input_dossier
            from genosapiens.track2_diagnosis_package import render_diagnosis_dossier
            bundle = load_analysis_input(args.bundle, input_origin=args.origin, require_review=not args.draft,
                                         allow_synthetic=args.allow_synthetic)
            from genosapiens.track2_storage import output_storage
            with output_storage(args.output, classification=bundle.data['classification'],
                                controlled_store=store, kind='track2-input-dossier') as target:
                with target.open("x") as stream:
                    stream.write(render_diagnosis_dossier(bundle) if args.origin == "diagnosis_guided_vcf"
                                 else render_input_dossier(bundle))
            print(json.dumps({"input_report_ready": bundle.input_report_ready,
                              "input_origin": args.origin, "scientific_readiness": "not_assessed"}))
    except ReportLanguageError as error:
        # This exception exposes only fixed operational messages, never remote bodies.
        parser.error(str(error))
    except (OSError, ValueError, TypeError, KeyError, RuntimeError, MemoryError, LLMGuardError) as error:
        if args.command == 'discover':
            from genosapiens.track2_discovery_outcome import safe_diagnostic
            print(json.dumps(safe_diagnostic('discovery','discovery-stage',error)),file=sys.stderr)
            return 3
        # Library exceptions may contain paths, identifiers or invalid source
        # values. Keep operational output safe for both input classifications.
        parser.error("Track 2 operation failed input or integrity validation")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
