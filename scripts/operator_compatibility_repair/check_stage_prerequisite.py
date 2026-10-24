"""Operator-only prerequisite diagnostics; never emit receipt contents."""
from pathlib import Path
import argparse
import json
import sys

PREPARATION_STAGES = {'input_authentication', 'source', 'source_audit', 'prefilter', 'prepared_case'}
COMPATIBILITY_SCHEMA = 'genosapiens.track1-approved-header-compatibility-execution'
BOUNDED_PREPARATION_SCHEMA = 'genosapiens.wgs-bounded-preparation-result'


PREPARATION_RUNTIME_PHASES = (
    'input_validation',
    'hpo_initialization',
    'source_ingestion',
    'source_accounting',
    'filter_publication',
    'filter_index_commit',
    'filter_closure',
    'working_set_publication',
    'working_set_closure',
    'preparation_terminal_publication',
)


def failure_reason(step, record, release, runtime=None, completed_source=None):
    if not isinstance(record, dict):
        return 'Prerequisite receipt must be a JSON object.'
    prerequisite = 'Compatibility' if step == 'preparation' else 'Preparation'
    expected_schema = (COMPATIBILITY_SCHEMA if step == 'preparation'
                       else BOUNDED_PREPARATION_SCHEMA)
    expected_version = 1 if step == 'preparation' else 3
    if (record.get('schema') != expected_schema
            or type(record.get('version')) is not int
            or record.get('version') != expected_version
            or record.get('classification') != 'controlled'):
        return prerequisite + ' receipt has the wrong schema; ' + step + ' remains blocked.'
    if step == 'analysis':
        from dataclasses import fields
        from genosapiens.wgs_preparation import PreparationRunReceipt
        try:
            if set(record) != {field.name for field in fields(PreparationRunReceipt)}:
                raise ValueError
            values = {**record, 'completed_journal_heads': tuple(record['completed_journal_heads'])}
            PreparationRunReceipt(**values)
        except (ValueError, TypeError, KeyError, RuntimeError):
            return 'Preparation receipt is invalid; analysis remains blocked.'
    if record.get('status') != 'passed':
        message = prerequisite + ' did not pass; ' + step + ' remains blocked.'
        if step == 'analysis':
            stage = record.get('terminal_stage')
            if isinstance(stage, str) and stage in PREPARATION_STAGES:
                message += ' Failed stage: ' + stage + '.'
        return message + ' Inspect the prerequisite status and logs locally; existing results are preserved.'
    if record.get('base_release_sha256') != release:
        return prerequisite + ' receipt belongs to a different release; ' + step + ' remains blocked.'
    if step == 'preparation':
        if record.get('complete_dictionary_compatible') is not True:
            return prerequisite + ' receipt has no confirmed compatible dictionary; preparation remains blocked.'
    elif (not isinstance(record.get('prepared_case_sha256'), str)
          or len(record['prepared_case_sha256']) != 64
          or any(value not in '0123456789abcdef' for value in record['prepared_case_sha256'])
          or not isinstance(record.get('early_snapshot_sha256'), str)
          or len(record['early_snapshot_sha256']) != 64
          or any(value not in '0123456789abcdef' for value in record['early_snapshot_sha256'])):
        return 'Preparation receipt has no authenticated prepared authority; analysis remains blocked.'
    if step == 'analysis':
        from genosapiens.wgs_runtime import (
            WgsRuntimeError,
            require_capacity_qualified_runtime,
        )
        try:
            require_capacity_qualified_runtime(
                runtime,
                completed_source=completed_source,
                required_phases=PREPARATION_RUNTIME_PHASES + (
                    ('source_audit',) if record.get('audit_requested') is True and record.get('source_origin') == 'produced_this_run' else ()
                ),
            )
        except WgsRuntimeError:
            return 'Preparation runtime evidence is incomplete or failed; analysis remains blocked.'
    return None


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('step', choices=('preparation', 'analysis'))
    parser.add_argument('receipt', type=Path)
    parser.add_argument('release')
    parser.add_argument('--runtime', type=Path)
    parser.add_argument('--runtime-attempt-id')
    parser.add_argument('--execution-config-sha256')
    args = parser.parse_args(argv)
    try:
        if args.receipt.is_symlink():
            raise ValueError
        with args.receipt.open('rb') as source:
            content = source.read(65537)
        if len(content) > 65536:
            raise ValueError
        record = json.loads(content)
    except (OSError, UnicodeDecodeError, ValueError):
        print('Prerequisite receipt is unavailable or invalid; ' + args.step + ' remains blocked. Inspect the prerequisite status locally.', file=sys.stderr)
        return 1
    runtime = None
    if args.step == 'analysis':
        try:
            if (
                args.runtime is None
                or args.runtime_attempt_id is None
                or args.execution_config_sha256 is None
                or args.runtime.is_symlink()
            ):
                raise ValueError
            from genosapiens.wgs_runtime import load_stage_runtime_evidence
            with args.runtime.open('rb') as source:
                runtime_content = source.read(65537)
            if len(runtime_content) > 65536:
                raise ValueError
            runtime = load_stage_runtime_evidence(
                runtime_content,
                expected_stage='preparation',
                expected_attempt_id=args.runtime_attempt_id,
                expected_execution_config_sha256=args.execution_config_sha256,
            )
        except (OSError, ValueError, RuntimeError):
            print('Preparation runtime evidence is unavailable or invalid; analysis remains blocked.', file=sys.stderr)
            return 1
    completed_source = None
    if args.step == "analysis" and record.get("source_origin") == "completed_input":
        from genosapiens.wgs_preparation import (
            load_preparation_source_binding, load_completed_source_input, resolve_completed_source,
        )
        try:
            binding = load_preparation_source_binding(args.receipt.parent,
                expected_sha256=record["completed_source_input_sha256"])
            descriptor = load_completed_source_input(binding.read())
            completed_source = resolve_completed_source(descriptor,
                expected_source=descriptor.producer.source, audit_requested=record["audit_requested"])
        except (OSError, ValueError, RuntimeError):
            print("Completed source dependency is unavailable; analysis remains blocked.", file=sys.stderr)
            return 1
    reason = failure_reason(args.step, record, args.release, runtime, completed_source)
    if reason is not None:
        print(reason, file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
