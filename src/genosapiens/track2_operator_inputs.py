"""Small operator intake and human review records. No scientific policy lives here."""

import argparse
import csv
from datetime import datetime, timezone
import io
import json
import os
from pathlib import Path
import stat
import sys

from genosapiens.bounded_json import read_bounded_bytes
from genosapiens.track2_published_genotype import canonical, sha256
from genosapiens.wgs_operator_invocation import positive_hours


class OperatorError(ValueError):
    """Static operator diagnostic. Never include case values or exception locals."""


def require(ok, message):
    if not ok:
        raise OperatorError(message)


def failure_diagnostic(error):
    """Code locations only: never exception messages, source lines or case locals."""
    lines=['The offline stage failed. Diagnostic code locations:']
    seen=set()
    while error is not None and id(error) not in seen and len(seen)<4:
        seen.add(id(error))
        label=type(error).__name__
        if type(getattr(error,'errno',None)) is int:
            label+=f' (errno {error.errno})'
        lines.append(label)
        locations=[];frame=error.__traceback__
        while frame is not None:
            module=frame.tb_frame.f_globals.get('__name__','')
            if module.startswith('genosapiens.'):
                locations.append(f'  {module}.{frame.tb_frame.f_code.co_name}:{frame.tb_lineno}')
            frame=frame.tb_next
        lines.extend(locations[-8:])
        error=error.__cause__ or error.__context__
    lines.append('Read the private current diagnostics before retrying.')
    return '\n'.join(lines)


def parser():
    result = argparse.ArgumentParser(prog='genosapiens track2 run',
        description='Build and review a diagnosis-guided report offline.')
    result.add_argument('--output', type=Path, required=True)
    result.add_argument('--vcf', type=Path)
    result.add_argument('--diagnosis')
    result.add_argument('--sample')
    for name in ('hpo', 'concordance', 'review', 'mechanism-review'):
        result.add_argument('--' + name, type=Path)
    result.add_argument('--resume', action='store_true')
    result.add_argument('--require-review', action='store_true',
        help='Pause for input and plan review before report generation. The default generates a draft report.')
    result.add_argument('--max-hours', type=positive_hours,
        help='Positive finite hours for this invocation. Omission means no overall deadline.')
    return result


def validate_command(args):
    require(args.mechanism_review is None or args.resume,'Use --resume with a completed mechanism worksheet.')
    require(args.mechanism_review is None or args.review is None,'Submit one review stage per invocation.')
    require(args.resume or args.vcf is not None and args.diagnosis is not None,
        'Supply --vcf and --diagnosis for a new run.')


def private_path(value, *, directory=False, missing=False):
    path = Path(value).absolute()
    require(path.resolve() == path and not any(p.is_symlink() for p in (path, *path.parents)),
            'Use a real private path without symbolic links.')
    if missing and not path.exists():
        return path
    info = path.stat()
    require(info.st_uid == os.getuid() and stat.S_IMODE(info.st_mode) & 0o077 == 0,
            'Private files must belong to you and deny access to other users.')
    require(stat.S_ISDIR(info.st_mode) if directory else stat.S_ISREG(info.st_mode),
            'The supplied path has the wrong file type.')
    return path


def small_json(path, *, maximum=32 * 1024 * 1024):
    try:
        return json.loads(read_bounded_bytes(path, maximum_bytes=maximum))
    except (OSError, ValueError):
        raise OperatorError('A required local record is missing or invalid.') from None


def write_new(path, value):
    """Publish a new small record; never overwrite a prior attempt or review."""
    path = Path(path)
    path.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
    raw = value if isinstance(value, bytes) else canonical(value)
    temporary = path.with_name('.' + path.name + '.' + os.urandom(8).hex())
    try:
        descriptor = os.open(temporary, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
        with os.fdopen(descriptor, 'wb') as handle:
            handle.write(raw)
            handle.flush()
            os.fsync(handle.fileno())
        os.link(temporary, path)  # Atomic publication, with no replacement.
    finally:
        temporary.unlink(missing_ok=True)
    return sha256(raw)


def intake(vcf, sample=None, hpo=None, *, private_output=None):
    """Bounded syntax checks only; callers authorize access before mounting inputs."""
    from genosapiens.vcf_qc import read_vcf_samples, VcfQcError
    from genosapiens.vcf_contracts import DEFAULT_MAX_HEADER_BYTES
    from genosapiens.phenotype_ingestion import _hpo_lines, MAX_HPO_LIST_BYTES, PhenotypeIngestionError
    try:
        samples = read_vcf_samples(vcf, max_header_bytes=DEFAULT_MAX_HEADER_BYTES)
    except VcfQcError as error:
        if 'capacity' in str(error):
            raise OperatorError('VCF intake is incomplete: header capacity exceeded.') from error
        raise OperatorError('Invalid VCF header. Supply unique sample names in a supported VCF.') from None
    except (ValueError, OSError, VcfQcError):
        raise OperatorError('Invalid VCF header. Supply unique sample names in a supported VCF.') from None
    if (sample is None and len(samples) != 1) or (sample is not None and sample not in samples):
        print('Private sample choices: ' + json.dumps(samples), file=private_output or sys.stderr)
        raise OperatorError('Choose an exact sample with --sample NAME from the private list.')
    selected = samples[0] if sample is None else sample
    hpo_raw = None
    if hpo is not None:
        try:
            hpo_raw = read_bounded_bytes(hpo, maximum_bytes=MAX_HPO_LIST_BYTES)
            _hpo_lines(hpo_raw)
        except (OSError, ValueError, PhenotypeIngestionError):
            raise OperatorError('Use UTF-8 HPO text with one HP: identifier and seven digits per line.') from None
    return selected, hpo_raw


REVIEW_FIELDS = ('object_id', 'object_sha256', 'decision', 'reviewer',
                 'independence', 'rationale', 'limitations')


def review_worksheet(kind, object_sha256):
    require(kind in {'input', 'plan'}, 'Invalid review object type.')
    buffer = io.StringIO(newline='')
    writer = csv.DictWriter(buffer, fieldnames=REVIEW_FIELDS, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    writer.writerow(dict(object_id=kind, object_sha256=object_sha256,
                         decision='pending_human_review'))
    return buffer.getvalue().encode('utf-8')


def read_review(path, *, kind, object_sha256, governance, host_uid):
    """Only the authenticated local reviewer can enter a human decision."""
    from genosapiens.participant_governance import authorize_controlled_access
    try:
        raw = read_bounded_bytes(private_path(path), maximum_bytes=64 * 1024)
        reader = csv.DictReader(io.StringIO(raw.decode('utf-8'), newline=''), delimiter='\t')
        require(tuple(reader.fieldnames or ()) == REVIEW_FIELDS, 'Use the generated review worksheet columns.')
        rows = list(reader)
        require(len(rows) == 1 and set(rows[0]) == set(REVIEW_FIELDS),
                'The review worksheet must contain exactly one complete object row.')
        row = rows[0]
        require(row['object_id'] == kind and row['object_sha256'] == object_sha256,
                'Review bindings changed. Use the worksheet for the current object.')
        require(row['decision'] in {'accepted', 'rejected', 'pending_human_review'},
                'Choose accepted, rejected or pending_human_review.')
        require(all(isinstance(row[k], str) and row[k].strip() for k in REVIEW_FIELDS),
                'Complete every human review field before continuing.')
        require(row['independence'] in {'independent', 'self_review_disclosed'},
                'Declare independent review or self_review_disclosed.')
        authorize_controlled_access(governance, row['reviewer'], host_uid)
        return {**row, 'reviewed_at': datetime.now(timezone.utc).isoformat(),
                'worksheet_sha256': sha256(raw)}
    except (OSError, UnicodeError, csv.Error):
        raise OperatorError('The review worksheet is unreadable. Use the generated UTF-8 TSV format.') from None


def input_review_record(bundle, row):
    return {'schema': 'genosapiens.track2-diagnosis-review/v1',
        'reviewer': row['reviewer'], 'reviewer_type': 'human',
        'independence': row['independence'], 'reviewed_at': row['reviewed_at'],
        'decision': row['decision'], 'rationale': row['rationale'],
        'limitations': [row['limitations']], 'input_sha256': bundle.input_sha256,
        'source_hashes': bundle.source_hashes,
        'scenario_ids': [s['scenario_id'] for s in bundle.components['scenarios']]}


def plan_review_record(plan_path, bundle, row):
    from genosapiens.track2_synthesis import assessment_bindings
    plan=small_json(plan_path)
    return {'schema': 'genosapiens.track2-plan-review/v1',
        'reviewer': row['reviewer'], 'reviewer_type': 'human',
        'independence': row['independence'], 'reviewed_at': row['reviewed_at'],
        'decision': row['decision'],
        'rationale': row['rationale'] + '\nLimitations: ' + row['limitations'],
        'plan_sha256': sha256(Path(plan_path).read_bytes()),
        'source_hashes': {**bundle.source_hashes, **assessment_bindings(plan),
            **{s['record']['import_id']:s['record']['document_checksum'] for s in plan['evidence_imports']}}}
