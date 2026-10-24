"""Interactive first-use registration, using the existing signed access contract.

Only the operator supplies identity and attestations. This module never reads a
case, performs scientific review, submits a result, or resets incident state.
"""

import argparse
import base64
from datetime import datetime, timedelta, timezone
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from uuid import uuid4

from genosapiens.track2_operator_inputs import OperatorError, private_path, require, small_json, write_new
from genosapiens.track2_operator_install import configure, validate_package
from genosapiens.track2_published_genotype import canonical, sha256


ACCEPTANCES = (
    ('individual_registration', 'I have individually registered for the Hackathon.'),
    ('age_eligible', 'I meet its age eligibility requirement.'),
    ('huggingface_terms_accepted', 'I have accepted the applicable Hugging Face terms.'),
    ('hackathon_terms_accepted', 'I have read and accepted the applicable Hackathon terms.'),
    ('gated_data_authorized', 'I have authorization to access the gated challenge data.'),
    ('no_recontact_accepted', 'I accept the no-recontact obligation.'),
    ('no_resharing_accepted', 'I accept the no-resharing obligation.'),
    ('applicable_law_accepted', 'I accept the applicable-law obligations.'),
    ('professional_obligations_accepted', 'I accept the applicable professional obligations.'),
)


def validate_registration(settings_path, reference_root, rules=None):
    """Authenticate retained access before preparation or registration changes."""
    from genosapiens.track2_operator import authorize
    from genosapiens.bounded_json import read_bounded_bytes
    settings = small_json(private_path(settings_path))
    require(settings['reference_root'] == str(Path(reference_root).absolute()),
            'The registration uses a different controlled storage root.')
    authorize(settings)
    if rules is not None:
        supplied = read_bounded_bytes(rules, maximum_bytes=1024*1024)
        require(sha256(supplied) == settings['rules_sha256'],
                'Explicit rules conflict with the retained registration. Preserve its original rules binding.')
    return settings


def setup(package, rules, reference_root, state, image, *, ask=input, tell=print):
    """Create a solo registration only after explicit local confirmation."""
    require(os.getuid() != 0, 'Run setup as your ordinary non-root account.')
    state = private_path(state, directory=True)
    package, reference_root = Path(package).absolute(), Path(reference_root).absolute()
    require(state.is_relative_to(reference_root/'controlled-work'),
            'Keep setup records within the installed controlled storage root.')
    package_sha = sha256((package/'package.json').read_bytes())
    manifest = validate_package(package, package_sha)
    require(manifest.get('image',manifest.get('images',{}).get('diagnosis')) == image, 'The package and runtime pins disagree.')
    settings_path = state/'command/track2-settings.json'
    if settings_path.exists():
        validate_registration(settings_path, reference_root, rules)
        tell('Reusing the existing registration and incident state.')
        return settings_path
    require(not any(p.name != '.setup.lock' for p in state.iterdir()),
            'Setup is incomplete. Preserve its records for repair; do not replace its incident state.')
    from genosapiens.bounded_json import read_bounded_bytes
    rules_bytes = read_bounded_bytes(rules, maximum_bytes=1024*1024)
    require(bool(rules_bytes.strip()), 'The local rules file is empty.')
    tell('One-time setup for a solo participant. No case inputs will be read.')
    tell('This records your statements locally. It does not grant challenge access.')
    identity = ask('Your registered name or username: ').strip()
    require(bool(identity) and len(identity) <= 1024 and identity.isprintable(),
            'Enter your actual registered name or username.')
    tell('Confirm each statement below before continuing:')
    tell('- I am the only team member and controlled-data reviewer in this registration.')
    for _, statement in ACCEPTANCES:
        tell('- '+statement)
    tell('The local record will expire in 30 days. Submission and scientific review remain separate.')
    answer = ask('Type I CONFIRM to record these statements, or press Enter to cancel: ')
    require(answer == 'I CONFIRM', 'Setup cancelled. No authorization was created.')
    from genosapiens.participant_governance import (
        Participant, ParticipantGovernanceManifest, PARTICIPANT_GOVERNANCE_ATTESTATION_SCHEMA,
    )
    from genosapiens.incident_response import register_environment
    now = datetime.now(timezone.utc)
    participant_id = 'participant:'+uuid4().hex
    participant = Participant(participant_id=participant_id, registered_identity=identity,
        host_uid=os.getuid(), roles=('controlled_reviewer', 'team_member'),
        **{key: True for key, _ in ACCEPTANCES}, submission_material_cc_by_accepted=False,
        organizer_rerun_accepted=False, public_name_attribution_accepted=False)
    governance = ParticipantGovernanceManifest(rules_sha256=sha256(rules_bytes),
        issued_at=now.isoformat(), expires_at=(now+timedelta(days=30)).isoformat(),
        roster_complete=True, participants=(participant,), team_member_ids=(participant_id,),
        controlled_review_roster=(participant_id,), package_review=None)
    old_umask = os.umask(0o077)
    try:
        def openssl(*args):
            return subprocess.run(['openssl', *map(str, args)], check=True,
                                  capture_output=True, timeout=30).stdout
        signing, verification = state/'signing.pem', state/'verification.pem'
        openssl('genpkey', '-algorithm', 'ED25519', '-out', signing)
        openssl('pkey', '-in', signing, '-pubout', '-out', verification)
        with tempfile.NamedTemporaryFile() as payload:
            payload.write(canonical(governance.as_dict())); payload.flush()
            signature = openssl('pkeyutl', '-sign', '-rawin', '-inkey', signing, '-in', payload.name)
        signed = state/'governance.json'
        write_new(signed, {'schema':PARTICIPANT_GOVERNANCE_ATTESTATION_SCHEMA,
            'algorithm':'Ed25519', 'manifest':governance.as_dict(),
            'signature':base64.b64encode(signature).decode('ascii')})
        local_rules = state/'rules.txt'; write_new(local_rules, rules_bytes)
        incident = state/'incident'; incident.mkdir(mode=0o700)
        environment_id = register_environment(incident)
        write_new(state/'environment-id', (environment_id+'\n').encode())
        inventory = state/'inventory'; inventory.mkdir(mode=0o700)
        launcher = configure(package, reference_root=reference_root, governance=signed,
            verification_key=verification, rules=local_rules, participant_id=participant_id,
            incident_root=incident, environment_id=environment_id,
            inventory=inventory/'case.json', case_id='track2-'+uuid4().hex,
            destination=state/'command')
    finally:
        os.umask(old_umask)
    tell('Your local registration is complete. No scientific review or submission was recorded.')
    return launcher


def main(argv=None):
    parser = argparse.ArgumentParser(description='Complete first-time local Track 2 setup.')
    for key in ('package', 'rules', 'reference-root', 'state'):
        parser.add_argument('--'+key, required=True, type=Path)
    parser.add_argument('--image', required=True)
    args = parser.parse_args(argv)
    try:
        require(sys.stdin.isatty() and sys.stdout.isatty(),
                'Run ./track2 setup in your terminal. Piped acceptance is not supported.')
        setup(**vars(args))
    except (EOFError, KeyboardInterrupt):
        print('Setup interrupted. Existing records were preserved.', file=sys.stderr)
        return 130
    except OperatorError as error:
        print(str(error), file=sys.stderr)
        return 1
    except (OSError, ValueError, subprocess.SubprocessError):
        print('Setup could not complete. Preserve the local records for repair.', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
