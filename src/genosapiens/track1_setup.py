"""Automatic public-only preparation for the ordinary Docker launcher."""
from dataclasses import replace
import json
import os
from pathlib import Path
import sys
import tempfile
from datetime import datetime, timezone

from genosapiens.bounded_json import load_bounded_json, read_bounded_bytes
from genosapiens.track1_installation import canonical, digest
from genosapiens.installed_public_index import validate_index_pin


def atomic_json(path, value):
    path = Path(path)
    descriptor, temporary = tempfile.mkstemp(prefix='.' + path.name, dir=path.parent)
    try:
        with os.fdopen(descriptor, 'wb') as stream:
            stream.write(canonical(value))
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        Path(temporary).unlink(missing_ok=True)


def current_index_request(repository):
    arguments, _ = load_bounded_json(repository / 'config/public-resources/track1-retention-index.json')
    # Code belongs to the execution package. Published data stays reusable when
    # its source versions and supported evidence format remain compatible.
    return arguments, digest(canonical(dict(
        schema='genosapiens.track1-public-index-request/v2', arguments=arguments)))


def check_completed_index(reference, pin, *, source_bindings=None):
    """Reuse completed production and read only metadata plus indexed coverage."""
    from genosapiens.common_allele_public_indexes import (
        verify_common_allele_public_indexes, CommonAllelePublicIndexSession)
    from genosapiens.bounded_wgs_config import PublicQueryLimits
    validate_index_pin(pin)
    verified = verify_common_allele_public_indexes(reference_root=reference,
        manifest_path=pin['manifest_path'], expected_manifest_sha256=pin['manifest_sha256'],
        verify_sources=False, reuse_qualified_reconciliation=True)
    with CommonAllelePublicIndexSession(verified, limits=PublicQueryLimits(64, 10, 16*1024**2, 60)) as session:
        coverage = session.retention_coverage()
    if source_bindings is not None and coverage['source_bindings'] != source_bindings:
        raise ValueError('Published public index source bindings differ.')


def completed_index_sources(reference, pin, arguments):
    """Match pinned source versions through the small published manifest only."""
    manifest, content = load_bounded_json(reference / pin['manifest_path'],
        reference_root=reference, maximum_bytes=1024*1024)
    if digest(content) != pin['manifest_sha256']:
        raise ValueError('Published public index manifest checksum differs.')
    closure = manifest.get('identity', {}).get('closure', {})
    contracts = {name: arguments[name + '_contract_sha256']
                 for name in ('clinvar', 'open_targets', 'reference')}
    manifests = {name: arguments[name + '_manifest_sha256']
                 for name in ('clinvar', 'open_targets', 'reference')}
    manifests['exceptions'] = arguments['exception_manifest_sha256']
    curated = closure.get('curated_source_binding', {})
    accessions = [row for row in closure.get('source_artifacts', [])
                  if row.get('role') == 'curated_accessions']
    if (closure.get('contracts') != contracts or closure.get('manifests') != manifests
            or curated.get('bundle_manifest_sha256') != arguments['curated_bundle_sha256']
            or len(accessions) != 1
            or accessions[0].get('path') != arguments['accession_report_path']
            or accessions[0].get('sha256') != arguments['accession_report_sha256']):
        return None
    gene_output = curated.get('outputs', {}).get('gene_assertion', {})
    return dict(clinvar=manifests['clinvar'], open_targets=manifests['open_targets'],
        reference=manifests['reference'], curated=gene_output.get('sha256'),
        accessions=arguments['accession_report_sha256'])


def prepare_public_index(repository, reference, store):
    from genosapiens.common_allele_public_indexes import build_common_allele_public_indexes
    arguments, request = current_index_request(repository)
    saved = store / 'public-index.json'
    if saved.exists():
        record, _ = load_bounded_json(saved, reference_root=store)
        pin = validate_index_pin(record['index'])
        source_bindings = completed_index_sources(reference, pin, arguments)
        if source_bindings is not None:
            check_completed_index(reference, pin, source_bindings=source_bindings)
            if record.get('request_sha256') != request:
                atomic_json(saved, dict(request_sha256=request, index=pin))
            atomic_json(store / 'preparation-progress.json', dict(status='reused', index=pin,
                request_sha256=request, updated_at=datetime.now(timezone.utc).isoformat()))
            print('Using the published public selection index.', file=sys.stderr, flush=True)
            return pin
    print('Preparing the required public selection index. This first build can take time.', file=sys.stderr, flush=True)
    # This producer receives only frozen public inputs. No case data is mounted.
    def progress(stage):
        print('Public preparation: ' + stage + '.', file=sys.stderr, flush=True)
        atomic_json(store / 'preparation-progress.json', dict(status='running', stage=stage,
            request_sha256=request, updated_at=datetime.now(timezone.utc).isoformat()))
    try:
        built = build_common_allele_public_indexes(repository_root=repository,
            reference_root=reference, progress=progress, **arguments)
    except Exception as error:
        atomic_json(store / 'preparation-progress.json', dict(status='incomplete', error=str(error),
            request_sha256=request, updated_at=datetime.now(timezone.utc).isoformat()))
        raise
    pin = validate_index_pin(dict(manifest_path=built.manifest_path, manifest_sha256=built.manifest_sha256))
    atomic_json(saved, dict(request_sha256=request, index=pin))
    atomic_json(store / 'preparation-progress.json', dict(status='complete', index=pin,
        request_sha256=request, updated_at=datetime.now(timezone.utc).isoformat()))
    return pin


def prepare_early_scope(repository, reference, store, index_pin):
    """Run the public producer inside the existing installation owner and lock."""
    from genosapiens.track1_early_scope import build_early_scope
    print('Preparing or reusing the public early scope.', file=sys.stderr, flush=True)
    try:
        pin = build_early_scope(repository, reference, index_pin, store / 'early-scope.json')
        atomic_json(store / 'early-scope-progress.json', dict(status='complete', scope=pin))
        return pin
    except (Exception, KeyboardInterrupt):
        atomic_json(store / 'early-scope-progress.json', dict(status='incomplete'))
        raise


def public_versions(repository, reference, pin, scope_pin):
    """Derive versions from current public metadata without a previous package."""
    from genosapiens import track1_public_validation as public
    from genosapiens.reference_concordance import load_reference_identity
    from genosapiens.variant_retention import POLICY_SHA256
    descriptor, descriptor_bytes = load_bounded_json(repository / public.MODEL_LOCUS_CONTRACT_PATH)
    reference_bytes = read_bounded_bytes(reference / descriptor['reference']['manifest_path'],
                                         reference_root=reference)
    if digest(reference_bytes) != public.REFERENCE_MANIFEST_SHA256:
        raise ValueError('The installed reference manifest differs from the current public version.')
    loci = public.load_public_model_loci(repository, reference)
    policies = public.load_public_policy_family(repository, model_loci=loci)
    if len(policies.runtimes) != 1:
        raise ValueError('Installation requires one current scientific policy.')
    bindings = replace(policies.runtimes[0].candidate.bindings,
                       common_allele_public_index_manifest_sha256=pin['manifest_sha256'])
    identity_path = repository / public.REFERENCE_IDENTITY_PATH
    identity = load_reference_identity(identity_path)
    handoff, handoff_bytes = load_bounded_json(repository / 'config/track1/handoff-contract.json')
    reference_versions = dict(
        reference_identity_id=identity.reference_role('public_annotation_benchmark').identity_id,
        reference_manifest_sha256=public.REFERENCE_MANIFEST_SHA256,
        reference_bindings_sha256=digest(canonical(dict(identity_contract_sha256=digest(identity_path.read_bytes()),
            reference_manifest_sha256=digest(reference_bytes), model_loci_descriptor_sha256=digest(descriptor_bytes)))),
        assessment_scope_sha256=digest(canonical(dict(mode='base-supplied-call-vcf-first',
            analysis_profile='focused', retention_policy_sha256=POLICY_SHA256))),
        evaluator_pin_sha256=digest((repository / 'config/track1/challenge-evaluator.json').read_bytes()),
        canonical_handoff_path_sha256=digest(canonical(handoff['canonical_artifact_path'])),
        path_selection_policy_sha256=digest(handoff_bytes))
    from genosapiens.track1_early_scope import load_completed_early_scope
    scope = load_completed_early_scope(reference, scope_pin)
    resource_versions = dict(scope_policy_sha256=scope.document["request"]["policy_sha256"],
        scope_manifest_sha256=scope.sha256, population_manifest_sha256=public.GNOMAD_JOINT_MANIFEST_SHA256,
        clinical_manifest_sha256=public.CLINVAR_MANIFEST_SHA256,
        exception_manifest_sha256=bindings.common_allele_exception_manifest_sha256,
        precision_manifest_sha256=bindings.gnomad_projection_precision_sha256,
        common_allele_public_index_manifest_sha256=pin['manifest_sha256'],
        ontology_manifest_sha256=public.HPO_MANIFEST_SHA256)
    return reference_versions, resource_versions, bindings


def prepare_installation(repository, reference, store):
    from genosapiens.track1_public_validation_coordinator import executable_stack_bytes
    from genosapiens.wgs_installation import build_package, authenticate_package
    repository, reference, store = (Path(p).resolve(strict=True) for p in (repository, reference, store))
    pin = prepare_public_index(repository, reference, store)
    scope_pin = prepare_early_scope(repository, reference, store, pin)
    stack = executable_stack_bytes(repository)
    request = digest(canonical(dict(stack_sha256=digest(stack), public_index=pin, early_scope=scope_pin)))
    saved = store / 'current-package.json'
    if saved.exists():
        record, _ = load_bounded_json(saved, reference_root=store)
        if record.get('request_sha256') == request:
            name = record['package_name']
            if not isinstance(name, str) or Path(name).name != name or not name.startswith('package-'):
                raise ValueError('Saved package location is invalid.')
            authenticate_package(store / name, record['package_sha256'])
            return record
    versions, resources, bindings = public_versions(repository, reference, pin, scope_pin)
    # Each build owns a fresh staging directory. Failed builds remain incomplete.
    staging = Path(tempfile.mkdtemp(prefix='.build-', dir=store))
    package = staging / 'package'
    sha = build_package(repository_root=repository, destination=package,
        reference_versions=versions, resource_versions=resources, policy_bindings=bindings,
        public_index_pin=pin, scope_pin=scope_pin)
    name = 'package-' + sha
    destination = store / name
    if destination.exists():
        authenticate_package(destination, sha)
        import shutil
        shutil.rmtree(staging)
    else:
        package.rename(destination)
        staging.rmdir()
    record = dict(request_sha256=request, package_name=name, package_sha256=sha,
                  runtime_image_digest=json.loads(stack)['runtime_image'].split('@')[-1])
    atomic_json(saved, record)
    print('Track 1 package is ready.', file=sys.stderr, flush=True)
    return record
