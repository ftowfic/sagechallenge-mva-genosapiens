"""Bounded host metadata checks for reference-root selection. Standard library only."""
import hashlib
import json
from pathlib import Path


def real_directory(value):
    path = Path(value).absolute()
    if path.resolve(strict=True) != path or not path.is_dir():
        raise ValueError('The reference root must be a real directory without symbolic links.')
    return path


def reference_selection(supplied, saved=None, *, resume=False):
    if saved is None and resume and supplied is None:
        raise ValueError('This run has no saved reference root. Supply --reference-root PATH for compatibility checks.')
    if saved is not None:
        if not isinstance(saved, str) or not Path(saved).is_absolute():
            raise ValueError('The saved reference root is invalid.')
        if supplied is not None and Path(supplied).absolute() != Path(saved):
            raise ValueError('Explicit reference root conflicts with the saved run. Restore its original root.')
    return real_directory(saved if saved is not None else supplied or '/mnt/reference')


def bound_json(path, expected=None, *, maximum=2*1024**2):
    path = Path(path).absolute()
    if path.resolve(strict=True) != path:
        raise ValueError('Bound resource metadata must not contain symbolic links.')
    with path.open('rb') as stream:
        content = stream.read(maximum+1)
    if len(content) > maximum or (expected is not None and hashlib.sha256(content).hexdigest() != expected):
        raise ValueError('The retained resource metadata binding differs. Restore the original metadata.')
    return json.loads(content)


def track1_legacy_reference(data, reference, repository):
    """Compare six small manifests to retained run contracts before preparation."""
    import sys
    sys.path.insert(0, str(Path(repository)/'src'))
    from genosapiens import track1_public_validation as pins
    paths = {
        'reference': bound_json(Path(repository)/'config/public-resources/grch38-model-loci.json')['reference']['manifest_path'],
        'population': bound_json(Path(repository)/'config/public-resources/gnomad-joint-publication.json')['manifest_path'],
        'clinical': pins.CLINVAR_MANIFEST_PATH,
        'exceptions': pins.COMMON_ALLELE_EXCEPTION_MANIFEST_PATH,
        'precision': pins.GNOMAD_PRECISION_REPORT_PATH,
        'ontology': 'manifests/phenotype/hpo-ic-060cf7c49671f9a98f9a8578c35f0fb8eeae62166d79cb40c2bcc81855abcc89.json',
    }
    bindings = data['contracts']['preparation']
    for role, relative in paths.items():
        expected = bindings[role]
        if not isinstance(expected, str) or len(expected) != 64:
            raise ValueError('The retained resource binding is invalid.')
        bound_json(reference/relative, expected)


def track2_legacy_reference(record, repository):
    """Recover the root only from a package authenticated by the saved run."""
    expected = record.get('source_versions', {}).get('package_sha256')
    if not isinstance(expected, str) or len(expected) != 64:
        return None
    for path in sorted((Path(repository)/'.track2-packages').glob('package-*/package.json')):
        try:
            package = bound_json(path, expected, maximum=32*1024**2)
        except ValueError:
            continue
        resources = bound_json(path.parent/'resources.json', package['files']['resources.json'], maximum=32*1024**2)
        return resources['reference']['root']
    return None
