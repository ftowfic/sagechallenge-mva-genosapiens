"""Resolve the public index selected by an immutable installed package."""
from pathlib import Path
import re
from genosapiens.bounded_json import load_bounded_json

PIN_PATH = 'config/track1/installed-public-index.json'


def validate_index_pin(value):
    if (type(value) is not dict or set(value) != {'manifest_path', 'manifest_sha256'}
            or type(value['manifest_sha256']) is not str
            or re.fullmatch(r'[0-9a-f]{64}', value['manifest_sha256']) is None
            or type(value['manifest_path']) is not str
            or re.fullmatch(r'manifests/common-allele-public-indexes/[0-9a-f]{64}\.json',
                            value['manifest_path']) is None):
        raise ValueError('Installed public index metadata is invalid.')
    return value


def load_installed_index_pin(repository_root):
    root = Path(repository_root).resolve()
    path = root / PIN_PATH
    if not path.exists() and not path.is_symlink():
        return None  # Public validation fixtures supply their own fixed index.
    value, _ = load_bounded_json(path, reference_root=root, maximum_bytes=65536)
    return validate_index_pin(value)
