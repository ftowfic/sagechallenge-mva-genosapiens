#!/usr/bin/env python3
"""Read private host paths from the current manifest before Docker mounts."""
import json
import os
from pathlib import Path
import stat
import sys
from portable_paths import reference_selection, track1_legacy_reference


def main():
    try:
        output=Path(sys.argv[1]).absolute()
        path=output/'.operator/run.json'
        if any(p.is_symlink() for p in (path,*path.parents)):
            raise ValueError
        metadata=path.stat()
        if (not stat.S_ISREG(metadata.st_mode) or metadata.st_uid != os.getuid()
                or stat.S_IMODE(metadata.st_mode) not in (0o400,0o600)
                or not 0 < metadata.st_size <= 2*1024**2):
            raise ValueError
        with path.open('rb') as stream:
            raw=stream.read(2*1024**2+1)
        data=json.loads(raw)
        if data['schema'] != 'genosapiens.track1-current-run' or data['version'] != 1:
            raise ValueError
        if data.get('preparation_contract') != 'track1-early-scope/v1':
            print('This run predates the early scope contract. Use a fresh output directory with this package. '
                  'Preserve the original run and package for separate recovery.', file=sys.stderr)
            return 2
        values=[]
        for name,supplied in zip(('vcf','hpo','concordance','sample'),sys.argv[2:6],strict=True):
            value=data['inputs'].get(name if name=='sample' else 'host_'+name)
            if supplied:
                current=supplied if name=='sample' else str(Path(supplied).resolve(strict=True))
                if current != value:
                    raise ValueError('conflict')
            if value is not None and (type(value) is not str or '\0' in value):
                raise ValueError
            if name != 'concordance' and not value:
                raise ValueError
            values.append(value or '')
        saved_root = data['inputs'].get('host_reference_root')
        reference = reference_selection(sys.argv[6] or None, saved_root, resume=True)
        if saved_root is None:
            track1_legacy_reference(data, reference, Path(__file__).resolve().parents[1])
        values.append(str(reference))
        sys.stdout.buffer.write(('\0'.join(values)+'\0').encode())
        return 0
    except (OSError,ValueError,KeyError,TypeError):
        print('Resume input or reference binding is missing, invalid, or conflicts with an explicit option. '
              'For a run without a saved root, supply --reference-root PATH with its original resource manifests.',file=sys.stderr)
        return 2

if __name__=='__main__':
    sys.exit(main())
