#!/usr/bin/env python3
"""Print bounded metadata for a user-selected Track 2 report transition."""
import argparse
import hashlib
import json
from pathlib import Path
import stat

LIMIT=2*1024**2
FILES=('current-routing.json','completions/report.json','synthesis/completion.json',
    'synthesis/production.json','retained-report-routing.json','synthesis/controlled-work/result.json',
    'synthesis/controlled-work/report/manifest.json','preparation/snapshot.json')
SCHEMAS={'genosapiens.track2-report-production/v1','genosapiens.track2-scientific-completion/v1',
    'genosapiens.track2-research-manifest/v1','genosapiens.track2-challenge-package/v1',
    'genosapiens.track2-challenge-package/v2','genosapiens.track2-private-input-snapshot/v2'}


def diagnose(root):
    root=Path(root).absolute()
    if root.resolve()!=root:raise ValueError('Source path links are unsupported.')
    records={}
    for name in FILES:
        path=root/name;row={};records[name]=row
        try:
            if path.resolve()!=path:raise ValueError('linked')
            info=path.stat()
            if not stat.S_ISREG(info.st_mode):raise ValueError('not_regular')
            if info.st_size>LIMIT:
                row['state']='exceeds_diagnostic_bound';continue
            with path.open('rb') as stream:raw=stream.read(LIMIT+1)
            if len(raw)>LIMIT:
                row['state']='exceeds_diagnostic_bound';continue
            value=json.loads(raw)
            if not isinstance(value,dict):raise ValueError('not_object')
            row.update(state='read',size=len(raw),sha256=hashlib.sha256(raw).hexdigest())
            schema=value.get('schema')
            row['schema']=schema if schema in SCHEMAS else 'absent_or_unsupported'
            row['classification']=value.get('classification') if value.get('classification') in {'public','controlled'} else 'unavailable'
            if name=='current-routing.json':
                row['report_complete']=value.get('report_complete') is True
                row['synthesis_complete']=value.get('synthesis_complete') is True
            code=value.get('binding',{}).get('code') if isinstance(value.get('binding'),dict) else None
            if isinstance(code,dict):
                row['producer_code_sha256']=hashlib.sha256(json.dumps(code,sort_keys=True,separators=(',',':'),ensure_ascii=True).encode()).hexdigest()
        except FileNotFoundError:row['state']='missing'
        except (OSError,ValueError,TypeError):row['state']='unreadable_or_invalid'
    return dict(schema='genosapiens.track2-report-transition-diagnostic/v1',diagnostic_only=True,
        safe_resume_established=False,per_file_byte_limit=LIMIT,records=records)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--run',type=Path,required=True)
    args=parser.parse_args(argv);print(json.dumps(diagnose(args.run),sort_keys=True,indent=2))


if __name__=='__main__':main()
