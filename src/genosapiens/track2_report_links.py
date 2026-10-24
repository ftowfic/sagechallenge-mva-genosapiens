"""Check report destinations with one private JSON index per referenced file."""
from contextlib import ExitStack
from pathlib import Path
import re

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_report_model import relative_destination
from genosapiens.track2_scientific_reader import read_scientific
from genosapiens.track2_scientific_storage import file_identity


def check_links(root, text):
    root=Path(root)
    identities={};documents={};sections={}
    with ExitStack() as owners:
        for link in re.findall(r'\]\(([^)]+)\)',text):
            relative,_,fragment=relative_destination(link).partition('#')
            path=root/relative
            require(path.is_file(), 'A challenge report link is missing.')
            if not fragment:continue
            current=file_identity(path)
            require(identities.setdefault(path,current)==current,
                    'A linked report artifact changed during validation.')
            if fragment.startswith('/'):
                if path not in documents:
                    value=read_scientific(path)
                    document=getattr(value,'document',None)
                    if document is not None:owners.callback(document.close)
                    documents[path]=value
                value=documents[path]
                try:
                    for part in fragment[1:].split('/'):
                        part=part.replace('~1','/').replace('~0','~')
                        value=value[int(part)] if not isinstance(value,dict) else value[part]
                except (KeyError,IndexError,ValueError,TypeError):
                    raise ValueError('A challenge report record reference is missing.') from None
            else:
                if path not in sections:
                    contents=path.read_text()
                    headings={re.sub(r'[^\w\- ]','',line.lstrip('#').strip().lower()).replace(' ','-')
                        for line in contents.splitlines() if line.startswith('#')}
                    sections[path]=contents,headings
                contents,headings=sections[path]
                require(fragment in headings or 'id="'+fragment+'"' in contents,
                        'A challenge report section reference is missing.')
        for path,expected in identities.items():
            require(file_identity(path)==expected, 'A linked report artifact changed during validation.')
