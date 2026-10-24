"""Invocation-local page proofs bound to source bytes and generator dependencies."""
import itertools
import json
from pathlib import Path

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_scientific_storage import file_identity, json_identity, publish_chunks, json_chunks


def generator_identity():
    from genosapiens.track2_construction_dependencies import python_dependencies
    root=Path(__file__).resolve().parents[2]
    return python_dependencies(root,[('genosapiens.track2_figure_stream','figure_pages'),
        ('genosapiens.track2_tissue_figures','tissue_condition_pages')])


class CheckedPageStreams:
    """Cache only digests. Retained receipts never seed a semantic proof."""
    def __init__(self):
        self._proofs={};self.closed=False;self.generations=0

    def close(self):
        self._proofs.clear();self.closed=True

    def _key(self,query,result,conditions):
        require(not self.closed,'The checked page invocation has ended.')
        return json_identity(dict(query=query,result=result,conditions=file_identity(conditions),
            generator=generator_identity()))['sha256']

    def publish(self,folder,query,result,conditions):
        from genosapiens.track2_figure_stream import figure_pages
        before=self._key(query,result,conditions)
        path=Path(folder)/'figure-data.jsonl'
        publish_chunks(path,(chunk for page in figure_pages(query,result,conditions) for chunk in json_chunks(page)))
        require(self._key(query,result,conditions)==before,'Figure sources changed during page generation.')
        self._proofs[before]=file_identity(path)['sha256'];self.generations+=1

    def check(self,folder,query,result,conditions):
        from genosapiens.track2_figure_stream import figure_pages
        before=self._key(query,result,conditions)
        folder=Path(folder);path=folder/'figure-data.jsonl';digest=file_identity(path)['sha256']
        if self._proofs.get(before)!=digest:
            with path.open() as stream:
                for expected,line in itertools.zip_longest(figure_pages(query,result,conditions),stream):
                    require(expected is not None and line is not None and json.loads(line)==expected,
                        'Native figure data differs from its scientific evidence.')
            require(self._key(query,result,conditions)==before and file_identity(path)['sha256']==digest,
                    'Figure evidence changed during page admission.')
            self._proofs[before]=digest;self.generations+=1
        # Assets require fresh checks even when page semantics are reusable.
        with path.open() as stream:
            for line in stream:
                page=json.loads(line)
                identifier=page['id']
                require(isinstance(identifier,str) and identifier not in {'.','..'}
                    and '/' not in identifier and '\\' not in identifier,'A figure identifier is not local.')
                for suffix in ('.png','.svg'):
                    asset=folder/(identifier+suffix)
                    require(asset.is_file() and file_identity(asset)['size']>0,'A required native figure is missing.')
