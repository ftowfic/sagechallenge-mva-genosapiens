"""Owned document views for one report helper chain."""
from pathlib import Path

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_scientific_reader import read_scientific
from genosapiens.track2_scientific_storage import json_identity


class _LocalModel:
    """Own opened records until the dependent helpers finish, including failures.

    Borrowers use the same evidence root and manifest. They never close views.
    Path checks reject replacement even when an open descriptor remains valid.
    This owner provides no authority across copies or independent invocations.
    """
    def __init__(self, root, manifest):
        self.root=Path(root).absolute()
        self.manifest=json_identity(manifest)
        self.documents=[]
        self.closed=False
        from genosapiens.track2_report_evidence import scientific_model
        try:
            self.model=scientific_model(self.root,manifest,_read=self._read)
            self.borrow(root,manifest)
        except BaseException:
            self.close()
            raise

    def _read(self, path):
        value=read_scientific(path)
        if hasattr(value,'document'):
            self.documents.append((Path(path).absolute(),value.document))
        return value

    def borrow(self, root, manifest):
        require(not self.closed and Path(root).absolute()==self.root and json_identity(manifest)==self.manifest,
                'The borrowed scientific model has another owner or evidence binding.')
        for path,document in self.documents:
            require(path.resolve()==path, 'A borrowed scientific document path changed.')
            info=path.stat()
            signature=(info.st_dev,info.st_ino,info.st_size,info.st_mtime_ns,info.st_ctime_ns)
            require(signature==document.signature==document._signature(),
                    'A borrowed scientific document changed.')
        return self.model

    def __enter__(self):return self

    def __exit__(self,*error):self.close()

    def close(self):
        self.closed=True
        for _,document in self.documents:document.close()
        self.documents.clear()
