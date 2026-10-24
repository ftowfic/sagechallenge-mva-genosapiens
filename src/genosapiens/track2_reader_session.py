"""Framed local R reader with bounded wall-clock deadlines and explicit failure."""
import json
import os
from pathlib import Path
import selectors
import subprocess
import time

from genosapiens.track2_operator_inputs import require
from genosapiens.track2_scientific_storage import publish_json

PROTOCOL = 'genosapiens.track2-reader-session/v1'
MAX_FRAME = 2 * 1024**2


class ReaderSession:
    def __init__(self, initialization, directory, *, timeout=900, counters=None, command=None):
        require(0 < timeout <= 900, 'Invalid reader request deadline.')
        self.timeout = timeout
        self.counters = counters if counters is not None else {}
        self.process = None
        self.selector = selectors.DefaultSelector()
        self.buffer = bytearray()
        self.serial = 0
        self.closed = False
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True, mode=0o700)
        self.initialization = dict(initialization, protocol=PROTOCOL)
        request = directory/'initialization.json'
        publish_json(request, self.initialization)
        self.stderr = (directory/'reader.stderr').open('wb')
        reader = Path(__file__).resolve().parents[2]/'scripts/read_track2_lincs.R'
        try:
            self.process = subprocess.Popen(command or ['Rscript',str(reader),'--session',str(request)],
                stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=self.stderr, start_new_session=True)
            os.set_blocking(self.process.stdin.fileno(), False)
            self.counters['reader_starts'] = self.counters.get('reader_starts',0)+1
            self.selector.register(self.process.stdout, selectors.EVENT_READ)
            response = self._receive(time.monotonic()+timeout)
            require(response['type']=='ready' and response['request_id']=='initialize' and
                response['matrix_identity']==initialization['matrix_identity'] and response['identifier_loads']==1,
                'The reader initialization contract changed.')
            self.counters['identifier_loads'] = self.counters.get('identifier_loads',0)+1
        except BaseException:
            self.abort()
            raise

    def _send(self, message, deadline):
        frame = json.dumps(dict(message, protocol=PROTOCOL), allow_nan=False).encode()+b'\n'
        require(len(frame)<=MAX_FRAME and not self.closed, 'Invalid reader frame or closed session.')
        with selectors.DefaultSelector() as writable:
            writable.register(self.process.stdin, selectors.EVENT_WRITE)
            pending = memoryview(frame)
            while pending:
                remaining = deadline-time.monotonic()
                if remaining <= 0 or not writable.select(remaining):
                    raise TimeoutError('The matrix reader request exceeded its deadline.')
                try:
                    pending = pending[os.write(self.process.stdin.fileno(),pending):]
                except BlockingIOError:
                    continue

    def _receive(self, deadline):
        while b'\n' not in self.buffer:
            remaining = deadline-time.monotonic()
            if remaining <= 0 or not self.selector.select(remaining):
                raise TimeoutError('The matrix reader request exceeded its deadline.')
            chunk = os.read(self.process.stdout.fileno(), 65536)
            if not chunk:
                raise RuntimeError('The matrix reader exited before request completion.')
            self.buffer.extend(chunk)
            require(len(self.buffer)<=MAX_FRAME, 'The reader response exceeds its frame bound.')
        line, _, rest = self.buffer.partition(b'\n')
        self.buffer = bytearray(rest)
        response = json.loads(line)
        require(isinstance(response,dict) and response.get('protocol')==PROTOCOL and
            isinstance(response.get('request_id'),str), 'Malformed reader response.')
        if response.get('type')=='error':
            raise RuntimeError('The matrix reader reported a terminal error: '+str(response.get('message','')))
        return response

    def extract(self, request):
        self.serial += 1
        request_id = str(self.serial)
        deadline = time.monotonic()+self.timeout
        blocks = []
        try:
            self._send(dict(type='extract', request_id=request_id, timeout_seconds=self.timeout, request=request),deadline)
            while True:
                response = self._receive(deadline)
                require(response['request_id']==request_id, 'The reader response belongs to another request.')
                if response['type']=='complete':
                    require(response['blocks']==len(blocks), 'The reader completion inventory changed.')
                    self.counters['extraction_blocks'] = self.counters.get('extraction_blocks',0)+len(blocks)
                    return blocks
                require(response['type']=='block', 'Unexpected reader response.')
                blocks.append(response['block'])
        except BaseException:
            self.abort()
            raise

    def abort(self):
        if self.closed:
            return
        self.closed = True
        if self.process is not None:
            if self.process.poll() is None:
                self.process.terminate()
                try:
                    self.process.wait(timeout=2)
                except subprocess.TimeoutExpired:
                    self.process.kill()
                    self.process.wait(timeout=2)
            self.process.stdin.close()
            self.process.stdout.close()
        self.selector.close()
        self.stderr.close()

    def close(self):
        if self.closed:
            return
        try:
            deadline = time.monotonic()+min(self.timeout,5)
            self._send(dict(type='shutdown', request_id='shutdown'),deadline)
            response = self._receive(deadline)
            require(response['type']=='closed' and response['request_id']=='shutdown', 'Reader shutdown was not acknowledged.')
            self.process.wait(timeout=2)
            require(self.process.returncode==0, 'The matrix reader failed during shutdown.')
        finally:
            self.abort()

    def __enter__(self):
        return self

    def __exit__(self, kind, value, traceback):
        self.abort() if kind else self.close()
