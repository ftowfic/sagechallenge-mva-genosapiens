FROM python:3.12.11-slim-bookworm@sha256:519591d6871b7bc437060736b9f7456b8731f1499a57e22e6c285135ae657bf7

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PYTHONPATH=/workspace/src

WORKDIR /workspace

RUN apt-get update \
    && apt-get install --yes --no-install-recommends git=1:2.39.5-0+deb12u3 \
    && rm -rf /var/lib/apt/lists/*

COPY requirements/runtime.txt requirements/runtime.txt
RUN python -m pip install --no-cache-dir --require-hashes --only-binary=:all: \
    --requirement requirements/runtime.txt

COPY --chown=65532:65532 . .
RUN chown 65532:65532 /workspace

USER 65532:65532

CMD ["python", "-m", "genosapiens", "--help"]
