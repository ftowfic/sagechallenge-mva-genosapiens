FROM python:3.12.11-slim-bookworm@sha256:519591d6871b7bc437060736b9f7456b8731f1499a57e22e6c285135ae657bf7

ARG VCS_REF=unknown
LABEL org.opencontainers.image.title="Genosapiens Track 2 public research" \
      org.opencontainers.image.revision="${VCS_REF}" \
      org.opencontainers.image.licenses="CC-BY-4.0" \
      org.genosapiens.role="track2-pipeline"
ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1 PYTHONPATH=/opt/genosapiens/src
WORKDIR /opt/genosapiens
COPY requirements/runtime.txt requirements/runtime.txt
RUN python -m pip install --no-cache-dir --require-hashes --only-binary=:all: -r requirements/runtime.txt
COPY requirements/supply-chain-policy.json requirements/supply-chain-policy.json
COPY requirements/track2-pipeline.lock requirements/runtime.lock
COPY src/ src/
COPY config/track2/ config/track2/
COPY schemas/track2/ schemas/track2/
RUN python -m genosapiens.supply_chain validate \
    --lock requirements/runtime.lock --policy requirements/supply-chain-policy.json \
    --requirements requirements/runtime.txt --role track2-pipeline
USER 65532:65532
ENTRYPOINT ["python", "-m", "genosapiens.track2_cli"]
CMD ["--help"]
