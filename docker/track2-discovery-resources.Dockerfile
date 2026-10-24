FROM genosapiens/track1-pipeline@sha256:640ba2fa699cc8ed677ed19f57aa39d256f05ff8377cf745f8f12131b95b0b92

LABEL org.opencontainers.image.title="Genosapiens Track 2 public discovery builder" \
      org.opencontainers.image.source="https://github.com/ftowfic/sagechallenge-raredisease-genosapiens" \
      org.opencontainers.image.licenses="CC-BY-4.0" \
      org.genosapiens.role="track2-discovery-resources" \
      org.genosapiens.base.digest="sha256:640ba2fa699cc8ed677ed19f57aa39d256f05ff8377cf745f8f12131b95b0b92"

USER 0:0
COPY requirements/track2-discovery-resources.txt requirements/track2-discovery-resources.txt
COPY requirements/track2-discovery-resources.lock requirements/track2-discovery-resources.lock
COPY requirements/track2-discovery-resources.provenance.json requirements/track2-discovery-resources.provenance.json
RUN python -m pip install --no-cache-dir --require-hashes --only-binary=:all: \
        --requirement requirements/track2-discovery-resources.txt && \
    python -m genosapiens.supply_chain validate \
        --lock requirements/track2-discovery-resources.lock \
        --policy requirements/supply-chain-policy.json \
        --requirements requirements/track2-discovery-resources.txt \
        --role track2-discovery-resources

USER 65532:65532
ENTRYPOINT ["python"]
CMD ["--version"]
