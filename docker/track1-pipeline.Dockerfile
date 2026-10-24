FROM python:3.12.11-slim-bookworm@sha256:519591d6871b7bc437060736b9f7456b8731f1499a57e22e6c285135ae657bf7 AS python-runtime

FROM python-runtime AS bcftools-packages
ADD --checksum=sha256:905a0cc21cb47f93e2d96de4eec80799a67b5bb6032362940f97cd353168058e \
    https://archive.ubuntu.com/ubuntu/pool/main/c/curl/libcurl3-gnutls_7.81.0-1ubuntu1.27_amd64.deb \
    /packages/libcurl3-gnutls_7.81.0-1ubuntu1.27_amd64.deb
ADD --checksum=sha256:a8fbeb21db6a42dbcc4402235c6aff4a0281ed9fd3004b336d7259318fefe745 \
    https://archive.ubuntu.com/ubuntu/pool/universe/h/htscodecs/libhtscodecs2_1.1.1-3_amd64.deb \
    /packages/libhtscodecs2_1.1.1-3_amd64.deb
ADD --checksum=sha256:8e9d1437f218af8c80f851e63de1e19bfbd96741fd78ab442418cf4b0b997420 \
    https://archive.ubuntu.com/ubuntu/pool/universe/h/htslib/libhts3_1.13%2bds-2build1_amd64.deb \
    /packages/libhts3_1.13+ds-2build1_amd64.deb
ADD --checksum=sha256:5ab4909fc958af5a4871b309ac3affd4bd82e7d47a048bc575e93f991c148004 \
    https://archive.ubuntu.com/ubuntu/pool/universe/b/bcftools/bcftools_1.13-1_amd64.deb \
    /packages/bcftools_1.13-1_amd64.deb

FROM ensemblorg/ensembl-vep@sha256:f22d1f89b28c6c4caf043ddb2269fdfff8940364da589cadbbabc6212a436949

ARG IMAGE_VERSION=0.1.0
ARG VCS_REF=unknown
ARG DEPENDENCY_LOCK=requirements/track1-pipeline.lock

LABEL org.opencontainers.image.title="Genosapiens Track 1 pipeline" \
      org.opencontainers.image.description="Offline supplied-call annotation, ranking, and preflight" \
      org.opencontainers.image.version="${IMAGE_VERSION}" \
      org.opencontainers.image.revision="${VCS_REF}" \
      org.opencontainers.image.source="https://github.com/ftowfic/sagechallenge-raredisease-genosapiens" \
      org.opencontainers.image.licenses="CC-BY-4.0" \
      org.genosapiens.role="track1-pipeline" \
      org.genosapiens.base.digest="sha256:f22d1f89b28c6c4caf043ddb2269fdfff8940364da589cadbbabc6212a436949" \
      org.genosapiens.python.base.digest="sha256:519591d6871b7bc437060736b9f7456b8731f1499a57e22e6c285135ae657bf7" \
      org.genosapiens.vep.version="116.2" \
      org.genosapiens.bcftools.version="1.13" \
      org.genosapiens.htslib.version="1.13+ds"

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PYTHONPATH=/opt/genosapiens/src

WORKDIR /opt/genosapiens

USER 0:0

COPY --from=python-runtime /usr/local/ /usr/local/
COPY --from=bcftools-packages /packages/ /tmp/bcftools-packages/

RUN dpkg --install \
    /tmp/bcftools-packages/libcurl3-gnutls_7.81.0-1ubuntu1.27_amd64.deb \
    /tmp/bcftools-packages/libhtscodecs2_1.1.1-3_amd64.deb \
    /tmp/bcftools-packages/libhts3_1.13+ds-2build1_amd64.deb \
    /tmp/bcftools-packages/bcftools_1.13-1_amd64.deb && \
    test "$(bcftools --version-only)" = "1.13+htslib-1.13+ds" && \
    tabix --version && bgzip --version && \
    rm -rf /tmp/bcftools-packages

COPY requirements/runtime.txt requirements/runtime.txt
RUN python -m pip install --no-cache-dir --require-hashes --only-binary=:all: \
    --requirement requirements/runtime.txt

COPY requirements/supply-chain-policy.json requirements/supply-chain-policy.json
COPY ${DEPENDENCY_LOCK} requirements/runtime.lock
COPY config/ config/
COPY src/ src/

RUN python -m genosapiens.supply_chain validate \
    --lock requirements/runtime.lock \
    --policy requirements/supply-chain-policy.json \
    --requirements requirements/runtime.txt \
    --role track1-pipeline

USER 65532:65532

ENTRYPOINT ["python", "-m", "genosapiens", "track1"]
CMD ["--help"]
