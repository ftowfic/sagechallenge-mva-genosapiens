FROM python:3.12.11-slim-bookworm@sha256:519591d6871b7bc437060736b9f7456b8731f1499a57e22e6c285135ae657bf7 AS python-runtime

FROM python-runtime AS gcs-connector
ADD --checksum=sha256:340ee851d2abf55e428f4553b51882a2917cef14dd0fea397079d2c937b4e3c6 \
    https://repo1.maven.org/maven2/com/google/cloud/bigdataoss/gcs-connector/hadoop3-2.2.24/gcs-connector-hadoop3-2.2.24-shaded.jar \
    /gcs-connector-hadoop3-2.2.24-shaded.jar

FROM python-runtime AS bwa-mem2-stock
ADD --checksum=sha256:b4cfdbce8cc07cdf3f6a920facabc29c976cf77dd53573369508111d6d1c555b \
    https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 \
    /tmp/bwa-mem2.tar.bz2
RUN python -c 'import tarfile; tarfile.open("/tmp/bwa-mem2.tar.bz2").extractall("/opt", filter="data")' && \
    mv /opt/bwa-mem2-2.2.1_x64-linux /opt/bwa-mem2 && \
    /opt/bwa-mem2/bwa-mem2 version

FROM ensemblorg/ensembl-vep@sha256:f22d1f89b28c6c4caf043ddb2269fdfff8940364da589cadbbabc6212a436949 AS bwa-mem2
COPY --from=python-runtime /usr/local/ /usr/local/
USER 0:0
ADD --checksum=sha256:4ddaeb704f0e0acbb57c8b1565ddd10abff64c28f83a6ceee8601628b5ae392e \
    https://github.com/bwa-mem2/bwa-mem2/archive/3c433d5cb6bac5c87e2b20279c090d0739da7210.tar.gz \
    /tmp/bwa-mem2-source.tar.gz
ADD --checksum=sha256:6dea4f4accdc14eea706c88bb11725e8c2cd69b0d8c5f5ee19a0669d1883a901 \
    https://github.com/intel/safestringlib/archive/245c4b8cff1d2e7338b7f3a82828fc8e72b29549.tar.gz \
    /tmp/safestringlib-source.tar.gz
RUN python -c 'import tarfile; [tarfile.open(path).extractall("/opt", filter="data") for path in ("/tmp/bwa-mem2-source.tar.gz", "/tmp/safestringlib-source.tar.gz")]' && \
    rmdir /opt/bwa-mem2-3c433d5cb6bac5c87e2b20279c090d0739da7210/ext/safestringlib && \
    mv /opt/safestringlib-245c4b8cff1d2e7338b7f3a82828fc8e72b29549 \
       /opt/bwa-mem2-3c433d5cb6bac5c87e2b20279c090d0739da7210/ext/safestringlib && \
    mv /opt/bwa-mem2-3c433d5cb6bac5c87e2b20279c090d0739da7210 /opt/bwa-mem2-source
RUN make --directory=/opt/bwa-mem2-source CXX=g++ portable=1 multi && \
    mkdir /opt/bwa-mem2 && \
    install --mode=0755 \
      /opt/bwa-mem2-source/bwa-mem2 \
      /opt/bwa-mem2-source/bwa-mem2.sse41 \
      /opt/bwa-mem2-source/bwa-mem2.sse42 \
      /opt/bwa-mem2-source/bwa-mem2.avx \
      /opt/bwa-mem2-source/bwa-mem2.avx2 \
      /opt/bwa-mem2-source/bwa-mem2.avx512bw \
      /opt/bwa-mem2/ && \
    install --mode=0444 /opt/bwa-mem2-source/LICENSE /opt/bwa-mem2/LICENSE && \
    /opt/bwa-mem2/bwa-mem2 version
COPY --from=bwa-mem2-stock /opt/bwa-mem2/ /opt/bwa-mem2-stock/
COPY tests/container/verify-bwa-mem2.py /tmp/verify-bwa-mem2.py
RUN python /tmp/verify-bwa-mem2.py \
    --stock /opt/bwa-mem2-stock/bwa-mem2 \
    --candidate /opt/bwa-mem2/bwa-mem2

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

FROM eclipse-temurin@sha256:fc451894669bc656f81082ced9daddfd8df0d8e8fea659dc26a232f0e0124837 AS java-runtime

FROM ensemblorg/ensembl-vep@sha256:f22d1f89b28c6c4caf043ddb2269fdfff8940364da589cadbbabc6212a436949

ARG IMAGE_VERSION=0.1.0
ARG VCS_REF=unknown
ARG DEPENDENCY_LOCK=requirements/reference-builder.lock

LABEL org.opencontainers.image.title="Genosapiens reference builder" \
      org.opencontainers.image.description="Public reference acquisition and validation" \
      org.opencontainers.image.version="${IMAGE_VERSION}" \
      org.opencontainers.image.revision="${VCS_REF}" \
      org.opencontainers.image.source="https://github.com/ftowfic/sagechallenge-raredisease-genosapiens" \
      org.opencontainers.image.licenses="CC-BY-4.0" \
      org.genosapiens.role="reference-builder" \
      org.genosapiens.base.digest="sha256:f22d1f89b28c6c4caf043ddb2269fdfff8940364da589cadbbabc6212a436949" \
      org.genosapiens.python.base.digest="sha256:519591d6871b7bc437060736b9f7456b8731f1499a57e22e6c285135ae657bf7" \
      org.genosapiens.java.base.digest="sha256:fc451894669bc656f81082ced9daddfd8df0d8e8fea659dc26a232f0e0124837" \
      org.genosapiens.hail.version="0.2.139" \
      org.genosapiens.spark.version="3.5.7" \
      org.genosapiens.gcs-connector.version="hadoop3-2.2.24" \
      org.genosapiens.bcftools.version="1.13" \
      org.genosapiens.htslib.version="1.13+ds" \
      org.genosapiens.bwa-mem2.version="2.2.1" \
      org.genosapiens.bwa-mem2.revision="3c433d5cb6bac5c87e2b20279c090d0739da7210" \
      org.genosapiens.bwa-mem2.source.sha256="4ddaeb704f0e0acbb57c8b1565ddd10abff64c28f83a6ceee8601628b5ae392e" \
      org.genosapiens.bwa-mem2.safestringlib.revision="245c4b8cff1d2e7338b7f3a82828fc8e72b29549" \
      org.genosapiens.java.version="11.0.28+6" \
      org.genosapiens.vep.version="116.2"

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PYTHONPATH=/opt/genosapiens/src \
    JAVA_HOME=/opt/java/openjdk \
    PATH=/opt/bwa-mem2:/opt/java/openjdk/bin:${PATH}

WORKDIR /opt/genosapiens

USER 0:0

COPY --from=python-runtime /usr/local/ /usr/local/
COPY --from=java-runtime /opt/java/openjdk/ /opt/java/openjdk/
COPY --from=bwa-mem2 /opt/bwa-mem2/ /opt/bwa-mem2/
COPY --from=bcftools-packages /packages/ /tmp/bcftools-packages/

RUN groupadd --gid 65532 genosapiens && \
    useradd --uid 65532 --gid 65532 --home-dir /nonexistent \
    --shell /usr/sbin/nologin --no-create-home genosapiens

RUN dpkg --install \
    /tmp/bcftools-packages/libcurl3-gnutls_7.81.0-1ubuntu1.27_amd64.deb \
    /tmp/bcftools-packages/libhtscodecs2_1.1.1-3_amd64.deb \
    /tmp/bcftools-packages/libhts3_1.13+ds-2build1_amd64.deb \
    /tmp/bcftools-packages/bcftools_1.13-1_amd64.deb && \
    test "$(bcftools --version-only)" = "1.13+htslib-1.13+ds" && \
    rm -rf /tmp/bcftools-packages

COPY requirements/reference-builder-build.txt requirements/reference-builder-build.txt
COPY requirements/reference-builder.txt requirements/reference-builder.txt
RUN python -m pip install --no-cache-dir --require-hashes --only-binary=:all: \
    --requirement requirements/reference-builder-build.txt && \
    python -m pip install --no-cache-dir --require-hashes \
    --only-binary=:all: --no-binary=pyspark --no-build-isolation \
    --requirement requirements/reference-builder.txt
COPY --chmod=0444 --from=gcs-connector /gcs-connector-hadoop3-2.2.24-shaded.jar \
    /usr/local/lib/python3.12/site-packages/pyspark/jars/gcs-connector-hadoop3-2.2.24-shaded.jar

# Compose defaults to the host UID/GID so mounted reference files retain host ownership.
RUN groupadd --gid 1000 genosapiens-host && \
    useradd --uid 1000 --gid 1000 --home-dir /nonexistent \
    --shell /usr/sbin/nologin --no-create-home genosapiens-host

COPY requirements/supply-chain-policy.json requirements/supply-chain-policy.json
COPY ${DEPENDENCY_LOCK} requirements/runtime.lock
COPY config/ config/
COPY src/ src/

RUN python -m genosapiens.supply_chain validate \
    --lock requirements/runtime.lock \
    --policy requirements/supply-chain-policy.json \
    --requirements requirements/reference-builder.txt \
    --build-requirements requirements/reference-builder-build.txt \
    --role reference-builder

USER 65532:65532

ENTRYPOINT ["python", "-m", "genosapiens", "reference"]
CMD ["--help"]
