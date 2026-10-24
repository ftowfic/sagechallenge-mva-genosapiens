# The installation declaration records availability of this exact base.
ARG LINCS_IMAGE=genosapiens/track2-signaturesearch@sha256:98ae55f3af076e7d392beb763cf6dabee06680d11dfd8fa82549cf1bd0fe2a40
FROM ${LINCS_IMAGE}
ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1 PYTHONPATH=/opt/genosapiens/src
WORKDIR /opt/genosapiens
COPY src/ src/
COPY config/track2/ config/track2/
COPY schemas/track2/ schemas/track2/
COPY scripts/read_track2_lincs.R scripts/render_track2_lincs.R scripts/render_track2_overviews.R scripts/
USER 1000:1000
ENTRYPOINT ["python3", "-m", "genosapiens.track2_cli"]
CMD ["--help"]
