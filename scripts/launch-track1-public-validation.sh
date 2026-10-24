#!/bin/sh
# Outer evaluator owns the Docker socket; scientific case containers never do.
set -eu
project_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
public_reference_root=/mnt/reference
public_work_root=${TRACK1_PUBLIC_SCRATCH_ROOT:-}
if [ -z "$public_work_root" ]; then
    public_work_root=$(mktemp -d /mnt/reference/scratch/public-build/track1-public-20260915.XXXXXXXX)
fi
case "$public_work_root" in
    /mnt/reference/scratch/public-build/*) ;;
    *) echo 'invalid public validation scratch root' >&2; exit 1 ;;
esac
[ -d "$public_work_root" ]
[ ! -L "$public_work_root" ]
public_release_root=/mnt/reference/manifests/track1/base-releases
mkdir -p "$public_release_root"
[ ! -L "$public_release_root" ]
[ "$#" -eq 0 ] || { echo 'use TRACK1_PUBLIC_SCRATCH_ROOT and TRACK1_PUBLIC_RESUME_OUTCOME_SHA256 for recovery' >&2; exit 1; }
# Sources and targets have the same absolute paths, so the evaluator can submit
# host bind mounts without a Docker socket or path-translation layer in workers.
mount_arguments_file=$(mktemp)
trap 'rm -f "$mount_arguments_file"' EXIT HUP INT TERM
# Use newline arguments through a shell function, avoiding eval/string commands.
launch_evaluator() {
    set -- "$@" \
        sha256:640ba2fa699cc8ed677ed19f57aa39d256f05ff8377cf745f8f12131b95b0b92 \
        "$project_root/scripts/run_track1_public_validation.py" \
        --repository-root "$project_root" --reference-root "$public_reference_root" \
        --scratch-root "$public_work_root" --output "$public_work_root/validation.json"
    if [ -n "${TRACK1_PUBLIC_SMOKE_CASE_ID:-}" ]; then
        [ -z "${TRACK1_PUBLIC_RESUME_OUTCOME_SHA256:-}" ] || { echo 'preparation smoke cannot resume an outcome' >&2; exit 1; }
        set -- "$@" --case-id "$TRACK1_PUBLIC_SMOKE_CASE_ID" --workers 1
    else
        set -- "$@" --workers "${TRACK1_PUBLIC_WORKERS:-2}"
    fi
    if [ -n "${TRACK1_PUBLIC_RESUME_OUTCOME_SHA256:-}" ]; then
        set -- "$@" --resume-outcome-sha256 "$TRACK1_PUBLIC_RESUME_OUTCOME_SHA256"
    fi
    if [ -n "${TRACK1_PUBLIC_CONTINUATION:-}" ]; then
        [ -z "${TRACK1_PUBLIC_RESUME_OUTCOME_SHA256:-}" ]
        [ -z "${TRACK1_PUBLIC_SMOKE_CASE_ID:-}" ]
        set -- "$@" --continuation "$TRACK1_PUBLIC_CONTINUATION"
    fi
    docker run --rm --network none --read-only --entrypoint python \
        --user "$(id -u):$(id -g)" \
        --group-add "$(stat -c %g /var/run/docker.sock)" \
        --cap-drop ALL --security-opt no-new-privileges \
        --memory 8g --memory-swap 8g --cpus 4 --pids-limit 128 \
        --tmpfs /tmp:rw,nosuid,nodev,size=64m \
        --env "PYTHONPATH=$project_root/src" \
        --env PYTHONDONTWRITEBYTECODE=1 --env PYTHONUNBUFFERED=1 \
        --volume /var/run/docker.sock:/var/run/docker.sock \
        --volume "$project_root:$project_root:ro" \
        --volume "$public_work_root:$public_work_root" \
        --volume "$public_release_root:$public_release_root" \
        "$@"
}
# Public-only resource trees; synchronized with PUBLIC_REFERENCE_TREES in the
# worker module. No whole reference root, controlled tree or remote credential.
public_mounts='releases/model-loci/grc-ensembl-nhgri/grch38-ensembl-116-nhgri-retrieved-2026-09-16
derived/model-loci/22b69fc15bbbb97dfeac76209e57be0962270ca6a1eaa8011151c9daa2b9fcaf
manifests/model-loci
releases/clinvar/2026-08-22
derived/gnomad/4.1/joint-frequency/gnomad-joint-v4.1-frequency-projection-v1
derived/gnomad/4.1/all-sites/gnomad-v4.1-all-sites-an-v1
manifests/clinvar
manifests/gnomad/4.1
releases/grch38/GCA_000001405.15-masked-v2-plus-hs38d1-current
releases/ensembl/vep/116_GRCh38
releases/open-targets/26.06
releases/orphadata/2026-06-23
releases/gencc/2026-08-30
releases/ontologies/core-ontology-references-v1
releases/repeats/broad-str-analysis/d41c85c97da3cd8c8ac9d7d47883570c95069874
releases/challenge-evaluator/SageBio-rare-disease-real-kid-mva-hackathon-2026/1c761cc23d90aebe6a011fd5b0b99517df42408c
derived/phenotype/hpo-ic/v1/060cf7c49671f9a98f9a8578c35f0fb8eeae62166d79cb40c2bcc81855abcc89
derived/rare-disease-evidence/orphadata-2026-06-23_gencc-2026-08-30_ot-26.06_v6
manifests/grch38
manifests/ensembl/vep
manifests/open-targets
manifests/rare-disease-evidence
manifests/ontologies
manifests/phenotype
manifests/variant-class-resources
manifests/common-allele-public-indexes
derived/common-allele-public-indexes
manifests/track1-common-allele-exceptions
derived/track1-common-allele-exceptions'
# Positional parameters contain only generated volume arguments.
set --
if [ -n "${TRACK1_PUBLIC_CONTINUATION:-}" ]; then
    case "${TRACK1_PUBLIC_PRIOR_SCRATCH_ROOT:-}" in
        /mnt/reference/scratch/public-build/*/acceptance) ;;
        *) echo 'invalid prior public validation scratch root' >&2; exit 1 ;;
    esac
    [ -d "$TRACK1_PUBLIC_PRIOR_SCRATCH_ROOT" ]
    [ ! -L "$TRACK1_PUBLIC_PRIOR_SCRATCH_ROOT" ]
    set -- "$@" --volume "$TRACK1_PUBLIC_PRIOR_SCRATCH_ROOT:$TRACK1_PUBLIC_PRIOR_SCRATCH_ROOT:ro"
fi
while IFS= read -r relative; do
    public_path=$public_reference_root/$relative
    [ -d "$public_path" ]
    [ ! -L "$public_path" ]
    printf '%s\n' "$public_path" >> "$mount_arguments_file"
done <<EOF
$public_mounts
EOF
while IFS= read -r public_path; do
    if [ "$public_path" = "$public_reference_root/derived/common-allele-public-indexes" ]; then
        set -- "$@" --volume \
            "$public_reference_root/derived/common-allele-public-indexes-cache-bounded-v1:$public_path:ro"
    else
        set -- "$@" --volume "$public_path:$public_path:ro"
    fi
done < "$mount_arguments_file"
# The evaluator must see the host source path used by its isolated worker spec.
set -- "$@" --volume \
    "$public_reference_root/derived/common-allele-public-indexes-cache-bounded-v1:$public_reference_root/derived/common-allele-public-indexes-cache-bounded-v1:ro"
launch_evaluator "$@"
