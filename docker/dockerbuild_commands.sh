#!/usr/bin/env bash
#
# Build, tag, and (optionally) push the per-module mg-proc Docker images.
#
# Usage (run from the repository root so the build context '.' is correct):
#   bash docker/dockerbuild_commands.sh                  # build + tag :latest locally
#   VERSION=v1.0.0 bash docker/dockerbuild_commands.sh   # also tag :v1.0.0
#   PUSH=1 bash docker/dockerbuild_commands.sh           # build, tag, push :latest to ghcr.io
#   PUSH=1 VERSION=v1.0.0 bash docker/dockerbuild_commands.sh  # push :latest and :v1.0.0
#
# Pushing requires being logged in to the registry first, e.g.:
#   echo "$GHCR_PAT" | docker login ghcr.io -u <github-user> --password-stdin
# Once pushed (and the packages made public), other machines pull the images
# automatically on the first `nextflow run` — no local rebuild needed.

set -euo pipefail

SCRIPT_DIR=$(dirname "$(readlink -f "$0")")

REGISTRY="ghcr.io/pereiramemo/mg-proc"
VERSION="${VERSION:-}"   # optional extra tag (e.g. v1.0.0); empty = only :latest
PUSH="${PUSH:-0}"        # set to 1 to push images after building

# Module base names — must match docker/<module>.Dockerfile and the image
# names referenced in modules/<module>.nf.
MODULES=(
  1.1-quality-check
  1.2-quality-check
  2-preprocess
  3-assembly-and-map
)

for module in "${MODULES[@]}"; do
  image="${REGISTRY}/${module}"

  echo ">>> Building ${image}:latest"
  docker build --network=host \
    -f "${SCRIPT_DIR}/${module}.Dockerfile" \
    -t "${image}:latest" \
    .

  if [[ -n "${VERSION}" ]]; then
    echo ">>> Tagging ${image}:${VERSION}"
    docker tag "${image}:latest" "${image}:${VERSION}"
  fi

  if [[ "${PUSH}" == "1" ]]; then
    echo ">>> Pushing ${image}:latest"
    docker push "${image}:latest"
    if [[ -n "${VERSION}" ]]; then
      echo ">>> Pushing ${image}:${VERSION}"
      docker push "${image}:${VERSION}"
    fi
  fi
done

echo "Done."
