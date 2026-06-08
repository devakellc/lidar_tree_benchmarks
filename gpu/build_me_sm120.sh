#!/usr/bin/env bash
# Build the MinkowskiEngine sm_120 base image (prerequisite for SAT / FF3D ME layers).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
TAG="${1:-me-sm120-test}"
docker build -t "$TAG" "$HERE/minkowski-sm120"
echo "Built $TAG — smoke: docker run --rm --gpus all $TAG"
