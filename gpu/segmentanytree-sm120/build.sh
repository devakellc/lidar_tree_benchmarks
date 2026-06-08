#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ME_BASE="${ME_BASE:-me-sm120-test}"
TAG="${1:-sat-sm120-test}"
docker build -t "$TAG" --build-arg ME_BASE="$ME_BASE" "$HERE"
echo "Built $TAG — smoke: docker run --rm --gpus all $TAG"
