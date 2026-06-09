#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ME_BASE="${ME_BASE:-me-sm120-test}"
TAG="${1:-sat-sm120-test}"
# Pin the SAT clone to gpu/SegmentAnyTree.commit unless SAT_REF is set; the
# Dockerfile defaults to `main` if neither is given.
PIN="$(cd "$HERE/.." && pwd)/SegmentAnyTree.commit"
SAT_REF="${SAT_REF:-$(grep -m1 -vE '^[[:space:]]*#' "$PIN" 2>/dev/null | tr -d '[:space:]')}"
ARGS=(--build-arg ME_BASE="$ME_BASE")
[ -n "$SAT_REF" ] && ARGS+=(--build-arg SAT_REF="$SAT_REF")
docker build -t "$TAG" "${ARGS[@]}" "$HERE"
echo "Built $TAG (SAT_REF=${SAT_REF:-main}) — smoke: docker run --rm --gpus all $TAG"
