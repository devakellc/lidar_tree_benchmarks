#!/bin/bash
# Container entry for the FF3D arm. run_docker_arm passes NO env into the
# container, so every path arrives positionally (input + out first, then the
# `extra` args) and is identity-mounted. Invoked as:
#   bash ff3d_entry.sh <in_dir> <out_laz> <ckpt> <repo> <patch> <driver>
# cd into the mounted repo, apply ff3d_repo.patch idempotently, run the driver
# with just (in_dir, out_laz, ckpt).
set -euo pipefail
in_dir="$1"; out_laz="$2"; ckpt="$3"; repo="$4"; patch="$5"; driver="$6"
cd "$repo"
if git apply --check "$patch" 2>/dev/null; then
  git apply "$patch"
fi
exec python "$driver" "$in_dir" "$out_laz" "$ckpt"
