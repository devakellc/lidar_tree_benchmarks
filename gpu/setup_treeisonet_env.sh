#!/usr/bin/env bash
# #I3a: headless TreeisoNet env. Pinned upstream commit + pinned cu128 torch
# (Blackwell sm_120). Pure PyTorch + timm (no compiled ops) -> the wheel swap is
# the only Blackwell change. Run from the repo root: bash gpu/setup_treeisonet_env.sh
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
COMMITF="$ROOT/TreeAIBox.commit"            # versioned, reviewed pin (see #I5/#I3)
PIN="$(cat "$COMMITF")"
if [ ! -d "$ROOT/TreeAIBox" ]; then
  git clone https://github.com/NRCan/TreeAIBox "$ROOT/TreeAIBox"
fi
git -C "$ROOT/TreeAIBox" fetch --depth 1 origin "$PIN" 2>/dev/null || true
git -C "$ROOT/TreeAIBox" checkout -q "$PIN"
echo "TreeAIBox pinned at $PIN"
# Pinned cu128 (sm_120) torch + the headless deps. timm<1.0 keeps models.layers;
# pyproj makes laspy's CRS preflight reliable. uv creates the venv without the
# system python3-venv/ensurepip package (and installs fast).
uv venv "$ROOT/.venv" --python 3.12
. "$ROOT/.venv/bin/activate"
uv pip install torch==2.7.0 torchvision==0.22.0 --index-url https://download.pytorch.org/whl/cu128
uv pip install "timm==0.9.16" "laspy[lazrs]" numpy numpy_indexed numpy_groupies \
               scipy scikit-image commentjson einops pyproj
# Copy the ALS treeloc/treeoff config JSONs beside the weights (glob = robust to
# the upstream paren-vs-underscore naming).
cp "$ROOT"/TreeAIBox/modules/treeisonet/*reclamation*treeloc*.json "$ROOT/store/treeaibox/" 2>/dev/null || true
cp "$ROOT"/TreeAIBox/modules/treeisonet/*reclamation*treeoff*.json "$ROOT/store/treeaibox/" 2>/dev/null || true
python - <<'PY'
import torch
print("torch", torch.__version__, "cuda", torch.version.cuda, "avail", torch.cuda.is_available())
print("device", torch.cuda.get_device_name(0) if torch.cuda.is_available() else "CPU")
print("arch", torch.cuda.get_arch_list() if torch.cuda.is_available() else [])
PY
echo "--- configs copied ---"; ls "$ROOT/store/treeaibox/"*.json 2>/dev/null || echo "NO configs copied — check repo layout"
