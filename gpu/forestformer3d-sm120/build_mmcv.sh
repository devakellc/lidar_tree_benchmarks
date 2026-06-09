#!/usr/bin/env bash
# Clean mmcv 2.0.0 _ext build for torch 2.7 / sm_120. FRESH clone (the earlier
# clone was corrupted by a buggy colon-split sed that turned `s/c10::optional/.../`
# into `s/c10/optional/g`, mangling every `#include <c10/...>`).
# The ONLY real fix needed: torch 2.7 headers require C++17; mmcv 2.0.0 hardcodes
# -std=c++14. c10::optional is still a valid alias in torch 2.7 — leave it alone.
set +e
NJ="$(nproc)"
echo "### MAX_JOBS=$NJ (all cores) ###"
export MMCV_WITH_OPS=1 FORCE_CUDA=1 TORCH_CUDA_ARCH_LIST=12.0 MAX_JOBS="$NJ"
pip uninstall -y mmcv mmcv-lite >/dev/null 2>&1
cd /opt && rm -rf mmcv
git clone --depth 1 --branch v2.0.0 https://github.com/open-mmlab/mmcv.git || { echo "CLONE_FAIL"; exit 1; }
cd mmcv
sed -i 's/c++14/c++17/g' setup.py
echo "remaining c++14 in setup.py: $(grep -c 'c++14' setup.py)"
pip install --no-deps -v .
echo "MMCV_INSTALL_EXIT=$?"
python -c "import mmcv; from mmcv.ops import box_iou_rotated, points_in_boxes_all; print('mmcv', mmcv.__version__, '_ext OK')"
echo "MMCV_IMPORT_EXIT=$?"
