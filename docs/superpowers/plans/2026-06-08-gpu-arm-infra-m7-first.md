# GPU-Arm Infra (#I5/#I3/#I4), TreeisoNet-first Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> superpowers:subagent-driven-development (recommended) or
> superpowers:executing-plans to implement this plan task-by-task. Steps use
> checkbox (`- [ ]`) syntax for tracking.

**Goal:** Stand up the GPU-arm infrastructure — weights/image mirror (#I5),
model-runner contract (#I3), and the I/O + CRS bridge (#I4) — and prove it
end-to-end with **TreeisoNet (#M7) as the first arm**, because TreeisoNet is the
only model that runs Blackwell-clean (pure PyTorch, no MinkowskiEngine/spconv).
The acceptance gate is one TreeisoNet inference on a real SOAP clip, reduced to
apex detections, **converted to height-above-ground via the frozen clip's DTM**,
and scored by the existing harness.

**Architecture:** TreeisoNet runs **headless in a Python venv** (cu128 torch),
not Docker — its `modules/treeisonet/{treeLoc,treeOff}.py` already run from
`laspy`+`torch`. An R `system2()` wrapper invokes a Python arm against a frozen
clip in a shared work dir and ingests a detections CSV; a no-op echo arm proves
the R↔Python↔scorer handshake before any model. For apex-detection scoring the
driver stops at `postPeakExtraction` (tree-tops `x,y,z`), skipping the offset
net. The Docker-backed runner and LAZ→PLY conversion that #M6 SegmentAnyTree /

## M8 ForestFormer3D need are **deferred** to those arms (they hit the sm_120

MinkowskiEngine/spconv rebuild wall; see `model-benchmark-plan.md` and the
`gpu-arm-blackwell-sm120` memo).

**Tech Stack:** Python 3.10+ venv (`torch>=2.7 cu128`, `timm<1.0`,
`laspy[lazrs]`, `numpy`, `numpy_indexed`, `numpy_groupies`, `scipy`,
`scikit-image`, `commentjson`, `einops`); R (`lidR`, `terra`, `data.table`) +
the existing `model_bench_lib.R` bridge / `sweep_lib.R` scorer; `curl` +
`sha256sum` for the mirror. One local RTX 5090 (32 GB, Blackwell sm_120),
Docker 29 (only for the deferred arms).

---

### Runtime reality this plan is built on (from the GPU-arm research)

- **TreeisoNet headless = yes, no lifting.** `treeLoc(config, pcd, model_path,
  use_cuda, if_stem=False, custom_resolution=...)` returns per-point
  `[x,y,z,conf,radius]`; `postPeakExtraction(preds[conf>thr])` returns discrete
  tops `[x,y,z,radius]` — the apex list. `treeOff(...)` (optional) returns
  per-point instance ids 1..N. All in `NRCan/TreeAIBox/modules/treeisonet/`.
- **Input:** LAS/LAZ via `laspy`; **XYZ only**; raw coords **min-subtracted**
  (`pcd[:,:3] -= pcd.min(0)`), metric **UTM** (never EPSG:3857). Feed the
  frozen clip's **raw-with-ground** variant (TreeisoNet does its own
  normalization, so the tops come out in **absolute elevation** — see the DTM
  step below).
- **Voxel:** checkpoint is `128³` window at `[0.1,0.1,0.2]` m = 12.8×12.8×25.6 m,
  **hard-edged non-overlapping tiling** (seams). The native baseline is the
  config's own `[0.1,0.1,0.2]` (pass **no** override). `custom_resolution=[r,r,r]`
  (r>0) overrides it to an isotropic grid: 0.8 m → 102.4 m window (whole plot,
  one block); untested off-distribution → **sweep {native, 0.8, 2.0}**.
- **Output:** labeled LAZ extra-dim `pred_itc` int32, **0 = unassigned, trees
  1..N**. Apex reduction = max-Z per id (drop 0). Tops CSV alternative from
  `postPeakExtraction`.
- **Height convention (load-bearing for scoring):** `score_plot`'s height gate
  (`sweep_lib.R` `greedy_match`) compares detection `z` to **field tree height
  (above ground)**. TreeisoNet apex `z` is **absolute elevation** (SOAP ground
  ~1000–2000 m), so it MUST be converted to height-above-ground by subtracting
  the frozen clip's persisted `ground_dtm.tif` before scoring — exactly the
  "raw-ground apex elevations leaking into the height gate" gap that the
  gap→issue map assigns to the `ground_dtm.tif` transform. This is `det_to_agl`
  in #I4.
- **Env:** pure PyTorch + `timm`, **no compiled CUDA ops** → cu128 wheel swap is
  the only Blackwell change. Pass a **local `model_path`** to skip the
  downloader; the matching `{model_name}.json` config must sit beside the
  `.pth`. The config filename uses upstream's own punctuation (parens vs
  underscores vary), so **discover it by glob in the cloned repo**, never
  hardcode the name. Pin `timm<1.0` (newer timm moved `timm.models.layers` →
  `timm.layers`).
- **Accuracy risk (not env):** ALS auto mIoU ~0.59; sparse-NEON zero-shot is the
  real unknown — hence the voxel sweep and an honest caveat in the eventual doc.

### #M7-first design decisions (call these out in review)

1. **venv, not Docker, for #M7** — #I3's first-class backend is a Python
   subprocess; the Docker backend is stubbed but deferred to #M6.
2. **Apex-only (stop at `postPeakExtraction`)** is the default arm; the full
   `treeOff` instance path is an optional variant (kept for parity but not on
   the critical path). Note `treeOff` labels are 1..N (a 0 class exists only for
   out-of-block points), so the apex-only path does not rely on a `treeOff`
   "non-tree 0".
3. **Weights live in a gitignored store** `gpu/store/`, but the
   **`gpu/CHECKSUMS.sha256` manifest is versioned** (outside the store) and
   verified with `sha256sum -c` on every rerun.
4. **Three steps need the GPU** (Task 2 smoke, Task 6) — the card is ~94%
   occupied now; those steps block until it frees. Everything else is GPU-free.

### File structure

- `gpu/mirror_weights.sh` — #I5 download + checksum-verify into `gpu/store/`,
  manifest at `gpu/CHECKSUMS.sha256`.
- `gpu/setup_treeisonet_env.sh` — #I3 venv + pinned `TreeAIBox` clone + config
  copy.
- `gpu/run_echo.py` — #I3 no-op echo arm (handshake proof, no model).
- `gpu/run_treeisonet.py` — #M7 headless driver (treeLoc→postPeakExtraction).
- `scripts/model_runner.R` — #I3 R `system2()` wrapper → detections `data.frame`.
- `scripts/io_bridge.R` — #I4 instance-cloud→`(x,y,z)`, `det_to_agl` DTM
  transform, CRS round-trip.
- `tests/testthat/test-model-runner.R`, `tests/testthat/test-io-bridge.R`.
- `.gitignore` — add `gpu/store/`, `gpu/.venv/`, `gpu/TreeAIBox/` (but NOT
  `gpu/CHECKSUMS.sha256` or `gpu/TreeAIBox.commit`).

---

### Task 1: #I5 — mirror weights/images + versioned checksum manifest

**Files:** Create `gpu/mirror_weights.sh`; Modify `.gitignore`.

No real GPU needed; downloads + hashing. The manifest lives **outside** the
ignored store so it is versioned; reruns verify against it.

- [ ] **Step 1: Add the store + clone paths to `.gitignore`**

```text
gpu/store/
gpu/.venv/
gpu/TreeAIBox/
```

(Do NOT ignore `gpu/CHECKSUMS.sha256` or `gpu/TreeAIBox.commit` — they are
versioned.)

- [ ] **Step 2: Write `gpu/mirror_weights.sh`**

```bash
#!/usr/bin/env bash
# #I5 mirror: fetch model weights/images into a local store, verify SHA256
# against a VERSIONED manifest (gpu/CHECKSUMS.sha256). First run records the
# full hash; later runs hard-fail on any mismatch.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
STORE="${1:-$HERE/store}"; MAN="$HERE/CHECKSUMS.sha256"
mkdir -p "$STORE/treeaibox" "$STORE/segmentanytree" "$STORE/forestformer3d"
touch "$MAN"
fetch() {  # url dest sha_prefix(sanity; "" to skip)
  local url="$1" dest="$2" pre="$3" base got exp
  base="$(basename "$dest")"
  [ -f "$dest" ] || { echo "GET $base"; curl -fL --retry 3 "$url" -o "$dest"; }
  got="$(sha256sum "$dest" | cut -d' ' -f1)"
  if [ -n "$pre" ]; then case "$got" in "$pre"*) : ;;
    *) echo "FAIL prefix $base: $got !~ $pre"; exit 1;; esac; fi
  if grep -q "  $base\$" "$MAN"; then
    exp="$(grep "  $base\$" "$MAN" | cut -d' ' -f1)"
    [ "$got" = "$exp" ] || { echo "FAIL manifest $base: $got != $exp"; exit 1; }
  else echo "$got  $base" >> "$MAN"; fi
  echo "OK $base ($got)"
}
B=https://github.com/NRCan/TreeAIBox/releases/download/v1.0
# --- TreeisoNet ALS (FIRST; #M7) ---
fetch "$B/treeisonet_als_reclamation_treeloc_esegformer3D_128_10cm_GPU4GB.pth" \
      "$STORE/treeaibox/als_treeloc.pth" c6e15a21
fetch "$B/treeisonet_als_reclamation_treeoff_esegformer3D_128_10cm_GPU4GB.pth" \
      "$STORE/treeaibox/als_treeoff.pth" 6cf5d890
# --- SegmentAnyTree (#M6 later): Git-LFS media blob ---
fetch "https://media.githubusercontent.com/media/SmartForest-no/SegmentAnyTree/main/model_file/PointGroup-PAPER.pt" \
      "$STORE/segmentanytree/PointGroup-PAPER.pt" 0b4d74b4
# --- ForestFormer3D (#M8 later): Zenodo zip (md5 on Zenodo; sha recorded here) ---
fetch "https://zenodo.org/api/records/16742708/files/clean_forestformer.zip/content" \
      "$STORE/forestformer3d/clean_forestformer.zip" ""
echo "manifest -> $MAN"; cat "$MAN"
# NOTE: the TreeisoNet config JSONs are copied from the cloned repo in Task 2
# (their upstream filenames vary in punctuation; glob, don't hardcode a URL).
```

- [ ] **Step 3: Run it**

Run: `bash gpu/mirror_weights.sh`
Expected: `OK als_treeloc.pth (...)` etc. and a populated `gpu/CHECKSUMS.sha256`.
The two TreeisoNet `.pth` (~21 MB) verify against prefixes `c6e15a21` /
`6cf5d890`. SAT (~635 MB) and FF3D (~198 MB) are deferred-arm artifacts (slow
fetch is fine). Re-run once to confirm it verifies (no re-download, no mismatch).

- [ ] **Step 4: Commit (manifest is versioned; payloads are not)**

```bash
git add gpu/mirror_weights.sh .gitignore gpu/CHECKSUMS.sha256
git commit -m "feat(#I5): weights/image mirror + versioned checksum manifest"
```

---

### Task 2: #I3a — TreeisoNet venv (pinned) + config copy + sm_120 weight-load smoke

**Files:** Create `gpu/setup_treeisonet_env.sh`; commit `gpu/TreeAIBox.commit`.

Env setup is GPU-free; the final smoke loads the **real** checkpoint + config and
runs `treeLoc` on the 5090 — the critical Blackwell de-risk (proves cu128 torch
runs the esegformer3D forward AND that weight loading works).

- [ ] **Step 1: Write `gpu/setup_treeisonet_env.sh`**

```bash
#!/usr/bin/env bash
# #I3a: headless TreeisoNet env. Pinned repo + pinned cu128 torch (Blackwell).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
COMMITF="$ROOT/TreeAIBox.commit"
if [ ! -d "$ROOT/TreeAIBox" ]; then
  git clone https://github.com/NRCan/TreeAIBox "$ROOT/TreeAIBox"
  if [ -f "$COMMITF" ]; then git -C "$ROOT/TreeAIBox" checkout "$(cat "$COMMITF")"
  else git -C "$ROOT/TreeAIBox" rev-parse HEAD > "$COMMITF"; fi
fi
# Pin tested versions (cu128 = Blackwell sm_120; timm<1.0 keeps models.layers).
python3 -m venv "$ROOT/.venv"; . "$ROOT/.venv/bin/activate"; pip install --upgrade pip
pip install torch==2.7.0 torchvision==0.22.0 --index-url https://download.pytorch.org/whl/cu128
pip install "timm==0.9.16" "laspy[lazrs]" numpy numpy_indexed numpy_groupies \
            scipy scikit-image commentjson einops
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
ls "$ROOT/store/treeaibox/"*.json
```

- [ ] **Step 2: Run env setup; record the pin; confirm sm_120 is in the arch list**

Run: `bash gpu/setup_treeisonet_env.sh`
Expected: torch `2.7.0`, `cuda 12.8`, `avail True`, device `NVIDIA GeForce RTX
5090`, `sm_120` present in `arch`; at least one `*treeloc*.json` copied into the
store. If `sm_120` is absent or `avail False`, stop — wrong wheel/driver (need
R570+). Then `git add gpu/TreeAIBox.commit` (created on first clone).

- [ ] **Step 3: sm_120 weight-load smoke (GPU-gated; needs the card + Task 1)**

```bash
. gpu/.venv/bin/activate
LOC=gpu/store/treeaibox/als_treeloc.pth
CFG=$(ls gpu/store/treeaibox/*reclamation*treeloc*.json | head -1)
PYTHONPATH=gpu/TreeAIBox python - "$LOC" "$CFG" <<'PY'
import sys, numpy as np
from modules.treeisonet.treeLoc import treeLoc
loc, cfg = sys.argv[1], sys.argv[2]
rng = np.random.default_rng(0)                 # synthetic 20 m occupancy box
pcd = rng.uniform([0,0,0], [20,20,15], (5000,3)).astype(float)
out = np.asarray(treeLoc(cfg, pcd, loc, use_cuda=True, if_stem=False))
print("treeLoc OK; preds shape", out.shape)    # expect [N>=1, 5] = x,y,z,conf,radius
PY
```

Expected: `treeLoc OK; preds shape (..., 5)` with NO `no kernel image is
available` error — confirms the real esegformer3D weights run on sm_120. (If it
errors on the config path, check the copied JSON name; if on args, reconcile
against `treeLoc.py`.)

- [ ] **Step 4: Commit**

```bash
git add gpu/setup_treeisonet_env.sh gpu/TreeAIBox.commit
git commit -m "feat(#I3): pinned TreeisoNet venv (cu128) + sm_120 weight-load smoke"
```

---

### Task 3: #I3b — R model-runner wrapper + no-op echo arm

**Files:** Create `gpu/run_echo.py`, `scripts/model_runner.R`,
`tests/testthat/test-model-runner.R`.

Proves the R↔Python↔scorer handshake with **no model**. Honors crash discipline:
non-zero exit OR missing output → `NULL` (the equal-set guard drops the cell;
never a fake 0-row from a stale/partial file).

- [ ] **Step 1: Write the failing test** `tests/testthat/test-model-runner.R`

```r
source(file.path("..", "..", "scripts", "model_runner.R"), local = TRUE)

test_that("run_python_arm returns x,y,z from a CSV the arm writes", {
  d <- tempfile(); dir.create(d)
  py <- file.path(d, "arm.py"); out <- file.path(d, "det.csv")
  writeLines(c("import sys",
               "open(sys.argv[2],'w').write('x y z\\n1 2 3\\n4 5 6\\n')"), py)
  det <- run_python_arm("python3", py, input = "ignored.laz", out_csv = out)
  expect_identical(names(det), c("x", "y", "z"))
  expect_equal(nrow(det), 2L); expect_equal(det$z, c(3, 6))
})

test_that("run_python_arm returns NULL on non-zero exit even if a stale CSV exists", {
  d <- tempfile(); dir.create(d)
  out <- file.path(d, "stale.csv"); writeLines("x y z\\n9 9 9", out)  # pre-existing
  py <- file.path(d, "boom.py"); writeLines("import sys; sys.exit(3)", py)
  expect_null(run_python_arm("python3", py, input = "x.laz", out_csv = out))
})
```

- [ ] **Step 2: Run it — expect FAIL** (`run_python_arm` undefined)

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-model-runner.R")'`

- [ ] **Step 3: Write `scripts/model_runner.R`**

```r
#!/usr/bin/env Rscript
# #I3 model-runner contract. Invoke a GPU model arm as a SUBPROCESS against a
# frozen clip and read its detections CSV back as a contract-valid
# data.frame(x,y,z). Backend = a Python venv (TreeisoNet, #M7); a Docker backend
# is added with #M6. Non-zero exit OR missing output -> NULL so the equal-set
# guard drops the cell (never a fake 0-row from a stale/partial file).
.find <- function(rel) Find(file.exists, c(file.path("scripts", rel),
                                           file.path("..", "..", "scripts", rel),
                                           file.path(getwd(), "scripts", rel)))
source(.find("model_bench_lib.R"))

run_python_arm <- function(venv_python, script, input, out_csv,
                           extra = character(), timeout = 1800) {
  if (file.exists(out_csv)) unlink(out_csv)          # never read a stale file
  args <- c(script, input, out_csv, extra)
  out <- tryCatch(system2(venv_python, args, stdout = TRUE, stderr = TRUE,
                          timeout = timeout), error = function(e) NULL)
  st <- if (is.null(out)) 1L else attr(out, "status")     # non-NULL only if != 0
  if (!is.null(st) && st != 0) return(NULL)          # crash -> skip cell
  if (!file.exists(out_csv)) return(NULL)
  d <- tryCatch(read.table(out_csv, header = TRUE), error = function(e) NULL)
  if (is.null(d) || !all(c("x", "y", "z") %in% names(d)))
    return(data.frame(x = numeric(), y = numeric(), z = numeric()))  # ran-empty
  det <- data.frame(x = as.numeric(d$x), y = as.numeric(d$y), z = as.numeric(d$z))
  assert_detection_contract(det)
  det
}
```

- [ ] **Step 4: Run the test — expect PASS**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-model-runner.R")'`

- [ ] **Step 5: Write the no-op echo arm `gpu/run_echo.py`**

```python
#!/usr/bin/env python3
# #I3 handshake proof: read the input cloud, emit 3 fixed detections near its
# centroid. No model, no GPU -- proves R -> python -> CSV -> scorer end to end.
import sys, numpy as np, laspy
inp, out = sys.argv[1], sys.argv[2]
las = laspy.read(inp)
cx, cy = float(np.median(las.x)), float(np.median(las.y))
det = np.array([[cx, cy, 20.0], [cx + 5, cy, 18.0], [cx, cy + 5, 15.0]])
np.savetxt(out, det, header="x y z", comments="")
```

- [ ] **Step 6: Smoke the echo arm against a real frozen clip (GPU-free)**

```bash
. gpu/.venv/bin/activate
CLIP=work/neon/SOAP/frozen/SOAP/SOAP_001/native/clip_rawground.laz
python gpu/run_echo.py "$CLIP" /tmp/echo_det.csv && cat /tmp/echo_det.csv
Rscript -e 'source("scripts/model_runner.R"); print(run_python_arm("gpu/.venv/bin/python","gpu/run_echo.py","'"$CLIP"'","/tmp/echo_det2.csv"))'
```

Expected: a 3-row CSV and a 3-row R `data.frame` — handshake proven.

- [ ] **Step 7: Commit**

```bash
git add gpu/run_echo.py scripts/model_runner.R tests/testthat/test-model-runner.R
git commit -m "feat(#I3): R model-runner wrapper + no-op echo arm (handshake)"
```

---

### Task 4: #I4 — instance→apex bridge + DTM height transform

**Files:** Create `scripts/io_bridge.R`, `tests/testthat/test-io-bridge.R`.

Two model-output-side bridge pieces: (a) `instances_to_det` — labeled cloud
(instance-id field; 0 = unassigned) → apex `(x,y,z)` via `reduce_instances`;
(b) `det_to_agl` — convert apex **absolute elevation → height above ground** via
the frozen clip's `ground_dtm.tif`, so `z` matches the scorer's height gate.

- [ ] **Step 1: Write the failing tests** `tests/testthat/test-io-bridge.R`

```r
source(file.path("..", "..", "scripts", "io_bridge.R"), local = TRUE)
suppressMessages({ library(lidR); library(terra) })

test_that("instances_to_det reduces a labeled cloud to per-id apexes, dropping 0", {
  dt <- data.frame(X = c(10,10,40,40, 25), Y = c(10,10,12,12, 25),
                   Z = c(18, 9, 12, 5,  3),
                   pred_itc = c(1L,1L,2L,2L, 0L))   # id 0 = unassigned
  det <- instances_to_det(dt, id_field = "pred_itc")
  expect_identical(names(det), c("x", "y", "z"))
  expect_equal(nrow(det), 2L)                 # two trees; the 0-point dropped
  expect_setequal(det$z, c(18, 12))           # max-Z apex per id
})

test_that("instances_to_det is NA-safe when the id column already has NA", {
  dt <- data.frame(X = c(1,2), Y = c(1,2), Z = c(5,6),
                   pred_itc = c(NA_integer_, 0L))    # NA and 0 both unassigned
  expect_equal(nrow(instances_to_det(dt, id_field = "pred_itc")), 0L)
})

test_that("det_to_agl subtracts the DTM to give height above ground", {
  r <- terra::rast(nrows=10, ncols=10, xmin=0, xmax=100, ymin=0, ymax=100)
  terra::values(r) <- 50                              # flat ground at 50 m
  f <- tempfile(fileext=".tif"); terra::writeRaster(r, f)
  det <- data.frame(x = c(25, 75), y = c(25, 75), z = c(77, 62))  # absolute elev
  agl <- det_to_agl(det, f)
  expect_equal(agl$z, c(27, 12))                      # 77-50, 62-50
})
```

- [ ] **Step 2: Run it — expect FAIL** (functions undefined)

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-io-bridge.R")'`

- [ ] **Step 3: Write `scripts/io_bridge.R`**

```r
#!/usr/bin/env Rscript
# #I4 I/O bridge. instances_to_det: labeled point table/LAS (instance-id field;
# 0 OR NA = unassigned) -> apex (x,y,z) via the bridge's reduce_instances
# (max-Z per id). det_to_agl: apex absolute elevation -> height above ground via
# the frozen clip's ground_dtm.tif, so z matches score_plot's height gate.
# (LAZ->PLY for SAT/FF3D is added with #M6.)
.find <- function(rel) Find(file.exists, c(file.path("scripts", rel),
                                           file.path("..", "..", "scripts", rel),
                                           file.path(getwd(), "scripts", rel)))
source(.find("model_bench_lib.R"))
suppressMessages(library(terra))

instances_to_det <- function(pts, id_field = "pred_itc",
                             x = "X", y = "Y", z = "Z") {
  pts <- as.data.frame(pts)
  if (!id_field %in% names(pts)) return(data.frame(x=numeric(), y=numeric(), z=numeric()))
  ids <- pts[[id_field]]
  pts[[id_field]][!is.na(ids) & ids == 0] <- NA      # 0 or NA = unassigned -> dropped
  det <- reduce_instances(pts, id_col = id_field, x = x, y = y, z = z)
  assert_detection_contract(det)
  det
}

read_instances_laz <- function(path, id_field = "pred_itc") {
  las <- lidR::readLAS(path)
  if (is.null(las) || lidR::is.empty(las)) return(data.frame(x=numeric(), y=numeric(), z=numeric()))
  instances_to_det(las@data, id_field = id_field)
}

# Subtract ground elevation at each apex's (x,y); apexes off the DTM are dropped.
det_to_agl <- function(det, dtm_path) {
  if (!nrow(det)) return(det)
  g <- terra::extract(terra::rast(dtm_path), cbind(det$x, det$y))[, 1]
  ok <- !is.na(g)
  det <- det[ok, , drop = FALSE]
  det$z <- det$z - g[ok]
  rownames(det) <- NULL
  det
}
```

- [ ] **Step 4: Run the tests — expect PASS**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-io-bridge.R")'`

- [ ] **Step 5: Commit**

```bash
git add scripts/io_bridge.R tests/testthat/test-io-bridge.R
git commit -m "feat(#I4): instance->apex bridge + DTM height-above-ground transform"
```

---

### Task 5: #I4b — CRS / units round-trip acceptance test

**Files:** Modify `tests/testthat/test-io-bridge.R`.

A known stem in UTM survives write→read→reduce to within tolerance, and `z` is
the apex height. Guards the CRS/units handoff (the repo's 3857 gotcha). Codex
verified lidR round-trips the `pred_itc` extra dim without renaming.

- [ ] **Step 1: Add the failing round-trip test (apex fixture is self-consistent)**

```r
test_that("CRS/units round-trip: a known UTM stem reduces to its apex within tol", {
  skip_if_not_installed("lidR")
  e0 <- 320000; n0 <- 4100000              # UTM 11N; apex (max-Z) at (e0+0.1, n0+0.1)
  dt <- data.frame(X = c(e0+0.1, e0,   e0-0.1),
                   Y = c(n0+0.1, n0,   n0-0.1),
                   Z = c(27,     14,   6),
                   pred_itc = c(1L, 1L, 1L))
  las <- lidR::LAS(data.frame(X=dt$X, Y=dt$Y, Z=dt$Z)); sf::st_crs(las) <- 32611L
  las <- lidR::add_lasattribute(las, dt$pred_itc, "pred_itc", "instance id")
  f <- tempfile(fileext = ".laz"); lidR::writeLAS(las, f)
  det <- read_instances_laz(f, id_field = "pred_itc")
  expect_equal(nrow(det), 1L)
  expect_lt(sqrt((det$x - (e0+0.1))^2 + (det$y - (n0+0.1))^2), 0.5)  # apex xy in UTM
  expect_equal(det$z, 27)                                            # apex elevation
})
```

- [ ] **Step 2: Run it — PASS** (if `pred_itc` does not survive `writeLAS`, read
the extra dim by its on-disk name in `read_instances_laz` and re-run)

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-io-bridge.R")'`

- [ ] **Step 3: Commit**

```bash
git add tests/testthat/test-io-bridge.R
git commit -m "test(#I4): CRS/units round-trip acceptance for instance->apex"
```

---

### Task 6: #M7 walking skeleton — TreeisoNet on one SOAP clip (GPU-gated)

**Files:** Create `gpu/run_treeisonet.py`.

The acceptance gate: real TreeisoNet inference on one frozen clip → tops → AGL →
scored. Needs the GPU free and Tasks 1–4 done.

- [ ] **Step 1: Write `gpu/run_treeisonet.py`** (apex-only; CRS preflight; empty-guard)

```python
#!/usr/bin/env python3
# #M7 headless driver (apex-only). treeLoc -> postPeakExtraction -> tops CSV in
# UTM with ABSOLUTE Z (R converts to height-above-ground via ground_dtm.tif).
# Usage: run_treeisonet.py <input.laz> <out.csv> <loc.pth> <loc.json> [voxel] [conf]
#   voxel <= 0 -> use the checkpoint's native [0.1,0.1,0.2]; >0 -> isotropic override.
import os, sys, numpy as np, laspy
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "TreeAIBox"))
from modules.treeisonet.treeLoc import treeLoc, postPeakExtraction
inp, out, loc_pth, loc_cfg = sys.argv[1:5]
voxel = float(sys.argv[5]) if len(sys.argv) > 5 else 0.0
conf  = float(sys.argv[6]) if len(sys.argv) > 6 else 0.5
las = laspy.read(inp)
try: epsg = las.header.parse_crs().to_epsg()
except Exception: epsg = None
if epsg == 3857:
    sys.exit("ERROR: input is EPSG:3857 (web-mercator); reproject to metric UTM")
pcd = np.transpose([np.asarray(las.x), np.asarray(las.y), np.asarray(las.z)]).astype(float)
pmin = pcd.min(0); pcd[:, :3] -= pmin
cr = np.array([voxel, voxel, voxel]) if voxel > 0 else np.zeros(3)
preds = np.atleast_2d(treeLoc(loc_cfg, pcd, loc_pth, use_cuda=True,
                              if_stem=False, custom_resolution=cr))
sel = preds[preds[:, 3] > conf] if preds.size else preds
if sel.shape[0] < 2:                      # too few peaks -> 0 tops (legit recall 0)
    np.savetxt(out, np.empty((0, 3)), header="x y z", comments=""); sys.exit(0)
tops = np.atleast_2d(postPeakExtraction(sel))[:, :3] + pmin[:3]   # back to UTM
np.savetxt(out, tops, header="x y z", comments="")
print(f"wrote {len(tops)} tops -> {out}")
```

- [ ] **Step 2: Run on one SOAP native clip, voxel sweep {native, 0.8, 2.0} (GPU-gated)**

```bash
. gpu/.venv/bin/activate
CLIP=work/neon/SOAP/frozen/SOAP/SOAP_001/native/clip_rawground.laz
LOC=gpu/store/treeaibox/als_treeloc.pth
CFG=$(ls gpu/store/treeaibox/*reclamation*treeloc*.json | head -1)
for V in 0 0.8 2.0; do   # 0 = native [0.1,0.1,0.2]
  PYTHONPATH=gpu/TreeAIBox python gpu/run_treeisonet.py "$CLIP" "/tmp/tops_$V.csv" "$LOC" "$CFG" "$V" 0.5
done
wc -l /tmp/tops_*.csv
```

Expected: each writes a tops CSV; counts are plausible (a SOAP plot has ~10–40
in-core stems; expect tens to low-hundreds of raw tops). Compare counts across
voxel sizes — the empirical sweep the research flagged.

- [ ] **Step 3: Score one clip through the existing harness (GPU-free, R)**

```bash
Rscript -e '
source("scripts/sweep_lib.R"); source("scripts/model_runner.R"); source("scripts/io_bridge.R")
gt <- read.csv("work/neon/SOAP/ground_truth_stems.csv"); pc <- read.csv("work/neon/SOAP/plot_centroids.csv")
ci <- pc[pc$plotID=="SOAP_001",][1,]; ph <- plot_half(ci$plotType)
stems <- gt[gt$live & gt$is_tree & gt$plotID=="SOAP_001" &
            abs(gt$E-ci$easting)<=ph & abs(gt$N-ci$northing)<=ph,]
cfg <- Sys.glob("gpu/store/treeaibox/*reclamation*treeloc*.json")[1]
det <- run_python_arm("gpu/.venv/bin/python","gpu/run_treeisonet.py",
        "work/neon/SOAP/frozen/SOAP/SOAP_001/native/clip_rawground.laz","/tmp/m7_det.csv",
        extra=c("gpu/store/treeaibox/als_treeloc.pth", cfg, "0.8", "0.5"))
stopifnot(!is.null(det))                              # acceptance must produce a table
det <- det_to_agl(det, "work/neon/SOAP/frozen/SOAP/SOAP_001/native/ground_dtm.tif")
print(score_plot(stems, det, tol_xy=4, core_cx=ci$easting, core_cy=ci$northing, core_half=ph))'
```

Expected: a one-row score (recall/precision/F1 + per-class), recall ∈ [0,1] — the
GPU arm is now a peer of every other detector in the harness.

- [ ] **Step 4: Commit**

```bash
git add gpu/run_treeisonet.py
git commit -m "feat(#M7): TreeisoNet headless apex driver + one-clip skeleton"
```

---

### Deferred (explicitly out of this plan)

- **Docker runner backend** + **LAZ→PLY converter** — needed by #M6
  SegmentAnyTree / #M8 ForestFormer3D, which also need the sm_120 image rebuild.
- **Full #M7 ladder run** (18 plots × 5 rungs, chosen voxel size, optional
  `treeOff` instance variant) and folding TreeisoNet into
  `analyze_model_benchmark.R` + the result doc — the next plan, once this
  skeleton scores clean.
- **#M6 / #M8 arms** themselves.

### Self-review notes

- **GPU-gated steps:** Task 2 Step 3, Task 6 Steps 2–3. The card is busy now;
  these wait. Everything else (mirror, env install, R wrappers, bridge, tests)
  is GPU-free and can land first.
- **The DTM height transform (`det_to_agl`) is load-bearing** — without it apex
  `z` is absolute elevation and the height gate rejects every match.
- **Config name is discovered by glob**, never hardcoded (upstream punctuation
  varies). The pinned commit (`gpu/TreeAIBox.commit`) keeps it reproducible.
- **Contract discipline kept:** non-zero exit / missing output → NULL; ran-empty
  → 0-row; apex = max-Z per id via the existing `reduce_instances`.
- **Two unverified-until-run specifics:** the `treeLoc` return columns at Task 2
  Step 3, and `pred_itc` surviving `writeLAS` at Task 5 — both have an in-step
  "adjust if it errors" instruction; neither is a placeholder.
- **#M7-first is load-bearing:** no Docker on the critical path; if review
  prefers SAT-first, Tasks 2/3/6 change substantially (Docker backend + image
  rebuild move onto the critical path).

## Implementation findings (empirical pivots during the skeleton)

The walking skeleton (SOAP_001) surfaced four input/reduction choices that
deviate from the as-planned Task 6 and are baked into `gpu/run_treeisonet.py`.

- **Feed the NORMALIZED clip, not raw-with-ground.** Empirically the normalized
  variant gives ~2x higher TreeLoc confidence (max 0.58 vs 0.33) and a clean
  Z range (53 m vs a 414 m raw clip corrupted by height outliers). Because
  normalized Z is already height-above-ground, `det_to_agl` is NOT needed for
  M7 (Codex critical 1 dissolves via the input choice). `det_to_agl` stays in
  `io_bridge.R` for the raw-output Docker arms (SAT).
- **Native voxel only.** Coarsening to 0.8 m collapsed confidence to ~0.10
  (off-distribution); the `custom_resolution` sweep is not useful here.
- **Snap apex z to the local canopy max.** TreeLoc returns a tree *location*
  (z ~ stem base, 0–2 m), not an apex. The driver snaps each top's z to the max
  normalized height within 2 m — the same canopy-surface height the CHM arms
  use, so z is comparable across arms (fair, not inflation).
- **Confidence threshold is sensitive and must be calibrated once** (zero-shot,
  not per-plot). conf 0.15 → recall 0.30 on SOAP_001; 0.2/0.3 → 0 (postPeak
  clustering is threshold-sensitive). The full-ladder run picks one fixed conf.

Skeleton result: TreeisoNet zero-shot on SOAP_001 ~ recall 0.30 / precision 0.08
at conf 0.15 — weak, consistent with the ALS auto-mIoU ~0.59 caveat. The full
ladder + fixed-conf run + folding into `analyze_model_benchmark.R` is the next
plan.
