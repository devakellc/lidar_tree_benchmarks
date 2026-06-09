# ForestFormer3D Benchmark Arm (#M8) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> superpowers:subagent-driven-development (recommended) or
> superpowers:executing-plans to implement this plan task-by-task. Steps use
> checkbox (`- [ ]`) syntax for tracking.

**Goal:** Wire ForestFormer3D (#M8) into the NEON density-ladder benchmark as a
zero-shot apex detector on SOAP (native + 8 pts/m²), scored by the existing
harness, and fold the result into `results/model-benchmark-results.md`.

**Architecture:** R driver tiles each plot into 16 m-radius cylinders → a Docker
container (`ff3d-sm120`) runs FF3D once per plot over all its cylinders → one
merged UTM LAZ → pure-R `ff3d_collapse` (cross-block apex-cluster dedup → reduce)
→ `agl_guard` (z→AGL with a wholesale-off-DTM backstop) → `score_plot`. Mirrors
`detect_segmentanytree_sweep.R`. The analyzer gains an additive native+8 FF3D
comparison that cannot shrink the existing full ladder.

**Tech Stack:** R (lidR, terra, data.table, testthat), Python (torch 2.7/cu128,
mmdet3d/oneformer3d, spconv-cu128, laspy) in the `ff3d-sm120` Docker image, the
design at `docs/superpowers/specs/2026-06-09-forestformer3d-benchmark-arm-design.md`.

**Working dir:** the worktree `.claude/worktrees/gpu-ff3d-arm` on branch
`worktree-gpu-ff3d-arm`. The FF3D repo + weights are at
`gpu/store/forestformer3d/ForestFormer3D/` (gitignored runtime data, already
present). Run tests with `Rscript tests/run_tests.R` from the worktree root.

---

## Task 1: `dedup_blocks()` cross-block dedup helper (pure R, TDD)

**Files:**

- Modify: `tests/testthat/helper-synth.R` (add `synth_block_points()`)
- Create: `tests/testthat/test-dedup-blocks.R`
- Modify: `scripts/model_bench_lib.R` (add `dedup_blocks()` after
  `reduce_instances`, ~line 24)
- [ ] **Step 1: Add the synthetic fixture to `helper-synth.R`**

Append to `tests/testthat/helper-synth.R`:

```r
# Block-labelled points for cross-block dedup. Columns block, inst, X, Y, Z.
# A (b0,10,10,18) & B (b0,11,10,15): SAME block, 1 m apart -> must stay distinct.
# C (b0,40,40,12) & C' (b1,40.1,40,11.8): DIFFERENT blocks, ~0.1 m -> must merge
# (merged apex = max-Z = 12). D (b1,60,60,10): isolated. inst 0 = unassigned.
synth_block_points <- function() {
  rows <- function(block, inst, x, y, ztop, zmid)
    data.table(block = as.integer(block), inst = as.integer(inst),
               X = c(x, x + 0.1, x), Y = c(y, y, y + 0.1),
               Z = c(ztop, zmid, zmid - 3))
  dt <- rbind(
    rows(0, 1, 10,   10,   18, 12),
    rows(0, 2, 11,   10,   15, 10),
    rows(0, 3, 40,   40,   12,  8),
    rows(1, 1, 40.1, 40.0, 11.8, 7),
    rows(1, 2, 60,   60,   10,  6))
  noise <- data.table(block = 1L, inst = 0L, X = 25, Y = 25, Z = 1)
  rbind(dt, noise)
}
```

- [ ] **Step 2: Write the failing test**

Create `tests/testthat/test-dedup-blocks.R`:

```r
source(file.path("..", "..", "scripts", "model_bench_lib.R"), local = TRUE)

test_that("dedup_blocks merges cross-block dups, keeps same-block neighbors distinct", {
  rel <- data.table::as.data.table(dedup_blocks(synth_block_points(), merge_tol = 2.0))
  gid <- function(bx, by) unique(rel[abs(X - bx) < 0.5 & abs(Y - by) < 0.5]$global_id)
  a <- gid(10, 10); b <- gid(11, 10)
  expect_length(a, 1L); expect_length(b, 1L)
  expect_false(a == b)                       # same block within tol -> distinct
  expect_equal(gid(40, 40), gid(40.1, 40))   # different blocks within tol -> merged
  expect_equal(length(unique(rel$global_id)), 4L)
})

test_that("dedup_blocks + reduce_instances yields the max-Z apex across merged blocks", {
  rel <- dedup_blocks(synth_block_points(), merge_tol = 2.0)
  det <- reduce_instances(rel, id_col = "global_id")
  expect_equal(nrow(det), 4L)
  c_row <- det[abs(det$x - 40) < 1 & abs(det$y - 40) < 1, ]
  expect_equal(nrow(c_row), 1L); expect_equal(c_row$z, 12)   # 12 > 11.8
})

test_that("dedup_blocks drops unassigned (0/NA) and is 0-row safe", {
  empty <- data.frame(block = integer(), inst = integer(),
                      X = numeric(), Y = numeric(), Z = numeric())
  rel <- dedup_blocks(empty)
  expect_equal(nrow(rel), 0L)
  det <- reduce_instances(rel, id_col = "global_id")
  expect_identical(names(det), c("x", "y", "z")); expect_equal(nrow(det), 0L)
})
```

- [ ] **Step 3: Run it; verify it fails**

Run: `Rscript tests/run_tests.R 2>&1 | grep -A2 dedup-blocks`
Expected: FAIL — `could not find function "dedup_blocks"`.

- [ ] **Step 4: Implement `dedup_blocks()` in `model_bench_lib.R`**

Insert after `reduce_instances` (after line 24, before the crown-diameter block):

```r
## ---- cross-block apex-cluster dedup (#M8) --------------------------------
# Stacked per-cylinder labelled points (block, inst, X, Y, Z; inst 0/NA =
# unassigned) -> the SAME table relabelled with a globally-consistent integer
# `global_id`. Per (block, inst) apex (max-Z); union-find over apexes restricted
# to pairs in DIFFERENT blocks within horizontal `merge_tol` (so the model's own
# within-cylinder over-segmentation is NEVER laundered). The driver then calls
# reduce_instances(id_col = "global_id"). Returns a 0-row frame (with global_id)
# when nothing is assigned.
dedup_blocks <- function(pts, merge_tol = 2.0, block = "block", id = "inst",
                         x = "X", y = "Y", z = "Z") {
  empty <- data.frame(block = integer(), inst = integer(), X = numeric(),
                      Y = numeric(), Z = numeric(), global_id = integer())
  dt <- as.data.table(pts)
  if (!nrow(dt) || !all(c(block, id, x, y, z) %in% names(dt))) return(empty)
  ids <- dt[[id]]
  dt <- dt[!is.na(ids) & ids != 0, c(block, id, x, y, z), with = FALSE]
  if (!nrow(dt)) return(empty)
  setnames(dt, c(block, id, x, y, z), c("block", "inst", "X", "Y", "Z"))
  dt[, key := .GRP, by = .(block, inst)]
  ap <- dt[, .(blk = block[1L], cx = X[which.max(Z)], cy = Y[which.max(Z)]),
           by = key][order(key)]
  n <- nrow(ap)
  parent <- seq_len(n)
  find <- function(i) { r <- i; while (parent[r] != r) r <- parent[r]
                        while (parent[i] != r) { nx <- parent[i]; parent[i] <<- r; i <- nx }
                        r }
  if (n > 1L) for (i in seq_len(n - 1L)) for (j in (i + 1L):n)
    if (ap$blk[i] != ap$blk[j] &&
        (ap$cx[i] - ap$cx[j])^2 + (ap$cy[i] - ap$cy[j])^2 <= merge_tol^2) {
      ri <- find(i); rj <- find(j); if (ri != rj) parent[rj] <- ri
    }
  roots <- vapply(seq_len(n), find, integer(1))
  ap[, global_id := match(roots, sort(unique(roots)))]
  dt <- merge(dt, ap[, .(key, global_id)], by = "key")
  dt[, key := NULL]
  as.data.frame(dt)
}
```

- [ ] **Step 5: Run tests; verify pass**

Run: `Rscript tests/run_tests.R 2>&1 | tail -20`
Expected: `dedup-blocks: ...` all pass; the whole suite is green.

- [ ] **Step 6: Commit**

```bash
git add scripts/model_bench_lib.R tests/testthat/helper-synth.R tests/testthat/test-dedup-blocks.R
git commit -m "feat(#18): dedup_blocks cross-block apex-cluster dedup

Different-block-only union-find on per-(block,instance) apexes within merge_tol;
same-block instances stay distinct so model over-segmentation is scored honestly.
Returns globally-consistent ids for the canonical reduce_instances.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 2: `ff3d_collapse()` + `agl_guard()` reader/guard helpers (pure R, TDD)

These are the FF3D-specific halves of the reader closure + score guard, factored
out so they are unit-testable without a GPU (the driver in Task 5 just composes
them, exactly as SAT composes `read_instances_laz` + `det_to_agl`).

**Files:**

- Create: `tests/testthat/test-forestformer3d-sweep.R`
- Modify: `scripts/io_bridge.R` (add both after `det_to_agl`, ~line 46)
- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-forestformer3d-sweep.R`:

```r
suppressMessages({ library(lidR); library(terra) })
source(file.path("..", "..", "scripts", "model_bench_lib.R"), local = TRUE)
source(file.path("..", "..", "scripts", "io_bridge.R"), local = TRUE)

# Write a merged labelled LAZ (UTM): UserData = cylinder/block, PointSourceID =
# per-cylinder instance id. Same geometry as synth_block_points().
write_ff3d_laz <- function(path) {
  p <- synth_block_points()
  las <- LAS(data.frame(X = p$X, Y = p$Y, Z = p$Z,
                        UserData = as.integer(p$block),
                        PointSourceID = as.integer(p$inst)))
  lidR::writeLAS(las, path)
}

test_that("ff3d_collapse reads UserData/PointSourceID, dedups across blocks, reduces", {
  f <- tempfile(fileext = ".laz"); write_ff3d_laz(f)
  det <- ff3d_collapse(f, merge_tol = 2.0)
  expect_identical(names(det), c("x", "y", "z"))
  expect_equal(nrow(det), 4L)                          # A, B, merged-C, D
  c_row <- det[abs(det$x - 40) < 1 & abs(det$y - 40) < 1, ]
  expect_equal(nrow(c_row), 1L); expect_equal(c_row$z, 12)
})

test_that("ff3d_collapse returns NULL on an unreadable file (schema failure -> skip)", {
  expect_null(ff3d_collapse(tempfile(fileext = ".laz")))
})

test_that("agl_guard: empty in -> 0-row; partial off-DTM -> AGL; all off-DTM -> NULL", {
  dtm <- tempfile(fileext = ".tif")
  r <- terra::rast(xmin = 0, xmax = 100, ymin = 0, ymax = 100,
                   resolution = 1, vals = 5)            # flat ground at z=5
  terra::writeRaster(r, dtm, overwrite = TRUE)
  # empty in -> empty out (legit ran-but-empty)
  e <- data.frame(x = numeric(), y = numeric(), z = numeric())
  expect_equal(nrow(agl_guard(e, dtm)), 0L)
  # on-DTM apex -> AGL (z - 5)
  on <- data.frame(x = 50, y = 50, z = 25)
  g <- agl_guard(on, dtm); expect_equal(g$z, 20)
  # all apexes off the raster -> wholesale drop -> NULL (skip the cell)
  off <- data.frame(x = c(1e6, 1e6), y = c(1e6, 1e6), z = c(25, 30))
  expect_null(agl_guard(off, dtm))
})
```

- [ ] **Step 2: Run it; verify it fails**

Run: `Rscript tests/run_tests.R 2>&1 | grep -A2 forestformer3d-sweep`
Expected: FAIL — `could not find function "ff3d_collapse"`.

- [ ] **Step 3: Implement both helpers in `io_bridge.R`**

Insert after `det_to_agl` (after line 46):

```r
## ---- ForestFormer3D arm helpers (#M8) ------------------------------------
# Reader for the merged per-plot LAZ ff3d_arm.py writes: UserData = cylinder
# (block), PointSourceID = per-cylinder instance id, coords UTM. Cross-block
# dedup -> canonical reduce -> apex det(x,y,z). NULL on an unreadable file
# (schema failure -> run_docker_arm skips the cell); a valid empty cloud -> a
# legit 0-row frame.
ff3d_collapse <- function(out_laz, merge_tol = 2.0) {
  empty <- data.frame(x = numeric(), y = numeric(), z = numeric())
  las <- tryCatch(lidR::readLAS(out_laz), error = function(e) NULL)
  if (is.null(las)) return(NULL)
  if (lidR::is.empty(las)) return(empty)
  d <- las@data
  if (!all(c("UserData", "PointSourceID") %in% names(d))) return(NULL)
  pts <- data.frame(block = as.integer(d$UserData),
                    inst  = as.integer(d$PointSourceID),
                    X = d$X, Y = d$Y, Z = d$Z)
  reduce_instances(dedup_blocks(pts, merge_tol = merge_tol), id_col = "global_id")
}

# z -> AGL with a wholesale-off-DTM backstop. The FF3D centering-offset restore
# is the fragile step; a frame bug puts every apex off the DTM. Empty in -> empty
# out (legit ran-but-empty, recall 0). Non-empty in but ALL apexes dropped
# off-DTM -> NULL so the driver SKIPS the cell (a frame bug must not masquerade as
# a valid 0-recall row). Partial drops (edge apexes) are legitimate.
agl_guard <- function(det_abs, dtm_path) {
  if (!nrow(det_abs)) return(det_abs)
  det <- det_to_agl(det_abs, dtm_path)
  if (!nrow(det) && isTRUE(attr(det, "n_dropped") > 0)) return(NULL)
  det
}
```

- [ ] **Step 4: Run tests; verify pass**

Run: `Rscript tests/run_tests.R 2>&1 | tail -20`
Expected: `forestformer3d-sweep:` all pass; suite green.

- [ ] **Step 5: Commit**

```bash
git add scripts/io_bridge.R tests/testthat/test-forestformer3d-sweep.R
git commit -m "feat(#18): ff3d_collapse + agl_guard io_bridge helpers

ff3d_collapse: read the merged labelled LAZ (UserData=block, PointSourceID=inst)
-> dedup_blocks -> reduce_instances -> apex det. agl_guard: det_to_agl with a
wholesale-off-DTM backstop so a coordinate-frame bug skips the cell instead of
scoring a fake 0-recall row.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 3: Vendor the mm-stack swaps + bake them into the image (build)

The three `replace_mmdetection_files` swaps modify *site-packages* (in the
image), so they must be baked at build. Vendor copies under the build context so
the build is reproducible from the branch (not the gitignored repo).

**Files:**

- Create: `gpu/forestformer3d-sm120/replace_mmdetection_files/{loops.py,base_model.py,transforms_3d.py}`
- Modify: `gpu/forestformer3d-sm120/Dockerfile` (add COPY + cp before `WORKDIR`)
- [ ] **Step 1: Vendor the three files**

```bash
cd /home/alex/projects/lidar_tree_benchmarks/.claude/worktrees/gpu-ff3d-arm
mkdir -p gpu/forestformer3d-sm120/replace_mmdetection_files
REPO=/home/alex/projects/lidar_tree_benchmarks/gpu/store/forestformer3d/ForestFormer3D
cp "$REPO"/replace_mmdetection_files/{loops.py,base_model.py,transforms_3d.py} \
   gpu/forestformer3d-sm120/replace_mmdetection_files/
ls -1 gpu/forestformer3d-sm120/replace_mmdetection_files/
```

Expected: the three files listed.

- [ ] **Step 2: Add the bake step to the Dockerfile**

In `gpu/forestformer3d-sm120/Dockerfile`, immediately before the `WORKDIR
/workspace` line (currently line ~62), insert:

```dockerfile
# Bake the mmengine/mmdet3d site-packages swaps (FF3D's replace_mmdetection_files)
# into the image so they are present for every --rm run. SITE is conda's
# site-packages. ff3d_repo.patch (which edits the mounted REPO, not site-packages)
# is applied at container start by ff3d_entry.sh instead.
COPY replace_mmdetection_files/ /opt/replace_mmdetection_files/
RUN SITE="$(python -c 'import site; print(site.getsitepackages()[0])')" && \
    cp /opt/replace_mmdetection_files/loops.py        "$SITE/mmengine/runner/loops.py" && \
    cp /opt/replace_mmdetection_files/base_model.py   "$SITE/mmengine/model/base_model/base_model.py" && \
    cp /opt/replace_mmdetection_files/transforms_3d.py "$SITE/mmdet3d/datasets/transforms/transforms_3d.py"
```

- [ ] **Step 3: Rebuild and smoke-test**

```bash
cd /home/alex/projects/lidar_tree_benchmarks/.claude/worktrees/gpu-ff3d-arm
docker build -t ff3d-sm120 -f gpu/forestformer3d-sm120/Dockerfile gpu/forestformer3d-sm120/ \
  > work/build_logs/ff3d_build_swaps.log 2>&1 && echo BUILD_OK
docker run --rm --gpus all ff3d-sm120 python -c \
  "import mmengine.runner.loops, mmengine.model.base_model.base_model, mmdet3d.datasets.transforms.transforms_3d; print('SWAPS_IMPORT_OK')" \
  2>&1 | grep -E "SWAPS_IMPORT_OK|Error"
```

Expected: `BUILD_OK` then `SWAPS_IMPORT_OK`.

- [ ] **Step 4: Commit**

```bash
git add gpu/forestformer3d-sm120/Dockerfile gpu/forestformer3d-sm120/replace_mmdetection_files/
git commit -m "build(#18): bake FF3D mm-stack site-packages swaps into the image

Vendor the three replace_mmdetection_files (loops/base_model/transforms_3d) and
COPY+cp them into mmengine/mmdet3d site-packages at build, so every --rm run has
them. ff3d_repo.patch stays a runtime step (it edits the mounted repo).

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 4: Container CLI `ff3d_arm.py` + `ff3d_entry.sh` (build)

**Files:**

- Create: `gpu/forestformer3d-sm120/ff3d_arm.py`
- Create: `gpu/forestformer3d-sm120/ff3d_entry.sh`
- [ ] **Step 1: Write `ff3d_arm.py`**

Create `gpu/forestformer3d-sm120/ff3d_arm.py` (refactor of `run_infer_save.py`;
runs from the repo root as CWD):

```python
#!/usr/bin/env python
"""FF3D benchmark arm: stage every cylinder clip in <in_dir> as a ForAINetV2
scene, run FF3D once over all of them, and write ONE merged UTM LAZ <out_laz>
(point_source_id = per-cylinder instance id, user_data = cylinder/block index,
instance score in an extra dim). Run from the FF3D repo root.

Usage: python ff3d_arm.py <in_dir> <out_laz> <ckpt>
"""
import os, sys, glob, shutil
import numpy as np
import laspy

_orig = __import__("torch").load
def _load(*a, **k):
    k.setdefault("weights_only", False); return _orig(*a, **k)
__import__("torch").load = _load

IN_DIR, OUT_LAZ, CKPT = sys.argv[1], sys.argv[2], sys.argv[3]
ROOT = "data/ForAINetV2"                       # relative to repo CWD (ephemeral)
INST = os.path.join(ROOT, "forainetv2_instance_data")
META = os.path.join(ROOT, "meta_data")
for d in (INST, META):
    shutil.rmtree(d, ignore_errors=True); os.makedirs(d, exist_ok=True)

# 1. stage each cylinder as a scene; remember its centering offset for UTM restore
clips = sorted(glob.glob(os.path.join(IN_DIR, "cyl_*.laz")))
offsets, scans = {}, []
for clip in clips:
    scan = os.path.splitext(os.path.basename(clip))[0]        # cyl_000
    las = laspy.read(clip)
    xyz = np.vstack([las.x, las.y, las.z]).T.astype(np.float64)
    if xyz.shape[0] < 50:                                     # too sparse -> skip
        continue
    off = np.array([xyz[:, 0].mean(), xyz[:, 1].mean(), xyz[:, 2].min()])
    offsets[scan] = off
    pts = (xyz - off).astype(np.float32)
    n = pts.shape[0]
    np.save(f"{INST}/{scan}_vert.npy", pts)
    np.save(f"{INST}/{scan}_sem_label.npy", np.zeros(n, np.int64))
    np.save(f"{INST}/{scan}_ins_label.npy", np.zeros(n, np.int64))
    np.save(f"{INST}/{scan}_unaligned_bbox.npy", np.zeros((0, 7), np.float32))
    np.save(f"{INST}/{scan}_aligned_bbox.npy", np.zeros((0, 7), np.float32))
    np.save(f"{INST}/{scan}_axis_align_matrix.npy", np.eye(4))
    scans.append(scan)
with open(f"{META}/test_list.txt", "w") as f:
    f.write("\n".join(scans) + "\n")
open(f"{META}/train_list.txt", "w").close(); open(f"{META}/val_list.txt", "w").close()
if not scans:
    print("no usable cylinders"); laspy.LasData(laspy.LasHeader(point_format=3,
        version="1.2")).write(OUT_LAZ); print("ARM_DONE"); sys.exit(0)

# 2. build points/*.bin + the test pkl
os.system(f"python tools/create_data_forainetv2.py forainetv2 --root-path {ROOT}")

# 3. one runner; iterate all scenes; collect UTM-restored labelled points
from mmengine.config import Config, ConfigDict
from mmengine.runner import Runner
import oneformer3d  # noqa: F401
cfg = Config.fromfile("configs/oneformer3d_qs_radius16_qp300_2many.py")
cfg.work_dir = "./work_dirs/arm"; cfg.load_from = CKPT
if cfg.model.get("test_cfg") is None:
    cfg.model.test_cfg = ConfigDict()
cfg.model.test_cfg["output_dir"] = cfg.work_dir
runner = Runner.from_cfg(cfg); runner.load_or_resume(); model = runner.model.eval()
import torch
def npy(x): return x.detach().cpu().numpy() if torch.is_tensor(x) else np.asarray(x)

XS, YS, ZS, INST_ID, BLK, SCORE = [], [], [], [], [], []
for block_i, data in enumerate(runner.test_dataloader):
    scan = data["data_samples"][0].lidar_path if hasattr(
        data["data_samples"][0], "lidar_path") else scans[block_i]
    scan = scans[block_i]                                     # dataloader is ordered
    with torch.no_grad():
        res = model.test_step(data, epoch=0)
    seg = res[0].pred_pts_seg
    per_point = npy(seg.pts_instance_mask[1]).astype(np.int64)  # (N,)
    pts = npy(data["inputs"]["points"][0] if isinstance(
        data["inputs"]["points"], list) else data["inputs"]["points"])[:, :3]
    n = min(len(pts), len(per_point))
    off = offsets[scan]
    XS.append(pts[:n, 0] + off[0]); YS.append(pts[:n, 1] + off[1])
    ZS.append(pts[:n, 2] + off[2]); INST_ID.append(per_point[:n])
    BLK.append(np.full(n, block_i, np.int64))
    # per-instance score broadcast to points (diagnostic only)
    scores = npy(seg.instance_scores).astype(np.float32)
    sc = np.zeros(n, np.float32)
    for k, s in enumerate(scores):
        sc[per_point[:n] == (k + 1)] = s                      # ids are 1-based
    SCORE.append(sc)
    print(f"{scan}: {n} pts, {len(np.unique(per_point[per_point>0]))} instances")

X = np.concatenate(XS); Y = np.concatenate(YS); Z = np.concatenate(ZS)
inst = np.concatenate(INST_ID); blk = np.concatenate(BLK); score = np.concatenate(SCORE)

# 4. write the merged UTM LAZ
h = laspy.LasHeader(point_format=3, version="1.2")
h.offsets = np.array([X.min(), Y.min(), Z.min()]); h.scales = [0.001, 0.001, 0.001]
las = laspy.LasData(h)
las.x, las.y, las.z = X, Y, Z
las.point_source_id = np.clip(inst, 0, 65535).astype(np.uint16)
las.user_data = np.clip(blk, 0, 255).astype(np.uint8)
las.add_extra_dim(laspy.ExtraBytesParams(name="ff3d_score", type=np.float32))
las.ff3d_score = score
las.write(OUT_LAZ)
print(f"wrote {OUT_LAZ}: {len(X)} pts from {len(scans)} cylinders"); print("ARM_DONE")
```

- [ ] **Step 2: Write `ff3d_entry.sh`**

Create `gpu/forestformer3d-sm120/ff3d_entry.sh`:

```bash
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
```

- [ ] **Step 3: Make the entry script executable + sanity-check Python syntax**

```bash
cd /home/alex/projects/lidar_tree_benchmarks/.claude/worktrees/gpu-ff3d-arm
chmod +x gpu/forestformer3d-sm120/ff3d_entry.sh
docker run --rm -v "$PWD/gpu/forestformer3d-sm120:$PWD/gpu/forestformer3d-sm120" \
  ff3d-sm120 python -m py_compile "$PWD/gpu/forestformer3d-sm120/ff3d_arm.py" \
  && echo PY_COMPILE_OK
```

Expected: `PY_COMPILE_OK` (no syntax errors). Full functional validation is the
live smoke in Task 7.

- [ ] **Step 4: Commit**

```bash
git add gpu/forestformer3d-sm120/ff3d_arm.py gpu/forestformer3d-sm120/ff3d_entry.sh
git commit -m "feat(#18): ff3d_arm.py multi-cylinder driver + ff3d_entry.sh

ff3d_arm.py stages every cylinder clip as a scene, runs FF3D once, and writes one
merged UTM LAZ (point_source_id=instance, user_data=block, ff3d_score extra dim),
restoring each cylinder's centering offset. ff3d_entry.sh cds to the mounted repo,
applies ff3d_repo.patch idempotently, and execs the driver.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 5: R driver `detect_forestformer3d_sweep.R` (composition)

**Files:**

- Create: `scripts/detect_forestformer3d_sweep.R`

No new unit test (mirrors SAT/treeisonet, which have none — their logic is in the
Task 1/2 helpers and `run_docker_arm`, all already tested). Validated by the live
run (Task 7).

- [ ] **Step 1: Write the driver**

Create `scripts/detect_forestformer3d_sweep.R`:

```r
#!/usr/bin/env Rscript
# ForestFormer3D (#M8) density-ladder arm. Tiles each plot core into 16 m-radius
# cylinders, runs FF3D zero-shot ONCE per plot over all cylinders in the
# ff3d-sm120 container (RAW-WITH-GROUND frozen clip), cross-block dedups the
# instances, reduces to apex detections, and scores against field stems. SOAP,
# native + 8 only (the heaviest arm). Serial -- one GPU.
#
#   rawground.laz --clip_circle(16) x N--> in_dir/cyl_*.laz
#     --run_docker_arm(ff3d-sm120, ff3d_entry.sh)--> merged.laz
#     --ff3d_collapse(dedup_blocks + reduce)--> apex(x,y,z UTM)
#     --agl_guard(ground_dtm.tif)--> apex(z AGL) --score_plot.
#
# Usage:
#   Rscript scripts/detect_forestformer3d_sweep.R [SITE=SOAP] [PLOTS=ALL]
#     [SPACING=24] [MERGE_TOL=2.0] [TOL=4] [IMAGE=ff3d-sm120]
#     [REPO=<abs>] [CKPT=<abs>] [TIMEOUT=3600]
# Output: $CLAUDE_JOB_DIR/neon/<SITE>/forestformer3d_results.csv (row per plot x rung).
suppressMessages({ library(lidR); library(data.table) })
options(lidR.progress = FALSE)
d <- Sys.getenv("CLAUDE_JOB_DIR", file.path(getwd(), "work"))
bs <- Find(file.exists, c(
  file.path("scripts", "bootstrap.R"),
  file.path("..", "..", "scripts", "bootstrap.R"),
  file.path(getwd(), "scripts", "bootstrap.R")))
if (!length(bs)) stop("bootstrap.R not found", call. = FALSE)
source(bs[1]); rm(bs)
source(.find("sweep_lib.R")); source(.find("model_bench_lib.R"))
source(.find("model_runner.R")); source(.find("io_bridge.R"))

args     <- strsplit(commandArgs(TRUE), "=")
A        <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
SITE     <- if (is.null(A$SITE))  "SOAP" else A$SITE
PLOTS    <- if (is.null(A$PLOTS) || A$PLOTS == "ALL") NULL else strsplit(A$PLOTS, ",")[[1]]
SPACING  <- as.numeric(if (is.null(A$SPACING)) 24 else A$SPACING)
MERGE_TOL<- as.numeric(if (is.null(A$MERGE_TOL)) 2.0 else A$MERGE_TOL)
TOL      <- as.numeric(if (is.null(A$TOL)) 4.0 else A$TOL)
IMAGE    <- if (is.null(A$IMAGE)) "ff3d-sm120" else A$IMAGE
TIMEOUT  <- as.numeric(if (is.null(A$TIMEOUT)) 3600 else A$TIMEOUT)
RADIUS   <- 16; RUNGS <- c(8); MINTREES <- 6
REPO  <- if (is.null(A$REPO)) file.path(.ROOT, "gpu/store/forestformer3d/ForestFormer3D") else A$REPO
CKPT  <- if (is.null(A$CKPT)) file.path(REPO, "work_dirs/clean_forestformer/epoch_3000_fix.pth") else A$CKPT
ENTRY <- file.path(.ROOT, "gpu/forestformer3d-sm120/ff3d_entry.sh")
PATCH <- file.path(REPO, "ff3d_repo.patch")          # the #27 patch, shipped in the repo
DRIVER<- file.path(.ROOT, "gpu/forestformer3d-sm120/ff3d_arm.py")

# Cylinder centers on a square grid covering [-ph, ph]^2 about the plot center,
# spacing SPACING; each cylinder processes radius RADIUS (overlap = 2*RADIUS-SPACING).
cyl_centers <- function(cx, cy, ph, spacing) {
  k <- max(1L, ceiling((2 * ph) / spacing) + 1L)
  off <- seq(-ph, ph, length.out = k)
  g <- expand.grid(dx = off, dy = off)
  data.frame(cx = cx + g$dx, cy = cy + g$dy)
}

run_main <- function() {
  stopifnot(file.exists(ENTRY), file.exists(DRIVER), file.exists(CKPT), file.exists(REPO))
  nd  <- file.path(d, "neon", SITE)
  gt  <- read.csv(file.path(nd, "ground_truth_stems.csv"), stringsAsFactors = FALSE)
  pc  <- read.csv(file.path(nd, "plot_centroids.csv"),     stringsAsFactors = FALSE)
  gt  <- gt[gt$live & gt$is_tree & !is.na(gt$E), ]
  laz <- list.files(file.path(nd, "lidar"), pattern = "\\.laz$",
                    recursive = TRUE, full.names = TRUE)
  ctg <- readLAScatalog(laz, progress = FALSE)
  counts <- table(gt$plotID); keep <- names(counts)[counts >= MINTREES]
  if (!is.null(PLOTS)) keep <- intersect(keep, PLOTS)
  keep <- intersect(keep, pc$plotID)
  cat(sprintf("[%s] forestformer3d plots: %d (image=%s spacing=%g merge_tol=%g)\n",
              SITE, length(keep), IMAGE, SPACING, MERGE_TOL))

  out <- list()
  for (pid in keep) {                       # SERIAL -- one GPU
    ci <- pc[pc$plotID == pid, ][1, ]
    cx <- ci$easting; cy <- ci$northing; ph <- plot_half(ci$plotType)
    stems <- gt[gt$plotID == pid & abs(gt$E - cx) <= ph & abs(gt$N - cy) <= ph, ]
    if (nrow(stems) < 1) next
    native_pdens <- NA_real_; ncell <- 0L
    for (rung in c(NA, RUNGS)) {
      prep <- tryCatch(frozen_clip(ctg, SITE, pid, rung, cx, cy, ph,
                                   out_root = file.path(nd, "frozen")),
                       error = function(e) NULL)
      if (is.null(prep)) next
      pdens <- prep$pdens; frdens <- prep$frdens
      if (is.na(rung)) native_pdens <- pdens
      else if (is.na(native_pdens) || rung >= native_pdens) next
      tag <- ifelse(is.na(rung), "native", as.character(rung))
      # tile the raw clip into cylinders
      raw <- tryCatch(lidR::readLAS(prep$rawground), error = function(e) NULL)
      if (is.null(raw) || lidR::is.empty(raw)) next
      in_dir <- file.path(tempdir(), sprintf("ff3d_%s_%s", pid, tag))
      unlink(in_dir, recursive = TRUE); dir.create(in_dir, recursive = TRUE)
      cc <- cyl_centers(cx, cy, ph, SPACING); n_cyl <- 0L
      for (i in seq_len(nrow(cc))) {
        cyl <- lidR::clip_circle(raw, cc$cx[i], cc$cy[i], RADIUS)
        if (lidR::is.empty(cyl) || lidR::npoints(cyl) < 50) next
        lidR::writeLAS(cyl, file.path(in_dir, sprintf("cyl_%03d.laz", n_cyl)))
        n_cyl <- n_cyl + 1L
      }
      if (n_cyl == 0L) next
      out_laz <- file.path(tempdir(), sprintf("ff3d_%s_%s.laz", pid, tag))
      # run_docker_arm passes no env, so REPO/PATCH/DRIVER ride in `extra` as
      # positional args (entry.sh reads $3..$6); all are identity-mounted.
      det_abs <- run_docker_arm(IMAGE, in_dir, out_laz,
                   cmd    = c("bash", ENTRY),
                   extra  = c(CKPT, REPO, PATCH, DRIVER),
                   mounts = c(REPO, dirname(CKPT), dirname(ENTRY)),
                   reader = function(p) ff3d_collapse(p, merge_tol = MERGE_TOL),
                   gpus = "all", timeout = TIMEOUT,
                   label = sprintf("%s/%s", pid, tag))
      if (is.null(det_abs)) next            # container crash/schema -> skip cell
      det <- agl_guard(det_abs, prep$dtm)
      if (is.null(det)) next                # wholesale off-DTM (frame bug) -> skip
      sc <- tryCatch(score_plot(stems, det, tol_xy = TOL, core_cx = cx,
                                core_cy = cy, core_half = ph),
                     error = function(e) NULL)
      if (is.null(sc)) next
      out[[length(out) + 1]] <- cbind(data.frame(site = SITE, plot = pid,
        plotType = ci$plotType, detector = "forestformer3d", rung = tag,
        pdens = round(pdens, 2), frdens = round(frdens, 2),
        n_cyl = n_cyl, n_apex = nrow(det)), sc)
      ncell <- ncell + 1L
    }
    cat(sprintf("  %s: %d cells\n", pid, ncell))
  }
  results <- do.call(rbind, out)
  if (is.null(results) || !nrow(results)) {
    cat("no forestformer3d results\n"); return(invisible())
  }
  results$tp_core <- round(results$precision * results$n_det)
  write.csv(results, file.path(nd, "forestformer3d_results.csv"), row.names = FALSE)
  cat(sprintf("[%s] forestformer3d DONE: %d rows -> %s\n", SITE, nrow(results),
              file.path(nd, "forestformer3d_results.csv")))
}

if (sys.nframe() == 0L) run_main()
```

- [ ] **Step 2: Lint-load the script (parse check, no run)**

Run: `Rscript -e 'invisible(parse("scripts/detect_forestformer3d_sweep.R")); cat("PARSE_OK\n")'`
Expected: `PARSE_OK`.

- [ ] **Step 3: Confirm the patch ships in the repo**

Run:

```bash
ls "$REPO"/ff3d_repo.patch && echo PATCH_PRESENT   # $REPO = the FF3D checkout
```

Expected: `PATCH_PRESENT`. If absent, copy
`gpu/forestformer3d-sm120/ff3d_repo.patch` into the repo (the patch is
repo-relative) so `ff3d_entry.sh` can apply it.

- [ ] **Step 4: Commit**

```bash
git add scripts/detect_forestformer3d_sweep.R
git commit -m "feat(#18): detect_forestformer3d_sweep.R density-ladder arm

Tiles each SOAP plot core into 16 m-radius cylinders, runs FF3D once per plot in
the ff3d-sm120 container, then ff3d_collapse (cross-block dedup + reduce) +
agl_guard + score_plot. Native + 8 only; mirrors detect_segmentanytree_sweep.R.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 6: Analyzer integration — additive native+8 FF3D comparison (R, TDD)

**Files:**

- Modify: `scripts/analyze_model_benchmark.R` (register CSV at line ~97; add
  `FF3D_CMP_ARMS` + `pool_ff3d()` near line 71; emit in `run_main`)
- Modify: `tests/testthat/test-model-benchmark-analysis.R` (add the no-shrink test)
- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-model-benchmark-analysis.R`:

```r
test_that("pool_ff3d compares native+8 without shrinking the full ladder", {
  # full-ladder arm at all 5 rungs on 2 plots; FF3D at native+8 only on plot p1
  mk <- function(det, plot, rung) data.frame(
    site = "SOAP", plot = plot, rung = rung, detector = det,
    n_ref = 10, n_det = 9, TP = 8, precision = 0.8, frdens = 5,
    recall = 0.8, F1 = 0.8, tp_core = 7,
    n_dominant = 4, rec_dominant = 0.9, n_codominant = 3, rec_codominant = 0.8,
    n_intermediate = 2, rec_intermediate = 0.5, n_suppressed = 1, rec_suppressed = 0.3)
  rows <- list()
  for (rg in c("native","8","4","2","1")) for (p in c("p1","p2"))
    rows[[length(rows)+1]] <- mk("chm_vwf", p, rg)
  rows[[length(rows)+1]] <- mk("forestformer3d", "p1", "native")
  rows[[length(rows)+1]] <- mk("forestformer3d", "p1", "8")
  u <- harmonize_union(rows)
  # the existing ladder keeps all 5 rungs for chm_vwf (FF3D not in LADDER_ARMS)
  ladder <- pool_arms(u, arms = intersect(LADDER_ARMS, unique(u$detector)))
  expect_setequal(as.character(unique(ladder$rung)), c("native","8","4","2","1"))
  # FF3D comparison exists at native+8 only
  ff <- pool_ff3d(u)
  expect_true("forestformer3d" %in% ff$detector)
  expect_setequal(as.character(unique(ff$rung)), c("native","8"))
})
```

- [ ] **Step 2: Run it; verify it fails**

Run: `Rscript tests/run_tests.R 2>&1 | grep -A2 model-benchmark-analysis`
Expected: FAIL — `could not find function "pool_ff3d"`.

- [ ] **Step 3: Add `FF3D_CMP_ARMS` + `pool_ff3d()`**

In `analyze_model_benchmark.R`, after the `NATIVE_ARMS` line (line 71), add:

```r
# FF3D (#M8) is native+8-only: compare it against the baselines in its OWN pool so
# it never enters LADDER_ARMS' equal-set guard (which would drop rungs 4/2/1 for
# every arm). Self-contained -> FF3D's plot coverage affects only this table.
FF3D_CMP_ARMS <- c("forestformer3d", "chm_vwf", "treeisonet")
pool_ff3d <- function(u) {
  arms <- intersect(FF3D_CMP_ARMS, unique(u$detector))
  if (!"forestformer3d" %in% arms) return(NULL)
  pool_arms(u, arms = arms, rungs = c("native", "8"))
}
```

- [ ] **Step 4: Register the CSV + emit the table in `run_main`**

In `run_main`, change `arm_files` (line 97-98) to include FF3D:

```r
  arm_files <- c(ams3d = "ams3d_results.csv", lidrplugins = "lidrplugins_results.csv",
                 li2012 = "li2012_results.csv", treeisonet = "treeisonet_results.csv",
                 forestformer3d = "forestformer3d_results.csv")
```

After the `ladder`/`native`/`dl` block (after line 110), add:

```r
  ff3d <- pool_ff3d(u)
  if (!is.null(ff3d)) {
    write.csv(ff3d, file.path(nd, "model_bench_ff3d.csv"), row.names = FALSE)
    cat("forestformer3d native+8 comparison written\n")
  }
```

And append its fragment to `frag` (inside the `c(...)` that builds the markdown,
before `writeLines`):

```r
            , if (!is.null(ff3d)) c("",
                "#### ForestFormer3D (#M8) native + 8, vs baselines", "",
                .md_table(ff3d[order(ff3d$detector, ff3d$rung), ], cols_l)) else NULL
```

- [ ] **Step 5: Run tests; verify pass**

Run: `Rscript tests/run_tests.R 2>&1 | tail -20`
Expected: the new analysis test passes; suite green.

- [ ] **Step 6: Commit**

```bash
git add scripts/analyze_model_benchmark.R tests/testthat/test-model-benchmark-analysis.R
git commit -m "feat(#18): additive FF3D native+8 comparison in the analyzer

Register forestformer3d_results.csv and add pool_ff3d() (FF3D vs chm_vwf/
treeisonet at native+8 in its own equal-set pool). LADDER_ARMS/NATIVE_ARMS are
untouched, so the full 5-rung ladder is not shrunk by FF3D's native+8-only rows.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 7: Live smoke + the SOAP benchmark run (GPU)

**Files:** none modified (produces `work/neon/SOAP/forestformer3d_results.csv`).

- [ ] **Step 1: One-cylinder live smoke (de-risk the offset round-trip)**

```bash
cd /home/alex/projects/lidar_tree_benchmarks/.claude/worktrees/gpu-ff3d-arm
export CLAUDE_JOB_DIR=$PWD/work
Rscript scripts/detect_forestformer3d_sweep.R SITE=SOAP PLOTS=$(\
  Rscript -e 'd=read.csv("work/neon/SOAP/plot_centroids.csv");cat(d$plotID[1])') \
  2>&1 | tail -30
```

Expected: `forestformer3d DONE: N rows` with N ≥ 1, and a `forestformer3d_results.csv`.
Inspect that detection coords are UTM (x ≈ plot easting ~3e5, y ≈ northing ~4.1e6),
not centered near 0 — proving the offset restore. If coords are near 0, the
container did not re-add the offset — fix `ff3d_arm.py` step 3 first.

- [ ] **Step 2: Run the full SOAP benchmark (native + 8, all plots)**

```bash
cd /home/alex/projects/lidar_tree_benchmarks/.claude/worktrees/gpu-ff3d-arm
export CLAUDE_JOB_DIR=$PWD/work
nohup Rscript scripts/detect_forestformer3d_sweep.R SITE=SOAP \
  > work/neon/SOAP/ff3d_fullrun.log 2>&1 &
echo "PID $!  — tail -f work/neon/SOAP/ff3d_fullrun.log"
```

Wait for the log to print `forestformer3d DONE: N rows` (and the CSV path).
This is the heaviest arm; budget time accordingly — one model load plus a
cylinder loop per plot×rung.

- [ ] **Step 3: Tune `SPACING`/`MERGE_TOL` if needed**

Inspect `forestformer3d_results.csv`: if recall is low because tall trees fall in
cylinder gaps, lower `SPACING` (more overlap) and re-run; if precision is low from
seam double-counts surviving, raise `MERGE_TOL` slightly. Record the final values.

- [ ] **Step 4: Run the analyzer**

```bash
cd /home/alex/projects/lidar_tree_benchmarks/.claude/worktrees/gpu-ff3d-arm
export CLAUDE_JOB_DIR=$PWD/work
Rscript scripts/analyze_model_benchmark.R SITE=SOAP 2>&1 | tail -10
cat work/neon/SOAP/model_bench_ff3d.csv
```

Expected: `model_bench_ff3d.csv` with `forestformer3d` rows at native + 8, and the
existing `model_bench_ladder.csv` still carrying all 5 rungs.

- [ ] **Step 5: Commit the run log (results CSVs are gitignored)**

```bash
git add -f work/neon/SOAP/ff3d_fullrun.log 2>/dev/null || true
git commit -m "chore(#18): FF3D SOAP run log (native+8)" --allow-empty || true
```

(Generated CSVs stay gitignored per the repo convention — only the doc carries
the numbers.)

---

## Task 8: Results doc (un-defer #M8) + README + memory + final suite

**Files:**

- Modify: `results/model-benchmark-results.md` (replace the #M8-deferred lines)
- Modify: `gpu/forestformer3d-sm120/README.md` (tick "benchmark arm" done)
- Update memory `gpu-arm-blackwell-sm120` (FF3D now runs + scored on sm_120)
- [ ] **Step 1: Replace the #M8-deferred passages**

In `results/model-benchmark-results.md`, find the lines that defer #M8 (around
lines 10, 254, 268, 283 — `grep -n "#M8\|ForestFormer3D" results/model-benchmark-results.md`)
and replace with a ForestFormer3D subsection carrying the pooled native + 8 table
from `model_bench_ff3d.csv` (recall / precision / F1 + per-crown-class recall vs
chm_vwf/treeisonet), the tuned `SPACING`/`MERGE_TOL` + cylinder counts, and the
zero-shot domain-gap caveat (FF3D trains on dense ULS/TLS; SOAP is sparse ALS —
consistent with the `gpu-arm-blackwell-sm120` note). Keep prose ≤ 80 cols.

- [ ] **Step 2: Lint the doc**

Run: `rumdl check results/model-benchmark-results.md`
Expected: `Success: No issues found`.

- [ ] **Step 3: Tick the README next-step**

In `gpu/forestformer3d-sm120/README.md`, change the "a full NEON benchmark run …
is follow-up" line to state the benchmark arm now exists
(`scripts/detect_forestformer3d_sweep.R`, SOAP native + 8).

- [ ] **Step 4: Update the memory note**

Update `gpu-arm-blackwell-sm120` memory: FF3D (#M8) now builds + runs + is scored
zero-shot on sm_120 (RTX 5090) via `ff3d-sm120` + `detect_forestformer3d_sweep.R`;
the Dockerfile build fixes (`--index-url`, transitive backfill, `nvidia-arch`).

- [ ] **Step 5: Full suite + commit**

```bash
cd /home/alex/projects/lidar_tree_benchmarks/.claude/worktrees/gpu-ff3d-arm
Rscript tests/run_tests.R 2>&1 | tail -8
git add results/model-benchmark-results.md gpu/forestformer3d-sm120/README.md
git commit -m "docs(#18): un-defer #M8 — ForestFormer3D SOAP native+8 results

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

Expected: suite green; doc committed.

---

## Self-Review notes (for the executor)

- **Spec coverage:** Task 1 ↔ §3 dedup; Task 2 ↔ §4 reader + §1 off-DTM guard;
  Tasks 3–5 ↔ §2 container + §1 packaging/checkpoint; Task 6 ↔ §7 analyzer;
  Tasks 7–8 ↔ §7 live run + results. All §-requirements have a task.
- **Type consistency:** `dedup_blocks()` returns `block,inst,X,Y,Z,global_id`;
  `reduce_instances(id_col = "global_id")` consumes `global_id`; `ff3d_collapse`
  feeds it `X/Y/Z` (capitalized) from `las@data`. `agl_guard` returns `det` or
  `NULL`. The container writes `point_source_id`/`user_data` (laspy); the reader
  reads `PointSourceID`/`UserData` (lidR) — intentional, per side.
- **No placeholders:** every code step is complete. The only deferred value is
  the live `SPACING`/`MERGE_TOL` tuning (Task 7 step 3), which the spec
  explicitly leaves empirical.

```text
```
