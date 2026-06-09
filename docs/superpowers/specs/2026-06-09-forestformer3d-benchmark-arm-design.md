# ForestFormer3D Benchmark Arm (#M8) — Design

**Status.** Approved design (brainstorming output, revised after a design
review). Ready for implementation. Last updated 2026-06-09. Implements
[issue #18](https://github.com/agrigoriev/lidar_tree_benchmarks/issues/18)
(#M8). Extends
[PR #27](https://github.com/agrigoriev/lidar_tree_benchmarks/pull/27) (the
sm_120 container recipe) and builds on the merged #19 Docker backend +
LAZ↔PLY bridge.

**Scope.** Turn the #27 *recipe* (FF3D runs end-to-end on one staged cloud)
into a *benchmark arm*: tile each SOAP plot into 16 m-radius cylinders, run
FF3D zero-shot one cylinder at a time, dedup instances across cylinders, reduce
to apex detections, and score against NEON field stems with the existing
harness. **SOAP only, native + 8 pts/m² rungs only** (per the issue — FF3D is
the heaviest arm). Deliverable is the arm + an analyzer integration + unit
tests + a **live SOAP run** folded into `results/model-benchmark-results.md`
(un-deferring #M8).

This spec assumes the working `ff3d-sm120` image from the repaired #27
Dockerfile (three build fixes already committed on this branch: `--index-url`
for the rathaROG cumm/spconv wheels, the `--no-deps` transitive backfill
including `pybind11`/`rich`, and `nvidia-arch==7.1.0` for spconv).

The driver mirrors **`detect_segmentanytree_sweep.R`** (the existing Docker arm:
`laz`→clip→`run_docker_arm`→labelled cloud→`reduce`→`det_to_agl`→`score_plot`,
sourcing `io_bridge.R`), **not** the venv TreeisoNet arm.

---

## 1. Locked decisions

Settled during brainstorming and a follow-up design review:

| Decision | Choice | Consequence |
|---|---|---|
| Cross-block dedup | **Apex-cluster union-find, DIFFERENT blocks only** | Per-(block,instance) apex → union pairs **from different cylinders** within `merge_tol` → one global tree id. Merges only the cross-tile seam artifact; **never** merges two instances of the *same* cylinder (that would launder the model's genuine within-tile over-segmentation and inflate precision). |
| Reducer ordering | **Dedup BEFORE the reducer** | `dedup_blocks()` assigns globally-consistent ids; the canonical `reduce_instances()` then collapses to apex **once**. A tree split across a seam becomes one id → its apex = the true (highest) treetop. |
| Tiling unit | **R clips deterministic 16 m-radius cylinders** over the plot core | FF3D's test pipeline has **no `CylinderCrop`** (config L148–163) — it ingests the whole staged scene, so *we* control extent. lidR `clip_circle` at known UTM centers is testable; centers also key the dedup. |
| Tiling granularity | **Tune `SPACING`/`merge_tol` empirically** | Grid spacing < 2·radius gives overlap (seam context). `SPACING=24` (8 m overlap at radius 16) is the starting default; the chosen values go in the result doc, not pinned here. |
| Model-load cost | **One container call per plot; model loaded once** | The container stages *all* of a plot's cylinders as scenes and loops the dataloader (~40 calls, not ~240 per-cylinder reloads). |
| Output shape | **One merged LAZ per plot** (block id in `user_data`) | Keeps `run_docker_arm`'s single-output-file contract intact. Overlap-duplicated points are fine — `dedup_blocks` + max-Z reduce absorb them. |
| Coordinate frame | **Container restores per-cylinder centering offset → emits UTM** | `prep_test_data` centers by (mean-x, mean-y, min-z); the driver re-adds the saved per-scene offset so output is UTM. Matches the #19 "emit in the frame you received" contract. The driver adds a wholesale-off-DTM guard (§4) so a frame bug becomes a *skipped* cell, not a fake 0-recall row. |
| Checkpoint path | **Explicit absolute arg + identity-mounted weights dir** | `run_infer_save.py`'s relative `load_from` won't resolve under an identity mount. `ff3d_arm.py` takes the checkpoint as a 3rd positional arg (`cfg.load_from = argv[3]`); the driver passes it via `run_docker_arm(extra=)` and mounts its dir via `mounts=`. |
| Input variant | **Raw-with-ground clip** (`clip_rawground.laz`) | FF3D has its own `ground` semantic class; feed the raw clip, then z→AGL via `det_to_agl` + the frozen clip's DTM (exactly as SAT does). |
| Image packaging | **Mount FF3D repo + weights at runtime; bake only the 3 mm-stack site-packages swaps; mount driver + entry script** | The FF3D repo is gitignored runtime data — baking it couples the image to a data dir. The 3 `replace_mmdetection_files` swaps modify *site-packages* (in the image), so they're baked at build (vendored under `gpu/forestformer3d-sm120/`). `ff3d_arm.py` + `ff3d_entry.sh` are version-controlled there and mounted; the entry script applies `ff3d_repo.patch` to the mounted repo **idempotently** (`git apply --check` guard) so repeated `--rm` calls are safe. Staging writes to the ephemeral layer. |
| Analyzer integration | **Register the CSV; report FF3D in a dedicated native+8 comparison; DO NOT add it to `LADDER_ARMS`/`NATIVE_ARMS`** | The full-ladder `equal_set_guard` requires every arm at every rung; a native+8-only arm would silently drop rungs 4/2/1 for all arms. FF3D gets its own self-contained pooled table. |
| Results | **Un-defer #M8 inside `results/model-benchmark-results.md`** | One arm = consistent with the other model arms already in that doc; no standalone result file. |

---

## 2. Component: container CLI `ff3d_arm.py`

Refactor of `run_infer_save.py` into a headless, contract-shaped driver, baked
at `/workspace/ff3d_arm.py` (source under `gpu/forestformer3d-sm120/`).

**Interface:** `python ff3d_arm.py <in_dir> <out_laz> <ckpt>`

- `<in_dir>` — directory of cylinder clips `cyl_NNN.laz` (UTM), from the driver.
- `<out_laz>` — a single merged LAZ: every cylinder's labelled points in **UTM**
  (input frame restored), the model's per-point instance id in `point_source_id`,
  the **cylinder index in `user_data`** (the "block" dedup keys on), and the FF3D
  instance score in an extra dim (diagnostic only — not thresholded; we take the
  model's own per-point assignment, as `run_infer_save` did).
- `<ckpt>` — absolute path to `epoch_3000_fix.pth`; sets `cfg.load_from`.

**Flow (model loaded once):**

1. For each `cyl_NNN.laz`: stage it as a ForAINetV2 scene (the `prep_test_data`
   bypass — dummy sem/ins/bbox), saving the per-scene centering offset
   `(mean_x, mean_y, min_z)` to a sidecar. All scenes listed in `test_list.txt`.
2. `tools/create_data_forainetv2.py` builds `points/*.bin` + the **test** pkl
   for all scenes (only the test pkl — `update_pkl_infos` crashes on the empty
   train/val pkls, per #27).
3. Build the runner once (`Runner.from_cfg`, `load_or_resume`); iterate the test
   dataloader. For each scene's `pred_pts_seg`: take per-point instance ids,
   re-add that scene's offset, tag the cylinder index as `user_data`; accumulate
   all scenes and write a single merged `<out_laz>` (UTM).
4. Print a per-scene one-line summary + a final `ARM_DONE` sentinel.

**Crash discipline:** a scene yielding no points contributes no rows (its
cylinder is absent from the merge). A total failure (no output file, non-zero
exit) is caught by `run_docker_arm` → cell skipped.

The #27 correctness patches split by what they touch: the two mmengine/mmdet3d
**site-packages** swaps are baked at **build time**; `ff3d_repo.patch` (the
`tools/test.py` spconv-`permute` disable + tensorboard-hook drop) is applied to
the **mounted repo** at container start by `ff3d_entry.sh`, idempotently
(`git apply --check` guard).

## 3. Component: cross-block dedup `dedup_blocks()` (in `model_bench_lib.R`)

The load-bearing new R helper. **Pure, GPU-free, unit-tested.**

```text
dedup_blocks(pts, merge_tol = 2.0,
             block = "block", id = "inst", x = "X", y = "Y", z = "Z")
  -> data.frame with a globally-consistent integer "global_id" column
```

- Input `pts`: stacked per-cylinder labelled points — `block` (cylinder index),
  `id` (per-cylinder instance id; `0`/`NA` = unassigned), `X,Y,Z` (UTM).
- Step 1 — per `(block, id)` apex (max-Z point).
- Step 2 — **union-find over apexes, restricted to pairs in DIFFERENT blocks**
  within horizontal distance ≤ `merge_tol`. Same-block instances are the model's
  own output and are **left distinct** (honest over-segmentation scoring). A
  cross-block link can still transitively bind two same-block pieces *via* a
  neighbour's single detection — acceptable, since that link is positive
  evidence the pieces are one tree.
- Step 3 — assign each cluster a `global_id`; map every point's `(block,id)` to
  its `global_id`; drop unassigned points.
- **Returns** the relabelled point table. The driver then calls
  `reduce_instances(.., id_col = "global_id")` → apex per global tree (max-Z over
  all merged blocks' points, robust to overlap-duplicated points).

Post-condition (assertable): no two distinct `global_id`s **from different
blocks** have apexes within `merge_tol`.

## 4. Component: R driver `detect_forestformer3d_sweep.R`

Mirrors `detect_segmentanytree_sweep.R` (serial — one GPU). **Sources
`sweep_lib.R`, `model_bench_lib.R`, `model_runner.R`, and `io_bridge.R`** (the
last for `det_to_agl`). Per plot × rung in `{native, 8}`:

1. `frozen_clip(... rung ...)` → `clip_rawground.laz` + `ground_dtm.tif` (reuse
   the cached clip; no-upsampling guard `rung >= native_pdens → skip`).
2. **Tile**: cylinder centers on a grid (spacing `SPACING`) covering the plot
   core (`plot_half(plotType)`); `lidR::clip_circle(radius = 16)` each from the
   raw clip; write `cyl_NNN.laz` into a per-cell temp `in_dir`.
3. **Infer**: `run_docker_arm(image = IMAGE, input = in_dir, out_csv = out_laz,
   cmd = c("bash", ENTRY_SH), extra = c(CKPT_ABS),
   mounts = c(FF3D_REPO, dirname(CKPT_ABS), dirname(ENTRY_SH)),
   reader = <closure>, gpus = "all", timeout)`. `ENTRY_SH` cds to `FF3D_REPO`,
   applies `ff3d_repo.patch` idempotently, then `exec python ff3d_arm.py "$@"`
   (so the container receives `<in_dir> <out_laz> <ckpt>`). The **reader
   closure** `lidR::readLAS(out_laz)` → `(X, Y, Z, inst = PointSourceID,
   block = UserData)` → `dedup_blocks()` → `reduce_instances(id_col =
   "global_id")` → UTM apexes `det(x, y, z)`. `run_docker_arm` keeps its
   single-file contract and asserts `det`.
4. **Score** (with the off-DTM guard):
   - `det_abs <- run_docker_arm(...)`; if `is.null(det_abs)` → skip (crash/schema).
   - `nrow(det_abs) == 0` → genuine ran-but-empty → score a **0-row** (legit
     recall 0).
   - else `det <- det_to_agl(det_abs, dtm)`; if `nrow(det) == 0 &&
     attr(det, "n_dropped") > 0` → **wholesale off-DTM drop = a frame bug → skip
     the cell** (not a fake 0-recall row). Otherwise `score_plot(stems, det,
     tol_xy = TOL, core_cx, core_cy, core_half)`.
5. Emit one row per (plot, rung): `site, plot, plotType,
   detector = "forestformer3d", rung, pdens, frdens, n_cyl, n_apex` + the
   `score_plot` columns; `tp_core = round(precision * n_det)`. Write
   `work/neon/SOAP/forestformer3d_results.csv`.

Args (KEY=VALUE): `SITE=SOAP PLOTS=ALL MERGE_TOL=2.0 SPACING=24 TOL=4.0
IMAGE=ff3d-sm120 REPO=<abs FF3D repo> CKPT=<abs .pth> TIMEOUT=3600`. `REPO`/`CKPT`
default to the `gpu/store/forestformer3d/ForestFormer3D` checkout + its
`work_dirs/clean_forestformer/epoch_3000_fix.pth`. `RUNGS <- c(8);
MINTREES <- 6`; loop `c(NA, RUNGS)`.

## 5. Coordinate & contract chain (end to end)

`rawground.laz (UTM)` → R `clip_circle` (UTM) → container centers per cylinder,
saves offset → model (local) → container re-adds offset → merged `out.laz (UTM)`
→ reader → `dedup_blocks` (UTM) → `reduce_instances` (UTM apex, abs-Z) →
`det_to_agl` (z = AGL) + wholesale-off-DTM guard → `score_plot` vs stems (UTM
x/y, AGL z). The fragile step — the centering-offset restore — is owned by the
container, unit-tested via a round-trip, and backstopped by the §4 guard.

## 6. Tests (pure-R testthat, plus one gated GPU smoke)

`tests/testthat/test-dedup-blocks.R` (per-helper, matching the repo's layout):

- One tree in the overlap detected in **both** cylinders (apexes within
  `merge_tol`, different blocks) → **one** `global_id`.
- **Two SAME-block instances within `merge_tol`, with no cross-block bridge →
  stay DISTINCT** (the over-segmentation guard from finding 2).
- Two genuinely distinct trees → two ids; `0`/`NA` dropped; empty → 0-row.
- Post-condition: after dedup + `reduce_instances`, no two different-block
  apexes within `merge_tol`; the merged apex is the **max-Z** across blocks.

`tests/testthat/test-forestformer3d-sweep.R` (driver, container stubbed):

- A fake `docker` (reuse the #19 fake-docker pattern) whose "output" is a
  prebuilt merged labelled `out.laz` (two cylinders, `UserData` block ids) → the
  reader closure dedups + reduces → asserts a known apex set; verifies the
  `run_docker_arm` argv (identity mounts incl. the checkpoint dir, `cmd`, the
  `CKPT` `extra`, `in_dir`/`out_laz`).
- **Off-DTM guard**: a `det_abs` whose apexes all fall off the DTM → cell
  **skipped** (no row); a genuinely empty `det_abs` → a scored **0-row**.
- No-upsampling guard skips `rung >= native_pdens`; a missing output file → cell
  skipped.

**Gated live smoke** (`skip_if` no GPU / no image): one real cylinder through
`ff3d_arm.py` → non-empty UTM labelled LAZ whose coords fall in the clip's bbox
(the offset round-trip is real, not just self-consistent).

## 7. Analyzer integration + results (un-defer #M8)

**`analyze_model_benchmark.R` changes** (the existing pool/guard/figures stay
untouched; we only *add* an FF3D path so it can't shrink the full ladder):

1. Add `forestformer3d = "forestformer3d_results.csv"` to `arm_files` so the CSV
   is discovered and unioned.
2. Add `FF3D_CMP_ARMS <- c("forestformer3d", "chm_vwf", "treeisonet")` and a
   dedicated `ff3d <- pool_arms(u, arms = intersect(FF3D_CMP_ARMS, present),
   rungs = c("native", "8"))`, written to `model_bench_ff3d.csv` + a `.md_table`
   fragment. Its `equal_set_guard` is over just these arms at native+8, so FF3D's
   plot coverage only affects **this** table.
3. **Do NOT** add `forestformer3d` to `LADDER_ARMS` or `NATIVE_ARMS` — both keep
   their current arm sets and equal-set populations exactly.

**Results doc.** Replace the "#M8 deferred behind sm_120" lines in
`results/model-benchmark-results.md` with the FF3D native+8 comparison: pooled
recall / precision / F1 and per-crown-class recall vs `chm_vwf`/`treeisonet`,
with the zero-shot domain-gap caveat (FF3D trains on dense ULS/TLS; SOAP is
sparse ALS — consistent with the `gpu-arm-blackwell-sm120` note). Record the
tuned `SPACING`/`merge_tol` and cylinder counts.

## 8. Scope boundary (non-goals)

- **Crown/diameter** metrics for FF3D — the #M8 follow-on, not this issue.
  Apex-only here.
- Other sites (SJER/TEAK) and other rungs (4/2/1) — issue says SOAP, native + 8.
- No change to `score_plot`, `pool`, `equal_set_guard`, or the existing
  `LADDER_ARMS`/`NATIVE_ARMS` ladder/native tables (FF3D is additive).
- Re-training / fine-tuning FF3D — strictly zero-shot.
