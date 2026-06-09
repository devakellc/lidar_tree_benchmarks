# ForestFormer3D on Blackwell (sm_120)

ForestFormer3D (#M8, ICCV 2025, `SmartForest-no/ForestFormer3D`) ported to run on
the **RTX 5090** (sm_120) under **torch 2.7.0 / CUDA 12.8** — upstream is torch
1.13 / cu116 and does not run on Blackwell. Verified end-to-end: the model
loads the official `epoch_3000_fix.pth`, runs inference on a point-cloud plot, and
emits per-tree instance masks.

## Result (zero-shot on sparse NEON ALS)

A NEON SJER plot (31,676 pts, ~3.5 pts/m²) through the pretrained model:

```text
raw instance candidates: 101   instance score range: 0.023 .. 0.273
model per-point assignment: 8 trees (18.4% of points)
score>=0.10: 35 candidate trees   score>=0.15: 14
```

Scores are low because FF3D is trained on **dense ULS/TLS/MLS**; zero-shot on
sparse airborne ALS is expected weak (see memory note `gpu-arm-blackwell-sm120`).
The port is correct — the domain gap is the limiter, not the implementation.

## The four hard problems (and fixes)

1. **MinkowskiEngine voxelizer** (sm_120 wall). FF3D's `collate()` is pure ME and
   `oneformer3d.py` imports ME at module top. Solved by the shared
   [`../minkowski-sm120`](../minkowski-sm120) build (ME 0.5.4 for arch 12.0).
2. **spconv backbone** (SpConvUNet). `spconv-cu128` 2.4.1 + `cumm-cu128` 0.9.1
   from the rathaROG index run on sm_120 (see issue #18). `cumm-cu128` must come
   from rathaROG, not PyPI; and keep `ccimport>=0.4.4` / `pccm>=0.4.16` (older
   pins break `cumm`'s `IsAppleSiliconMacOs`).
3. **mmcv `_ext` build on torch 2.7.** No cu128 wheel exists, so it compiles from
   source — and mmcv 2.0.0 hardcodes `-std=c++14` while torch 2.7 headers need
   **C++17** (`jit_type.h operator*`, `c10` hash errors). `sed c++14 -> c++17` in
   `setup.py` fixes it. See `build_mmcv.sh`.
4. **spconv checkpoint layout.** `epoch_3000_fix.pth` is already in spconv-cu128
   2.4.1 layout `(out,K,K,K,in)`; FF3D's in-memory `permute(1,2,3,4,0)` in
   `tools/test.py` **re-breaks** it. Disable the permute (load as-is).

## Dependency-port notes

- **numpy < 2 track** (`numpy==1.26.4`, `numba==0.59.1`): the old mm-stack needs
  it; `spconv-cu128` + ME both import fine under it.
- **mm-stack** pinned to FF3D's: `mmengine 0.7.3`, `mmdet 3.0.0`, `mmsegmentation
  1.0.0`, `mmdet3d @ 22aaa47` — all `--no-deps`. Exact env in `ff3d_freeze.txt`.
- **opencv-python-headless** (slim base lacks libGL); **open3d** needs apt
  `libgl1 libgomp1 libx11-6 libxext6 libxrender1 libsm6 libglib2.0-0 libusb-1.0-0`
  plus `plotly`/`dash`; `mmdet3d.evaluation` eagerly imports `lyft_dataset_sdk`
  `nuscenes-devkit`. The Karbo123 **segmentator is NOT needed** — it is mesh-only
  and the point-cloud test pipeline computes voxel-superpoints in the model.
- **torch-scatter 2.1.2 / torch-cluster 1.6.3 / torch-points-kernels 0.7.0** built
  for arch 12.0.

## Data prep (no batch_load)

`batch_load_ForAINetV2_data.py` pulls in segmentator/open3d/Delaunay and forces
labelled mode. `prep_test_data.py` bypasses it: writes the
`forainetv2_instance_data/*_vert/_sem/_ins/_bbox.npy` that `converter_forainetv2`
consumes (dummy sem=0/ins=0, empty `(0,7)` bboxes → `gt_num=0`). Then
`tools/create_data_forainetv2.py` builds `points/*.bin` + the test pkl. (Update
only the **test** pkl — `update_pkl_infos` crashes on the empty train/val pkls.)

## Run

```sh
# build the image (FROM the proven ME base) — see Dockerfile
docker build -t ff3d-sm120 .
# inside a container with the FF3D repo + weights mounted:
python prep_test_data.py                       # stage a plot (SRC=*.laz)
python tools/create_data_forainetv2.py forainetv2 --root-path ./data/ForAINetV2
python run_infer_save.py                        # -> work_dirs/infer/pred_instances.laz
```

`ff3d_repo.patch` carries the four `tools/test.py` + config edits. `ff3d_freeze.txt`
is the authoritative version lock (177 packages) for the working container.

## Benchmark arm (#M8)

The headless arm is
[`scripts/detect_forestformer3d_sweep.R`](../../scripts/detect_forestformer3d_sweep.R):
it tiles each NEON plot into 16 m-radius cylinders, runs FF3D once per plot via
[`ff3d_arm.py`](ff3d_arm.py) + [`ff3d_entry.sh`](ff3d_entry.sh) in this image
(`run_docker_arm`), cross-block-dedups the instances (apex-cluster, different
blocks only), and scores apexes against field stems. SOAP native + 8 results
(zero-shot) are in
[results/model-benchmark-results.md](../../results/model-benchmark-results.md).

```sh
Rscript scripts/detect_forestformer3d_sweep.R SITE=SOAP REPO=<FF3D repo>   # ~8 min, one GPU
```

The `Dockerfile` here builds clean: the three `replace_mmdetection_files`
site-packages swaps are baked at build, and `ff3d_repo.patch` is applied to the
mounted repo at container start by `ff3d_entry.sh` (idempotent). `ff3d_arm.py` is
the multi-cylinder, single-model-pass driver; `run_infer_save.py` remains the
single-plot demo.
