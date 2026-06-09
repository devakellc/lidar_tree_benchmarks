# SegmentAnyTree on Blackwell (sm_120) — #17 / #M6

Rebuild of the SegmentAnyTree Docker stack for the RTX 5090. The published
image `maciekwielgosz/segment-any-tree:latest` is torch 1.9 / CUDA 11.1 and
cannot run on sm_120.

## Status

**Inference runs end-to-end on the 5090** (verified 2026-06-09, SOAP-native
frozen clip): full pipeline → merged `.laz` with per-point `PredInstance`
(0=non-tree, 1..N) → host-side `read_instances_laz` → `reduce_instances` →
`det_to_agl` gave 35 tree apexes at AGL 3.9–53.5 m. The sm_120 MinkowskiEngine
build was always sound (zero `no kernel image` errors); the work was porting the
2023-era SAT/torch_points3d code to the modern torch-2.7 / cu128 / py3.11 stack
Blackwell forces — a set of dep pins + compat shims (see **Inference bring-up**).

| Layer | Status |
|---|---|
| MinkowskiEngine 0.5.4 + sm_120 | **Done** — `gpu/minkowski-sm120/`; kernels verified in the real model forward |
| PyG companions on torch 2.7 cu128 | **Done** — official wheels; `torch-geometric==2.3.1` pinned |
| torch-points-kernels (region_grow) on sm_120 | **Done** — patched numpy pin + source build |
| torch-points3d + SAT inference | **Done** — full forward on the 5090 |
| Docker backend + LAZ↔PLY bridge | **Done** — `run_docker_arm` + `laz_to_ply`/`read_instances_laz` (#19) |
| Headless driver + sweep harness | **Done** — `gpu/run_segmentanytree.py` + `scripts/detect_segmentanytree_sweep.R` |
| End-to-end inference → labeled cloud → apex det | **Done** — validated host-side on real output |

## Build

```sh
# 1) MinkowskiEngine base (~5 min compile)
bash gpu/build_me_sm120.sh

# 2) SegmentAnyTree layer (PyG 2.3.1 pin + hdbscan + pyarrow + SAT clone)
bash gpu/segmentanytree-sm120/build.sh

# 3) Import smoke on the 5090
docker run --rm --gpus all sat-sm120-test
```

Run an inference smoke (folder-in/out wrapped to single-file by the driver):

```sh
docker run --rm --gpus all --shm-size=8g --ipc=host \
  -v /abs/io:/abs/io -v "$PWD/gpu":"$PWD/gpu" sat-sm120-test \
  python3 "$PWD/gpu/run_segmentanytree.py" /abs/io/clip.ply /abs/io/out.las
```

## Inference bring-up (what it took)

A SOAP-native clip through `run_inference.sh` failed in turn on each place the
SAT stack assumed its pinned (old) deps. All five are fixed:

1. **Hydra 1.3** rejects the absolute `--config-name` upstream passes, and an
   extra `--config-dir` is outranked by the in-tree default config. The driver
   rewrites that line to copy the `modify_eval.py`-edited config over
   `conf/eval.yaml` and load it by bare name.
2. **omegaconf ≥ 2.2** raises on a missing attribute where 2.0 returned `None`
   (the `self._cfg.<group>` idiom, pervasive in torch_points3d). Shim restores it.
3. **torch ≥ 2.6** defaults `torch.load(weights_only=True)`, refusing the
   omegaconf objects pickled in the trusted checkpoint → set `False`.
4. **torch-geometric** `Data.keys` became a method (`set(data.keys)` →
   TypeError). Pinned `torch-geometric==2.3.1` (last property-`keys` release;
   imports fine on torch 2.7). Also added `hdbscan` (model import) and `pyarrow`
   (dask merge step) — both were simply missing.
5. **torch 2.7** strict cross-device indexing (`cpu_tensor[gpu_index]`) +
   **CUDA-after-`fork`** deadlock in the clustering `multiprocessing.Pool`. Fixed
   by a `Tensor.__getitem__` device shim and forcing the `spawn` start method;
   run with `--shm-size=8g --ipc=host`.

Fixes 2, 3, 5 live in `gpu/sat_compat/usercustomize.py` (auto-applied via
PYTHONPATH by the driver — upstream files stay unmodified); fix 1 + the output
glob are in the driver; fix 4 is in the Dockerfile.

## Remaining to ship the arm

- [ ] Rebuild the image so the pins are baked: `bash gpu/segmentanytree-sm120/build.sh`
  (the current local `sat-sm120-test` predates the Dockerfile pins; bring-up used
  runtime `pip install`).
- [ ] Run `scripts/detect_segmentanytree_sweep.R SITE=SOAP` (one cell ≈ 4 min).
- [ ] Add `"segmentanytree"` to `LADDER_ARMS` in `analyze_model_benchmark.R`
  *after* a results CSV exists (the equal-set guard drops all cells otherwise).
