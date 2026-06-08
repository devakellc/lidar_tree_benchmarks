# SegmentAnyTree on Blackwell (sm_120) — WIP port (#17 / #M6)

Rebuild of the SegmentAnyTree Docker stack for the RTX 5090. The published
image `maciekwielgosz/segment-any-tree:latest` is torch 1.9 / CUDA 11.1 and
cannot run on sm_120.

## Status

| Layer | Status |
|---|---|
| MinkowskiEngine 0.5.4 + sm_120 | **Done** — see `gpu/minkowski-sm120/` |
| PyG (scatter/sparse/cluster) on torch 2.7 cu128 | **Done** — official PyG wheels |
| torch-points-kernels (region_grow) on sm_120 | **Done** — patched numpy pin + source build |
| torch-points3d + SAT repo import | **Done** — `test_imports.py` smoke on 5090 |
| Full `run_inference.sh` pipeline | **Not yet** — needs LAZ→PLY bridge (#I4) |

## Build

```sh
# 1) MinkowskiEngine base (~5 min compile)
bash gpu/build_me_sm120.sh

# 2) SegmentAnyTree layer (PyG + upstream clone)
bash gpu/segmentanytree-sm120/build.sh

# 3) Import smoke on the 5090
docker run --rm --gpus all sat-sm120-test
```

## Next steps (same PR series)

- Pin SAT commit + record in `gpu/SegmentAnyTree.commit`
- Add `gpu/run_segmentanytree.py` headless driver (folder-in → apex CSV)
- Wire Docker backend in `scripts/model_runner.R`
- LAZ→PLY conversion for raw-with-ground UTM clips
