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
| Docker backend + LAZ↔PLY bridge | **Done** — `run_docker_arm` + `laz_to_ply`/`read_instances_ply` (#19) |
| Full `run_inference.sh` pipeline | **Not yet** — needs the headless driver (#M6) |

## Build

```sh
# 1) MinkowskiEngine base (~5 min compile)
bash gpu/build_me_sm120.sh

# 2) SegmentAnyTree layer (PyG + upstream clone)
bash gpu/segmentanytree-sm120/build.sh

# 3) Import smoke on the 5090
docker run --rm --gpus all sat-sm120-test
```

## Next steps (#M6)

Shared GPU-arm infra is now in place (#19):

- [x] Docker backend in `scripts/model_runner.R` (`run_docker_arm`)
- [x] LAZ→PLY conversion + instance-labeled output parsing (`scripts/io_bridge.R`:
  `laz_to_ply` / `read_instances_ply`)

Remaining, specific to the SAT arm:

- [ ] Pin SAT commit + record in `gpu/SegmentAnyTree.commit`
- [ ] Add `gpu/run_segmentanytree.py` headless driver (folder-in → labeled PLY/LAS)
- [ ] Emit the SAT-shaped input PLY via `laz_to_ply(props = …)` and feed one real
  Docker smoke
