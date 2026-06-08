# MinkowskiEngine on Blackwell (sm_120)

Proven recipe to build **MinkowskiEngine 0.5.4 for the RTX 5090** (compute
capability sm_120) on the torch 2.7.0 / CUDA 12.8 base. Both GPU arms that were
thought to be ME-blocked depend on this: SegmentAnyTree (#17) and ForestFormer3D
(#18 — voxelizer only; its conv backbone is spconv, see
[[../store/forestformer3d]] repro).

## Why this exists

Upstream `NVIDIA/MinkowskiEngine` ships no sm_120 build and is effectively
unmaintained ([ME#620](https://github.com/NVIDIA/MinkowskiEngine/issues/620)).
The fork `CiSong10/MinkowskiEngine@cuda12-installation` plus one libstdc++ patch
builds cleanly for arch `12.0`.

## Build and test

```sh
cd gpu/minkowski-sm120
docker build -t me-sm120-test .
docker run --rm --gpus all me-sm120-test
```

Expected: `SparseTensor ok` and `MinkowskiConvolution ok`, with the torch arch
list containing `sm_120` / `compute_120`.

## Recipe notes (the parts that bite)

- Base image `pytorch/pytorch:2.7.0-cuda12.8-cudnn9-devel`.
- `TORCH_CUDA_ARCH_LIST` / `CMAKE_CUDA_ARCHITECTURES` must include `12.0`.
- The `std::__to_address` sed fixes the `shared_ptr_base.h` ambiguity
  ([ME#596](https://github.com/NVIDIA/MinkowskiEngine/issues/596)); the gcc
  include path varies by base image, so the Dockerfile patches both 11 and 12.
- The CiSong10 branch **already disables NVTX** — do not re-apply an
  `NVTX_DISABLE` sed; it duplicates an existing `define_macros` and breaks the
  build.

Verified 2026-06-08 on driver 595.71.05 (image `me-sm120-test`, 17.5 GB):
ME 0.5.4, torch 2.7.0+cu128, MinkowskiConvolution forward on the 5090.
