# SAM2Point seed→refine arm on Blackwell (sm_120) — #P3

The seed→refine stage: the strongest overstory detector (CHM-VWF) drives a
zero-shot promptable 3-D segmenter (SAM2Point, Apache-2.0) for crown masks,
decoupling "where is the tree" (#P1 seeds) from "what is its crown" (the deep
refiner). SAM2Point voxelizes the LiDAR, slices the voxel grid into a "video",
and runs SAM2's video predictor from a 3-D point prompt — so each CHM-VWF apex
becomes a prompt and each returned mask a crown instance.

## Why this base

Built `FROM me-sm120-test:latest` — the MinkowskiEngine sm_120 base
(`pytorch/pytorch:2.7.0-cuda12.8-cudnn9-devel`, `TORCH_CUDA_ARCH_LIST=12.0`) the
SegmentAnyTree (#M6) and ForestFormer3D (#M8) arms already use. SAM2Point needs
torch ≥ 2.3.1; torch 2.7.0/cu128 satisfies that **and** carries Blackwell
sm_120, so no PyTorch-on-Blackwell rebuild is required. The SAM2
connected-components CUDA extension (`sam2/csrc/connected_components.cu`) is left
uncompiled — it is optional and falls back to CPU mask postprocessing.

## Build + smoke

```sh
# prerequisite: the me-sm120 base (see ../build_me_sm120.sh)
docker build -t sam2point-sm120:test gpu/sam2point-sm120/
# self-contained demo on bundled Objaverse data (point prompt):
docker run --rm --gpus all sam2point-sm120:test bash run.sh
```

The image bakes in the `sam2_hiera_large.pt` checkpoint (~900 MB, official FAIR
mirror). `scripts/detect_sam2point_sweep.R` drives it per (plot, rung): it
voxelizes the frozen normalized clip, feeds `detect_lasr` CHM-VWF apexes as 3-D
point prompts, and takes each per-prompt mask as a crown instance (compared to
SAM2Point automatic mode, bare CHM-VWF tops, and the #32 marker-watershed seeded
from the same tops).

Cite: Guo et al., *SAM2Point: Segment Any 3D as Videos in Zero-shot and
Promptable Manners* (2024); Ravi et al., *SAM 2* (2024).
