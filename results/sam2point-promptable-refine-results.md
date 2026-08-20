# Seed→refine: CHM-VWF tops → SAM2Point promptable 3-D segmenter (#P3)

All existing deep arms (#M6 SegmentAnyTree, #M8 ForestFormer3D, #M7 TreeisoNet)
run in **automatic** mode. This builds the **seed→refine** stage: the strongest
overstory detector (CHM-VWF) drives a zero-shot **promptable** 3-D segmenter
(SAM2Point, Apache-2.0) for crown masks, decoupling "where is the tree" (the
CHM-VWF seed) from "what is its crown" (the deep refiner). SAM2Point voxelizes
the clip, slices the voxel grid into a "video", and runs SAM2's video predictor
from each 3-D apex prompt; each returned mask is a crown instance.

## Blackwell (sm_120) bring-up — verified

The headline engineering result: **SAM2Point runs on the RTX 5090 (sm_120)**.
`gpu/sam2point-sm120/` builds `FROM me-sm120-test:latest` — the
`pytorch/pytorch:2.7.0-cuda12.8-cudnn9-devel` base the SegmentAnyTree (#M6) and
ForestFormer3D (#M8) arms already use — so no PyTorch-on-Blackwell rebuild was
needed (SAM2Point wants torch ≥ 2.3.1; torch 2.7/cu128 satisfies that *and*
carries sm_120). The bundled demo and the custom LiDAR runner both complete on
the GPU with no CUDA errors. `run_sam2point_arm.py` (added here) reads a LAZ,
normalizes to a unit cube, feeds each CHM-VWF apex as a 3-D prompt, and writes
per-point `sam2point` instance labels; `scripts/detect_sam2point_sweep.R` drives
it per plot and scores the seeded masks next to the bare seeds.

## Cost (why this is a bounded proof-of-concept)

Each prompt is a **full 3-axis SAM2 video segmentation** (~10 s/prompt on the
5090). A native SOAP clip yields ~200 CHM-VWF local maxima, so one plot is
~10–35 min and a full 3-site density ladder is **many GPU-hours to days**. This
study therefore caps prompts (core+tol apexes, `MAXPROMPTS`) and runs a 2-plot
proof-of-concept; a full sweep + a SAM2Point-automatic / marker-watershed
comparison is future GPU-time, not a code gap (the driver accepts `PLOTS=ALL`).

## Proof-of-concept (SOAP, native, 2 plots, ≤40 seeds/plot)

| arm | n_ref | recall | precision | F1 |
|---|--:|--:|--:|--:|
| chm_vwf_seeds (the seeds) | 42 | 0.381 | 0.410 | 0.395 |
| **sam2point_seeded** | 42 | 0.262 | **0.733** | 0.386 |

Per plot: SOAP_031 — seeded precision 0.70 (vs seeds 0.36), recall 0.21 (vs
0.35); SOAP_021 — seeded precision 0.80 (vs 0.67), recall 0.50 (= seeds).

## Readings

- **The seed→refine mechanism works, and the refiner buys precision, not
  recall.** Given CHM-VWF tops, SAM2Point returns clean crown masks that lift
  precision sharply (0.73 vs the seeds' 0.41) — the masks that survive are
  usually real trees. But recall drops (0.26 vs 0.38): the refiner *refines*, it
  does not *detect more* than its seeds, and on these 2 plots 40 prompts collapse
  to ~10 distinct instances (conservative, overlapping masks merge under the
  first-claim assembly). F1 is a wash (0.386 vs 0.395).
- **This is the expected division of labour.** A promptable segmenter's job is
  "what is the crown of THIS seed", so it cannot exceed the seed set's recall —
  its value is mask quality/precision, which is exactly what it delivers. To gain
  recall the seed set itself must improve (e.g. fuse seeds from #P1) before
  refining.
- **Blackwell is no longer the blocker for promptable 3-D models in this repo.**
  The torch-2.7/cu128 `me-sm120` base carries SAM2 video inference unmodified, so
  future promptable arms can reuse it.

## Caveats

- **2-plot proof-of-concept, default parameters.** `voxel_size`, the
  height-ranked `MAXPROMPTS=40` cap, and the first-claim per-point assembly are
  untuned; the recall drop and prompt-collapse are partly artifacts of a coarse
  voxelization of sparse ALS crowns. A tuned, full-ladder run is future work.
- **LiDAR has no RGB**, so the SAM2Point voxel feature is the normalized height
  replicated to 3 channels (a pseudo-intensity) — a reasonable but unvalidated
  adaptation of an RGB-trained video model.
- **Scoring is apex-proximity here**; the per-point masks are persisted to
  `sam2point_instances/`, so the #V1 IoU/PQ harness and crown-diameter scoring can
  be run on them once a full sweep exists.
- **Cost, not correctness, bounds this arm** — the driver is `PLOTS=ALL`-ready.
