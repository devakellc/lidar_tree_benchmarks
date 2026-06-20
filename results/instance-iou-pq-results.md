# Point-set IoU / Coverage / Panoptic-Quality scorer (#V1)

A mask-aware evaluation of the deep instance-segmentation arms on the NEON
frozen clips, added **alongside** (never replacing) the `greedy_match`
apex-distance scoring that grades the rest of the benchmark. The apex scorer
reduces every model to one point per tree, so it is blind to point-level
instance IoU and crown-shape quality — the limitation
[`model-benchmark-results.md`](model-benchmark-results.md) already flags in its
Caveats, and the exact failure mode the meta-pipeline's fusion arm (#P1) is
meant to fix. FGI-EMIT and FOR-instanceV2 grade with point-set IoU≥0.5 +
Coverage + Panoptic Quality; adopting that suite here makes the repo's arms
comparable to those external leaderboards and gives fusion a mask-aware
objective. This is the evaluation backbone the consensus arm (#P1) and the
seed→refine arm (#P3) reference.

Regenerate:

```sh
export CLAUDE_JOB_DIR=$(pwd)/work
Rscript scripts/score_instances_iou.R SITES=SOAP,SJER,TEAK CORES=4
# -> work/neon/<SITE>/instance_iou_pq.csv (one row per plot x rung x model)
```

## What this is

For each plot the scorer builds **one fixed evaluation cloud** — the cached
`frozen_clip` normalized canopy points inside the plot core — and scores two
partitions of it:

- **Reference** instances: each evaluation point is assigned to the nearest
  field stem within that stem's crown radius (`maxCrownDiameter / 2`). This is
  an explicit **Voronoi-on-stems proxy** for true field crown masks — the best
  reference obtainable from stem points plus a measured crown width, not a
  hand-delineated truth. It is reported as a proxy throughout.
- **Predicted** instances: the per-point labels the deep arms already persist,
  projected onto the evaluation cloud (see below). No container is ever re-run.

It then computes, with the pure helpers in
[`model_bench_lib.R`](../scripts/model_bench_lib.R) (all unit-tested with no
NEON data, in `tests/testthat/test-instance-iou-pq.R`):

- **P / R / F1 at IoU≥0.5** — `iou_match()` is the IoU sibling of
  `greedy_match`: greedy best-IoU-first, 1:1, accepting only pairs that clear
  the 0.5 gate. True positives are the matched pairs.
- **Coverage** — `coverage_table()`: the mean over reference instances of the
  maximum IoU against any prediction (the **ungated** companion to recall@0.5).
- **Panoptic Quality** — `panoptic_quality()`: PQ = SQ · RQ, where SQ
  (segmentation quality) is the mean IoU over matched pairs and RQ (recognition
  quality) is TP / (TP + ½FP + ½FN) — exactly the F1 over matched pairs.

Everything is stratified by `crown_class` and **pooled by SUMMING the panoptic
accumulators** (sum matched-IoU, sum TP/FP/FN, sum per-class counts), then
recomputing the rates from the sums — mirroring the `sum(TP)/sum(n_ref)` rule
in `pool()`. A small plot never dominates a site rate.

### Arms scored

- **SegmentAnyTree** (#M6) — `segmentanytree_instances/<plot>_<rung>.laz`, the
  `PredInstance` extra dim (0 = non-tree → dropped).
- **ForestFormer3D** (#M8) — `forestformer3d_instances/<plot>_<rung>.laz`, with
  `UserData` (block) + `PointSourceID` (per-cylinder instance) run through
  `dedup_blocks()` to a globally consistent label per point (the block id is
  required, so `read_instance_points_laz` alone is insufficient).
- **TreeisoNet** — deferred until `detect_treeisonet_crowns.R` persists its
  per-point `treeOff` labels (today it writes only
  `treeisonet_crown_metrics.csv`, and `detect_treeisonet_sweep.R` is apex-only).

### The label projection

Standard instance-IoU is computed on a single fixed evaluation cloud. The
persisted model clouds are **not** the frozen normalized cloud point-for-point —
different Z datum (absolute UTM vs AGL), container voxel downsampling, and FF3D
cylinder overlap all differ — so each model's labels are projected onto the
normalized cloud by **nearest-neighbour in XY** within `XFER_TOL` (0.5 m),
implemented as a dependency-free grid hash (`transfer_labels()`; no RANN/FNN is
installed). On SegmentAnyTree, whose raw-with-ground clip and the normalized
clip derive from the *same* decimated set, the projection is effectively exact
(spot-checked transfer hit rate ≈ 0.85 of canopy points; the remainder are
points the model itself left unlabelled). The substrate is restricted to canopy
points (Z ≥ 2 m AGL, the detector `min_height`) inside the plot core
(`plot_half` by plot type), mirroring `score_plot`'s core precision denominator
so buffer trees never count as false commissions.

## Generated tables

Native density, all three sites. Reference is the Voronoi-on-stems proxy; IoU
gate 0.5. 699 mapped stems fall in the 46 plot cores; **n_ref = 686** is the
recall denominator — the stems that captured ≥1 canopy substrate point, since a
stem with no point-mask cannot be detected by a point-set metric (13 stems,
~2 %, captured none and are excluded). 7.3 % of stems lacked a measured crown
width and used the 2 m fallback radius.

### Per site × model (pooled by SUM)

| site | model | cells | n_ref | P | R | F1 | Cov | SQ | RQ | PQ |
|---|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| SOAP | segmentanytree | 18 | 231 | 0.104 | 0.216 | 0.140 | 0.342 | 0.718 | 0.140 | 0.101 |
| SOAP | forestformer3d | 18 | 231 | 0.028 | 0.061 | 0.038 | 0.171 | 0.638 | 0.038 | 0.024 |
| SJER | segmentanytree | 8 | 69 | 0.052 | 0.246 | 0.085 | 0.357 | 0.718 | 0.085 | 0.061 |
| SJER | forestformer3d | 8 | 69 | 0.025 | 0.101 | 0.040 | 0.191 | 0.677 | 0.040 | 0.027 |
| TEAK | segmentanytree | 20 | 386 | 0.096 | 0.202 | 0.130 | 0.313 | 0.684 | 0.130 | 0.089 |
| TEAK | forestformer3d | 20 | 386 | 0.054 | 0.054 | 0.054 | 0.115 | 0.668 | 0.054 | 0.036 |

### All sites combined (pooled by SUM)

| model | cells | n_ref | P | R | F1 | Cov | SQ | RQ | PQ |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| segmentanytree | 46 | 686 | 0.089 | 0.211 | 0.126 | 0.327 | 0.700 | 0.126 | 0.088 |
| forestformer3d | 46 | 686 | 0.036 | 0.061 | 0.045 | 0.141 | 0.660 | 0.045 | 0.030 |

### Per crown class, all sites combined

recall@IoU0.5 / Coverage / SQ (SQ over matched pairs only; "--" = no match).

| model | class | n_ref | recall@0.5 | Coverage | SQ |
|---|---|--:|--:|--:|--:|
| segmentanytree | dominant | 178 | 0.326 | 0.408 | 0.738 |
| segmentanytree | codominant | 359 | 0.209 | 0.324 | 0.685 |
| segmentanytree | intermediate | 95 | 0.095 | 0.240 | 0.632 |
| segmentanytree | suppressed | 9 | 0.222 | 0.317 | 0.550 |
| segmentanytree | understory | 104 | 0.106 | 0.247 | 0.617 |
| forestformer3d | dominant | 178 | 0.118 | 0.190 | 0.680 |
| forestformer3d | codominant | 359 | 0.056 | 0.142 | 0.639 |
| forestformer3d | intermediate | 95 | 0.011 | 0.076 | 0.646 |
| forestformer3d | suppressed | 9 | 0.000 | 0.090 | -- |
| forestformer3d | understory | 104 | 0.010 | 0.077 | 0.646 |

## Readings

- **Apex matching overstates instance quality.** SegmentAnyTree's pooled
  recall falls from the ~0.5–0.7 it scores under apex-distance matching in
  [`model-benchmark-results.md`](model-benchmark-results.md) to **0.211** under
  point-set IoU≥0.5. Finding a tree's apex is not the same as recovering its
  crown: roughly two of every three apex matches do not carry a mask that
  overlaps the reference crown by half. This is the quantitative case for the
  meta-pipeline's mask-aware fusion.
- **The ranking survives, the magnitudes do not.** SegmentAnyTree beats
  ForestFormer3D on every metric, site, and crown class — the same ordering the
  apex scorer gives — but both sit far below their apex scores, so the suite
  adds discrimination the apex metric hides rather than re-deriving it.
- **Coverage > recall@0.5, and SQ is healthy.** Mean max-IoU (0.327 for
  SegmentAnyTree) exceeds recall@0.5 (0.211) while matched-pair SQ is ~0.70:
  predictions usually overlap the right crown, and when a match clears the gate
  its mask IoU is good — the 0.5 threshold, not gross mislabelling, is the
  binding constraint. Coverage is the gentler signal fusion can optimise first.
- **The canopy-occlusion gradient is explicit.** recall@0.5 and Coverage fall
  monotonically dominant → codominant → intermediate/understory for both arms,
  the point-mask analogue of the per-class recall gradient the apex scorer
  reports.

## Caveats

- **The reference is a proxy, not truth.** Voronoi-on-stems assigns each canopy
  point to the nearest mapped stem within `maxCrownDiameter / 2`, so reference
  crowns are circular and non-overlapping. Absolute IoU magnitudes are
  proxy-relative; the cross-arm ranking and the apex-vs-IoU gap are the robust
  results, not the exact PQ values. A future hand-delineated crown set would
  raise every number but is not required to read the comparison.
- **Precision and PQ are conservative lower bounds.** The reference covers only
  field-mapped stems, so any prediction over an unmapped tree — regeneration,
  an unmeasured neighbour — is counted as a false commission. This is the same
  reason apex precision in the benchmark is a lower bound; the
  recall / Coverage / SQ side is the more interpretable signal here.
- **Zero-shot on sparse ALS.** As elsewhere, these are zero-shot transfers of
  models trained on denser ULS/UAS/TLS to NEON discrete-return airborne LiDAR;
  the numbers measure that transfer, not ceiling performance.
- **Native density only.** The persisted instance clouds exist at native
  density for both arms (SegmentAnyTree also has the sparse rungs); the scorer
  accepts a `RUNGS=` list and extends to the sparse ladder for free once those
  clouds are present.
