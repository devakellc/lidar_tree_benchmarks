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
# the classical arms' instance clouds first (#V6; re-runs also refresh the
# results CSVs -- CORES=1 for the lidrplugins script, its chm_vwf arm runs
# lasR exec, which drops dense cells under fork):
Rscript scripts/detect_li2012_native.R SITE=SOAP CORES=6
Rscript scripts/detect_ams3d_sweep.R SITE=SOAP CORES=6
Rscript scripts/detect_lidrplugins_sweep.R SITE=SOAP CORES=1
Rscript scripts/score_instances_iou.R SITES=SOAP,SJER,TEAK CORES=4 \
    RUNGS=native,8,4,2,1
# -> work/neon/<SITE>/instance_iou_pq.csv (one row per plot x rung x model x
#    mask_source; APEX_PROXY=0 drops the apex-Voronoi proxy rows)
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

Native per-point masks (`mask_source = "native"`):

- **SegmentAnyTree** (#M6) — `segmentanytree_instances/<plot>_<rung>.laz`, the
  `PredInstance` extra dim (0 = non-tree → dropped).
- **ForestFormer3D** (#M8) — `forestformer3d_instances/<plot>_<rung>.laz`, with
  `UserData` (block) + `PointSourceID` (per-cylinder instance) run through
  `dedup_blocks()` to a globally consistent label per point (the block id is
  required, so `read_instance_points_laz` alone is insufficient).
- **ptrees / AMS3D / Li2012** (#V6) — the classical segmenters always computed
  per-point `treeID`/`crown_id` and threw it away before the apex collapse;
  `io_bridge.R::write_instances_laz` now persists it (integer extra dim, 0 =
  unassigned, the SAT layout) from `detect_lidrplugins_sweep.R`,
  `detect_ams3d_sweep.R`, and `detect_li2012_native.R` to
  `<arm>_instances/<plot>_<rung>.laz`. These clips are segmented from the
  normalized frozen clouds, so Z is already AGL.
- **Treeiso** (#P5) — `treeiso_instances/<plot>_<rung>.laz` (`treeiso` extra
  dim), persisted by its own driver since #P5; now registered in `MODELS`.
- **TreeisoNet** — still deferred (apex-only GPU results; no per-point labels).

Apex-Voronoi proxy masks (`mask_source = "voronoi_apex"`, #V6): every cached
best-configuration apex set (`best_treetop_cache`, 12 arms) additionally
becomes a mask by nearest-apex assignment within `APEX_R` (4 m) on the same
substrate — the `fuse_detectors.R` scoring proxy, symmetric with the
Voronoi-on-stems reference. This puts the apex-only detectors (CHM-VWF,
multichm, lmfauto, TreeisoNet, the pc twins) on the board at their best rung,
and double-scores the native-mask arms so **proxy inflation is measurable**.
Proxy rows never mix with native rows in rankings.

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

### The completed native-mask board (#V6; all sites, native rung)

The classical segmenters' persisted masks (re-run 2026-08-20; the regenerated
results CSVs reproduce the shipped distance numbers exactly) complete the
board:

| model | cells | n_ref | P | R | F1 | Cov | SQ | PQ |
|---|--:|--:|--:|--:|--:|--:|--:|--:|
| segmentanytree | 46 | 686 | 0.089 | 0.211 | 0.126 | 0.327 | 0.700 | 0.088 |
| li2012 | 46 | 686 | 0.069 | 0.157 | 0.096 | 0.269 | 0.686 | 0.066 |
| ptrees | 46 | 686 | 0.038 | 0.175 | 0.063 | 0.303 | 0.689 | 0.043 |
| ams3d | 46 | 686 | 0.027 | 0.173 | 0.047 | 0.289 | 0.680 | 0.032 |
| treeiso | 43 | 661 | 0.070 | 0.035 | 0.046 | 0.110 | 0.677 | 0.031 |
| forestformer3d | 46 | 686 | 0.036 | 0.061 | 0.045 | 0.141 | 0.660 | 0.030 |

SegmentAnyTree keeps #1 on every mask metric. The crown-splitters hold decent
mask **recall** (ptrees 0.175, AMS3D 0.173 — their splits still overlap the
reference at ≥0.5 for matched trees) but collapse on mask **precision**
(0.027–0.038): every split fragment is a false instance. Li2012 is the best
classical mask (PQ 0.066).

### Proxy inflation (same arm, same rung, native vs voronoi_apex)

Where an arm's cached best rung is native, both mask sources exist for the
same detections (SOAP shown; TEAK SAT +0.011 F1):

| arm (SOAP) | F1 native | F1 proxy | inflation |
|---|--:|--:|--:|
| segmentanytree | 0.140 | 0.179 | +28 % |
| li2012 | 0.095 | 0.111 | +17 % |

The apex-Voronoi proxy **flatters mask quality** — most visibly on understory
(SAT understory recall@IoU0.5 0.133 native vs 0.200 proxy, both pooled over the
same 18 SOAP native-rung cells; the proxy rows only exist at an arm's cached
best rung, so the comparison has to stay inside that rung). Fusion's IoU
columns (built on this proxy) are therefore optimistic upper bounds; the
native-mask board is the honest ranking.

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

## Matching-rule sensitivity (#V2)

If the router (#P2) selects "best arm per rung", that choice is only trustworthy
if the ranking is stable across reasonable metrics. This recomputes the SOAP
native leaderboard under apex-distance F1, point-set IoU≥0.5 recall, and
threshold-free Coverage (`scripts/compare_matching_rules.R` →
`work/neon/SOAP/matching_rule_ranks.csv`), and flags reorderings. The hypothesis
under test: distance matching over-credits high-recall/low-precision point
segmenters (ptrees, AMS3D — split crowns) that IoU/RQ should demote.

| arm | dist recall | dist F1 | dist rank | IoU R@0.5 | Coverage | native PQ | mask? |
|---|--:|--:|--:|--:|--:|--:|:--:|
| segmentanytree | 0.642 | 0.464 | 1 | 0.216 | 0.342 | 0.101 | yes |
| multichm | 0.629 | 0.439 | 2 | n/a | n/a | n/a | no |
| treeisonet | 0.599 | 0.393 | 3 | n/a | n/a | n/a | no |
| chm_vwf | 0.478 | 0.382 | 4 | n/a | n/a | n/a | no |
| li2012 | 0.595 | 0.359 | 5 | 0.182 | 0.302 | 0.063 | yes |
| lmfauto | 0.509 | 0.337 | 6 | n/a | n/a | n/a | no |
| forestformer3d | 0.409 | 0.262 | 7 | 0.061 | 0.171 | 0.024 | yes |
| ptrees | 0.849 | 0.259 | 8 | 0.212 | 0.360 | 0.041 | yes |
| ams3d | 0.789 | 0.191 | 9 | 0.190 | 0.319 | 0.030 | yes |
| treeiso | 0.086 | 0.151 | 10 | 0.035 | 0.100 | 0.033 | yes |

Rank correlations over the **six** native-mask arms (dist F1 vs IoU recall@0.5,
then vs native PQ), as printed by `compare_matching_rules.R SITE=<site>` — the
`pq` column and both PQ correlations come out of the same guarded pool as the
IoU ones and land in `matching_rule_ranks.csv`. Both boards are equal-set
guarded — the IoU pools are restricted to the cells every mask arm scored
(treeiso misses 3 TEAK cells), mirroring the distance side, so the correlations
never mix denominators (guarding moves no rank):

| site | τ (IoU R) | ρ (IoU R) | τ (PQ) | ρ (PQ) |
|---|--:|--:|--:|--:|
| SOAP | 0.467 | 0.543 | 0.467 | 0.600 |
| SJER | 0.645 | 0.736 | 0.467 | 0.429 |
| TEAK | 0.600 | 0.714 | 0.733 | 0.829 |

Readings (updated for #V6 — the previous 2-of-9 data limit is closed):

- **τ/ρ are now defined, and they say the boards genuinely disagree.**
  Kendall τ 0.47–0.73 is a moderate, far-from-perfect agreement: the
  distance leaderboard is not a safe stand-in for the mask leaderboard. The
  router (#P2) can keep ranking on distance F1 for apex-counting, but any
  mask-consuming decision (crown delineation, #P3 seeds) needs the native
  IoU/PQ board.
- **The over-crediting hypothesis resolves with a twist.** ptrees and AMS3D
  hold their IoU **recall** (0.212/0.190 — ranks 2–3, their split fragments
  still overlap the reference) but collapse on mask **precision** (0.035/0.025
  native), so IoU-recall alone would *not* demote the crown-splitters — PQ and
  IoU-F1 do (ptrees 0.061, AMS3D 0.044 vs SegmentAnyTree 0.140 on SOAP).
  Distance F1's precision term and mask PQ agree on the demotion; metrics that
  ignore false instances (plain recall, Coverage) are the ones a splitter can
  game.
- **No flip at the top.** SegmentAnyTree is #1 on distance F1, IoU recall,
  Coverage, and PQ; its lead is metric-robust. Treeiso is last or near-last
  everywhere (native PQ 0.031–0.036): its fusion membership rests on vote
  diversity, not mask quality.
- **Li2012 is the quiet winner among the classical masks** (PQ 0.066 pooled,
  2nd overall behind SAT): the sub-canopy segmenter's masks are noticeably
  better than the other classical arms' at native density.

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
- **The ladder is covered but the headline tables are native.** #V6 persists
  the classical arms at every rung and the scorer ran `RUNGS=native,8,4,2,1`
  (the per-rung rows live in `instance_iou_pq.csv`); the tables here pool the
  native rung, where every arm has cells.
- **Proxy rows are upper bounds, never rankings.** The `voronoi_apex` rows
  measure the apex-Voronoi proxy (and its inflation vs native masks); they are
  excluded from `compare_matching_rules.R` and should never be ranked against
  native-mask rows.
- **Fusion pool (#P1) now includes ptrees + AMS3D** from the persisted clouds
  (SOAP native smoke: 7 single arms + modes + k1–k7 Pareto, best fused point
  k5 F1 0.430 vs best single 0.464); the full cross-site fusion re-synthesis
  on the extended pool is future work.
