# Matcher robustness: scaled tolerance + Hungarian re-score (#V4)

The whole benchmark grades detections with `greedy_match`: a greedy
nearest-distance 1:1 assignment within a **flat** 4 m radius, gated by a hard
height band `bz in [0.5*az, az + 8]`. Two known weaknesses motivate this arm:
(1) the apex↔stem-base horizontal offset grows with height, crown size, lean,
and slope (SOAP/TEAK are steep), so a flat 4 m is too tight for tall dominants
and too loose for dense understory; (2) greedy-by-distance is globally
suboptimal in dense clusters. This re-scores the CHM-VWF detector on the **same**
frozen clips with hardened matchers and reports whether any pooled rate moves.

The matchers are added to [`sweep_lib.R`](../scripts/sweep_lib.R) and unit-tested
(`tests/testthat/test-matcher-robustness.R`):

- `match_tol()` — per-stem `tol_i = max(base_tol, k·maxCrownDiameter/2, pos_unc)`
  (both fields already carried in `ground_truth_stems.csv`); the flat-4 m path is
  the back-compat default.
- `optimal_match()` — optimal 1:1 assignment (Hungarian, `clue::solve_LSAP`) on
  a finite-sentinel-gated cost matrix with dummy padding for non-square sets and
  unmatched stems; a drop-in for `greedy_match`.
- a soft scaled-3D cost `d = sqrt(dxy² + (λ·dz)²)` that replaces the hard height
  band, gated by a 3-D radius (`d ≤ tol`, so it still rejects height-impossible
  pairs), threaded through `optimal_match`/`score_plot`.
- `fp_structure()` — splits core false positives into *near a matched stem*
  (over-segmentation) vs *isolated* (real understory / field-map gap).

Regenerate:

```sh
export CLAUDE_JOB_DIR=$(pwd)/work
# CORES=1: detect_lasr uses lasR exec, which can drop dense-native cells under fork
Rscript scripts/matcher_robustness.R SITES=SOAP,SJER,TEAK CORES=1
# -> work/neon/<SITE>/matcher_robustness.csv (one row per plot x rung x config)
```

## What this is

For each plot the CHM-VWF detector is regenerated **deterministically** from the
cached frozen normalized clip (`detect_lasr` at the canonical density-derived
`chm_res`, `vwf_a = 0.10`), so every matcher scores identical apexes. Native
density, three sites, 699 pooled field stems. Pooled with the canonical `pool()`
(sum counts, never average rates). The baseline row reproduces the benchmark's
CHM-VWF native recall/F1 (0.40/0.37 pooled), confirming the harness is unchanged.

## Generated tables

### Matcher configs, all sites combined (pooled by SUM)

`hRMSE` = TP-weighted pooled apex-height RMSE over matched pairs — a matcher that
"recovers" matches by admitting height-implausible pairs inflates it.

| config | recall | precision | F1 | ΔF1 | hRMSE (m) | rec_understory |
|---|--:|--:|--:|--:|--:|--:|
| baseline | 0.401 | 0.336 | 0.365 | +0.000 | 3.17 | 0.200 |
| scaled | 0.411 | 0.343 | 0.374 | +0.009 | 3.20 | 0.200 |
| optimal | 0.413 | 0.342 | 0.374 | +0.009 | 3.17 | 0.200 |
| **optimal_scaled** | **0.425** | **0.351** | **0.385** | **+0.019** | 3.20 | 0.200 |
| soft3d | 0.382 | 0.318 | 0.347 | −0.018 | **2.35** | 0.162 |

### Per site: ΔF1 vs baseline (optimal_scaled and soft3d)

| site | n_ref | base F1 | optimal_scaled ΔF1 | soft3d ΔF1 | soft3d hRMSE | base hRMSE |
|---|--:|--:|--:|--:|--:|--:|
| SOAP | 232 | 0.382 | +0.028 | −0.036 | 2.34 | 3.30 |
| SJER | 71 | 0.319 | +0.019 | +0.009 | 1.91 | 2.04 |
| TEAK | 396 | 0.367 | +0.013 | −0.012 | 2.47 | 3.29 |

### tol_xy × tol_z_up sensitivity, combined (greedy F1)

| tol_xy | tz=5 | tz=8 | tz=12 |
|--:|--:|--:|--:|
| 2 | 0.214 | 0.227 | 0.247 |
| 3 | 0.288 | 0.308 | 0.330 |
| 4 | 0.337 | 0.365 | 0.384 |
| 5 | 0.373 | 0.399 | 0.418 |

### False-positive error structure (baseline, core FPs, combined)

near-a-matched-stem (over-segmentation) = **29 (5.7 %)**; isolated (real
understory / field-map gap) = **478 (94.3 %)**.

| class (nearest stem) | near (over-seg) | isolated |
|---|--:|--:|
| dominant | 12 | 171 |
| codominant | 16 | 214 |
| intermediate | 1 | 30 |
| suppressed | 0 | 0 |

## Readings

- **Modest but consistent gains — close to the "greedy is adequate" the issue
  anticipated.** The best principled variant is **optimal_scaled** (Hungarian
  optimal assignment + per-stem size/uncertainty-scaled tolerance): +0.019 F1
  pooled, and it is the only variant positive at every site (+0.028 SOAP, +0.019
  SJER, +0.013 TEAK) with no height-RMSE cost (3.20 vs 3.17 m). Scaled tolerance
  and optimal assignment each contribute ~+0.009 F1 and compose to +0.019. None
  of the increments is dramatic: the flat-4 m greedy baseline is already close to
  adequate, and the hardened matcher is worth adopting for the steady gain and
  the dropped magic numbers, not for a step change.
- **The soft 3-D cost buys height fidelity, not recall.** Replacing the hard
  `[0.5·az, az+8]` band with `sqrt(dxy² + (λ·dz)²)` gated by a 3-D radius
  (`d3 ≤ tol`, which still rejects height-impossible pairs — capping |dz| at
  tol/λ) slightly *lowers* recall/F1 (−0.018) but gives the **best matched-pair
  height RMSE of any variant, 2.35 m vs the baseline's 3.17 m**. It declines the
  marginal, height-implausible matches the hard band let through, trading a hair
  of recall for cleaner geometry — useful where matched-apex height quality
  matters (e.g. feeding #V3 error bars), not as a recall booster.
  (An earlier revision of this arm reported soft3d as a large F1 winner; that was
  a bug — the soft path had dropped the height gate entirely and was scoring
  height-impossible matches as true positives, which a #V4 self-review caught.
  The 3-D-radius gate is the fix.)
- **The flat-4 m / gate-8 baseline sits on a rising slope.** F1 climbs
  monotonically across the whole tol_xy {2→5} × tol_z_up {5→12} grid, so the
  benchmark's defaults are on the tight side. The principled fix is the per-stem
  scaled tolerance (loosen for big crowns, not globally, and capped at 12 m so a
  corrupt field record — SOAP carries a 344 m `maxCrownDiameter` — cannot blow up
  the radius), since a globally larger flat radius would eventually inflate true
  positives with spurious matches.
- **Low precision is a ground-truth coverage gap, not over-segmentation.**
  94.3 % of core false positives are *isolated* (no matched stem within 4 m);
  only 5.7 % sit beside a matched tree. The benchmark's modest precision is
  dominated by detections of real-but-unmapped trees (regeneration, unmeasured
  neighbours), not by the detector splitting one crown into many. This is the
  signal #P2 (router) and #P1 (fusion) need: isolated detections should be
  treated as probably-real, not suppressed as commission.

## Caveats

- **Scoring is apex-proximity**, which #V1 shows overstates instance quality for
  the mask-capable arms; for SegmentAnyTree/ForestFormer3D prefer the point-set
  IoU/PQ scorer (`results/instance-iou-pq-results.md`, #V1).
  This arm hardens the proximity matcher the CHM/apex detectors still rely on,
  and the same `match_tol`/height-gate logic feeds the #P1 fusion dedup.
- **The per-class FP attribution is by nearest stem**, so it locates each false
  positive in a crown-class neighbourhood rather than labelling the (unmapped)
  detection itself; the near/isolated split is the primary signal.
- **Native density only**; the driver accepts a `RUNGS=` list and extends to the
  sparse ladder once those detections are wanted. Soft-3D `λ` defaults to 0.5
  (tunable via `LAMBDA=`), and the scaled tolerance is capped at 12 m
  (`TOL_CAP=`); the sensitivity grid above uses the greedy hard-gate path.
