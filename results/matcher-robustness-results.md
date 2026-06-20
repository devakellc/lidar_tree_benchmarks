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
  band, threaded through `optimal_match`/`score_plot`.
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

| config | recall | precision | F1 | ΔRecall | ΔF1 | rec_understory |
|---|--:|--:|--:|--:|--:|--:|
| baseline | 0.401 | 0.336 | 0.365 | +0.000 | +0.000 | 0.200 |
| scaled | 0.412 | 0.343 | 0.375 | +0.011 | +0.009 | 0.200 |
| optimal | 0.413 | 0.342 | 0.374 | +0.013 | +0.009 | 0.200 |
| optimal_scaled | 0.426 | 0.351 | 0.385 | +0.026 | +0.020 | 0.200 |
| soft3d | 0.485 | 0.384 | 0.429 | +0.084 | +0.063 | 0.267 |

### Per site: baseline vs soft3d (the best variant)

| site | n_ref | base recall | base F1 | soft3d recall | soft3d F1 | ΔF1 |
|---|--:|--:|--:|--:|--:|--:|
| SOAP | 232 | 0.478 | 0.382 | 0.552 | 0.430 | +0.048 |
| SJER | 71 | 0.493 | 0.319 | 0.507 | 0.329 | +0.009 |
| TEAK | 396 | 0.338 | 0.367 | 0.442 | 0.462 | +0.095 |

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

- **Not a null result: the matcher is conservative and worth hardening.** Every
  variant lifts pooled F1, and the soft 3-D cost — which drops the magic
  `[0.5·az, az+8]` band for `sqrt(dxy² + (λ·dz)²)` — is the clear winner: +0.063
  F1 and +0.084 recall pooled, recovering matches the hard band rejected.
  Per-stem scaled tolerance and optimal assignment each add smaller, independent
  increments (+0.009 to +0.020 F1), so they compose.
- **The gain tracks terrain, as predicted.** soft3d helps most on steep **TEAK**
  (+0.095 F1) and least on flat, open **SJER** (+0.009). The base→apex offset
  that a flat radius + hard height band mishandle is largest exactly where slope
  and tall dominants are, which is where the hardened matcher pays off.
- **The flat-4 m / gate-8 baseline sits on a rising slope.** F1 climbs
  monotonically across the whole tol_xy {2→5} × tol_z_up {5→12} grid, so the
  benchmark's defaults are on the tight side. The principled fix is the per-stem
  scaled tolerance (loosen for big crowns, not globally), since a globally larger
  flat radius would eventually inflate true positives with spurious matches.
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
  (tunable via `LAMBDA=`); the sensitivity above fixes the hard-gate path.
