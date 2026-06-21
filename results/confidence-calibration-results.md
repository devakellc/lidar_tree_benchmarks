# Per-detection confidence calibration (#P4)

Weighted consensus in the fusion arm (#P1) needs each detection to carry a
**trustworthy probability** so a high-confidence SegmentAnyTree mask can outvote
a low-confidence apex — but raw scores are not comparable across arms (a CNN's
0.6 ≠ a watershed prominence of 0.6), and the repo keeps every detection
unweighted. This arm makes scores cross-arm-comparable: it labels each arm's
detections TP/FP against field stems, builds a per-arm reliability diagram +
Expected Calibration Error (ECE), fits a post-hoc **isotonic** calibrator
(`stats::isoreg`) mapping raw score → empirical precision, and measures the
ensemble payoff — precision at fixed recall when the pooled multi-arm detections
are ranked by calibrated vs raw score.

Regenerate:

```sh
export CLAUDE_JOB_DIR=$(pwd)/work
Rscript scripts/calibrate_confidence.R SITES=SOAP,SJER,TEAK
# -> work/neon/<SITE>/confidence_calibration.csv (one row per labelled detection)
#    work/neon/<SITE>/confidence_lookup.csv      (per-arm isotonic knots for #P1)
```

## What this is

Per-arm raw confidence, materialized on the same native frozen cells as #P1 and
labelled TP/FP with `greedy_match` restricted to the plot core (mirroring
`score_plot`'s precision denominator):

- **forestformer3d** — the **native** per-instance mask score: mean `ff3d_score`
  over the instance's points. `ff3d_arm.py` writes this extra dim and the R
  loaders ignored it; this exposes it (the issue's explicit ask).
- **segmentanytree** — crown point count (`run_segmentanytree.py` writes no score
  field): a size proxy.
- **chm_vwf / multichm / li2012** — apex CHM height (AGL): the classical-arm
  proxy (taller apices are more often real dominant trees).

Each arm's raw score is min-max normalized per arm to a predicted probability in
[0,1] (the uncalibrated *p*); the helpers `reliability_table`,
`expected_calibration_error`, `isotonic_calibrate`, and `precision_at_recall` (in
[`model_bench_lib.R`](../scripts/model_bench_lib.R), unit-tested in
`tests/testthat/test-confidence-calibration.R`) do the rest. **ECE_cal is 5-fold
cross-validated** (fit on 4 folds, scored on the held-out fold) — isotonic is
in-sample-perfect, so only a held-out number is honest.

## Generated tables

Native density, three sites, 5340 labelled detections.

### Per-arm calibration

| arm | n | base precision | ECE_raw | ECE_cal (CV) | ΔECE |
|---|--:|--:|--:|--:|--:|
| chm_vwf | 763 | 0.336 | 0.162 | 0.017 | −0.145 |
| multichm | 1023 | 0.333 | 0.168 | 0.033 | −0.135 |
| li2012 | 1189 | 0.266 | 0.140 | 0.024 | −0.116 |
| segmentanytree | 1245 | 0.334 | 0.256 | 0.026 | −0.230 |
| forestformer3d | 1120 | 0.224 | 0.126 | 0.004 | −0.121 |

### Ensemble precision @ fixed recall (pooled multi-arm, held-out)

Ranking the pooled 5-arm detection set by calibrated vs raw score:

| recall | P_raw | P_calib | gain |
|--:|--:|--:|--:|
| 0.50 | 0.317 | 0.366 | +0.049 |
| 0.60 | 0.322 | 0.351 | +0.029 |
| 0.70 | 0.318 | 0.351 | +0.033 |
| 0.80 | 0.319 | 0.347 | +0.028 |
| 0.90 | 0.313 | 0.326 | +0.014 |

### Per-site ensemble gain @ recall 0.70 (held-out)

| site | n_det | P_raw | P_calib | gain |
|---|--:|--:|--:|--:|
| SOAP | 2037 | 0.347 | 0.383 | +0.036 |
| SJER | 1203 | 0.216 | 0.222 | +0.006 |
| TEAK | 2100 | 0.378 | 0.401 | +0.023 |

### SegmentAnyTree raw reliability (the worst-ECE arm, 5 bins)

| confidence bin | n | mean score | observed precision |
|--:|--:|--:|--:|
| 1 (0.0–0.2) | 1130 | 0.060 | 0.327 |
| 2 (0.2–0.4) | 88 | 0.271 | 0.375 |
| 3 (0.4–0.6) | 17 | 0.505 | 0.588 |
| 4 (0.6–0.8) | 5 | 0.689 | 0.600 |
| 5 (0.8–1.0) | 5 | 0.940 | 0.000 |

## Readings

- **Raw scores are meaningfully miscalibrated; calibration fixes it
  out-of-sample.** Per-arm ECE falls from 0.13–0.26 raw to **0.004–0.033
  cross-validated** — the calibrator generalizes (a held-out ECE near zero, not
  the in-sample 0.000 isotonic trivially achieves). The worst raw arm is
  SegmentAnyTree's crown-point-count proxy (ECE 0.256); the best is
  ForestFormer3D's **native** `ff3d_score` (ECE 0.126 → 0.004), which makes
  sense — a learned mask confidence is a better raw signal than a geometric proxy,
  and it calibrates almost perfectly.
- **"Bigger mask = more confident" is false at the extreme.** SegmentAnyTree's
  reliability diagram is non-monotone at the top: its five largest-point-count
  masks (mean score 0.94) have **observed precision 0.000** — the biggest masks
  are usually over-grown / merged over-segmentation, not the most certain trees.
  Isotonic calibration caps that high end to empirical precision, which is exactly
  the correction a weighted vote needs.
- **Calibration buys real ensemble precision.** Ranking the pooled five-arm
  detection set by calibrated (comparable) score beats ranking by raw
  (incomparable) score at every recall, **+0.014 to +0.049 precision**, held-out
  and consistent across sites (+0.036 SOAP, +0.023 TEAK, +0.006 on the
  smaller/easier SJER set). This is the mechanism #P1 weighted fusion and the #P2
  router consume: a high-confidence detection from any arm now means the same
  thing.
- **The shipped lookup is fit on all data.** `confidence_lookup.csv` carries each
  arm's isotonic knots fit on the full per-arm sample (max data for downstream
  use); the ECE/precision numbers above are the honest held-out estimates of how
  well that mapping generalizes.

## Caveats

- **Labels are apex-proximity TP/FP** (`greedy_match`, 4 m + height gate), so a
  detection of a real-but-unmapped tree is labelled FP — the same field-map
  coverage gap #V4 quantified (≈94 % of FPs are isolated real trees). Calibration
  targets *this* precision; absolute precision is a lower bound, but the
  raw-vs-calibrated comparison is internally consistent (same labels both ways).
- **Proxies, not native confidences, for four of five arms.** Only ForestFormer3D
  exposes a learned score; SegmentAnyTree's `run_segmentanytree.py` writes none
  (point count is a stand-in until a mask score is persisted on the next run), and
  the CHM/classical arms use apex height. Calibration makes them *comparable*, not
  equally *informative* — a better raw signal (e.g. local-maximum prominence) would
  raise the ceiling.
- **Native density only**; the calibrators are fit at native and would need
  re-fitting per rung for a density-stratified weighting.
