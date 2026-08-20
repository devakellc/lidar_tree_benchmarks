# Coverage-gap crediting: de-biased precision/F1 (#V5)

## V4 (`matcher_robustness.R`) found that ~94% of the benchmark's core false

positives are **isolated** — not over-segmentation of a matched crown but
detections with no mapped stem anywhere near them, i.e. almost certainly real
trees the NEON woody-veg map never recorded (the field protocol maps a subset
of stems, not every tree). That makes **every precision, F1, and PQ number in
the benchmark a biased lower bound**, and the bias is not uniform across arms:
an arm that finds more real-but-unmapped trees is punished harder. This study
makes the bias measurable and re-grades the leaderboard on de-biased metrics.

The crediting rule (issue #93, option a — co-detection): an isolated core FP
of a target arm is reclassified **probable-real** when arms from ≥ `MIN_FAM`
**other modality families** also leave an isolated FP within `CRED_R` m of it
— several independent systems all see a tree exactly where the map has
nothing. Families (`coverage_lib.R::FAMILY_MAP`): `chm` (CHM-VWF, lmfauto,
multichm), `pc` (ptrees, AMS3D, Li2012 ×2, lmf-pc, lasR-pc), `deep`
(TreeisoNet, SegmentAnyTree, ForestFormer3D), `rgb` (DeepForest, Detectree2;
SOAP only). Errors correlate within a family (all CHM maxima share surface
artefacts; both RGB arms share the imagery), so independence is counted at the
family level and the target's own family never testifies. Two guards keep the
rule honest: **one credit per probable tree** (credited FPs within `CRED_R` of
an already-counted credit are over-segmentation of the same unmapped tree and
stay in the denominator), and witnesses must themselves be *isolated* FPs (a
near-FP is over-seg of a mapped tree, not evidence of an unmapped one).
Credited FPs leave the **precision denominator** (`pool()`'s
`precision_cred`/`F1_cred`); raw precision stays alongside so the bias is
visible, never hidden. Recall and TP are untouched — an unmapped tree can
never become a TP, it just stops counting against the arm.

New, unit-tested pieces (`tests/testthat/test-coverage-gap.R`):

- `sweep_lib.R::fp_points()` — coordinate-level companion to `score_plot`'s
  `fp_near`/`fp_isolated` counts (identical match, region, and core masks).
- `sweep_lib.R::co_detect_credit()` — the ≥ `min_fam`-distinct-families-within-
  `r` test.
- `coverage_lib.R` — `FAMILY_MAP`/`arm_family`, `best_treetop_cache` readers
  (pinned to `best_treetop_selection.csv`, glob fallback, multi-variant
  warning), optical-box → detection conversion (apex z from the native frozen
  CHM, floor 2 m, mirroring `detect_deepforest_sweep.R`), `credit_isolated`
  (family strike + one-credit-per-tree dedup), `read_selection`.
- `model_bench_lib.R::pool()` — pools `fp_credited` by SUM and emits
  `precision_cred`/`F1_cred` (denominator floored at pooled TP).

Regenerate:

```sh
export CLAUDE_JOB_DIR=$(pwd)/work
# CORES=1: the LADDER path runs lasR exec, which drops dense cells under fork
Rscript scripts/coverage_gap.R SITE=SOAP CORES=1
Rscript scripts/coverage_gap.R SITE=SJER CORES=1
Rscript scripts/coverage_gap.R SITE=TEAK CORES=1
# -> work/neon/<SITE>/coverage_gap.csv (pool() -> precision_cred / F1_cred)
```

### What this is

Each arm testifies **at its best tested operating point**: detections come
from the `best_treetop_cache` written by `export_best_treetops_geojson.R`
(each arm at its selected best rung per site — exactly the cells behind the
leaderboard, pinned by `best_treetop_selection.csv`), plus the persisted RGB
crown boxes on SOAP (a plot whose in-plot box set is empty is scored as
ran-and-found-nothing, keeping its stems in the recall denominator). Witness
evidence is rung-independent: a tree's physical existence does not depend on
the decimation rung, so an arm's best-rung isolated FPs testify for every
target cell. `LADDER=1` additionally regenerates the canonical CHM-VWF
detector on every frozen rung (`detector = "chm_vwf_ladder"`) so the bias
curve is measured against density on one arm. Scoring is the benchmark
baseline (`score_plot`, greedy flat-4 m, height gate); pooling is the
canonical sum-of-counts `pool()`. Because arms cover different plot subsets,
the report re-derives both orderings on the **common plot set**
(`equal_set_guard` keyed by site × plot) alongside the per-arm pools.
Defaults: `CRED_R = 2.0` m (the fusion apex-cluster radius), `MIN_FAM = 2`.

### Generated tables

`prec'`/`F1'` = corrected (credited) precision/F1; `cred/iso` = credited over
isolated core FPs.

#### SOAP (18 plots, n_ref 232; rgb family present)

| arm | rung | recall | prec | F1 | prec' | F1' | ΔF1 | cred/iso |
|---|---|--:|--:|--:|--:|--:|--:|--:|
| chm_vwf | 4 | 0.366 | 0.482 | 0.416 | 0.853 | 0.513 | +0.096 | 73/81 |
| lmfauto | 4 | 0.664 | 0.254 | 0.367 | 0.440 | 0.529 | +0.162 | 232/342 |
| multichm | 2 | 0.647 | 0.367 | 0.468 | 0.606 | 0.626 | +0.157 | 147/199 |
| ptrees | 4 | 0.504 | 0.339 | 0.405 | 0.601 | 0.548 | +0.143 | 138/181 |
| ams3d | 2 | 0.677 | 0.353 | 0.464 | 0.528 | 0.593 | +0.129 | 134/190 |
| li2012 | native | 0.595 | 0.257 | 0.359 | 0.477 | 0.530 | +0.171 | 227/321 |
| lidr_li2012 | native | 0.595 | 0.257 | 0.359 | 0.477 | 0.530 | +0.171 | 227/321 |
| lidr_lmf_pc | 8 | 0.487 | 0.271 | 0.349 | 0.510 | 0.498 | +0.150 | 181/258 |
| lasr_lmax_pc | 8 | 0.487 | 0.271 | 0.349 | 0.510 | 0.498 | +0.150 | 181/258 |
| treeisonet | 2 | 0.599 | 0.312 | 0.410 | 0.582 | 0.591 | +0.181 | 185/219 |
| segmentanytree | native | 0.642 | 0.364 | 0.464 | 0.584 | 0.612 | +0.147 | 137/182 |
| forestformer3d | 8 | 0.448 | 0.248 | 0.319 | 0.349 | 0.393 | +0.074 | 119/247 |
| deepforest | rgb | 0.552 | 0.389 | 0.457 | 0.663 | 0.602 | +0.146 | 125/156 |
| detectree2 | rgb | 0.168 | 0.369 | 0.231 | 0.691 | 0.270 | +0.039 | 48/61 |

Credited ordering (identical on the 18-plot equal set): **multichm (0.626) >
segmentanytree (0.612) > deepforest (0.602) > ams3d (0.593) > treeisonet
(0.591)** — DeepForest climbs past ams3d, treeisonet passes ptrees, and
chm_vwf falls from 5th to 10th: its precision "lead" was under-detection (top
corrected precision 0.85 but the lowest recall).

#### SJER (n_ref 71; no rgb family)

| arm | rung | recall | prec | F1 | prec' | F1' | ΔF1 | cred/iso |
|---|---|--:|--:|--:|--:|--:|--:|--:|
| chm_vwf | 4 | 0.437 | 0.308 | 0.361 | 0.757 | 0.554 | +0.193 | 54/58 |
| lmfauto | native | 0.606 | 0.210 | 0.312 | 0.477 | 0.534 | +0.222 | 112/147 |
| multichm | 4 | 0.718 | 0.236 | 0.355 | 0.441 | 0.547 | +0.192 | 97/151 |
| ptrees | 4 | 0.563 | 0.253 | 0.349 | 0.534 | 0.548 | +0.199 | 81/108 |
| ams3d | 1 | 0.451 | 0.250 | 0.322 | 0.345 | 0.391 | +0.069 | 33/72 |
| lidr_li2012 | native | 0.732 | 0.168 | 0.274 | 0.283 | 0.409 | +0.135 | 123/215 |
| treeisonet | native | 0.789 | 0.178 | 0.291 | 0.321 | 0.457 | +0.166 | 135/207 |
| segmentanytree | 2 | 0.507 | 0.298 | 0.376 | 0.486 | 0.496 | +0.121 | 44/69 |
| lidr_lmf_pc | 8 | 0.662 | 0.198 | 0.305 | 0.400 | 0.499 | +0.194 | 105/147 |
| lasr_lmax_pc | 8 | 0.662 | 0.198 | 0.305 | 0.400 | 0.499 | +0.194 | 105/147 |

On the 7-plot equal set the credited ordering is **multichm > ptrees >
chm_vwf > lmfauto > segmentanytree … > ams3d (last)**. In the open savanna the
crediting *discriminates*: AMS3D's FPs are mostly **not** co-detected (33/72),
i.e. genuinely low precision, while chm_vwf's 54/58 are.

#### TEAK (n_ref 358–396; no rgb family)

| arm | rung | recall | prec | F1 | prec' | F1' | ΔF1 | cred/iso |
|---|---|--:|--:|--:|--:|--:|--:|--:|
| chm_vwf | native | 0.394 | 0.421 | 0.407 | 0.746 | 0.515 | +0.108 | 130/162 |
| lmfauto | 2 | 0.571 | 0.347 | 0.432 | 0.514 | 0.541 | +0.110 | 200/328 |
| multichm | 8 | 0.523 | 0.378 | 0.439 | 0.593 | 0.555 | +0.116 | 168/253 |
| ptrees | 8 | 0.336 | 0.340 | 0.338 | 0.626 | 0.437 | +0.099 | 157/207 |
| ams3d | 4 | 0.535 | 0.337 | 0.414 | 0.492 | 0.513 | +0.099 | 181/306 |
| lidr_li2012 | 8 | 0.422 | 0.358 | 0.387 | 0.630 | 0.505 | +0.118 | 181/235 |
| lidr_lmf_pc | 8 | 0.336 | 0.367 | 0.351 | 0.717 | 0.457 | +0.107 | 158/193 |
| lasr_lmax_pc | 8 | 0.336 | 0.367 | 0.351 | 0.717 | 0.457 | +0.107 | 158/193 |
| treeisonet | native | 0.359 | 0.332 | 0.345 | 0.604 | 0.450 | +0.105 | 169/222 |
| segmentanytree | native | 0.654 | 0.383 | 0.483 | 0.567 | 0.608 | +0.124 | 200/302 |

SegmentAnyTree leads raw **and** credited (F1' 0.608, equal-set confirmed);
chm_vwf passes ams3d under crediting.

#### CHM-VWF ladder: how the bias moves with density

Pooled ΔF1 = F1' − F1 for the regenerated canonical CHM-VWF per rung:

| rung | SOAP ΔF1 | SJER ΔF1 | TEAK ΔF1 |
|---|--:|--:|--:|
| native | +0.167 | +0.225 | +0.100 |
| 8 | +0.092 | +0.199 | +0.070 |
| 4 | +0.093 | +0.192 | +0.067 |
| 2 | +0.091 | +0.171 | +0.063 |
| 1 | +0.089 | +0.151 | +0.056 |

The bias **grows with density** on every site: denser data finds more of the
real-but-unmapped trees, so raw precision punishes native hardest. Density
ladders scored on raw F1 therefore **under-state the native advantage**.

#### Crediting-rule sensitivity (pooled ΔF1, leaderboard cells)

| site | r=1.5 f=1 | r=1.5 f=2 | r=2 f=1 | r=2 f=2 | r=3 f=1 | r=3 f=2 |
|---|--:|--:|--:|--:|--:|--:|
| SOAP | +0.179 | +0.127 | +0.184 | +0.141 | +0.162 | +0.139 |
| SJER | +0.255 | +0.143 | +0.269 | +0.167 | +0.245 | +0.189 |
| TEAK | +0.141 | +0.102 | +0.145 | +0.110 | +0.131 | +0.109 |

Smooth in both knobs, and no longer monotone in `r`: past ~2 m the
one-credit-per-tree dedup suppresses more than the wider radius admits — the
credit is a cluster count, not a rubber stamp. The default (r=2, f=2) is the
conservative middle.

### Readings

- **The coverage gap is the dominant precision error.** 46–93% of isolated
  core FPs are co-detected by ≥2 independent modality families. Corrected
  pooled precision rises by +0.1 to +0.4 absolute; every arm's F1 was a lower
  bound, exactly as #V4 predicted.
- **The leaderboard reorders on every site** (equal-set confirmed). SOAP:
  multichm holds #1 but the gap to SegmentAnyTree/DeepForest collapses to
  0.01–0.02, and chm_vwf drops 5th → 10th. SJER: multichm/ptrees lead and
  AMS3D falls to last — its savanna FPs are genuinely uncorroborated. TEAK:
  SegmentAnyTree's #1 is confirmed and widens.
- **The crediting discriminates rather than inflates**: credit rates range
  from 33/72 (AMS3D SJER) and 119/247 (ForestFormer3D SOAP) to 54/58
  (chm_vwf SJER) — arms with genuinely noisy detections keep low precision,
  and the dedup guard stops over-segmentation of unmapped trees from
  laundering (it cost lmfauto 20 and ForestFormer3D 28 credits on SOAP).
- **Optical co-detection matters on SOAP**: DeepForest's credited F1' 0.602
  puts an RGB-only detector level with the best LiDAR arms, and the rgb family
  testifies for LiDAR arms where chm/pc/deep thin out — the #P7
  fusion-membership argument. Detectree2, honestly scored over all 18 plots
  (zero-box plots kept), sits at recall 0.168 — fine-tuning (#92 checklist)
  remains its gate.
- **Router/fusion consequences (#P2/#P1/#92)**: k-of-N and router thresholds
  were tuned on raw F1; with corrected metrics the precision cost of the union
  mode shrinks substantially (most union-only detections are credited), so the
  fusion Pareto and the router's per-cell argmax should be re-derived on F1'.

### Caveats

- **Crediting is evidence-based reclassification, not ground truth.** A
  correlated cross-family artefact (e.g. a boulder tall enough for the CHM and
  textured enough for RGB) can still slip through; the rigorous close-outs
  remain hand-delineated crowns (#93 option b) and an externally-complete
  benchmark (FGI-EMIT, #D1). Treat `F1'` as a de-biased estimate bracketed by
  the sensitivity grid, not a replacement truth.
- The Li2012 twins (`li2012`/`lidr_li2012`) and the LM twins
  (`lidr_lmf_pc`/`lasr_lmax_pc`) are near-duplicate detectors; the family rule
  already collapses them for crediting, but their table rows are not
  independent evidence.
- Per-arm pools sit on each arm's own cell set (TEAK `chm_vwf` 17 plots, the
  SJER pc_ladder arms 6); the equal-set orderings are the cross-arm claim.
- Raw chm_vwf here (SOAP F1 0.416) differs from the historical sweep number
  (0.426) by regeneration variance: the cache regenerates from the frozen
  clips, the original sweep drew its own decimation.
- **PQ correction is deferred to #V6/#94**: this study credits at the apex
  level; the instance-level analogue (drop credited FP instances from RQ's
  denominator via a `pool_pq` hook mirroring `pool()`) becomes worthwhile once
  the classical arms persist masks and the IoU/PQ board covers all arms.
- SJER/TEAK lack the rgb family, so `MIN_FAM=2` there means chm/pc/deep
  agreement only; SOAP's numbers are the strongest-evidence configuration.
