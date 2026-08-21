# Coverage-gap crediting: de-biased precision/F1 (#V5)

The #V4 study (`matcher_robustness.R`) found that ~94% of the benchmark's core
false positives are **isolated** — not over-segmentation of a matched crown but
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
family level and the target's own family never testifies. Three guards keep the
rule honest: **one credit per probable tree** (credited FPs within `CRED_R` of
an already-counted credit are over-segmentation of the same unmapped tree and
stay in the denominator); witnesses must themselves be eligible FPs (a near-FP
is over-seg of a mapped tree, not evidence of an unmapped one); and
**eligibility is measured against every mapped stem of the plot, matched or
not** (`fp_points`' `credit_eligible`, fed the plot's full stem set via
`elig_stems`) — an FP sitting on a stem the arm merely failed to match, or on
one of the 15-47 stems per site whose mapped position falls just outside the
scored core, is already explained by the field map, so it never credits.
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
  (rung **and** `chm_res`/`vwf_a` pinned from `best_treetop_selection.csv` so a
  stale parameter variant can never shadow the selected leaderboard cell; glob
  fallback warns), optical-box → detection conversion (apex z from the native
  frozen CHM, floor 2 m, mirroring `detect_deepforest_sweep.R`),
  `credit_isolated` (family strike + eligibility + one-credit-per-tree dedup),
  `read_selection`.
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

## What this is

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

`prec'`/`F1'` = corrected (credited) precision/F1; `cred/elig` = credited
over credit-ELIGIBLE core FPs (isolated AND not near any mapped stem of the
plot, in-core or not), the denominator `coverage_gap.R` prints.

#### SOAP (18 plots, n_ref 232; rgb family present)

| arm | rung | recall | prec | F1 | prec' | F1' | ΔF1 | cred/elig |
|---|---|--:|--:|--:|--:|--:|--:|--:|
| chm_vwf | 4 | 0.366 | 0.482 | 0.416 | 0.794 | 0.501 | +0.085 | 66/69 |
| lmfauto | 4 | 0.664 | 0.254 | 0.367 | 0.426 | 0.519 | +0.152 | 222/328 |
| multichm | 2 | 0.647 | 0.367 | 0.468 | 0.573 | 0.608 | +0.139 | 134/178 |
| ptrees | 4 | 0.504 | 0.339 | 0.405 | 0.578 | 0.539 | +0.134 | 131/166 |
| ams3d | 2 | 0.677 | 0.353 | 0.464 | 0.516 | 0.586 | +0.122 | 128/184 |
| li2012 | native | 0.595 | 0.257 | 0.359 | 0.453 | 0.514 | +0.156 | 213/301 |
| lidr_li2012 | native | 0.595 | 0.257 | 0.359 | 0.453 | 0.514 | +0.156 | 213/301 |
| lidr_lmf_pc | 8 | 0.487 | 0.271 | 0.349 | 0.488 | 0.488 | +0.139 | 172/244 |
| lasr_lmax_pc | 8 | 0.487 | 0.271 | 0.349 | 0.488 | 0.488 | +0.139 | 172/244 |
| treeisonet | 2 | 0.599 | 0.312 | 0.410 | 0.554 | 0.575 | +0.166 | 174/205 |
| segmentanytree | native | 0.642 | 0.364 | 0.464 | 0.562 | 0.599 | +0.135 | 128/165 |
| forestformer3d | 8 | 0.448 | 0.248 | 0.319 | 0.344 | 0.389 | +0.070 | 114/230 |
| deepforest | rgb | 0.552 | 0.389 | 0.457 | 0.634 | 0.590 | +0.134 | 117/142 |
| detectree2 | rgb | 0.168 | 0.369 | 0.231 | 0.667 | 0.269 | +0.038 | 46/58 |

Credited ordering (identical on the 18-plot equal set): **multichm (0.608) >
segmentanytree (0.599) > deepforest (0.590) > ams3d (0.586) > treeisonet
(0.575)** — DeepForest climbs past ams3d, treeisonet passes ptrees, and
chm_vwf falls from 5th to 10th: its precision "lead" was under-detection (top
corrected precision 0.79 but the lowest recall).

#### SJER (n_ref 71; no rgb family)

| arm | rung | recall | prec | F1 | prec' | F1' | ΔF1 | cred/elig |
|---|---|--:|--:|--:|--:|--:|--:|--:|
| chm_vwf | 4 | 0.437 | 0.308 | 0.361 | 0.757 | 0.554 | +0.193 | 54/58 |
| lmfauto | native | 0.606 | 0.210 | 0.312 | 0.472 | 0.530 | +0.219 | 111/146 |
| multichm | 4 | 0.718 | 0.236 | 0.355 | 0.438 | 0.544 | +0.189 | 96/150 |
| ptrees | 4 | 0.563 | 0.253 | 0.349 | 0.527 | 0.545 | +0.195 | 80/107 |
| ams3d | 1 | 0.451 | 0.250 | 0.322 | 0.341 | 0.388 | +0.067 | 32/70 |
| lidr_li2012 | native | 0.732 | 0.168 | 0.274 | 0.282 | 0.407 | +0.133 | 122/213 |
| treeisonet | native | 0.789 | 0.178 | 0.291 | 0.320 | 0.455 | +0.164 | 134/204 |
| segmentanytree | 2 | 0.507 | 0.298 | 0.376 | 0.479 | 0.493 | +0.117 | 43/67 |
| lidr_lmf_pc | 8 | 0.662 | 0.198 | 0.305 | 0.392 | 0.492 | +0.188 | 100/145 |
| lasr_lmax_pc | 8 | 0.662 | 0.198 | 0.305 | 0.392 | 0.492 | +0.188 | 100/145 |

On the 7-plot equal set the credited ordering is **multichm > ptrees >
chm_vwf > lmfauto > segmentanytree … > ams3d (last)**. In the open savanna the
crediting *discriminates*: AMS3D's FPs are mostly **not** co-detected (32/70),
i.e. genuinely low precision, while chm_vwf's 54/58 are.

#### TEAK (n_ref 358–396; no rgb family)

| arm | rung | recall | prec | F1 | prec' | F1' | ΔF1 | cred/elig |
|---|---|--:|--:|--:|--:|--:|--:|--:|
| chm_vwf | native | 0.394 | 0.421 | 0.407 | 0.689 | 0.501 | +0.094 | 116/143 |
| lmfauto | 2 | 0.571 | 0.347 | 0.432 | 0.500 | 0.533 | +0.102 | 188/310 |
| multichm | 8 | 0.523 | 0.378 | 0.439 | 0.564 | 0.543 | +0.104 | 153/226 |
| ptrees | 8 | 0.336 | 0.340 | 0.338 | 0.585 | 0.427 | +0.089 | 144/183 |
| ams3d | 4 | 0.535 | 0.337 | 0.414 | 0.479 | 0.506 | +0.092 | 170/293 |
| lidr_li2012 | 8 | 0.422 | 0.358 | 0.387 | 0.584 | 0.490 | +0.102 | 162/205 |
| lidr_lmf_pc | 8 | 0.336 | 0.367 | 0.351 | 0.643 | 0.441 | +0.090 | 139/165 |
| lasr_lmax_pc | 8 | 0.336 | 0.367 | 0.351 | 0.643 | 0.441 | +0.090 | 139/165 |
| treeisonet | native | 0.359 | 0.332 | 0.345 | 0.551 | 0.434 | +0.089 | 149/186 |
| segmentanytree | native | 0.654 | 0.383 | 0.483 | 0.555 | 0.601 | +0.117 | 191/281 |

SegmentAnyTree leads raw **and** credited (F1' 0.601). TEAK is the one site
whose per-arm ordering does **not** move under crediting; on the 17-plot equal
set the only move is TreeisoNet rising past the two LM twins
(`lidr_lmf_pc`/`lasr_lmax_pc`), 9th to 7th.

#### CHM-VWF ladder: how the bias moves with density

Pooled ΔF1 = F1' − F1 for the regenerated canonical CHM-VWF per rung:

| rung | SOAP ΔF1 | SJER ΔF1 | TEAK ΔF1 |
|---|--:|--:|--:|
| native | +0.151 | +0.225 | +0.084 |
| 8 | +0.081 | +0.199 | +0.061 |
| 4 | +0.082 | +0.192 | +0.059 |
| 2 | +0.084 | +0.166 | +0.052 |
| 1 | +0.081 | +0.146 | +0.049 |

The bias **grows with density** on every site: denser data finds more of the
real-but-unmapped trees, so raw precision punishes native hardest. Density
ladders scored on raw F1 therefore **under-state the native advantage**.

#### Crediting-rule sensitivity (pooled ΔF1, leaderboard cells)

| site | r=1.5 f=1 | r=1.5 f=2 | r=2 f=1 | r=2 f=2 | r=3 f=1 | r=3 f=2 |
|---|--:|--:|--:|--:|--:|--:|
| SOAP | +0.164 | +0.117 | +0.168 | +0.130 | +0.148 | +0.128 |
| SJER | +0.250 | +0.141 | +0.263 | +0.164 | +0.240 | +0.185 |
| TEAK | +0.122 | +0.091 | +0.126 | +0.098 | +0.113 | +0.096 |

Smooth in both knobs, and no longer monotone in `r` in five of the six columns
(SJER `f=2` still rises): past ~2 m the one-credit-per-tree dedup suppresses
more than the wider radius admits — the credit is a cluster count, not a rubber
stamp. The default (r=2, f=2) is the conservative middle.

### Readings

- **The coverage gap is the dominant precision error.** 46–96% of
  credit-eligible core FPs are co-detected by ≥2 independent modality
  families. Corrected pooled precision rises by +0.09 (AMS3D SJER) to +0.45
  (chm_vwf SJER) absolute; every arm's F1 was a lower bound, exactly as #V4
  predicted.
- **The leaderboard reorders on two of three sites** (equal-set confirmed).
  SOAP: multichm holds #1 but the gap to SegmentAnyTree/DeepForest collapses
  to 0.009–0.018, and chm_vwf drops 5th → 10th. SJER: multichm/ptrees lead and
  AMS3D falls to last — its savanna FPs are genuinely uncorroborated. TEAK is
  stable: SegmentAnyTree's #1 is confirmed and only TreeisoNet moves (past the
  LM twins) on the equal set.
- **The crediting discriminates rather than inflates**: credit rates range
  from 32/70 (AMS3D SJER) and 114/230 (ForestFormer3D SOAP) to 66/69
  (chm_vwf SOAP) — arms with genuinely noisy detections keep low precision.
  Both honesty guards bite: dedup plus stem-adjacency eligibility cost ptrees
  13 credits and ForestFormer3D 33 on SOAP versus the naive rule.
- **Optical co-detection matters on SOAP**: DeepForest's credited F1' 0.590
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
- **Eligibility spans the plot's whole mapped stem set, scoring does not.**
  `score_plot`/`n_ref` stay on the core box (±10 m distributed, ±20 m tower),
  but a stem whose reconstructed position lands outside that box is still
  mapped, so it blocks crediting (`elig_stems`). Restricting eligibility to the
  core stems instead credits 1–7 more FPs per arm and lifts F1' by up to 0.008.
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
- **Cache pinning is exact only for CHM-VWF.** `best_treetop_selection.csv`
  records `chm_res`/`vwf_a` but not the GPU arms' `conf`/`voxel`/`image`/
  `spacing`/`merge` suffixes, so `treeisonet`/`segmentanytree`/`forestformer3d`
  cells still resolve by glob. That is exact whenever only one variant was ever
  cached (the committed state) and warns whenever more than one exists, but
  pinning those arms properly needs the selection writer to record their knobs.
  A pinned-but-missing CHM-VWF variant also warns and falls back to the glob
  rather than dropping the cell, so a re-parameterized re-run is loud, not
  silent.
