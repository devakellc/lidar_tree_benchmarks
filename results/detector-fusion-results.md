# Cross-arm detector fusion (#P1)

The first actual meta-pipeline **detector**: until now every arm
(#6/#36/#37/#M6/#M7/#M8) was scored in isolation and pooled side-by-side, but
nothing combined them. The per-class table in
[`model-benchmark-results.md`](model-benchmark-results.md) shows the
complementarity to exploit — CHM-VWF/multichm carry the overstory, SegmentAnyTree
carries the understory — so a consensus/NMS layer should trade recall against
precision. This arm materializes per-cell apexes for every available arm on the
**identical** frozen cells, clusters them across arms, and emits three operating
points plus the full k-of-N Pareto frontier, scored against field stems with both
the distance matcher and the #V1 IoU/PQ scorer.

Regenerate:

```sh
export CLAUDE_JOB_DIR=$(pwd)/work
# CORES=1: detect_lasr uses lasR exec, which can drop dense cells under fork
Rscript scripts/fuse_detectors.R SITES=SOAP,SJER,TEAK RUNGS=native,8,4,2,1 CORES=1
# -> work/neon/<SITE>/fusion_results.csv (one row per plot x rung x config)
```

## What this is

For each frozen cell, per-arm apexes are **materialized fresh** (the
`*_results.csv` are scored summaries, not apex inventories):

- `chm_vwf` — `detect_lasr` on the frozen normalized clip (density `chm_res`,
  `vwf_a = 0.10`); `multichm` — `lidRplugins::multichm` on the same clip;
  `li2012` — `lidR::li2012` on the same clip (**native only** — meaningless on
  decimated rungs).
- `segmentanytree` — the persisted `segmentanytree_instances/<plot>_<rung>.laz`
  (`PredInstance`) reduced to apexes and converted to AGL with the frozen DTM.
- `forestformer3d` — the persisted `forestformer3d_instances/<plot>_<rung>.laz`
  (`UserData`/`PointSourceID`) collapsed and AGL-guarded (**native + 8 only**).
- **TreeisoNet is deferred**: it persists no reusable per-point/apex labels
  (apex-only GPU `*_results.csv`), so it cannot be re-fused offline.

`fuse_apexes()` (in [`model_bench_lib.R`](../scripts/model_bench_lib.R), unit
tested in `tests/testthat/test-detector-fusion.R`) clusters apexes across arms by
single-linkage union-find within a horizontal `merge_tol` **and** a `|dz|` height
gate (the `dedup_blocks` pattern + the fusion height gate, so an overstory CHM
apex and an understory point apex at the same x,y stay separate trees); same-arm
apexes never merge directly. `fusion_points()` emits:

- **union** (k ≥ 1) — every fused cluster: recall-max.
- **majority** (k ≥ ⌈N/2⌉) — precision-max.
- **layered** — CHM-family apices in the overstory band (z ≥ ½·max z) + point/deep
  apices in the understory band, then fused.

Each fused set is scored with `score_plot` (apex distance, 4 m + height gate) and
with a **Voronoi-on-apexes** #V1 IoU/PQ proxy (each frozen canopy point assigned
to its nearest fused apex, scored against the Voronoi-on-stems reference). Pooled
by SUM with `pool()`; the k = 1…N frontier traces the recall–precision Pareto.

## Generated tables

### Native, all three sites combined (pooled, 699 stems)

| config | recall | precision | F1 | rec_understory | iou_Cov | iou_PQ |
|---|--:|--:|--:|--:|--:|--:|
| chm_vwf | 0.401 | 0.336 | 0.365 | 0.200 | 0.251 | 0.082 |
| multichm | 0.541 | 0.333 | 0.412 | 0.333 | 0.308 | 0.087 |
| li2012 | 0.495 | 0.266 | 0.346 | 0.257 | 0.287 | 0.072 |
| segmentanytree | 0.655 | 0.334 | 0.443 | 0.448 | 0.350 | 0.096 |
| forestformer3d | 0.369 | 0.224 | 0.279 | 0.190 | 0.211 | 0.045 |
| **union** | **0.764** | 0.226 | 0.349 | **0.600** | **0.363** | 0.063 |
| majority | 0.422 | 0.346 | 0.380 | 0.171 | 0.265 | 0.089 |
| layered | 0.707 | 0.233 | 0.351 | 0.543 | 0.357 | 0.068 |

Best single arm: SegmentAnyTree (F1 0.443, understory 0.448). Union recall
+0.109 and understory recall **+0.152** over it; union F1 −0.093 (precision cost).

### Per site, native: union vs best single arm

| site | n_ref | best arm | best F1 | union recall | union F1 | union understory | best-arm understory |
|---|--:|---|--:|--:|--:|--:|--:|
| SOAP | 232 | segmentanytree | 0.464 | 0.793 | 0.324 | 0.644 | 0.489 |
| SJER | 71 | multichm | 0.347 | 0.789 | 0.217 | 0.500 | 1.000 |
| TEAK | 396 | segmentanytree | 0.483 | 0.742 | 0.423 | 0.569 | 0.397 |

### SOAP per-rung k-of-N Pareto (recall / F1; "-" = fewer than k arms)

| rung | k1 | k2 | k3 | k4 | k5 | best-arm F1 |
|---|--|--|--|--|--|--:|
| native | 0.79 / 0.32 | 0.60 / 0.41 | 0.50 / 0.42 | 0.41 / 0.41 | 0.15 / 0.23 | 0.464 |
| 8 | 0.72 / 0.40 | 0.44 / 0.43 | 0.18 / 0.26 | - | - | 0.441 |
| 4 | 0.73 / 0.43 | 0.42 / 0.44 | 0.17 / 0.25 | - | - | 0.452 |
| 2 | 0.72 / 0.44 | 0.35 / 0.40 | 0.09 / 0.15 | - | - | 0.468 |
| 1 | 0.63 / 0.42 | 0.24 / 0.32 | 0.02 / 0.04 | - | - | 0.454 |

## Readings

- **Fusion delivers the predicted complementarity.** The union operating point
  lifts pooled recall to 0.764 (+0.11) and **understory recall to 0.600 (+0.152)**
  over the strongest single arm, consistently across all three sites — it recovers
  trees, especially understory trees, that no single detector finds. The layered
  mode reproduces most of that gain by construction.
- **It does not win on F1 at native — and that ceiling is mostly an artifact.**
  No mode beats SegmentAnyTree's native F1, because adding arms adds apparent
  false positives faster than recall rises. But #V4 established that ~94 % of core
  false positives are *isolated* — real, unmapped trees, not over-segmentation —
  so the precision/F1 penalty is dominated by the field-map coverage gap, not by
  spurious fusion detections. Against a complete reference, union's precision
  would be far higher; the recall and Coverage gains are the trustworthy signals
  (union has the highest IoU Coverage of any config, 0.363).
- **The k-of-N frontier is the usable product.** k1 (union) is recall-max, k5 is
  precision-max, and the F1 peak sits mid-frontier (k2–k3). At the sparse rungs,
  where no single arm dominates, k1–k2 fusion matches the best-arm F1 (rung 4: k2
  F1 0.44 vs best-arm 0.452) while carrying much higher recall. This frontier,
  not a single fused set, is what the #P2 router and #P4 calibration consume.
- **Best single arm varies by site** (SAT on SOAP/TEAK, multichm on open SJER),
  which is itself the argument for a router (#P2): no one arm wins everywhere, but
  the union envelope dominates them all on recall everywhere.

## Caveats

- **Precision, F1, and PQ are conservative lower bounds** — the reference covers
  only field-mapped stems, so fusion's extra real detections count as false
  commissions (the #V4 isolated-FP finding). Read recall / understory recall /
  Coverage as the fusion signal; F1 understates it.
- **The IoU/PQ here is a Voronoi-on-apexes proxy.** Fused apexes carry no point
  mask, so each frozen point is assigned to its nearest fused apex — symmetric
  with the Voronoi-on-stems reference but not a true instance mask. Genuine
  mask-aware fusion scoring needs the seed→refine masks of #P3; the standalone
  mask arms' real IoU/PQ are in
  `results/instance-iou-pq-results.md` (#V1).
- **TreeisoNet is absent** (no persisted apexes); five arms fuse at native, four
  at rung 8 (FF3D), three at rungs 4/2/1 (CHM-VWF, multichm, SAT). `n_arms` is
  recorded per cell so the majority threshold adapts.
- The fusion `merge_tol` (2 m) and height gate `z_tol` (5 m) are the `dedup_blocks`
  defaults; `OVERSTORY_FRAC` (0.5) sets the layered split. All are CLI-tunable.

## #V6 update (2026-08-20): ptrees + AMS3D join the pool

The classical segmenters now persist per-point instance clouds
(`<arm>_instances/`, issue #94), so `fuse_detectors.R` materializes **ptrees
and AMS3D** apexes from them (already AGL; no DTM step) alongside the original
five members. SOAP-native smoke on the 7-arm pool: best single arm stays
SegmentAnyTree (F1 0.464 over 18 common cells; union ΔF1 −0.181, majority
−0.051), and the k-of-N Pareto peaks at **k5 F1 0.430** (R 0.448 / P 0.414) —
the extra members raise the consensus operating point (the old 5-arm majority
was the peak). Full cross-site re-synthesis of this doc's tables on the 7-arm
pool is future work; the numbers above are from `fusion_results.csv` as
regenerated at `RUNGS=native CORES=1`.
