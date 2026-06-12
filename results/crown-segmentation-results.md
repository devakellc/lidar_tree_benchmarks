# Crown-Segmentation Benchmark — Results

*Closes the loop on [GitHub issue #7](https://github.com/agrigoriev/lidar_tree_benchmarks/issues/7):
the density-ladder sweep scores tree-top **detection** only; this benchmark
delineates **crowns** from the detected tops and scores their diameter against
NEON field crown diameter. Five segmenters are seeded from the **same** detected
tops, run per plot on a native-density pit-free CHM, matched back to field
stems, and scored RMSE/MAE/bias/R² by crown class. A SOAP-only TreeisoNet
`treeOff` crown arm and the deep instance segmenters SegmentAnyTree (#M6) and
ForestFormer3D (#M8) are unioned by
[`scripts/analyze_crown_metrics.R`](../scripts/analyze_crown_metrics.R) into the
deep-model sections below. Driver:
[`scripts/crown_metrics_sweep.R`](../scripts/crown_metrics_sweep.R). Sites: NEON
SJER (open oak savanna), SOAP (mixed conifer), TEAK (red-fir). Last run:
2026-06-05.*

---

## TL;DR

- **One detector, five crowns, on a shared CHM.** Tree-tops are detected once
  per plot (`locate_trees(lmf, ws=ws_factory(0.10))` on a pit-free CHM at 0.5 m)
  and reused as the seed set for every seeded segmenter. lasR `region_growing`
  cannot take an external point set as seeds, so it is seeded by injecting the
  *same* lidR pit-free CHM into the lasR pipeline (`load_raster`) and running
  `local_maximum_raster` with the *same* `ws` — yielding a near-identical seed
  set (81–94% of seeds within 0.5 m of a shared top), so the crown differences
  are a growing-rule effect on a shared surface and seed set, not a seed
  confound. 225 field stems matched across 40 plots (SJER 22, SOAP 87,
  TEAK 116).
- **lasR `region_growing` still tracks NEON crown diameter best, but the
  margin is narrow once the seed confound is removed.** Pooled over all sites it
  has the lowest RMSE/MAE/bias on both diameter definitions. Equivalent-circle
  diameter vs `ninetyCrownDiameter`: RMSE 2.62 m, R² +0.10 (dalponte2016 is a
  close second, 2.70 m / R² +0.05). Max-caliper vs `maxCrownDiameter`:
  RMSE 3.72 m, R² −0.14 — best of the five but **no method reaches a positive
  R² on the caliper definition**. The earlier "only method with positive R² on
  either definition" claim no longer holds: it rested on lasR's old independent
  seed set/CHM and shrinks to a modest lead on a fair shared-seed footing.
- **The random walker ran successfully** (sparse Dirichlet solve via `Matrix`,
  sub-second per plot, never timed out) but **underperforms**: like the
  marker-free watershed it over-grows crowns into a full Voronoi-style tiling of
  the canopy, inflating diameter (pooled bias +1.85 m / +3.15 m). Seed-less
  canopy islands are now dropped to background instead of being mislabelled into
  crown 1 (see below). It is an honest negative result, documented below.
- **Geometric caveat is large and systematic.** Equivalent-circle diameter
  `d_eq = 2√(area/π)` underestimates the widest-axis `maxCrownDiameter` because
  a real crown is not a disc: pooled `d_caliper`-vs-`maxCD` RMSE (3.72–5.61 m)
  is ~1.4–1.6× the `d_eq`-vs-`ninetyCD` RMSE (2.62–3.42 m) for the same crowns.
  Pick the diameter definition that matches the field column.

---

## Method

For each site, plots with ≥6 live mapped trees carrying a non-NA NEON field
crown diameter are processed at **native density** (no decimation;
`prepare_clip(rung=NA)` from [`sweep_lib.R`](../scripts/sweep_lib.R)):

1. **Shared seeds.** One pit-free CHM (lidR `pitfree`, mirrors
   [`segment_lidr.R`](../scripts/segment_lidr.R)) at `RES=0.5 m`; detect tops once
   via `locate_trees(chm, lmf(ws=ws_factory(0.10), hmin=2, shape="circular"))`.
2. **Five segmenters**, each producing a per-crown polygon → crown area:
   - lidR `dalponte2016(chm, ttops, th_seed=0.45, th_cr=0.55, max_cr=10/res px)`
     — seeded region-growing.
   - lidR `silva2016(chm, ttops, max_cr_factor=0.6, exclusion=0.3)` — seeded
     Voronoi-like.
   - lidR `watershed(chm, th_tree=2)` — **marker-free** (EBImage); crowns
     matched to stems by polygon containment of the seed point.
   - lasR `region_growing(pit_fill, seed, th_tree=2, th_seed=0.45, th_cr=0.55,
     max_cr=10)` (Dalponte growing rule, mirroring
     [`segment_lasr.R`](../scripts/segment_lasr.R)). lasR's `region_growing` takes
     a seed *stage*, not the shared `ttops` point set, so to grow from the same
     seeds the **same lidR pit-free CHM** is injected into the lasR pipeline via
     `load_raster` and `local_maximum_raster` is run on it with the **same**
     `ws`. The resulting seed set matches the shared `locate_trees` tops closely
     (81–94% of seeds within 0.5 m of a shared top across the three sites); the
     residual is the local-maximum implementation difference (lidR `lmf`
     circular vs lasR `local_maximum_raster`), not a different surface. So the
     lasR-vs-lidR crown differences are a growing-rule effect on a shared CHM
     and a shared seed set, not a seed-set confound.
   - **Random walker** (Grady 2006) on the CHM raster, seeded from the same
     detected tops, implemented with `Matrix`: 4-neighbour pixel graph over
     canopy pixels (Z ≥ 2 m), Gaussian edge weights
     `w_ij = exp(−β·(chmᵢ−chmⱼ)²)` on the normalized CHM (β = 1), one marker
     label per seed, solving the combinatorial Dirichlet problem
     `Lᵤ x = −B m` for marker probabilities and argmax-labelling each pixel.
     Pixels in a canopy component that no seed reaches get ~0 probability for
     every marker; those are assigned to background, not forced into crown 1 by
     argmax.
3. **Two diameter estimates per crown** (report both, see caveat):
   - `d_eq = 2√(area/π)` (equivalent-circle) → compared to `ninetyCrownDiameter`.
   - `d_caliper = max pairwise polygon-vertex distance` (max axis) → compared to
     `maxCrownDiameter` (NEON's widest-axis measurement).
4. **Matching.** The shared seeds are matched to field stems with
   `greedy_match` (global nearest-distance 1:1, position tol 4 m + height gate,
   from `sweep_lib.R`). Each matched stem's crown is the crown owning that seed
   (by `treeID` for the seeded lidR methods; by nearest seed-carrying crown for
   lasR/watershed/RW). The pair *(detected diameter, field diameter)* feeds the
   error stats.
5. **Scoring.** Pooled over matched trees (sum of squared errors / n, never a
   mean of per-plot rates): RMSE, MAE, bias (detected − field), R².

Field crown diameter is joined from the cached
`work/neon/<SITE>/vst/<site>_vst_allyears.rds`
(`vst_apparentindividual` → `maxCrownDiameter`, `ninetyCrownDiameter`), deduped
to the nearest-to-2021 measurement per `individualID`. The same two columns were
added to [`scripts/neon_ground_truth.R`](../scripts/neon_ground_truth.R) (additive,
not re-run) so future ground-truth regenerations carry them.

Field crown diameter (225 matched stems across 40 plots): `ninetyCrownDiameter`
mean 4.94 m (median 4.30, range 1.1–17.6); `maxCrownDiameter` mean 5.97 m
(median 5.20, range 1.3–26.1).

---

## Results — pooled over all sites

### Equivalent-circle `d_eq` vs `ninetyCrownDiameter`

| Algorithm | n | RMSE (m) | MAE (m) | bias (m) | R² |
|-----------|---:|---:|---:|---:|---:|
| **lasr_region_growing** | 225 | **2.62** | **2.03** | **+1.11** | **+0.102** |
| dalponte2016 | 225 | 2.70 | 2.06 | +1.16 | +0.046 |
| silva2016 | 225 | 2.79 | 2.22 | +1.34 | −0.021 |
| random_walker | 225 | 3.19 | 2.57 | +1.85 | −0.334 |
| watershed (marker-free) | 225 | 3.42 | 2.70 | +2.31 | −0.530 |

### Max-caliper `d_caliper` vs `maxCrownDiameter`

| Algorithm | n | RMSE (m) | MAE (m) | bias (m) | R² |
|-----------|---:|---:|---:|---:|---:|
| **lasr_region_growing** | 225 | **3.72** | **2.94** | **+2.04** | **−0.142** |
| silva2016 | 225 | 4.40 | 3.63 | +2.90 | −0.596 |
| dalponte2016 | 225 | 4.43 | 3.47 | +2.73 | −0.621 |
| random_walker | 225 | 4.69 | 3.85 | +3.15 | −0.817 |
| watershed (marker-free) | 225 | 5.61 | 4.54 | +4.21 | −1.595 |

**lasR `region_growing` has the lowest RMSE/MAE/bias on both definitions, but
the lead is narrow on a shared seed set.** On `d_eq` it edges dalponte2016
(2.62 vs 2.70 m; R² +0.10 vs +0.05) — both are positive-R², so lasR is no longer
the *only* method above zero there. On `d_caliper` it is best of the five
(3.72 m) but **every method, lasR included, has a negative R²**; the old
positive caliper R² was an artefact of lasR's previous independent seed set and
CHM. Its conservative growing rule (`max_cr=10` data-units, stopped at
`th_cr=0.55` of region mean height) still keeps crowns the most compact, whereas
the marker-free watershed and the random walker tile the canopy and run +2 to
+4 m high on diameter.

### RMSE by crown class — `d_eq` vs `ninetyCrownDiameter`

| Algorithm | dominant | codominant | intermediate | suppressed |
|-----------|---:|---:|---:|---:|
| | n=101 | n=109 | n=13 | n=2 |
| **lasr_region_growing** | **2.83** | **2.39** | 2.54 | **3.65** |
| dalponte2016 | 2.87 | 2.53 | **2.32** | 4.39 |
| silva2016 | 2.95 | 2.64 | 2.44 | 4.48 |
| random_walker | 3.40 | 3.03 | 2.48 | 4.41 |
| watershed (marker-free) | 3.39 | 3.40 | 3.56 | 4.77 |

### RMSE by crown class — `d_caliper` vs `maxCrownDiameter`

| Algorithm | dominant | codominant | intermediate | suppressed |
|-----------|---:|---:|---:|---:|
| | n=101 | n=109 | n=13 | n=2 |
| **lasr_region_growing** | **3.76** | **3.67** | 3.90 | **2.54** |
| silva2016 | 4.24 | 4.53 | 4.36 | 4.95 |
| dalponte2016 | 4.23 | 4.55 | 4.43 | 6.97 |
| random_walker | 4.99 | 4.52 | **3.60** | 4.24 |
| watershed (marker-free) | 5.40 | 5.78 | 5.81 | 4.74 |

lasR `region_growing` is best or tied-best in most crown classes on both
definitions (dalponte2016 edges it on intermediate `d_eq`, the random walker on
intermediate `d_caliper`). Accuracy degrades toward the suppressed class for the
over-growing methods because a sub-canopy stem's true crown is small but a
region-grower still claims a full canopy patch; the intermediate/suppressed bins
are also thin (13 and 2 trees) so treat those cells as indicative, not
definitive.

---

## Results — by site

The structure gradient SJER → SOAP → TEAK behaves as expected: the sparse oak
savanna (SJER) is hardest (few matched trees, all methods biased low into the
small open crowns), while the closed conifer canopies (SOAP, TEAK) bring
`region_growing` closest to neutral R² (SOAP `d_eq` R² ≈ 0). On the shared seed
set no method reaches a clearly positive per-site R², but `region_growing` stays
nearest zero and least biased across the gradient.

### `d_eq` vs `ninetyCrownDiameter` (RMSE m / bias m / R²)

| Algorithm | SJER (n=22) | SOAP (n=87) | TEAK (n=116) |
|-----------|---|---|---|
| lasr_region_growing | 3.54 / −1.92 / −0.53 | **2.78 / +1.15 / −0.00** | **2.26 / +1.65 / −0.11** |
| dalponte2016 | 3.54 / −1.75 / −0.52 | 3.00 / +1.36 / −0.17 | **2.23 / +1.56 / −0.08** |
| silva2016 | 3.48 / −1.74 / −0.48 | 3.00 / +1.33 / −0.17 | 2.46 / +1.93 / −0.31 |
| random_walker | 3.47 / −1.40 / −0.47 | 3.36 / +1.93 / −0.46 | 3.00 / +2.42 / −0.96 |
| watershed (m-free) | 4.65 / +2.11 / −1.64 | 3.64 / +2.41 / −0.72 | 2.93 / +2.27 / −0.86 |

### `d_caliper` vs `maxCrownDiameter` (RMSE m / bias m / R²)

| Algorithm | SJER (n=22) | SOAP (n=87) | TEAK (n=116) |
|-----------|---|---|---|
| lasr_region_growing | **3.71 / −0.95 / −0.10** | **4.02 / +1.79 / −0.12** | **3.48 / +2.80 / −0.93** |
| dalponte2016 | 4.15 / −0.08 / −0.38 | 5.02 / +2.91 / −0.75 | 3.99 / +3.13 / −1.53 |
| silva2016 | 3.86 / −0.55 / −0.19 | 4.60 / +2.62 / −0.47 | 4.34 / +3.77 / −1.99 |
| random_walker | 4.21 / −0.19 / −0.41 | 5.02 / +3.06 / −0.75 | 4.52 / +3.85 / −2.25 |
| watershed (m-free) | 7.08 / +4.62 / −3.00 | 6.32 / +4.36 / −1.77 | 4.64 / +4.01 / −2.42 |

---

## The random-walker outcome (honest notes)

**It worked, and it was not unstable or slow.** Per-plot CHMs are small
(~20k–33k pixels), so the sparse 4-neighbour Laplacian solve (`Matrix::solve`)
completes in well under a second; the whole 3-site run (all five segmenters over
40 plots) finished in 1.6 min and the random walker never hit its 120 s
timebox. A small ε on the Laplacian diagonal keeps the Dirichlet system
numerically solvable, **but solvable is not the same as correctly labelled**: a
canopy component that no seed reaches gets ~0 marker probability for every label,
and a plain `argmax` would force those pixels into label 1, inflating crown 1's
area. They are now dropped to **background** (NA) instead. The effect on the
pooled metrics is small (seed-less islands are rare in these dense-seed plots),
but the labelling is now correct rather than relying on "always solvable".

**But it is not the right growing rule for crown diameter here.** Like a
marker-controlled watershed, the random walker partitions the *entire* canopy
mask among the seeds (every canopy pixel is argmax-assigned to some marker),
producing a complete Voronoi-on-the-CHM tiling with no "stop growing" criterion.
On NEON's open-to-moderate canopies that systematically over-grows crowns into
the gaps between stems, inflating diameter (pooled bias +1.85 m on `d_eq`,
+3.15 m on `d_caliper`) and giving negative R². It lands close to the
marker-free watershed, which is the expected behaviour for a label-everything
partition. A height-percentile cutoff per crown (analogous to lidR's `th_cr`)
or a no-grow band of low-probability pixels would likely close most of the gap,
but as specified — pure argmax labelling — it is an over-segmenter for this
metric. Kept in the benchmark as an honest negative result.

---

## Geometric caveat — equivalent-circle vs max-axis

The two diameter definitions are **not interchangeable** and the gap is
systematic:

- `d_eq = 2√(area/π)` is the diameter of a circle with the crown's area. It is
  the natural match to `ninetyCrownDiameter` (a width-like measure) and is
  *insensitive to crown shape*.
- `d_caliper` is the crown polygon's widest chord (max pairwise vertex
  distance). It is the natural match to `maxCrownDiameter`, which NEON defines
  as the crown's **widest axis**, but it is *biased high* for any non-convex or
  elongated polygon and amplifies polygonization jaggedness.

Empirically, for the *same* crowns, `d_caliper`-vs-`maxCD` RMSE (3.72–5.61 m
pooled) runs ~1.4–1.6× the `d_eq`-vs-`ninetyCD` RMSE (2.62–3.42 m), and field
`maxCD` (mean 5.97 m) exceeds field `ninetyCD` (mean 4.94 m) by ~1 m — the same
"widest axis > equivalent width" relationship. **Compare like with like:** use
`d_eq` against `ninetyCrownDiameter` and `d_caliper` against `maxCrownDiameter`,
never cross them.

---

## Deep model: TreeisoNet treeOff crowns (SOAP, issues #M7 / #20)

The crown analogue of the TreeisoNet detection arm: its offset net (`treeOff`)
is run zero-shot on the same native frozen SOAP clips
(`gpu/run_treeisonet_crowns.py`), the per-point instances reduced to apexes +
convex-hull crown diameters (the bridge's `reduce_instances` +
`crown_diameter_table` on canopy points >= 2 m), matched to field stems, and
scored against field crown diameter exactly like the five CHM segmenters. The
comparison is **SOAP-only** (TreeisoNet's GPU clips are SOAP), so the five
classical arms are re-pooled on their SOAP rows here — their numbers differ
slightly from the multi-site tables above.

### Crown-diameter accuracy on SOAP (pooled matched trees)

**Equivalent-circle d_eq vs ninetyCrownDiameter**

| algo | n | rmse | mae | bias | r2 |
| --- | --- | --- | --- | --- | --- |
| lasr_region_growing | 87 | 2.78 | 2.08 | +1.15 | -0.004 |
| dalponte2016 | 87 | 3.00 | 2.27 | +1.36 | -0.167 |
| silva2016 | 87 | 3.00 | 2.31 | +1.33 | -0.169 |
| random_walker | 87 | 3.36 | 2.66 | +1.93 | -0.461 |
| watershed_markerfree | 87 | 3.64 | 2.95 | +2.41 | -0.720 |
| treeisonet | 96 | 8.97 | 4.72 | +2.36 | -10.197 |

**Max-caliper d_caliper vs maxCrownDiameter**

| algo | n | rmse | mae | bias | r2 |
| --- | --- | --- | --- | --- | --- |
| lasr_region_growing | 87 | 4.02 | 3.05 | +1.79 | -0.119 |
| silva2016 | 87 | 4.60 | 3.74 | +2.62 | -0.467 |
| dalponte2016 | 87 | 5.02 | 3.96 | +2.91 | -0.747 |
| random_walker | 87 | 5.02 | 4.06 | +3.06 | -0.749 |
| watershed_markerfree | 87 | 6.32 | 5.10 | +4.36 | -1.769 |
| treeisonet | 96 | 11.98 | 6.46 | +3.49 | -10.854 |

**TreeisoNet treeOff crowns are far worse than every classical segmenter** —
RMSE 8.97 m (`d_eq` vs ninetyCD) and 11.98 m (`d_caliper` vs maxCD), ~2.5–3x the
best CHM arm (2.78 / 4.02 m), with a strongly negative R² (about −10). But the
gap is a **heavy tail, not a uniform shift**: MAE (4.72 / 6.46 m) is far below
RMSE, and the median crown is reasonable (`d_eq` median ≈ 6 m vs field ≈ 5 m).
The tail comes from `treeOff`'s nearest-seed assignment lumping spatially distant
points into a few crowns zero-shot, inflating their hulls to tens of metres. So
even for the trees it detects, the deep model delineates crowns unreliably on
sparse NEON ALS — mirroring its detection collapse.

Caveats specific to this arm: (a) geometry is the **convex hull of canopy
points**, not the dissolved CHM polygon the classical arms use — not identical
estimators, though both target the same field diameter; (b) seeds are the
detection arm's `conf = 0.22` tops (few seeds → some under-segmentation, which
feeds the large-crown tail); (c) only TreeisoNet's *detected* trees contribute,
and its detection is poor, so n is limited to its matches.

---

## Deep instance segmenters: SegmentAnyTree, ForestFormer3D (issue #34)

The two deep instance-segmentation **detection** arms — SegmentAnyTree
([#M6](https://github.com/agrigoriev/lidar_tree_benchmarks/issues/17)) and
ForestFormer3D
([#M8](https://github.com/agrigoriev/lidar_tree_benchmarks/issues/18)) — already
produce per-point instance labels on the **same** frozen clips the detection
benchmark scored, but were never crown-scored. This section closes that gap by
piping each model's instance labelling through the **same** #7/#30 crown-scoring
harness (`crown_diameter_table` → `greedy_match` to field stems → `d_eq` /
`d_caliper` vs the NEON columns), so they sit head-to-head with the classical CHM
arms and the TreeisoNet negative result. Driver:
[`scripts/crown_metrics_deepmodel.R`](../scripts/crown_metrics_deepmodel.R).

The comparison **starts SOAP-native** — the site/rung where the GPU detection
arms run today (SAT covers the SOAP density ladder; FF3D runs SOAP native + 8).
The detection arms now persist their merged per-point instance cloud under
`work/neon/<SITE>/{segmentanytree,forestformer3d}_instances/<plot>_<rung>.laz`,
which is what this crown arm consumes; that persistence hook was added with this
arm, so the SOAP instance clouds must be **regenerated** by re-running the
detection arms before the crown CSVs can be produced (the earlier SOAP detection
runs wrote only the `*_results.csv` apex/score rows and discarded the per-point
cloud to tempdir). The arm extends to SJER + TEAK automatically once those sites'
GPU instance clouds exist; a requested site/model with no persisted instance
cloud is **skipped with a message, never fabricated**.

### Method

Both arms consume **persisted** per-point instance clouds (the script never runs
a GPU container) and reuse the #30 glue verbatim:

- **SegmentAnyTree** — the merged per-point LAS carries the instance label as
  the `PredInstance` extra dim (`0` = non-tree, `1..N`; see
  [`scripts/detect_segmentanytree_sweep.R`](../scripts/detect_segmentanytree_sweep.R)).
  The full labelled table is read via `read_instance_points_laz(PredInstance)`
  (the new full-table reader in
  [`scripts/io_bridge.R`](../scripts/io_bridge.R), `0` → NA so it is dropped),
  then `crown_diameter_table()` + `instance_apex()` come from that **same**
  labelling. `algo = "segmentanytree"`.
- **ForestFormer3D** — the merged per-cylinder labelled LAZ (`UserData` = block,
  `PointSourceID` = per-cylinder instance id; see
  [`scripts/detect_forestformer3d_sweep.R`](../scripts/detect_forestformer3d_sweep.R))
  is stacked into `(block, inst, X, Y, Z)` and passed through the new
  `ff3d_crown_table()` helper in
  [`scripts/model_bench_lib.R`](../scripts/model_bench_lib.R): it runs
  `dedup_blocks()` (the #M8 cross-block apex-cluster merge — which merges only
  duplicate detections across overlapping cylinders and **never** launders the
  model's within-cylinder over-segmentation), then derives
  `crown_diameter_table(id_col = "global_id")` + `instance_apex(id_col =
  "global_id")` so diameters are **per tree**, not per cylinder.
  `algo = "forestformer3d"`.

Both clouds keep **absolute UTM Z**, but the crown diameter (a horizontal X/Y
convex hull) is invariant to the Z datum; only the apex z feeds the matching
height gate, where the field heights are AGL. Each instance apex is therefore
converted to AGL via the cached frozen clip's `ground_dtm.tif` (`det_to_agl`,
the same transform the detection arms use) before `greedy_match`. Crown diameter
is scored at **one rung per plot** (native by default) so the pooled table is
one row per matched tree, exactly as #30 — mixing rungs would double-count a
stem. Pooling is by **summed** squared errors (RMSE/MAE/bias/R²), never a mean of
per-plot rates.

Estimator caveat (same as #30): this is a **convex hull of the instance's
points**, not the dissolved-CHM polygon the #7 classical arms use — not
identical estimators, though both target the same field column. Compare `d_eq`
only against `ninetyCrownDiameter` and `d_caliper` only against
`maxCrownDiameter`.

### Results — deep instance segmenters (SOAP-native)

*Results pending regeneration (run the command above on a data-equipped machine).*

The numbers require the persisted GPU instance clouds under
`work/neon/SOAP/{segmentanytree,forestformer3d}_instances/`. Those clouds are
written by the detection arms' persistence hook (added with this crown arm), so
they must be **regenerated** by re-running the SOAP detection arms on a
GPU-equipped machine — the prior SOAP runs predate the hook and left only the
`*_results.csv` rows. SJER/TEAK remain pending GPU clips. Once the clouds exist,
the crown CSVs are produced by:

```sh
# (re-run the detection arms first so the instance clouds are persisted)
Rscript scripts/detect_segmentanytree_sweep.R SITE=SOAP
Rscript scripts/detect_forestformer3d_sweep.R SITE=SOAP
# then crown-score the persisted clouds:
Rscript scripts/crown_metrics_deepmodel.R SITE=SOAP
Rscript scripts/analyze_crown_metrics.R   SITES=SOAP
```

The expectation to test on regeneration: like TreeisoNet, these deep arms are
zero-shot on sparse NEON ALS, so a heavy-tailed diameter error (RMSE ≫ MAE) and
a strongly negative R² versus the classical CHM arms (`lasr_region_growing`,
`dalponte2016`) would be the documented-negative-result pattern; the union into
[`scripts/analyze_crown_metrics.R`](../scripts/analyze_crown_metrics.R) re-pools
the classical and TreeisoNet arms on the same rows for the head-to-head. Do not
assume the outcome — record what regeneration produces.

---

## 3-D instance segmenters: Li 2012, ptrees, AMS3D (issue #30)

The five arms above all delineate crowns on a **CHM**. The model-benchmark
**detection** arms — lidR Li 2012, lidRplugins ptrees (Vega 2014), and AMS3D
(`crownsegmentr` adaptive mean shift) — already run on the same native frozen
normalized clips but were only ever scored for *detection* (apex recall), never
for crown diameter. This section closes that gap by crown-scoring those 3-D
point-instance segmenters head-to-head with the CHM controls.

### Method

For each site, the SAME plot set as #7 (plots with ≥6 live mapped trees carrying
a non-NA NEON field crown diameter) is processed at **native density** from the
cached frozen normalized clip (`frozen_clip(rung=NA)`; clips are reused, never
regenerated). Per plot, each segmenter is run **once** on the normalized clip,
reusing the verbatim invocation from its detection arm
([`detect_li2012_native.R`](../scripts/detect_li2012_native.R),
[`detect_lidrplugins_sweep.R`](../scripts/detect_lidrplugins_sweep.R),
[`detect_ams3d_sweep.R`](../scripts/detect_ams3d_sweep.R)). From the resulting
per-point instance labels (`seg@data` `treeID` / `crown_id`) two artefacts are
derived from the *same* labelling:

1. **Crown diameter** via `crown_diameter_table()` (`model_bench_lib.R`,
   `min_pts=5`): per instance the **2-D convex hull** of its points gives
   `d_eq = 2√(hull_area/π)` (→ `ninetyCrownDiameter`) and
   `d_caliper = max pairwise point distance` (→ `maxCrownDiameter`). Both are NA
   for instances with fewer than 5 points (a diameter over a few points is
   noise). `area` is the equivalent-circle area implied by `d_eq`.
2. **Instance apex** (`instance_apex()`): the max-Z point of each instance,
   keyed by instance id so the diameter row can be re-joined after matching.

Each instance apex is matched to a field stem with the **same** `greedy_match`
harness as everywhere else (global nearest-distance 1:1, position tol 4 m + the
height-consistency gate), then the matched stem's diameter row and field crown
diameter are joined. The match → diameter → field-join glue is the pure helper
`score_crowns_against_field()` in
[`model_bench_lib.R`](../scripts/model_bench_lib.R), unit-tested with synthetic
instances + stems. Pooling is by **summed** squared errors over matched trees
(RMSE/MAE/bias/R²), never a mean of per-plot rates, exactly as #7.

**Estimator caveat (important for the head-to-head).** This arm's diameter is a
**convex hull of the instance's points**, not the dissolved-CHM polygon the #7
classical arms use. The two are not identical estimators (a point hull is
sensitive to outlier returns; a CHM polygon to pixel resolution and pit-fill),
though both target the same field column. The same geometric caveat as #7
applies on top: `d_caliper` (widest axis) is biased high relative to `d_eq`
(equivalent-circle), so compare `d_eq` only against `ninetyCrownDiameter` and
`d_caliper` only against `maxCrownDiameter`.

### Results — 3-D segmenters (pooled, all sites)

*Results pending regeneration (run the command above on a data-equipped machine).*

### Head-to-head: 3-D vs the #7 CHM controls (same plot set)

The #7 CHM controls (`dalponte2016`, `lasr_region_growing`) are **re-pooled on
the same plot set** scored by the 3-D arms for a like-for-like comparison;
[`scripts/analyze_crown_metrics.R`](../scripts/analyze_crown_metrics.R) unions
`crown_metrics_results.csv` (CHM arms) with `crown_metrics_3d_results.csv` (3-D
arms) and pools per algorithm across the requested sites.

*Results pending regeneration (run the command above on a data-equipped machine).*

---

## Reproduce

```sh
export CLAUDE_JOB_DIR=$(pwd)/work
Rscript scripts/crown_metrics_sweep.R SITES=SJER,SOAP,TEAK CORES=4 \
        TOL=4 RES=0.5 A=0.10
# -> work/neon/<SITE>/crown_metrics_results.csv  (one row per matched tree:
#    site, plot, algo, crown_class, individualID, d_eq, d_caliper, area,
#    field_maxCD, field_ninetyCD), plus the pooled RMSE tables on stdout.

# 3-D instance-segmenter crown arm (Li 2012 / ptrees / AMS3D; issue #30).
# Native clips at SJER+SOAP+TEAK; reuses the cached frozen clips + ground truth
# (no LiDAR re-download). Writes the NEW crown_metrics_3d_results.csv:
Rscript scripts/crown_metrics_3d.R SITES=SJER,SOAP,TEAK CORES=8 TOL=4
# -> work/neon/<SITE>/crown_metrics_3d_results.csv
# Then union the CHM (+3-D, +TreeisoNet) arms and pool per algorithm:
Rscript scripts/analyze_crown_metrics.R SITES=SJER,SOAP,TEAK
# -> work/neon/crown_compare_tables.md

# TreeisoNet treeOff crown arm (SOAP, GPU; needs gpu/setup_treeisonet_env.sh):
Rscript scripts/detect_treeisonet_crowns.R SITE=SOAP PLOTS=ALL CONF=0.22
Rscript scripts/analyze_crown_metrics.R    SITE=SOAP   # union + SOAP RMSE table

# Deep instance-segmenter crown arm (SegmentAnyTree #M6 + ForestFormer3D #M8;
# issue #34). Reads PERSISTED per-point instance clouds (no GPU container is run
# here) under work/neon/<SITE>/{segmentanytree,forestformer3d}_instances/ and
# reuses the cached frozen DTMs. Those clouds are written by the detection arms'
# persistence hook, so re-run the SOAP detection arms first (they predate the
# hook); SJER/TEAK once GPU clips exist. Writes the NEW per-model
# segmentanytree_/forestformer3d_crown_metrics.csv:
Rscript scripts/detect_segmentanytree_sweep.R SITE=SOAP   # persists instance clouds
Rscript scripts/detect_forestformer3d_sweep.R SITE=SOAP   # persists instance clouds
Rscript scripts/crown_metrics_deepmodel.R SITE=SOAP
Rscript scripts/analyze_crown_metrics.R   SITES=SOAP  # union + SOAP RMSE table
```

Requires lidR + lasR (`pre-devel`), terra, sf, data.table, **Matrix** (random
walker), and EBImage (lidR `watershed`). All generated CSVs are gitignored;
regenerate with the command above.
