# Crown-Segmentation Benchmark — Results

*Closes the loop on [GitHub issue #7](https://github.com/agrigoriev/lidar_tree_benchmarks/issues/7):
the density-ladder sweep scores tree-top **detection** only; this benchmark
delineates **crowns** from the detected tops and scores their diameter against
NEON field crown diameter. Five segmenters are seeded from the **same** detected
tops, run per plot on a native-density pit-free CHM, matched back to field
stems, and scored RMSE/MAE/bias/R² by crown class. Driver:
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

## Reproduce

```sh
export CLAUDE_JOB_DIR=$(pwd)/work
Rscript scripts/crown_metrics_sweep.R SITES=SJER,SOAP,TEAK CORES=4 \
        TOL=4 RES=0.5 A=0.10
# -> work/neon/<SITE>/crown_metrics_results.csv  (one row per matched tree:
#    site, plot, algo, crown_class, individualID, d_eq, d_caliper, area,
#    field_maxCD, field_ninetyCD), plus the pooled RMSE tables on stdout.
```

Requires lidR + lasR (`pre-devel`), terra, sf, data.table, **Matrix** (random
walker), and EBImage (lidR `watershed`). All generated CSVs are gitignored;
regenerate with the command above.
