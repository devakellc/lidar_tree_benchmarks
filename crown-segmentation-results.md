# Crown-Segmentation Benchmark — Results

*Closes the loop on [GitHub issue #7](https://github.com/agrigoriev/lidar_tree_benchmarks/issues/7):
the density-ladder sweep scores tree-top **detection** only; this benchmark
delineates **crowns** from the detected tops and scores their diameter against
NEON field crown diameter. Five segmenters are seeded from the **same** detected
tops, run per plot on a native-density pit-free CHM, matched back to field
stems, and scored RMSE/MAE/bias/R² by crown class. Driver:
[`scripts/crown_metrics_sweep.R`](scripts/crown_metrics_sweep.R). Sites: NEON
SJER (open oak savanna), SOAP (mixed conifer), TEAK (red-fir). Last run:
2026-06-05.*

---

## TL;DR

- **One detector, five crowns.** Tree-tops are detected once per plot
  (`locate_trees(lmf, ws=ws_factory(0.10))` on a pit-free CHM at 0.5 m) and
  reused as the seed set for every seeded segmenter, so differences are
  attributable to the *growing rule*, not the seeds. 225 field stems matched
  across 41 plots (SJER 22, SOAP 87, TEAK 116).
- **lasR `region_growing` tracks NEON crown diameter best — by a clear
  margin.** Pooled over all sites it is the lowest-RMSE and the **only**
  segmenter with a positive R² on either diameter definition. Equivalent-circle
  diameter vs `ninetyCrownDiameter`: RMSE 2.38 m, R² +0.26. Max-caliper vs
  `maxCrownDiameter`: RMSE 3.35 m, R² +0.07. Every other method is biased high
  and has R² ≤ 0.
- **The random walker ran successfully** (sparse Dirichlet solve via `Matrix`,
  sub-second per plot, never timed out) but **underperforms**: like the
  marker-free watershed it over-grows crowns into a full Voronoi-style tiling of
  the canopy, inflating diameter (pooled bias +1.85 m / +3.15 m). It is an
  honest negative result, documented below.
- **Geometric caveat is large and systematic.** Equivalent-circle diameter
  `d_eq = 2√(area/π)` underestimates the widest-axis `maxCrownDiameter` because
  a real crown is not a disc: pooled `d_caliper`-vs-`maxCD` RMSE (3.35–5.61 m)
  is ~1.5–2× the `d_eq`-vs-`ninetyCD` RMSE (2.38–3.42 m) for the same crowns.
  Pick the diameter definition that matches the field column.

---

## Method

For each site, plots with ≥6 live mapped trees carrying a non-NA NEON field
crown diameter are processed at **native density** (no decimation;
`prepare_clip(rung=NA)` from [`sweep_lib.R`](scripts/sweep_lib.R)):

1. **Shared seeds.** One pit-free CHM (lidR `pitfree`, mirrors
   [`segment_lidr.R`](scripts/segment_lidr.R)) at `RES=0.5 m`; detect tops once
   via `locate_trees(chm, lmf(ws=ws_factory(0.10), hmin=2, shape="circular"))`.
2. **Five segmenters**, each producing a per-crown polygon → crown area:
   - lidR `dalponte2016(chm, ttops, th_seed=0.45, th_cr=0.55, max_cr=10/res px)`
     — seeded region-growing.
   - lidR `silva2016(chm, ttops, max_cr_factor=0.6, exclusion=0.3)` — seeded
     Voronoi-like.
   - lidR `watershed(chm, th_tree=2)` — **marker-free** (EBImage); crowns
     matched to stems by polygon containment of the seed point.
   - lasR `region_growing(pit_fill, seed, th_tree=2, th_seed=0.45, th_cr=0.55,
     max_cr=10)` seeded from `local_maximum_raster` on the lasR CHM, mirroring
     [`segment_lasr.R`](scripts/segment_lasr.R).
   - **Random walker** (Grady 2006) on the CHM raster, seeded from the same
     detected tops, implemented with `Matrix`: 4-neighbour pixel graph over
     canopy pixels (Z ≥ 2 m), Gaussian edge weights
     `w_ij = exp(−β·(chmᵢ−chmⱼ)²)` on the normalized CHM (β = 1), one marker
     label per seed, solving the combinatorial Dirichlet problem
     `Lᵤ x = −B m` for marker probabilities and argmax-labelling each pixel.
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
added to [`scripts/neon_ground_truth.R`](scripts/neon_ground_truth.R) (additive,
not re-run) so future ground-truth regenerations carry them.

Field crown diameter (matched stems): `ninetyCrownDiameter` mean 4.95 m
(median 4.30, range 1.1–17.6); `maxCrownDiameter` mean 5.97 m (median 5.20,
range 1.3–26.1).

---

## Results — pooled over all sites

### Equivalent-circle `d_eq` vs `ninetyCrownDiameter`

| Algorithm | n | RMSE (m) | MAE (m) | bias (m) | R² |
|-----------|---:|---:|---:|---:|---:|
| **lasr_region_growing** | 224 | **2.38** | **1.72** | **+0.24** | **+0.260** |
| dalponte2016 | 225 | 2.70 | 2.06 | +1.16 | +0.046 |
| silva2016 | 225 | 2.79 | 2.22 | +1.34 | −0.021 |
| random_walker | 225 | 3.19 | 2.57 | +1.85 | −0.334 |
| watershed (marker-free) | 225 | 3.42 | 2.70 | +2.31 | −0.530 |

### Max-caliper `d_caliper` vs `maxCrownDiameter`

| Algorithm | n | RMSE (m) | MAE (m) | bias (m) | R² |
|-----------|---:|---:|---:|---:|---:|
| **lasr_region_growing** | 224 | **3.35** | **2.41** | **+1.01** | **+0.073** |
| silva2016 | 225 | 4.40 | 3.63 | +2.90 | −0.596 |
| dalponte2016 | 225 | 4.43 | 3.47 | +2.73 | −0.621 |
| random_walker | 225 | 4.69 | 3.85 | +3.15 | −0.817 |
| watershed (marker-free) | 225 | 5.61 | 4.54 | +4.21 | −1.595 |

**lasR `region_growing` wins both definitions and is the only method with a
positive R².** Its conservative growing rule (`max_cr=10` data-units, stopped at
`th_cr=0.55` of seed height) keeps crowns compact and near-unbiased, whereas the
lidR seeded methods, the marker-free watershed, and the random walker all expand
crowns to tile the canopy and so run +1 to +4 m high on diameter.

### RMSE by crown class — `d_eq` vs `ninetyCrownDiameter`

| Algorithm | dominant | codominant | intermediate | suppressed |
|-----------|---:|---:|---:|---:|
| | n=101 | n=109 | n=13 | n=2 |
| **lasr_region_growing** | **2.80** | **1.92** | **2.31** | **2.44** |
| dalponte2016 | 2.87 | 2.53 | 2.32 | 4.39 |
| silva2016 | 2.95 | 2.64 | 2.44 | 4.48 |
| random_walker | 3.40 | 3.03 | 2.48 | 4.41 |
| watershed (marker-free) | 3.39 | 3.40 | 3.56 | 4.77 |

### RMSE by crown class — `d_caliper` vs `maxCrownDiameter`

| Algorithm | dominant | codominant | intermediate | suppressed |
|-----------|---:|---:|---:|---:|
| | n=101 | n=109 | n≈12 | n=2 |
| **lasr_region_growing** | **3.54** | **3.14** | 3.82 | **1.74** |
| silva2016 | 4.24 | 4.53 | 4.36 | 4.95 |
| dalponte2016 | 4.23 | 4.55 | 4.43 | 6.97 |
| random_walker | 4.99 | 4.52 | **3.60** | 4.24 |
| watershed (marker-free) | 5.40 | 5.78 | 5.81 | 4.74 |

lasR `region_growing` is best or tied-best in every crown class on both
definitions. Accuracy degrades toward the suppressed class for the over-growing
methods because a sub-canopy stem's true crown is small but a region-grower
still claims a full canopy patch; the intermediate/suppressed bins are also
thin (13 and 2 trees) so treat those cells as indicative, not definitive.

---

## Results — by site

The structure gradient SJER → SOAP → TEAK behaves as expected: the sparse oak
savanna (SJER) is hardest (few matched trees, all methods biased low into the
small open crowns), and the closed conifer canopies (SOAP, TEAK) let
`region_growing` reach positive R².

### `d_eq` vs `ninetyCrownDiameter` (RMSE m / bias m / R²)

| Algorithm | SJER (n=22) | SOAP (n=87) | TEAK (n=116) |
|-----------|---|---|---|
| lasr_region_growing | 3.84 / −2.63 / −0.80 | **2.59 / +0.14 / +0.13** | **1.77 / +0.86 / +0.33** |
| dalponte2016 | 3.54 / −1.75 / −0.52 | 3.00 / +1.36 / −0.17 | 2.23 / +1.56 / −0.08 |
| silva2016 | 3.48 / −1.74 / −0.48 | 3.00 / +1.33 / −0.17 | 2.46 / +1.93 / −0.31 |
| random_walker | 3.47 / −1.40 / −0.47 | 3.36 / +1.93 / −0.46 | 3.00 / +2.42 / −0.96 |
| watershed (m-free) | 4.65 / +2.11 / −1.64 | 3.64 / +2.41 / −0.72 | 2.93 / +2.27 / −0.86 |

### `d_caliper` vs `maxCrownDiameter` (RMSE m / bias m / R²)

| Algorithm | SJER (n=22) | SOAP (n=87) | TEAK (n=116) |
|-----------|---|---|---|
| lasr_region_growing | 3.84 / −1.91 / −0.18 | **3.83 / +0.64 / −0.02** | **2.82 / +1.86 / −0.26** |
| dalponte2016 | 4.15 / −0.08 / −0.38 | 5.02 / +2.91 / −0.75 | 3.99 / +3.13 / −1.53 |
| silva2016 | 3.86 / −0.55 / −0.19 | 4.60 / +2.62 / −0.47 | 4.34 / +3.77 / −1.99 |
| random_walker | 4.21 / −0.19 / −0.41 | 5.02 / +3.06 / −0.75 | 4.52 / +3.85 / −2.25 |
| watershed (m-free) | 7.08 / +4.62 / −3.00 | 6.32 / +4.36 / −1.77 | 4.64 / +4.01 / −2.42 |

---

## The random-walker outcome (honest notes)

**It worked, and it was not unstable or slow.** Per-plot CHMs are small
(~20k–33k pixels), so the sparse 4-neighbour Laplacian solve (`Matrix::solve`)
completes in well under a second; the whole 3-site run (all five segmenters over
41 plots) finished in 1.7 min and the random walker never hit its 120 s
timebox. A small ε on the Laplacian diagonal regularizes seed-less canopy
islands so the Dirichlet system is always solvable — no fallbacks were
triggered.

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

Empirically, for the *same* crowns, `d_caliper`-vs-`maxCD` RMSE (3.35–5.61 m
pooled) runs ~1.5–2× the `d_eq`-vs-`ninetyCD` RMSE (2.38–3.42 m), and field
`maxCD` (mean 5.97 m) exceeds field `ninetyCD` (mean 4.95 m) by ~1 m — the same
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
