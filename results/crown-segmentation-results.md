# Crown-Segmentation Benchmark — Results

*Closes the loop on [GitHub issue #7](https://github.com/agrigoriev/lidar_tree_benchmarks/issues/7):
the density-ladder sweep scores tree-top **detection** only; this benchmark
delineates **crowns** from the detected tops and scores their diameter against
NEON field crown diameter. Several segmenters are seeded from the **same**
detected tops, run per plot on a native-density pit-free CHM, matched back to
field stems, and scored RMSE/MAE/bias/R² by crown class (the original #7 tables
score five; issues #35 and #32 add the `random_walker_thcr` and `watershed_seeded`
arms, and issue #31 adds a multichm-seeded variant of the lidR segmenters, in
dedicated sections below). A SOAP-only TreeisoNet
`treeOff` crown arm is unioned by
[`scripts/analyze_crown_metrics.R`](../scripts/analyze_crown_metrics.R) into the
deep-model section below. Driver:
[`scripts/crown_metrics_sweep.R`](../scripts/crown_metrics_sweep.R). Sites: NEON
SJER (open oak savanna), SOAP (mixed conifer), TEAK (red-fir). Last run:
2026-06-05.*

---

## TL;DR

- **One detector, five crowns, on a shared CHM.** Tree-tops are detected once
  per plot (`locate_trees(lmf, ws=ws_factory(0.10))` on a pit-free CHM at 0.5 m)
  and reused as the seed set for every seeded segmenter. lasR `region_growing`
  cannot take an external point set as seeds, so it is seeded by injecting the
  _same_ lidR pit-free CHM into the lasR pipeline (`load_raster`) and running
  `local_maximum_raster` with the _same_ `ws` — yielding a near-identical seed
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
  sub-second per plot, never timed out) but the pure-argmax arm **underperforms**:
  like the marker-free watershed it over-grows crowns into a full Voronoi-style
  tiling of the canopy, inflating diameter (pooled bias +1.85 m / +3.15 m).
  Seed-less canopy islands are dropped to background instead of being mislabelled
  into crown 1 (see below). **Issue #35 retests with a `th_cr` stop rule**
  (`random_walker_thcr`): each crown is truncated at `TH_CR=0.55` of its own seed
  apex height so it can no longer tile the inter-crown gaps. Both arms run in one
  pass; the before/after RMSE table is below (numbers pending a data-equipped
  regeneration run).
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

   Seed-source tagging (issue #31): the two lidR segmenters that accept an
   external `treetops` sf — `dalponte2016` and `silva2016` — are now tagged in the
   CSV with their seed source (`dalponte2016_seedlmf` / `silva2016_seedlmf` for
   the lmf-seeded control, `_seedmultichm` for the multichm-seeded variant; see
   the seed-sensitivity section below). The `_seedlmf` arms are the **identical**
   pipeline these historical tables score under the bare `dalponte2016` /
   `silva2016` names; the `lasr_region_growing`, `watershed_*`, and
   `random_walker*` arms keep their names (lmf-seeded only — they cannot consume
   the multichm tops).
   - lidR `watershed(chm, th_tree=2)` — **marker-free** (EBImage); crowns
     matched to stems by polygon containment of the seed point.
   - **Marker-controlled (seeded) watershed** (`watershed_seeded`, issue #32) —
     the same shared `ttops` are rasterized into a seed-label image aligned to
     the CHM (`seeds_to_marker_raster`, label k = the seed's `treeID`) and a
     priority-flood watershed grows **exactly one basin per marker** over the
     canopy mask (Z ≥ 2 m), so crowns are keyed by `treeID` (no containment
     rematch, no spurious extra basins). `imager::watershed(img, seeds)` does
     this priority-flood from a labelled seed image but imager is not installed
     in this environment, so `priority_flood_watershed()`
     ([`sweep_lib.R`](../scripts/sweep_lib.R)) runs the flood directly on the CHM
     (descending-height flood; seed-less canopy islands → background, like the
     random walker). This is the marker-**controlled** path the approach doc
     (sec 6) recommends, distinct from the marker-**free** EBImage arm above.
   - lasR `region_growing(pit_fill, seed, th_tree=2, th_seed=0.45, th_cr=0.55,
     max_cr=10)` (Dalponte growing rule, mirroring
     [`segment_lasr.R`](../scripts/segment_lasr.R)). lasR's `region_growing` takes
     a seed _stage_, not the shared `ttops` point set, so to grow from the same
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
   - **Random walker + `th_cr` stop rule** (`random_walker_thcr`, issue #35) —
     the _same_ solve, with each crown then truncated to background below
     `TH_CR=0.55` of its own seed apex height (see "The `th_cr` stop rule"
     below). Added so the before/after comparison runs on one identical seed set.
3. **Two diameter estimates per crown** (report both, see caveat):
   - `d_eq = 2√(area/π)` (equivalent-circle) → compared to `ninetyCrownDiameter`.
   - `d_caliper = max pairwise polygon-vertex distance` (max axis) → compared to
     `maxCrownDiameter` (NEON's widest-axis measurement).
4. **Matching.** The shared seeds are matched to field stems with
   `greedy_match` (global nearest-distance 1:1, position tol 4 m + height gate,
   from `sweep_lib.R`). Each matched stem's crown is the crown owning that seed
   (by `treeID` for the seeded lidR methods; by nearest seed-carrying crown for
   lasR/watershed/RW). The pair _(detected diameter, field diameter)_ feeds the
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
the _only_ method above zero there. On `d_caliper` it is best of the five
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

**But pure argmax is not the right growing rule for crown diameter here.** Like
a marker-controlled watershed, the random walker partitions the _entire_ canopy
mask among the seeds (every canopy pixel is argmax-assigned to some marker),
producing a complete Voronoi-on-the-CHM tiling with no "stop growing" criterion.
On NEON's open-to-moderate canopies that systematically over-grows crowns into
the gaps between stems, inflating diameter (pooled bias +1.85 m on `d_eq`,
+3.15 m on `d_caliper`) and giving negative R². It lands close to the
marker-free watershed, which is the expected behaviour for a label-everything
partition. The original `#7` doc noted that a height-percentile cutoff per crown
(analogous to lidR's `th_cr`) would likely close most of the gap; issue #35
retests exactly that.

### The `th_cr` stop rule (issue #35)

We add a second arm, **`random_walker_thcr`**, that reuses the _same_ Dirichlet
solve (one solve per plot, no extra cost) and applies a per-crown height cutoff,
approach (b) of the issue: after argmax labelling, a pixel assigned to crown k is
dropped to **background** when its CHM height falls below `TH_CR * seed_height_k`,
with `TH_CR = 0.55` to match dalponte2016's `th_cr`. Truncating each crown at a
fraction of _its own_ seed apex height stops it creeping down into the inter-crown
gaps that inflate diameter, without changing the marker probabilities themselves.
The cutoff is a pure relabel (`apply_thcr_cutoff()` in
[`sweep_lib.R`](../scripts/sweep_lib.R), unit-tested in
`tests/testthat/test-rw-thcr.R`): labels with no seed apex height are left
untouched, and the existing seed-less-island rule (rmax ≤ 1e-8 → background) is
unchanged, so the two rules compose. A crown truncated out of existence simply
contributes no diameter (its matched stem is skipped), which is the correct
behaviour for a stem whose canopy signal sits entirely below the cutoff.

#### Before/after pooled RMSE — random walker

The two RW arms run in one pass so the comparison is on identical seeds, CHM, and
matched stems; `lasr_region_growing` and `dalponte2016` are repeated from the
tables above as the region-growing controls. Regenerate with:

```sh
Rscript scripts/crown_metrics_sweep.R SITES=SJER,SOAP,TEAK CORES=1
```

_Results pending regeneration (run the command above on a data-equipped machine)._

**Equivalent-circle `d_eq` vs `ninetyCrownDiameter`**

| Algorithm | n | RMSE (m) | MAE (m) | bias (m) | R² |
|-----------|---:|---:|---:|---:|---:|
| random_walker (argmax) | _pending_ | _pending_ | _pending_ | _pending_ | _pending_ |
| random_walker_thcr | _pending_ | _pending_ | _pending_ | _pending_ | _pending_ |
| lasr_region_growing (control) | _pending_ | _pending_ | _pending_ | _pending_ | _pending_ |
| dalponte2016 (control) | _pending_ | _pending_ | _pending_ | _pending_ | _pending_ |

**Max-caliper `d_caliper` vs `maxCrownDiameter`**

| Algorithm | n | RMSE (m) | MAE (m) | bias (m) | R² |
|-----------|---:|---:|---:|---:|---:|
| random_walker (argmax) | _pending_ | _pending_ | _pending_ | _pending_ | _pending_ |
| random_walker_thcr | _pending_ | _pending_ | _pending_ | _pending_ | _pending_ |
| lasr_region_growing (control) | _pending_ | _pending_ | _pending_ | _pending_ | _pending_ |
| dalponte2016 (control) | _pending_ | _pending_ | _pending_ | _pending_ | _pending_ |

The pure-argmax `random_walker` row stays in the benchmark unchanged as the
documented honest negative; `random_walker_thcr` is the with-stop-rule retest.
Whether the cutoff reaches the region-growing controls is to be read off the
regenerated table above — do not infer it here.

---

## Marker-controlled (seeded) watershed (issue #32)

Issue #7's watershed arm is **marker-free**: lidR/EBImage `watershed(th_tree=2)`
finds its own basins from CHM minima, which on these open-to-moderate NEON
canopies systematically over-grows crowns (the largest pooled bias of the five
arms, **+2.31 m on `d_eq`**, +4.21 m on `d_caliper`). The approach doc
([`treetop-detection-approach.md`](../docs/treetop-detection-approach.md) sec 6)
recommends instead a **marker-controlled** watershed — treetops as basin
markers, exactly one crown per seed — which #7 never scored against NEON crown
diameter. This arm (`watershed_seeded`) closes that gap.

**How it differs from the marker-free arm.** The shared `ttops` (the same seed
set every other arm uses) are rasterized into a seed-label image aligned to the
CHM (`seeds_to_marker_raster()`: the seed's CHM cell carries its `treeID`; a
shared cell goes to the first seed; off-extent seeds are dropped). A
priority-flood watershed then grows one basin per marker over the canopy mask
(Z ≥ 2 m), flooding in descending CHM height so basins meet along the canopy
valleys between crowns. Because every basin is anchored to a marker, crowns are
keyed by `treeID` (matched like dalponte2016, `by="treeID"`) — no containment
rematch, and no spurious extra basins from spurious CHM minima. A canopy island
that no marker can reach over the mask stays **background**, the same
"floodable ≠ force-labelled" rule the random walker uses for seed-less islands.

**Library note.** `imager::watershed(img, seeds)` implements exactly this
priority-flood from a labelled seed image and would be the natural choice, but
`imager` is **not installed** in this environment (verified with
`requireNamespace`); EBImage's `watershed` is the marker-**free** algorithm
already used by the #7 arm, so it cannot serve here. The flood is therefore
implemented directly on the CHM in `priority_flood_watershed()`
([`sweep_lib.R`](../scripts/sweep_lib.R), unit-tested in
`tests/testthat/test-watershed-seeded.R`) — a 4-neighbour, highest-CHM-first
frontier in base R, pure (no I/O), so it needs no `work/` data to test.

### Pooled RMSE — seeded vs marker-free watershed and the region-growing controls

All arms run in one pass so the comparison is on identical seeds, CHM, and
matched stems. `dalponte2016` and `lasr_region_growing` are repeated from the
tables above as the region-growing controls; `watershed_markerfree` is #7's
documented negative. Regenerate with:

```sh
Rscript scripts/crown_metrics_sweep.R SITES=SJER,SOAP,TEAK CORES=1
```

_Results pending regeneration (run the command above on a data-equipped machine)._

**Equivalent-circle `d_eq` vs `ninetyCrownDiameter`**

| Algorithm | n | RMSE (m) | MAE (m) | bias (m) | R² |
|-----------|---:|---:|---:|---:|---:|
| watershed_seeded | _pending_ | _pending_ | _pending_ | _pending_ | _pending_ |
| watershed_markerfree | _pending_ | _pending_ | _pending_ | _pending_ | _pending_ |
| lasr_region_growing (control) | _pending_ | _pending_ | _pending_ | _pending_ | _pending_ |
| dalponte2016 (control) | _pending_ | _pending_ | _pending_ | _pending_ | _pending_ |

**Max-caliper `d_caliper` vs `maxCrownDiameter`**

| Algorithm | n | RMSE (m) | MAE (m) | bias (m) | R² |
|-----------|---:|---:|---:|---:|---:|
| watershed_seeded | _pending_ | _pending_ | _pending_ | _pending_ | _pending_ |
| watershed_markerfree | _pending_ | _pending_ | _pending_ | _pending_ | _pending_ |
| lasr_region_growing (control) | _pending_ | _pending_ | _pending_ | _pending_ | _pending_ |
| dalponte2016 (control) | _pending_ | _pending_ | _pending_ | _pending_ | _pending_ |

The expectation under test is that anchoring one basin per marker curbs the
marker-free arm's over-grow (it can no longer spawn a basin in every CHM
minimum), narrowing the bias toward the region-growing controls. Whether it
reaches them — and whether one-basin-per-seed under-segments where the marker-free
arm over-segments — is to be read off the regenerated table above; do not infer
it here.

---

## Crown-diameter seed sensitivity: multichm vs lmf tops (issue #31)

Issue #7 seeds every crown segmenter from one shared **lmf** (CHM-VWF) top set.
The model benchmark, however, shows `multichm` (Eysn-style multi-layer CHM local
maxima, `lidRplugins::multichm`) is the best classical **detector** on SOAP
(F1 0.44 vs CHM-VWF 0.38). Detection and segmentation are decoupled, so better
tops may or may not translate into better crown width. This arm tests that
directly: it re-seeds the lidR segmenters from `multichm` tops and scores crown
diameter with the **same** #7 harness, side-by-side against the lmf-seeded
control.

**What is re-seeded, and what is not.** Only the two lidR segmenters that accept
an external `treetops` sf — `dalponte2016` and `silva2016` — are run from both
seed sets, on the **same** per-plot pit-free CHM. The multichm tops are detected
on the point cloud (`multichm_seed_tops()` in
[`sweep_lib.R`](../scripts/sweep_lib.R)) at a density-derived resolution (0.25 m
when first-return density ≥ 8 pts/m², else 0.5 m) and the **same** clamped
variable window `ws_factory(0.10)` as the lmf seeds, so the only thing that
changes vs the control is the detector. `multichm` geometry is 2-D, so each top
reads its apex height from the same pit-free CHM the segmenters grow on (the
sf `Z` is a fallback only for an off-extent top). Each seed set is matched to
the field stems independently (`greedy_match`, position tol 4 m + height gate),
so the multichm arm is scored on the stems **its own** tops found.

**lasR `region_growing` cannot be re-seeded from multichm and stays lmf-seeded
(honest limitation).** Its API takes a seed _stage_ (a `local_maximum` producer),
not an external point set, and there is no matching lasR `local_maximum` that
emits the multichm tops, so the multichm tops cannot be injected as a lasR seed
stage without a parallel lasR detector. Re-seeding it from multichm is therefore
out of scope here; it (and `watershed_*`, `random_walker*`) remain the
lmf-seeded control. The crown tags encode the seed source
(`dalponte2016_seedlmf` / `dalponte2016_seedmultichm`, likewise `silva2016`) so
the canonical CSV column schema is **unchanged** — no new `seed` column — which
keeps [`analyze_crown_metrics.R`](../scripts/analyze_crown_metrics.R) and the
issue #33 density ladder (stacked on this PR) reading the same ten columns.

All arms run in one pass on identical plots and the same CHM. Each seed set is
scored on the stems ITS OWN tops matched, so n differs (multichm finds more
tops: 225 lmf vs 330 multichm — see the n columns below). Regenerate with:

```sh
Rscript scripts/crown_metrics_sweep.R SITES=SJER,SOAP,TEAK CORES=1
```

_Regenerated 2026-06-12 — SJER+SOAP+TEAK, native density, CORES=1._

### Pooled — multichm seeds minus lmf seeds (Δ = multichm − lmf)

Negative ΔRMSE / Δbias means multichm seeds improve crown-diameter accuracy.
The `_seedlmf` rows reproduce the historical `dalponte2016` / `silva2016` numbers
from the pooled tables above (identical pipeline).

**Equivalent-circle `d_eq` vs `ninetyCrownDiameter`**

| Segmenter | seed | n | RMSE (m) | bias (m) | R² | ΔRMSE | Δbias | ΔR² |
|-----------|------|---:|---:|---:|---:|---:|---:|---:|
| dalponte2016 | lmf | 225 | 2.70 | +1.16 | +0.046 | — | — | — |
| dalponte2016 | multichm | 330 | 2.47 | +0.31 | +0.121 | −0.22 | −0.85 | +0.075 |
| silva2016 | lmf | 225 | 2.79 | +1.34 | −0.021 | — | — | — |
| silva2016 | multichm | 330 | 2.44 | +0.40 | +0.144 | −0.35 | −0.94 | +0.165 |

**Max-caliper `d_caliper` vs `maxCrownDiameter`**

| Segmenter | seed | n | RMSE (m) | bias (m) | R² | ΔRMSE | Δbias | ΔR² |
|-----------|------|---:|---:|---:|---:|---:|---:|---:|
| dalponte2016 | lmf | 225 | 4.43 | +2.73 | −0.621 | — | — | — |
| dalponte2016 | multichm | 330 | 3.61 | +1.45 | −0.155 | −0.82 | −1.28 | +0.465 |
| silva2016 | lmf | 225 | 4.40 | +2.90 | −0.596 | — | — | — |
| silva2016 | multichm | 330 | 3.50 | +1.40 | −0.085 | −0.90 | −1.50 | +0.511 |

### Δ by crown class and site (`d_eq` vs `ninetyCrownDiameter`)

| Segmenter | scope | n (lmf/mc) | ΔRMSE | Δbias | ΔR² |
|-----------|-------|---:|---:|---:|---:|
| dalponte2016 | dominant | 101/123 | +0.01 | −1.23 | −0.107 |
| dalponte2016 | codominant | 109/175 | −0.33 | −0.71 | +0.507 |
| dalponte2016 | intermediate | 13/28 | −0.18 | −0.45 | −0.066 |
| dalponte2016 | suppressed | 2/4 | −1.76 | −4.81 | n/a† |
| dalponte2016 | SJER | 22/25 | −0.04 | −0.21 | +0.149 |
| dalponte2016 | SOAP | 87/139 | −0.29 | −1.22 | +0.215 |
| dalponte2016 | TEAK | 116/166 | −0.19 | −0.78 | +0.020 |
| silva2016 | dominant | 101/123 | −0.11 | −1.45 | −0.032 |
| silva2016 | codominant | 109/175 | −0.45 | −0.74 | +0.647 |
| silva2016 | intermediate | 13/28 | −0.25 | −0.28 | +0.003 |
| silva2016 | suppressed | 2/4 | −3.51 | −3.90 | n/a† |
| silva2016 | SJER | 22/25 | +0.02 | −0.19 | +0.094 |
| silva2016 | SOAP | 87/139 | −0.36 | −1.15 | +0.261 |
| silva2016 | TEAK | 116/166 | −0.43 | −0.99 | +0.261 |

† suppressed n is 2 (lmf) / 4 (multichm) trees — too few for a meaningful R²; the
ΔRMSE/Δbias entries are shown but should be read as indicative only.

The script prints this same ΔRMSE/Δbias/ΔR² breakdown (pooled + per crown class

- per site, both diameter definitions) at the end of a run
(`print_seed_sensitivity`).

### Recommended per-density-rung seed choice (consumed by issue #33)

Issue #33 (the density ladder, stacked on top of this PR) needs a per-rung
decision: at each density rung, seed the crown segmenters from the detector that
delineates crown diameter best. The rule to apply once the table above is
regenerated:

- **Default to the seed source with the lower pooled ΔRMSE at that rung's
  density**, read from the regenerated tables. Detection and segmentation are
  decoupled, so a detector that wins on F1 (multichm on SOAP) need **not** win on
  crown-diameter RMSE — the seed choice for #33 must be made on the
  crown-diameter Δ, not on detection F1.
- **Tie-break toward `lmf`** when |ΔRMSE| is within noise (the lmf seeds are the
  established #7 control and `lasr_region_growing`, the overall RMSE leader, is
  lmf-seeded anyway), so the ladder stays on one consistent seed family unless
  multichm shows a clear crown-width gain.

At **native density** the table makes the call: multichm seeds improve both
segmenters on both diameter definitions (`d_eq` ΔRMSE −0.22 to −0.35 m, Δbias
−0.85 to −0.94 m; larger gains on `d_caliper`), so #33 seeds from multichm at
native. Whether that holds as density drops — multichm needs enough returns to
resolve sub-canopy maxima — is decided per rung in #33's ladder, on the
crown-diameter Δ at each rung, not on detection F1.

---

## Geometric caveat — equivalent-circle vs max-axis

The two diameter definitions are **not interchangeable** and the gap is
systematic:

- `d_eq = 2√(area/π)` is the diameter of a circle with the crown's area. It is
  the natural match to `ninetyCrownDiameter` (a width-like measure) and is
  _insensitive to crown shape_.
- `d_caliper` is the crown polygon's widest chord (max pairwise vertex
  distance). It is the natural match to `maxCrownDiameter`, which NEON defines
  as the crown's **widest axis**, but it is _biased high_ for any non-convex or
  elongated polygon and amplifies polygonization jaggedness.

Empirically, for the _same_ crowns, `d_caliper`-vs-`maxCD` RMSE (3.72–5.61 m
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
feeds the large-crown tail); (c) only TreeisoNet's _detected_ trees contribute,
and its detection is poor, so n is limited to its matches.

---

## Reproduce

```sh
export CLAUDE_JOB_DIR=$(pwd)/work
# Both random-walker arms (argmax + th_cr) and the four CHM controls run in one
# pass. Use CORES=1: lasR exec under mclapply transiently fork-drops dense-native
# cells, so CORES=1 is the reproducible-pooling setting for the lasR arm.
Rscript scripts/crown_metrics_sweep.R SITES=SJER,SOAP,TEAK CORES=1
# (parameter defaults TOL=4 RES=0.5 A=0.10; add CORES=4 only when not relying on
#  the lasR_region_growing arm for exact pooling)
# -> work/neon/<SITE>/crown_metrics_results.csv  (one row per matched tree:
#    site, plot, algo, crown_class, individualID, d_eq, d_caliper, area,
#    field_maxCD, field_ninetyCD), plus the pooled RMSE tables on stdout.
# The dalponte2016/silva2016 arms emit BOTH _seedlmf (the #7 control) and
# _seedmultichm rows (issue #31); the run also prints the lmf-vs-multichm
# seed-sensitivity ΔRMSE/Δbias/ΔR² breakdown at the end.

# TreeisoNet treeOff crown arm (SOAP, GPU; needs gpu/setup_treeisonet_env.sh):
Rscript scripts/detect_treeisonet_crowns.R SITE=SOAP PLOTS=ALL CONF=0.22
Rscript scripts/analyze_crown_metrics.R    SITE=SOAP   # union + SOAP RMSE table
```

Requires lidR + lasR (`pre-devel`), **lidRplugins** (multichm seed arm), terra,
sf, data.table, **Matrix** (random walker), and EBImage (lidR `watershed`). All
generated CSVs are gitignored; regenerate with the command above.
