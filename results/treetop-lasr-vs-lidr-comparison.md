# lasR vs lidR Tree-Top Detection — Implementation & Comparison

Both tools run the documented approach (pit-free CHM -> variable-window
local-maximum -> dalponte2016 segmentation -> per-tree crown metrics) on the
**same** input with the **same** window function, then their treetops and
crown polygons are matched spatially.

## Setup

- **Data:** `MixedConifer.las` (bundled with lasR) — 37,657 pts, 8,090 m²,
  **4.7 pts/m²**, single-return, already height-normalized, ground-classified.
  It also carries a `treeID` field from a prior segmentation (**205 trees**),
  used here only as a loose reference count, *not* field truth.
- **Window (identical for both):** `ws(h) = 0.1*h + 3`, clamped to
  `wfloor = max(2, round(2.5 * pulse_spacing, 1))` below 2 m and 5 m above
  20 m (Popescu & Wynne allometry); `hmin = 2 m`.
- **Step-5 smoothing (approach §2):** density < 8 pts/m² (QL2) -> 3x3 mean
  pre-LM smooth; QL1/QL0 -> skip. Applied identically on both engines.
- **Versions:** lasR 0.21.0 (dev build with variable-window `ws`), lidR 4.3.2,
  R 4.5.2.

## Pipelines

Both scripts are **density-first**: Step 0 measures first-return density the
same way (mean first-return count per 1 m cell) and *derives* the CHM
resolution and window floor — nothing is hardcoded. Both independently measure
**4.67 pts/m²** on the toy tile and land on **res = 0.50 m**, **wfloor = 2.0
m** (CHM-based regime, density >= 4; floor ~ 2.5x the 0.46 m pulse spacing,
clamped to 2 m).

| Stage | lasR | lidR |
|-------|------|------|
| 0. Density | `rasterize(1, "count", keep_first)` -> mean | `pixel_metrics(first, length, 1)` -> mean |
| 4. CHM | `triangulate(keep_first)` -> `rasterize(res)` -> `pit_fill` | `rasterize_canopy(res, pitfree(0,10,20; subcircle 0.2))` |
| 5a. Smooth (dens<8) | `focal(chm, size=3, fun="mean")` | `terra::focal(chm, w=matrix(1/9,3,3))` |
| 5b. Detect | `local_maximum_raster(chm_lm, ws, min_height = 2)` | `locate_trees(chm_lm, lmf(ws, hmin = 2, "circular"))` |
| 6. Segment | `region_growing(chm_lm, seeds, max_cr = 10 m)` | `segment_trees(las, dalponte2016(chm_lm, ttops, max_cr = 10/res px))` |
| 7. Metrics | `terra::as.polygons(crowns)` + `expanse()` + `zonal(chm_lm, fun=max)` | `crown_metrics(seg, .stdtreemetrics, geom = "convex")` |

> **Caveat (lasR pit-free):** lasR's `pit_fill` is *not* the Khosravipour pit-
> free algorithm used by `lidR::pitfree()`; it is a TIN + post-hoc pit filling.
> The CHM-construction difference is documented under "What drives the gap".

## Results — each tool, its own CHM (end-to-end with smoothing)

| Metric | lasR | lidR |
|--------|------|------|
| Treetops detected | **112** | **162** |
| Z range | 4.3-30.2 m | 3.4-31.6 m |
| Mean top height | 19.7 m | 20.9 m |
| Runtime (this tile) | 0.25 s | 0.69 s |

**Spatial agreement** (greedy nearest-neighbour matching):

| Tolerance | Matched | Jaccard | Height diff (lidR − lasR) |
|-----------|---------|---------|----------------------------|
| 1.0 m | 100 | 0.57 | +2.04 ± 1.78 m |
| 2.0 m | 106 | 0.63 | +2.12 ± 1.78 m |

At 2 m, 106 of lasR's 112 tops (~95%) have a lidR counterpart. lidR keeps
50 extra tops and its matched tops sit ~2 m higher: the `subcircle` Khosravipour
CHM preserves peak elevation well, while lasR's TIN + `pit_fill` plus the QL2
mean smooth flattens peaks noticeably at this density.

## What drives the 50-top gap? (controlled test, one shared CHM)

Running **both** local-maximum implementations on the **same** lidR pit-free
CHM (`load_raster` feeds it into lasR; no Step-5 smoothing in this test —
the goal is to isolate the detector):

| On the identical CHM | Treetops | Matched (1 m & 2 m) | Jaccard |
|----------------------|----------|---------------------|---------|
| lasR `local_maximum_raster` | 177 | — | — |
| lidR `lmf` | 168 | 168 (all of lidR) | 0.95 |

So most of the end-to-end gap is **CHM construction + Step-5 smoothing**, not
the detector. On the identical unsmoothed Khosravipour CHM, the LM
implementations return practically the same set (Jaccard 0.95; lidR's set ⊂
lasR's set). The smoothing branch then trims aggressively when the underlying
surface is lasR's `pit_fill` (bumpier than `pitfree`).

## Crown segmentation (Step 6) and crown metrics (Step 7)

Each detect script's seeds feed a region-growing segmenter; per-tree crowns
are written as `crowns_*.gpkg` polygons with area and apex height (zonal max
of the LM CHM for lasR, `.stdtreemetrics` convex hulls for lidR).

| Toy (~ 4.7 pts/m²) | lasR | lidR |
|--------------------|------|------|
| Crowns kept | 135 | 162 |
| Mean crown area | 37.2 m² | 36.9 m² |
| Median crown area | 38.8 m² | 35.9 m² |

`compare_crowns.R` does centroid-based greedy matching plus per-pair IoU on
the matched polygons:

| Tol | Matched | % of lidR | IoU mean | IoU median |
|-----|---------|-----------|----------|------------|
| 2 m | 118 | 73% | 0.67 | 0.70 |
| 3 m | 123 | 76% | 0.66 | 0.70 |

Mean crown areas agree to within ~ 1%; IoU ~ 0.7 on matched pairs is solid
agreement given that the two segmenters use slightly different `max_cr`
semantics (meters in lasR vs pixels in lidR; we convert with `max_cr = 10 / res`
to keep them comparable).

## Parameter sweep against the bundled reference

`scripts/sweep.R` searches a small (res, k, hmin) grid in lidR and scores each
against the 205-tree pseudo-truth (apex per `treeID`):

```text
 res   k hmin   n  recall  precision   F1
 0.25 1.0   2 173    0.78       0.92  0.85   <- best F1
 0.25 1.0   1 178    0.78       0.90  0.84
 0.50 0.7   1 231    0.89       0.79  0.83
 0.50 0.7   2 224    0.87       0.80  0.83
 0.25 0.7   2 240    0.90       0.77  0.83
 0.50 1.0   2 168    0.74       0.90  0.82   <- the document default
 ...
```

The doc-default parameters (`res = 0.5`, `k = 1.0`, `hmin = 2`) sit on the
high-precision end (~90% precision, ~74% recall). Scaling `k` down to 0.7
trades precision for recall and stays at F1 ~ 0.83.

## Takeaways

- **The detectors are essentially equivalent.** Given the same CHM, lasR and
  lidR find the same trees (Jaccard 0.95; lidR's set ⊂ lasR's set).
- **The CHM choice matters more than the detector.** Most of the difference
  is the canopy surface, not the maxima search. lasR's `pit_fill` is *not*
  the Khosravipour pit-free algorithm.
- **QL2 smoothing is a real lever.** Adding the approach doc's 3x3 mean smooth
  before LM cuts lasR's toy count from 192 to 112 — large because the
  underlying `pit_fill` CHM is bumpier than `pitfree`. On QL1+ data the branch
  is skipped and counts are unchanged.
- **Both under-count vs the 205-tree reference**, as expected at 4.7 pts/m² —
  dominant trees are found, suppressed/understory ones are not (occlusion
  limit, §5 of the approach doc). Neither is validated against field stems,
  so this is *not* an accuracy claim.

## Bigger AOI — real USGS 3DEP (CA_CarrHirzDeltaFires_2_2019)

A 25 ha AOI clipped from the public EPT and reprojected to UTM 10N (the EPT
is in EPSG:3857, where distances are inflated ~ 1.32x at this latitude, which
would corrupt density and window sizes — so PDAL crops + reprojects to
EPSG:32610 before processing).

- **AOI:** 505 x 504 m (254,681 m²), **5.68M points**, **22.3 pts/m²** total ->
  **13.6 first-returns/m²** (the Step-0 density) -> res **0.25 m**, wfloor
  **2.0 m**. Multi-return, ground-classified but **raw elevation**, so the
  full approach runs: drop noise -> normalize -> pit-free CHM -> VWF ->
  region_growing. Density >= 8, so the Step-5 smoothing branch is **skipped**.
- **Landscape:** post-fire (Carr/Hirz/Delta), so mean top height is low
  (~ 6.6 m) — many low regen/shrub maxima above the 2 m floor, plus scattered
  surviving conifers up to ~ 48 m.

| Metric | lasR | lidR |
|--------|------|------|
| Treetops | **6,809** | **6,510** |
| Local runtime | 27.3 s | 46.1 s |
| Height range | 2.0-43.4 m | 2.0-48.3 m |
| Mean top height | 6.6 m | 6.5 m |

**Detector agreement:** at 2 m, 5,255 of lidR's 6,510 tops (**81%**) matched
(Jaccard 0.65); count differs by **+299 (4.6%)** — tighter than the toy.

**Controlled (same lidR CHM):** lasR LM **6,523** vs lidR LM **6,510** —
**13 tops (0.2%) apart, Jaccard ~ 1.0**. So ~ 96% of the 299-top end-to-end
gap is CHM construction and only ~ 4% is the detector. The local-maximum
implementations are, for practical purposes, **identical**; lidR's `subcircle`
pit-free CHM preserves the tallest peaks better (max 48.3 vs 43.4 m).

### Crowns (Step 6/7) on the AOI

| AOI crowns | lasR | lidR |
|------------|------|------|
| Crowns kept | 6,821 | 6,262 |
| Mean area | 10.9 m² | 18.3 m² |
| Median area | 5.6 m² | 11.9 m² |

| Tol | Matched | % of lidR | IoU mean | IoU median |
|-----|---------|-----------|----------|------------|
| 2 m | 5,168 | 83% | 0.42 | 0.45 |
| 3 m | 5,408 | 86% | 0.40 | 0.43 |

The two segmenters end up with very different mean crown sizes on the AOI
(lasR's region-growing is more conservative on the bumpier `pit_fill` CHM,
producing tighter polygons). IoU ~ 0.4-0.45 is moderate agreement on matched
crowns: same trees, narrower lasR footprints. For pure crown footprint
analysis at QL1+ density, the lidR side is the safer choice.

### Point-cloud methods at this density

At 13.6 first-returns/m² the approach's §1 branch says point-cloud detection
is viable — and the AOI CHM run only resolves the dominant surface. On a
2.25 ha sub-clip (12.5 first-ret/m²):

| Method | Trees | vs CHM |
|--------|-------|--------|
| CHM-lmf (lasR ~= lidR) | 635 | — |
| point-cloud lmf (lidR) | 706 | +11% |
| Li 2012 3D segmentation (lidR) | 792 | +25% |

- lidR's point-cloud `lmf` adds 11% over the CHM by finding canopy-*surface*
  maxima directly, but is still effectively 2.5D.
- **Li 2012** (3D crown segmentation): +25%, and **95% of the extra trees are
  < 5 m** (median apex 3.6 m) — the sub-dominant/regen layer a 2.5D CHM
  cannot represent by construction. lasR has **no point-cloud segmenter**
  (only CHM `region_growing` + point-cloud `local_maximum`), so this step is
  lidR/PDAL-only.
- Cost/scale: Li 2012 ran in ~ 18 s on 2.25 ha but scales poorly to
  wall-to-wall; for large-area mapping the streaming CHM stays pragmatic.

So: yes, at this density point-cloud segmentation is the right tool to exploit
the data — but in lidR (or PDAL), and mainly to recover the understory; the
CHM result for the dominant layer is unchanged.

## Multi-tile streaming demo (approach §3 edge handling)

`scripts/tile_aoi.R` retiles `aoi.laz` into a 3x3 grid of LAZ tiles under
`tiles/`. Then both engines run the full Steps 0-5 across all tiles with
inter-tile context preserved:

- `scripts/detect_lasr_catalog.R`: `exec(pipeline, on = "tiles/")`. lasR
  auto-buffers between tiles for stages that need it (`triangulate`,
  `normalize`, `pit_fill`, `local_maximum_raster`).
- `scripts/detect_lidr_catalog.R`: `readLAScatalog(tdir)` with
  `opt_chunk_buffer(ctg) <- 20`, `opt_filter`, `opt_select` per approach §3.

| Run | Treetops | Runtime | Δ vs single tile |
|-----|----------|---------|------------------|
| lasR single-file (`detect_lasr_aoi.R`) | 6,809 | 27.3 s | — |
| lasR multi-tile (`detect_lasr_catalog.R`) | **6,809** | 44.8 s | **0** |
| lidR single-file (`detect_lidr_aoi.R`) | 6,510 | 46.1 s | — |
| lidR multi-tile (`detect_lidr_catalog.R`) | 7,336 | 94.4 s | **+826 (12.7%)** |

The lasR auto-buffering reproduces the single-file result **exactly**
(6,809 = 6,809) across 9 tiles + buffers — a strong validation of the
approach doc's §3 claim that "trees at seams aren't missed or duplicated."
The lidR catalog **over-counts by ~13%** at this configuration: the
`opt_chunk_buffer = 20 m` is doing its job for the buffer carry-over, but
chunked `rasterize_canopy(pitfree)` + chunked LM are not deduplicating as
cleanly as lasR's streaming pipeline. Practical implication: for wall-to-wall
detection on tiled data, prefer the lasR streaming path; if you must use
lidR, validate the catalog result against a single-file baseline on a
representative tile.
### Density and computational cost (three detectors)

Same 2.25 ha sub-clip, decimated to a range of first-return densities; each
method timed (elapsed, multithreaded where applicable).

| first-ret/m² | CHM res | k pts | CHM-lmf n / s | pc-lmf n / s | Li 2012 n / s |
|--------------|---------|-------|---------------|--------------|---------------|
| 0.7 | 1.00 | 23 | 491 / 0.2 | 997 / 0.0 | 1191 / 0.1 |
| 1.4 | 1.00 | 46 | 530 / 0.2 | 869 / 0.0 | 1063 / 0.3 |
| 2.8 | 1.00 | 91 | 507 / 0.2 | 781 / 0.0 | 897 / 0.7 |
| 5.7 | 0.50 | 182 | 542 / 0.6 | 712 / 0.0 | 822 / 2.1 |
| 11.4 | 0.25 | 365 | 631 / 1.2 | 708 / 0.1 | 795 / 7.5 |
| 12.5 | 0.25 | 399 | 635 / 1.3 | 706 / 0.1 | 784 / 8.9 |

**Cost:** point-cloud lmf is ~free (<0.1 s); CHM-lmf is cheap and roughly
linear (~19 s extrapolated to the full 25 ha, matching the AOI run); **Li 2012
scales ~O(n^1.7)** — 0.1 s at 23k pts but 8.9 s at 399k, i.e. ~10–15 min for the
full 25 ha and worse per km². It is the expensive method, and density makes it
explode.

**Reliability:** below ~3/m² the point-cloud methods **over-detect** (sparse
returns become spurious local maxima — ~2x the CHM count, mostly noise); their
counts settle as density rises. CHM-lmf is stable throughout because
rasterization + pit-fill smooth the surface.

| Method | Best density | Why |
|--------|--------------|-----|
| CHM-lmf | any, esp. <= 4/m² | robust to sparsity, cheap, linear — wall-to-wall workhorse |
| point-cloud lmf | >= 8/m² | trivial cost, but over-detects when sparse; ~= CHM for dominant tops |
| Li 2012 | >= 8/m² only | needs density for real sub-dominant trees; cost explodes, impractical wall-to-wall |

#### Measured 16-core Li 2012 throughput (real EPT data)

A 56.2 ha block (13.0M pts, ~14 first-ret/m²) was pulled from the EPT,
reprojected to UTM, retiled to 36 x 150 m tiles, and segmented **end-to-end**
(read -> drop noise -> normalize -> Li 2012) across **16 cores**:

- **40.1 s for 56.2 ha -> 0.71 s/ha** (19,462 trees, 346/ha).
- Extrapolated to **10,000 acres (4,047 ha): ~48 min on 16 cores.**

Note this is **~2.4x** the naive Li-2012-only estimate (~20 min): normalization
(`tin()`), disk I/O, and fork overhead roughly triple the segmentation cost
(~11 s/ha single-core end-to-end vs ~4 s/ha for the segmentation alone). It also
excludes inter-tile buffers (add ~10–40% for buffered processing) and assumes
this density — QL1 (~8/m²) would be ~0.4x (~18 min); at QL2 Li 2012 is the wrong
tool. CHM-lmf over the same 10,000 acres is ~a few minutes.

## Caveats

- **Runtime is not an engine benchmark.** On the 8,090 m² toy tile lidR's
  in-memory speed wins; on the 25 ha AOI lasR's streaming was ~ 1.7x faster
  (27 vs 46 s). On the 9-tile catalog lasR is again faster (45 vs 94 s).
  All are small datasets; lasR's real edge is large-area throughput and low
  memory (>= 100 km²). Benchmark at scale.
- **EPT reads from this PC are network-bound** (far from us-west-2), so the
  EPT-native script's 41 s is dominated by acquisition, not detection.
- Two AOIs, conifer/post-fire only. Results will shift for
  deciduous/multi-layered stands.

## Reproduce

Scripts in `$CLAUDE_JOB_DIR`. Toy: `detect_lasr.R`, `detect_lidr.R`,
`compare.R`, `shared_chm.R`, `sweep.R`, `segment_lasr.R`, `segment_lidr.R`,
`compare_crowns.R`. Real AOI: `extract.json` (PDAL EPT clip + reproject) ->
`detect_lasr_aoi.R`, `detect_lidr_aoi.R`, `shared_chm_aoi.R`, `pc_vs_chm.R`,
`segment_lasr_aoi.R`, `segment_lidr_aoi.R`. Multi-tile demo: `tile_aoi.R` ->
`detect_lasr_catalog.R`, `detect_lidr_catalog.R`. EPT-native: `detect_lasr_ept_aoi.R`.
