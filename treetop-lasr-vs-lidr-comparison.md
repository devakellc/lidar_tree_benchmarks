# lasR vs lidR Tree-Top Detection — Implementation & Comparison

Both tools run the documented approach (pit-free CHM -> variable-window
local-maximum) on the **same** input with the **same** window function, then
their treetops are matched spatially.

## Setup

- **Data:** `MixedConifer.las` (bundled with lasR) — 37,657 pts, 8,090 m²,
  **4.7 pts/m²**, single-return, already height-normalized, ground-classified.
  It also carries a `treeID` field from a prior segmentation (**205 trees**),
  used here only as a loose reference count, *not* field truth.
- **Window (identical for both):** `ws(h) = 0.1*h + 3`, clamped to 3 m below
  2 m and 5 m above 20 m (Popescu & Wynne allometry); `hmin = 2 m`.
- **Versions:** lasR 0.21.0 (dev build with variable-window `ws`), lidR from
  CRAN binaries; R 4.3.3.

## Pipelines

Both scripts are **density-first**: Step 0 measures first-return density the
same way (mean first-return count per 1 m cell) and *derives* the CHM resolution
and window floor — nothing is hardcoded. Both independently measure **4.67
pts/m²** and land on **res = 0.50 m**, **window floor = 2.0 m** (CHM-based
regime, since density >= 4; floor ~2.5x the 0.46 m pulse spacing).

| Stage | lasR | lidR |
|-------|------|------|
| 0. Density | `rasterize(1, "count", keep_first)` -> mean | `pixel_metrics(first, length, 1)` -> mean |
| 4. CHM | `triangulate(keep_first)` -> `rasterize(res)` -> `pit_fill` | `rasterize_canopy(res, pitfree(0,10,20; subcircle 0.2))` |
| 5. Detect | `local_maximum_raster(chm, ws, min_height = 2)` | `locate_trees(chm, lmf(ws, hmin = 2, "circular"))` |

## Results — each tool, its own CHM (end-to-end)

| Metric | lasR | lidR |
|--------|------|------|
| Treetops detected | **192** | **168** |
| Z range | 2.6–31.7 m | 2.1–32.0 m |
| Mean top height | 20.1 m | 21.2 m |
| Runtime (this tile) | 1.24 s | 0.28 s |

**Spatial agreement** (greedy nearest-neighbour matching):

| Tolerance | Matched | Jaccard | Height diff (lidR − lasR) |
|-----------|---------|---------|----------------------------|
| 1.0 m | 140 | 0.64 | +0.50 ± 0.78 m |
| 2.0 m | 162 | 0.82 | +0.69 ± 1.18 m |

At 2 m, **162 of lidR's 168 tops (96%) have a lasR counterpart** — the two
agree almost completely on dominant/co-dominant trees. lasR returns ~24 extra
tops and lidR's matched tops sit ~0.7 m higher (its `subcircle` pit-free CHM
preserves peak height slightly better).

## What drives the 24-top gap? (controlled test, one shared CHM)

Running **both** local-maximum implementations on the **same** lidR pit-free
CHM (`load_raster` feeds it into lasR):

| On the identical CHM | Treetops | Matched (1 m & 2 m) | Jaccard |
|----------------------|----------|---------------------|---------|
| lasR `local_maximum_raster` | 177 | — | — |
| lidR `lmf` | 168 | 168 (all of lidR) | 0.95 |

So the 24-top end-to-end gap splits into:

- **~15 tops from CHM construction.** lasR's TIN + `pit_fill` surface is
  slightly bumpier than lidR's Khosravipour pit-free CHM, yielding more local
  maxima (192 on its own CHM vs 177 on lidR's).
- **~9 tops from the detector.** On the identical surface lasR's
  `local_maximum_raster` is marginally more sensitive — it returns a strict
  **superset** of lidR's tops (every one of lidR's 168 is matched within 1 m,
  plus 9 more). The local-maximum algorithms themselves agree to **Jaccard
  0.95**.

## Takeaways

- **The detectors are essentially equivalent.** Given the same CHM, lasR and
  lidR find the same trees (Jaccard 0.95; lidR's set ⊂ lasR's set). Pick either
  for the local-maximum step.
- **The CHM choice matters more than the detector.** Most of the difference is
  the canopy surface, not the maxima search. If you need bit-for-bit parity,
  share one CHM between tools.
- **Both under-count vs the 205-tree reference**, as expected at 4.7 pts/m² —
  dominant trees are found, suppressed/understory ones are not (occlusion
  limit, §5 of the approach doc). lasR (192) lands closer to the reference
  count than lidR (168), but neither is validated against field stems, so this
  is *not* an accuracy claim.

## Bigger AOI — real USGS 3DEP (CA_CarrHirzDeltaFires_2_2019)

A 25 ha AOI clipped from the public EPT and reprojected to UTM (the EPT is in
EPSG:3857, where distances are inflated ~1.32x at this latitude, which would
corrupt density and window sizes — so PDAL crops + reprojects to EPSG:32610
before processing).

- **AOI:** 505 x 504 m (254,681 m²), **5.68M points**, **22.3 pts/m²** total ->
  **13.6 first-returns/m²** (the Step-0 density) -> res **0.25 m**, window floor
  2.0 m. Multi-return, ground-classified, but **raw elevation** (so the full
  approach runs: drop noise -> normalize -> pit-free CHM -> VWF).
- **Landscape:** post-fire (Carr/Hirz/Delta), so mean top height is low
  (~6.6 m) — many low regen/shrub maxima above the 2 m floor, plus scattered
  surviving conifers up to ~45 m.

| Metric | lasR | lidR |
|--------|------|------|
| Treetops | **6,813** | **6,512** |
| Local runtime | 9.5 s | 19.2 s |
| Height range | 2.0–43.4 m | 2.0–48.3 m |
| Mean top height | 6.6 m | 6.5 m |

**Agreement:** at 2 m, 5,263 of lidR's 6,512 tops (**81%**) matched
(Jaccard 0.65); count differs by only **+301 (4.6%)** — tighter than the toy.

**Controlled (same lidR CHM):** lasR LM **6,525** vs lidR LM **6,512** —
**13 tops (0.2%) apart, Jaccard ~1.0**. So ~96% of the 301-top end-to-end gap
is CHM construction and only ~4% is the detector. The local-maximum
implementations are, for practical purposes, **identical**; lidR's `subcircle`
pit-free CHM preserves the tallest peaks slightly better (max 48.3 vs 43.4 m).

This reproduces the toy-tile finding on real, dense (QL1+), un-normalized 3DEP
data: same trees, the CHM is the only meaningful difference.

### Point-cloud methods at this density

At 13.6 first-returns/m² the approach's §1 branch says point-cloud detection is
viable — and the AOI run used a CHM, so it only resolves the dominant surface.
On a 2.25 ha sub-clip (12.5 first-ret/m²):

| Method | Trees | vs CHM |
|--------|-------|--------|
| CHM-lmf (lasR ~= lidR) | 635 | — |
| point-cloud lmf (lasR 706, lidR 706) | 706 | +11% |
| Li 2012 3D segmentation (lidR) | 784 | +23% |

- lasR's `local_maximum()` runs on the **point cloud** (`local_maximum_raster`
  is the raster one) and matched lidR's point-cloud `lmf` **exactly**
  (706 = 706). But point-cloud lmf still finds only canopy-*surface* maxima, so
  it adds just ~11% over the CHM.
- The real gain is **Li 2012** (3D crown segmentation): +23%, and **95% of the
  extra trees are < 5 m** (median apex 3.6 m) — the sub-dominant/regen layer a
  2.5D CHM cannot represent by construction. lasR has **no point-cloud
  segmenter** (only CHM `region_growing` + point-cloud `local_maximum`), so this
  step is lidR/PDAL-only.
- Cost/scale: Li 2012 ran in 8 s on 2.25 ha but scales poorly to wall-to-wall;
  for large-area mapping the streaming CHM stays pragmatic. In this post-fire
  stand the extra low detections are regen/shrub clumps — useful signal or
  noise depending on the goal, and unvalidated here.

So: yes, at this density point-cloud segmentation is the right tool to exploit
the data — but in lidR (or PDAL), and mainly to recover the understory; the
CHM result for the dominant layer is unchanged.

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

## Caveats

- **Runtime is not an engine benchmark.** On the 8,090 m² toy tile lidR's
  in-memory speed wins; on the 25 ha AOI lasR's streaming was ~2x faster
  (9.5 vs 19.2 s) — but both are single small tiles. lasR's real edge is
  large-area throughput and low memory (>=100 km²); benchmark at scale.
- **EPT reads from this PC are network-bound** (far from us-west-2), so the
  21 s clip time is not representative either — benchmark acquisition on
  in-region compute (e.g. us-west-2), not locally.
- Two AOIs, conifer/post-fire only. Results will shift for
  deciduous/multi-layered stands.

## Reproduce

Scripts in `$CLAUDE_JOB_DIR`. Toy tile: `detect_lasr.R`, `detect_lidr.R`
(density-first end-to-end), `shared_chm.R` (controlled same-CHM test),
`compare.R` (matching), `sweep.R` (parameter sweep vs the reference). Real AOI:
`extract.json` (PDAL EPT clip + reproject), `detect_lasr_aoi.R`,
`detect_lidr_aoi.R`, `shared_chm_aoi.R`, `pc_vs_chm.R` (point-cloud vs CHM on a
sub-clip); compare with `compare.R tops_lasr_aoi.csv tops_lidr_aoi.csv`.
