# Best Approach: Tree-Top Detection & Tree Segmentation from Airborne LiDAR

*Focus: large tiled airborne LiDAR (e.g. USGS 3DEP), with point density as the
deciding variable. Primary toolchain: `lasR` / `lidR`, with alternatives noted.*

---

## TL;DR — the recommended pipeline

For wall-to-wall USGS 3DEP-class data (QL2, ~2 pulses/m²), the proven, scalable
approach is **CHM-based**, run tile-by-tile with buffers:

```text
measure density -> ground classify -> DTM -> height-normalize ->
pit-free CHM -> variable-window local-maxima (treetops) ->
marker-controlled / region-growing segmentation -> crown metrics
```

- **Detection winner:** variable-window local-maximum filter (VWF) on a
  pit-free CHM. It topped both major international benchmarks (Eysn 2015,
  Kaartinen 2012).
- **Segmentation winner:** `dalponte2016` (region growing seeded by the
  treetops) or marker-controlled watershed.
- **Engine:** `lasR` for >=100 km² (single streaming pass, low memory); `lidR`
  LAScatalog for smaller jobs or to mix detection/segmentation algorithms. The
  variable-window filter now runs natively in `lasR` too (dev branch — see §3).
- **Hard limit to accept up front:** at QL2 you reliably get **dominant +
  co-dominant** trees (~80–90% in open conifer). Suppressed/understory trees are
  largely invisible — this is an occlusion/physics limit, not a tuning problem.

---

## 1. Why density decides everything

Measure your actual **first-return canopy density** before choosing anything —
the nominal Quality Level understates it because forests return multiple points
per pulse.

| QL | Pulse density (ANPD) | Pulse spacing (ANPS) | Use |
|----|------|------|-----|
| QL0 | >= 8.0 /m² | <= 0.35 m | research-grade |
| QL1 | >= 8.0 /m² | <= 0.35 m | **practical threshold for individual-tree work** |
| **QL2** | **>= 2.0 /m²** | **<= 0.71 m** | **3DEP national baseline (most US coverage)** |
| QL3 | >= 0.5 /m² | <= 1.41 m | legacy/coarse |

> The 3DEP spec is written in **pulses**; `lidR::density()` / `lasinfo` report
> **points** (returns), which is higher in forest. Base your design on measured
> first-return density, not the QL label.

**The branch that follows from density:**

- **<~ 4 pts/m² (QL2, typical 3DEP) -> CHM/raster-based detection.** The cloud
  is too sparse for robust 3D clustering. You get the canopy *surface* reliably
  but not stems or understory.
- **>~ 8 pts/m² (QL1/QL0, UAV/drone) -> point-cloud detection becomes viable**
  (Li 2012, deep-learning instance segmentation), and you can start resolving
  sub-dominant trees. For *tops* specifically, though, a pit-free CHM still
  captures the dominant canopy well enough that point-based local maxima rarely
  beat raster LM until density is much higher — the point-cloud win is in
  *crowns* and sub-canopy stems, not apex detection.

**Density rules of thumb (from the literature):**

- CHM cell size ~ pulse spacing or slightly coarser. QL2 -> **0.5–1.0 m**;
  QL1/QL0 -> 0.25–0.5 m. Never go finer than the spacing — you just create pits.
- Smallest reliably detectable crown ~ **2–3x pulse spacing** (QL2 => ~1.5–2 m
  diameter).
- Detection accuracy **saturates around 8–10 pulses/m²** (Sparks et al. 2022: 8
  vs 22 pulses/m² gave F ~ 0.47 vs 0.50 — barely any gain). Below ~1 pulse/m²
  even basic metrics collapse (Jakubowski et al. 2013).
- Fixing understory detection would need **~170 pts/m²** (Hamraz et al. 2017) —
  far beyond any operational ALS.

---

## 2. The recommended pipeline, step by step

### Step 0 — Measure density & pick parameters

In `lasR`, `summarise()` now returns `$density` (points/m²) and `$pulse_density`
(first returns/m²); in `lidR` use `density(las)` or map it with
`rasterize_density()`. Check first-return density and spatial uniformity
(flight-line overlap inflates aggregate density). This sets CHM resolution and
detection window before anything else.

### Step 1 — Ground classification

Cloth Simulation Filter (CSF) or Progressive Morphological Filter (PMF). At QL2
under canopy, ground returns are sparse — this is your main source of *height*
error downstream.

### Step 2 — DTM

Triangulation (TIN) is parameter-free and fast (needs a tile buffer); IDW
(`knnidw`) is more robust at edges. Typical resolution 1 m.

### Step 3 — Height-normalize

Subtract terrain so Z becomes height-above-ground. Per-point TIN normalization
avoids raster discretization error. **Every later step assumes normalized
heights.**

### Step 4 — Build a **pit-free CHM** *(the step that most affects accuracy)*

Use the **pit-free algorithm (Khosravipour et al. 2014)** — a stack of partial
CHMs merged at successive height thresholds — rather than blanket Gaussian
smoothing.

- Pits (deep-penetrating pulses) create false local maxima and split crowns.
- Over-smoothing merges adjacent tops (under-segmentation) and biases
  heights/crown widths low.
- A `subcircle` (~0.15–0.2 m disc per point) fills gaps in sparse QL2 data.
- A *light* smooth before detection is still common.

### Step 5 — Detect treetops with a **Variable Window Filter (VWF)**

Local-maximum filter whose search window grows with canopy height (taller trees
-> wider crowns, the Popescu & Wynne 2004 allometry). Variable windows beat
fixed windows (~85% vs ~80% accuracy) and won both benchmarks.

- **Rule:** search window <= mean crown diameter of the stand.
- QL2 -> larger/clamped windows (~5 m) + smoothing to suppress noise peaks.
- QL1/QL0 -> smaller windows (~3 m), less smoothing, catches more small trees.
- Set `hmin` ~ 2 m (nothing below is a treetop).

### Step 6 — Segment crowns from the detected tops

- **`dalponte2016`** (CHM + treetops, region growing) — most popular,
  well-validated for ALS. Defaults: `th_seed=0.45, th_cr=0.55, max_cr=10`.
- **Marker-controlled watershed** — robust default for conifers; the treetops
  are the markers so you get one basin per tree.
- **`silva2016`** (Voronoi/nearest-top) — fast, good for open/plantation
  conifer.
- **`li2012`** (pure point-cloud, no CHM/seeds) — only worth it at high
  density; its region growing scales super-linearly with point count (the lidR
  implementation is unparallelized, worse than O(N²)), so it gets painful over
  large tiles at >=8 pts/m².
- **AMS3D (adaptive mean-shift 3D)** — the one practical point-cloud method
  that can edge out `dalponte2016` for crowns. It topped the most rigorous ALS
  broadleaf comparison (Dalponte2016 a close second), avoids CHM interpolation
  artifacts, and handles leaning/interlocking crowns better — at the cost of
  speed and no CHM seeds. Reserve it for plots where CHM pits or interlocking
  broadleaf crowns visibly hurt the raster pipeline; needs ~8+ pts/m² and the
  compute budget. No CRAN/lidR implementation (see deep-research-report.md).

### Step 7 — Crown geometry & metrics

Extract per-tree height, crown area/diameter, convex/concave hull polygons.

---

## 3. Tooling — lasR, lidR, and the alternatives

### lasR (recommended engine for large/tiled data — your toolchain)

Same author as lidR (Jean-Romain Roussel), rebuilt for **large-area,
low-memory** streaming. Lazy pipeline assembled with `+`, run in one read pass
per file, with **automatic inter-tile buffering**. Its docs explicitly say for
>=100 km²: "lidR will probably fail; use lasR."

```r
library(lasR)

# Step 0 — native density (dev branch): summarise() returns $density
# (points/m²) and $pulse_density (first returns / area).
s <- exec(summarise(), on = "usgs_tiles/")
s$density        # point density -> picks CHM res
s$pulse_density  # pulse density (ReturnNumber == 1) -> the USGS QL unit

# Detection + segmentation pipeline
del  <- triangulate(filter = keep_first())   # DSM TIN
chm  <- rasterize(0.5, del)                  # CHM (res from density)
chm2 <- pit_fill(chm)                        # pit-free CHM

# Variable window (dev branch): ws may be a function of height.
ws   <- function(h) { y <- 0.1*h + 3; y[h < 2] <- 3; y[h > 20] <- 5; y }
seed <- local_maximum_raster(chm2, ws)       # VWF treetops
tree <- region_growing(chm2, seed)           # crowns

pipeline <- del + chm + chm2 + seed + tree
ans  <- exec(pipeline, on = "usgs_tiles/")
```

As of the development (`pre-devel`) branch, both density and the variable window
are **native to lasR** — the workflow above needs no `lidR`:

- **Density:** `summarise()` returns `$density` (point density = npoints / area)
  and `$pulse_density` (pulse density = first returns / area). Following lidR's
  convention, pulses are counted as first returns (`ReturnNumber == 1`), not from
  `gpstime`. This is a full streaming pass, not header-only — for an instant
  header-level estimate, `lidR::density(ctg)` is still cheaper.
- **Variable window:** `local_maximum` and `local_maximum_raster` now accept `ws`
  as a *function* of the working attribute (height `Z` for the point cloud, the
  cell value for a CHM), so the VWF (Popescu & Wynne allometry) runs natively. A
  constant `ws` (a single number) still works as before.

Caveat: these two features live on the development branch and are not yet in a
released/CRAN version — pin the dev build if you depend on them. (Also verify
lasR's multi-threading API in the lasR3 vignette before scripting the parallel
config.)

### lidR (decoupled algorithms, or for smaller jobs)

The mature, function-per-stage package; **detection and segmentation are
decoupled** so you can mix algorithms.

```r
library(lidR); library(future)
plan(multisession)                              # parallel chunks

ctg <- readLAScatalog("usgs_tiles/")
opt_chunk_buffer(ctg) <- 20                     # 10-30 m: avoid edge effects
opt_filter(ctg)       <- "-drop_withheld"
opt_select(ctg)       <- "xyzrn"

# 1. Ground + DTM + normalize
ctg <- classify_ground(ctg, csf())
ctg_norm <- normalize_height(ctg, tin())

# 2. Pit-free CHM
chm <- rasterize_canopy(ctg_norm, res = 0.5,
          pitfree(thresholds = c(0,10,20), max_edge = c(0,1.5), subcircle = 0.2))

# 3. Variable-window treetop detection
f <- function(x) { y <- x*0.1 + 3; y[x<2] <- 3; y[x>20] <- 5; y }
ttops <- locate_trees(ctg_norm, lmf(ws = f, hmin = 2, shape = "circular"))

# 4. Segmentation
ctg_seg <- segment_trees(ctg_norm, dalponte2016(chm, ttops))

# 5. Crowns
crowns <- crown_metrics(ctg_seg, func = .stdtreemetrics, geom = "convex")
```

**Edge handling is critical at tile boundaries:** `opt_chunk_buffer` >= largest
crown radius (10–30 m). The engine loads neighbor points for context, then
strips the buffer and drops treetops whose seed falls in it — no
missed/duplicated trees at seams.

### Best non-R alternatives

| Tool | ITD step | Type | Cost | When to pick |
|------|----------|------|------|--------------|
| **ForestTools** (R) | `vwf()` + `mcws()` | CHM | Free | Cleanest open CHM detect+segment; pairs with lidR CHMs. Top open pick. |
| **FUSION** (USFS) | `CanopyMaxima` + `TreeSeg` | CHM | Free | USGS-lineage CLI, scales over thousands of tiles, USFS-aligned outputs. |
| **LAStools** | `lasground`/`lasheight`/`las2dem`/`blast2dem` | CHM prep | Mixed/commercial | Fastest engine to ground-classify & build CHMs at scale; feed into FUSION/ForestTools (no native treetop tool). |
| **PDAL** | `filters.smrf`->`hag_nn`->`filters.litree` | Point cloud (Li 2012) | Free | All-in-one open pipeline, **reads EPT** (matters for 3DEP-on-S3). Li 2012 weaker than CHM watershed at QL2. |
| **WhiteboxTools** | `IndividualTreeDetection` | Point cloud | Free/MIT | Fast Rust binary, point-based seeds, scriptable from R/Python; density-sensitive. |
| **LiDAR360** / **ArcGIS Pro** | ALS Forest Module / surface-analysis path | Both | Commercial | Turnkey GUI; ArcGIS "surface analysis" path is explicitly the low-density option. |

**Avoid at QL2** (all assume dense TLS/MLS/UAV with visible stems &
understory): TreeLS, 3DFin, treeX, treeseg, treeiso, and the deep-learning 3D
nets (TreeLearn, ForAINet, SegmentAnyTree, TreeAIBox).

---

## 4. Should you use deep learning?

**Not for QL2 airborne data.** The state-of-the-art 3D instance-segmentation
systems are genuinely strong (F1 ~ 0.85–0.97) **but only on dense clouds (>~50
pts/m²: TLS/MLS/UAV).** The leading sensor-agnostic model, **SegmentAnyTree**,
already degrades at ~10 pts/m² — five times denser than QL2. No leading 3D
system is validated at ~2 pts/m²; at that density it is out of distribution.

- **If you must use DL on sparse data**, the realistic route is **image-based**:
  the **DeepForest** Python package (RetinaNet on RGB/CHM, recall ~72% /
  precision ~64%) or **Segment Anything (SAM) on a CHM raster**. Rasterization
  tolerates low density. SAM out-of-the-box does *not* yet beat a trained
  Mask R-CNN, though.
- The 3D nets are also data-hungry on labels; the main labeled benchmark
  (FOR-instance) is high-density ULS — there is no comparable sparse-ALS
  instance-labeled dataset.
- Note: a tuned *unsupervised* classical method (treeX) and a point-based
  hierarchical clustering method matched or beat DL on some forests **with zero
  training labels** — DL's edge is real but contingent on dense data + matching
  training domain.

---

## 5. What accuracy to expect (and how to report it)

Accuracy is governed more by **forest structure** than by algorithm or density
(above the saturation point). The benchmarks converge on this gradient:

| Forest type | Expected detection | Source |
|-------------|--------------------|--------|
| Open / plantation / single-layer conifer | **~80–90%+** | Li 2012, Véga 2014, Eysn 2015 |
| Managed / natural conifer | high for dominants | Kaartinen 2012, Wang 2016 |
| Mixed / multi-layered | substantially lower | Eysn 2015 (overall matching ~47%) |
| Deciduous / broadleaf | lower (~81% even with DL) | Ayrey 2017, Vauhkonen 2012 |
| Tropical | sub-canopy "all but impossible" | Coomes 2017 |

- **Dominant/co-dominant trees** detect well everywhere; **intermediate and
  suppressed trees are largely missed** (omission ~40% to ~100% for suppressed;
  Wang 2016). This is the field's defining limitation.
- **Report accuracy stratified by crown class** — a single "detection rate"
  hides the overstory/understory split.
- For deciduous targets, **leaf-off acquisition is strongly preferred**;
  leaf-off + leaf-on fusion segments ~12% more trees correctly.

**At 3DEP QL2 specifically:** proven for CHMs (R² ~ 0.85 vs field), dominant
height, and site index (RMSE < 10%); feasible for dominant-tree detection in
open/plantation conifer. **QL1 (~8 pts/m²) is the recommended step-up for a
comprehensive individual-tree inventory.**

---

## 6. Pitfalls

- **Over-segmentation** (one tree -> several tops): window/smoothing too small,
  CHM pits, bumpy deciduous crowns, high density picking up large branches.
- **Under-segmentation** (trees merged / small ones missed): window/smoothing
  too large, touching crowns in closed canopy, suppressed trees below the
  surface.
- **CHM pits** are the classic CHM-method failure — always use a pit-free CHM.
- **Tile edges** — without a buffer, trees at seams are missed or duplicated.
- **Parameter sensitivity** — CHM-smoothing window and detection-window size
  matter as much as density; calibrate against field plots per forest type.
- **Conifer vs deciduous** — local-maxima/watershed assume one sharp apex;
  broad/multi-domed deciduous crowns break that assumption (use layer stacking
  or multi-scale PTrees for hard deciduous cases).

---

## Key references

**Methods:** Popescu & Wynne (2004, VWF); Li et al. (2012, point-cloud region
growing); Dalponte & Coomes (2016); Silva et al. (2016); Véga et al. (2014,
PTrees); Ayrey et al. (2017, layer stacking); Khosravipour et al. (2014,
pit-free CHM); Hyyppä et al. (2001, CHM seed+region-grow).

**Benchmarks/reviews:** Kaartinen et al. (2012, EuroSDR/ISPRS); Eysn et al.
(2015, Alpine); Wang et al. (2016, per-crown-class); Vauhkonen et al. (2012);
Lindberg & Holmgren (2017, review).

**Density/structure:** Jakubowski et al. (2013); Sparks et al. (2022); Hamraz et
al. (2017); Falkowski et al. (2008); Coomes et al. (2017); Pourshamsi et al.
(2024, 3DEP QL2 site index).

**Deep learning:** TreeLearn (Henrich 2024); SegmentAnyTree (Wielgosz 2024);
ForAINet (Xiang 2023); FOR-instance (Puliti 2023); DeepForest (Weinstein 2020).

**Software docs:** lidR book (r-lidar.github.io/lidRbook), lasR docs
(r-lidar.github.io/lasR), ForestTools, FUSION, PDAL `filters.litree`,
WhiteboxTools.
