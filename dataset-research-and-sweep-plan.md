# LiDAR Tree Detection — Ground-Truth Datasets & Parameter-Sweep Plan

**Scope.** Selecting ground-truth datasets for benchmarking individual-tree detection (ITD) and crown-segmentation algorithms on **airborne** LiDAR, with point density treated as the deciding variable, and a concrete plan for a parameter sweep of a CHM-based `lasR`/`lidR` pipeline. Geographic priority: California and Oregon conifer forests.

**Status.** Research synthesis + proposed experimental protocol. Last updated 2026-06-02.

---

## Contents

1. [Part 1 — Ground-truth datasets paired with high-density LiDAR](#part-1)
2. [Part 2 — Dataset selection for the parameter sweep](#part-2)
3. [Part 3 — Suggested protocol to proceed](#part-3)
4. [References & data access](#references)

---

<a name="part-1"></a>
## Part 1 — Ground-truth datasets paired with high-density LiDAR (US, CA/OR focus)

### Requirements used to screen datasets

- **Ground truth = manually-validated individual tree stem locations** (field-surveyed stem positions or hand-delineated crowns), suitable for scoring ITD/segmentation. Model-derived-only annotations are out of scope (usable for weak supervision only).
- **LiDAR density > 8 pts/m²** (QL1/QL0, research-grade airborne, or UAV/mobile laser). Any sensor source.
- **Preference for California / Oregon** conifer forest; other US/relevant datasets included with caveats.

### Why 8 pts/m² is a defensible cutoff

The threshold is not arbitrary. LaRue et al. (2022, *Ecosphere*) found the minimum density to reliably estimate structural-diversity metrics ranged from ~2 to ~7.5 pts/m², and recommended **at least ~8 pts/m²** for the most robust spatial/temporal comparisons. This aligns with the operational rule that QL2 (~2 pulses/m², the 3DEP national baseline) reliably yields the canopy *surface* and dominant/co-dominant trees, while QL1/QL0 (≥8 pts/m²) is the practical step-up for a comprehensive individual-tree inventory.

> **Pulses vs points.** The 3DEP spec is written in **pulses**; `lidR::density()` / `lasinfo` report **points** (returns), which is higher under canopy. Always design against *measured first-return density*, not the QL label.

### Headline conclusion

There is currently **no single openly-released benchmark that perfectly combines (a) field-surveyed stem coordinates and (b) co-located >8 pts/m² LiDAR in California or Oregon.** Each strong candidate carries a caveat. The two best practical routes are:

- **California:** NEON **SOAP 2021** acquisition (~20 pts/m², Optech Galaxy) + NEON field stem maps — the one openly-available CA pairing that clears 8 pts/m².
- **Oregon:** **HJ Andrews** field-surveyed ForestGEO/FOREG stem plots + the 2008/2011 Leica ALS50 LiDAR (~8 pts/m² native), assembled by manual co-registration.

### Datasets (most relevant CA/OR high-density first)

#### A. NEON Woody Vegetation Structure + NEON discrete-return LiDAR — California (SJER, SOAP, TEAK)
- **Access:** NEON Data Portal — field stems `DP1.10098.001` (Woody Plant Vegetation Structure); LiDAR `DP1.30003.001` (discrete-return point cloud). Open, no registration; also mirrored on OpenTopography.
- **Ground truth (manual):** Field-surveyed stems. In 20×20 m distributed and 40×40 m tower plots, stems >10 cm DBH are mapped by an offset technique — a TruPulse 360R rangefinder (≈30 cm accuracy) records distance + azimuth from a reference point with known coordinates; species, height, crown diameter recorded. Precise per-tree coordinates are computed with the `geoNEON` R package (`getLocTOS()`); raw download stores `stemDistance`/`stemAzimuth` + plot-level `coordinateUncertainty`.
- **LiDAR density — the key variable:**
  - **SOAP 2021 ≈ 20 pts/m² (Optech Galaxy Prime) — QUALIFIES.** Best openly-available CA high-density pairing.
  - SOAP 2019 / most SJER/TEAK years ≈ 4–6 pts/m² (Optech Gemini era) — **below threshold.**
  - Sensor timeline: Optech Gemini 2013–2020 (~4–6 pts/m²); Riegl Q780 from 2018 (up to 400 kHz); Optech Galaxy Prime from 2021 (up to 1000 kHz, 20+ pts/m²); Riegl VQ-780 II-S from 2025 (prospective, should comfortably exceed 8 pts/m²).
- **Forest type:** SJER = oak / foothill-pine open woodland; SOAP = mixed conifer/deciduous; TEAK = red-fir / subalpine mixed conifer. (Domain D17, Pacific Southwest.)
- **CRS:** UTM zone 11N (EPSG:32611) for D17; WGS84 lat/long also provided.
- **Caveat:** Density varies by site-year — select high-PRF (Galaxy/Riegl) acquisitions. Plot-level positional uncertainty limits exact stem-to-return matching.

#### B. NeonTreeEvaluation Benchmark
- **Access:** Zenodo `10.5281/zenodo.5914554` (data), `10.5281/zenodo.5912283` (training); GitHub `weecology/NeonTreeEvaluation` + companion R evaluation package; paper Weinstein et al. 2021, *PLOS Comp. Bio.* (`10.1371/journal.pcbi.1009180`).
- **Ground truth (manual):** (1) RGB-annotated crown bounding boxes for 22 NEON sites, (2) field-annotated crown polygons for 2 sites, (3) field stem points for 14 sites (from NEON WVS), plus hand-labeled "draped" LiDAR labels. ~30,975 tree annotations total.
- **Location/forest type:** 22 US sites incl. California SJER, SOAP, TEAK.
- **LiDAR density:** ~5 pts/m² (NEON discrete return). **BELOW 8 pts/m².** Value for this work is the **scoring harness**, not the density.
- **License:** Open.

#### C. IDTReeS 2020 Competition Data
- **Access:** Zenodo `10.5281/zenodo.3934932` (competition), `10.5281/zenodo.3700197` (crown polygons); idtrees.org.
- **Ground truth (manual):** Hand-delineated crown polygons + NEON field stem/species data. Tasks: crown delineation + species classification.
- **Location/forest type:** 3 NEON sites incl. California SJER (+ OSBS Florida, MLBS Virginia).
- **LiDAR density:** NEON discrete return ~5 pts/m². **BELOW 8 pts/m².**
- **License:** Open.

#### D. HJ Andrews Experimental Forest, Oregon — LiDAR + stem maps
- **Access:** Andrews Forest data catalog (andlter.forestry.oregonstate.edu); LiDAR datasets `GI010` (Aug 2008) and `GI011` (Oct 2011).
- **LiDAR density / sensor:** Leica ALS50 Phase II, designed for **average ~8 pts/m²** native density, up to 4 returns/pulse, ~100% overlap. **MEETS the threshold (borderline).**
- **Ground truth:** The publicly catalogued HJA stem map (`TV081`) is **LiDAR-model-derived (TreeVaW from the CHM)** — *not* field-validated. **However**, HJ Andrews also hosts genuine **field-surveyed stem-mapped plots**: the Andrews Forest Dynamics Plots / Smithsonian ForestGEO network (one 12-ha + six 3-ha plots, 450–1,350 m) and the FOREG plots (established 2019, total-station mapped, every stem >1 cm DBH). These hand-mapped stems can be paired with the 2008/2011 8 pts/m² LiDAR — a strong but manually-assembled Oregon pairing.
- **Forest type:** Tall Douglas-fir / western hemlock PNW conifer — among the tallest, highest-biomass forests on Earth (the hardest structural case in-domain).

#### E. TreeScope (Cheng et al., IEEE ICRA 2024)
- **Access:** treescope.org; GitHub `KumarRobotics/treescope`; arXiv `2310.02162`. License CC BY-NC-SA 4.0. ~2.2 TB.
- **Ground truth (manual):** ~1,860 manually annotated stem semantic labels + field DBH; per-tree point clouds; ground-truth map PNGs/JSON.
- **LiDAR sensor / density:** Mobile (MLS) Ouster OS0-128 and UAV (ULS) Ouster OS1-64 (Velodyne VLP-16 at one NJ site). **Explicit pts/m² not published — FLAG as not reported** (mobile/UAV laser of this type is typically orders of magnitude denser than 8 pts/m²).
- **Location/forest type:** Forests = Virginia (loblolly/Virginia/pitch/white pine, oak) + New Jersey (pitch pine, oak, Atlantic white cedar). Orchards = Merced County, California (almond, pistachio). **Not CA/OR forest; CA component is agricultural orchards.**

#### F. FOR-instance / FOR-instanceV2 (Puliti et al. 2023; ForestFormer3D 2025)
- **Access:** arXiv `2309.01279`; Zenodo `10.5281/zenodo.8287792`. FOR-instanceV2 + ForestFormer3D (arXiv `2506.16991`) described by authors as "to be released."
- **Ground truth (manual):** Hand-annotated semantic + instance segmentation of individual trees (stem / branches / terrain / low veg). ~1,100 trees, ~2.8 ha, 29 plots.
- **LiDAR:** UAV laser scanning (ULS), high density.
- **Location:** Norway, Czech Republic, Austria, New Zealand, Australia. **Not US.** Leading high-density 3D instance-segmentation benchmark; listed for completeness.
- **License:** Open (CC).

#### G. Open Forest Observatory — "Emerald Point," Lake Tahoe, California
- **Access:** docs.openforestobservatory.org / openforestobservatory.org/data; paper Young et al. 2022, *Methods Ecol. Evol.* (`10.1111/2041-210X.13860`). License CC BY 4.0.
- **Ground truth (manual):** 3.23-ha field stem map; **1,775 trees >5 m** (≈2,122 if all >7.5 cm DBH). Rangefinder distance + compass azimuth from reference grid.
- **3D-data caveat:** Co-located 3D is **drone RGB structure-from-motion photogrammetry, NOT LiDAR.** No co-located >8 pts/m² LiDAR for this site. CHM built by subtracting a USGS DEM from the photogrammetric DSM. OFO hosts >150 ground-mapped plots within drone footprints, but the 3D products are photogrammetric.
- **Forest type:** Mixed conifer (ponderosa/Jeffrey pine, incense cedar, white fir), Lake Tahoe, CA.

#### H. FIA (Forest Inventory & Analysis) — NOT usable for co-location
- Plot coordinates are deliberately perturbed: ~95% fuzzed within 0.8 km, ~5% within 1.6 km, and 1–20% of private-land plots swapped; true coordinates are protected by federal law (Food Security Act 1985) and require a Memorandum of Cooperation (the confidential-data request process has been suspended). **Unsuitable for precise stem-to-LiDAR pairing.**

#### Other relevant
- **Dubrovin et al. 2024** (*Sci. Reports* `10.1038/s41598-024-72669-5`): field inventory of 3,600 trees over 10 plots + co-located UAV LiDAR + RGB; CRS EPSG:32640. **Methodologically ideal template, but located in Perm Krai, Russia — not US.**
- **NEON Crowns** (Weinstein et al., *eLife* 2020), **SegmentAnyTree**, **TreeLearn**: model-derived predictions / non-CA-OR high-density training data — weak-supervision / pre-training only.

### Summary comparison

| Dataset | Location | Forest type | LiDAR source & density | Ground truth (manual?) | Trees / extent | Access / license | Meets >8 pts/m² + manual stems in CA/OR? |
|---|---|---|---|---|---|---|---|
| **NEON WVS + LiDAR — SOAP 2021** | Soaproot Saddle, CA | Mixed conifer/deciduous | Optech Galaxy, **~20 pts/m²** | Field stems ✅ | 10,000s network-wide | NEON Portal, open | **YES (best CA)** |
| NEON WVS + LiDAR — SJER/TEAK (most years) | Central CA | Oak woodland / red fir | NEON, **~4–6 pts/m²** | Field stems ✅ | 10,000s network-wide | NEON Portal, open | No (LiDAR too sparse) |
| NeonTreeEvaluation | 22 US sites incl. SJER/SOAP/TEAK | Multiple | NEON, **~5 pts/m²** | RGB boxes + field stems + draped labels ✅ | ~30,975 trees | Zenodo/GitHub, open | No (LiDAR <8) |
| IDTReeS 2020 | 3 NEON sites incl. SJER | Multiple | NEON, **~5 pts/m²** | Hand crown polygons + field stems ✅ | 3 sites | Zenodo, open | No (LiDAR <8) |
| **HJ Andrews — ForestGEO/FOREG + LiDAR** | HJ Andrews, OR | Douglas-fir / W. hemlock old growth | Leica ALS50, **~8 pts/m²** | Field total-station stems ✅ (public TV081 is model-derived ✗) | 12-ha + 6×3-ha; FOREG | Andrews catalog; partly request-based | **YES if field plots paired manually (best OR)** |
| TreeScope | VA & NJ forests; CA orchards | Pine; almond/pistachio | Ouster OS0/OS1, VLP-16; **density unstated** | Manual stem labels + DBH ✅ | ~1,860 trees | treescope.org, CC BY-NC-SA | Partial (not CA/OR forest) |
| FOR-instance / V2 | NO/CZ/AT/NZ/AU | Multiple | UAV laser, high density | Hand instance + semantic ✅ | ~1,100 trees; 2.8 ha | Zenodo, open | No (not US) |
| OFO "Emerald Point" | Lake Tahoe, CA | Mixed conifer | **Photogrammetry, NO LiDAR** | Field stems ✅ | 1,775 trees >5 m; 3.23 ha | OFO/UC Davis, CC BY 4.0 | No (no LiDAR) |
| Dubrovin et al. 2024 | Perm Krai, Russia | Mixed dense | UAV LiDAR + RGB | Field inventory stems ✅ | 3,600 trees; 10 plots | SciRep, open | No (not US) |
| FIA | US-wide | All | n/a | Field stems, **coords fuzzed/swapped** ✗ | National | Public (perturbed) | No (coords unusable) |

### Caveats carried forward

- The 8 pts/m² cutoff eliminates most NEON airborne data; only Riegl-Q780 (2018+) and Optech-Galaxy (2021+) site-years qualify. The popular NeonTreeEvaluation/IDTReeS benchmarks were built on ~5 pts/m².
- "Manual ground truth" must be checked per dataset: NeonTreeEvaluation boxes are mostly RGB hand-annotations; NEON Crowns and the HJA `TV081` map are model-derived; OFO stems are field-surveyed but its 3D data is photogrammetry.
- Coordinate precision limits matching: NEON stems carry plot-level uncertainty; FIA is fuzzed/swapped; TreeScope ships maps + per-tree clouds rather than a stem-coordinate table.
- The cleanest high-density + hand-validated benchmarks (FOR-instance, Dubrovin) are outside the US; the cleanest CA/OR field stems (NEON, OFO, HJA ForestGEO) are paired with sub-8 pts/m² LiDAR, photogrammetry, or model-derived stems.

---

<a name="part-2"></a>
## Part 2 — Dataset selection for the parameter sweep

**Reframe.** The pipeline under test is **airborne and CHM-based**, with density as the deciding variable. Many parameters being swept (CHM cell size, the variable-window-filter allometry, smoothing, `hmin`, the `dalponte2016`/watershed thresholds) are *meant to move with density*. So the dataset must let us vary density while holding ground truth and forest structure constant — and the experimental design matters as much as the dataset choice. This excludes the UAV/mobile-laser and photogrammetry datasets despite their density.

### Primary recommendation: NEON as a density ladder via decimation

1. **Anchor** on NEON **SOAP 2021** (~20 pts/m², Optech Galaxy) — the openly-available CA pairing that clears 8 pts/m² with field-mapped stems.
2. **Decimate** to a density ladder — e.g. **20 → 8 → 4 → 2 → 1 pts/m²** — using `lidR::decimate_points(homogenize())` or `random()`.
3. **Run the full parameter sweep at each rung against the same stem ground truth.**

This isolates density as the single variable while holding forest structure and ground truth constant — precisely the experiment the "density decides everything" thesis demands, and the approach used by Jakubowski et al. (2013) and Sparks et al. (2022).

> **Caveat.** Natively-acquired low-density data differs from decimated data (pulse characteristics, flight geometry). Decimation is the standard, most controlled way to map parameter sensitivity vs density, but it should be reported as a simulation, ideally cross-checked against a native QL2 acquisition (see 3DEP note below).

**Why NEON for this sweep specifically**
- It is **airborne** (matches the pipeline; UAV/mobile-laser sets do not).
- Field stems carry **height/DBH**, enabling accuracy **stratified by crown class** — essential, because a single detection rate hides the overstory/understory split.
- Multiple sites span the **forest-structure gradient** that governs accuracy more than density: **SJER** (open oak / foothill pine) → **SOAP** (mixed conifer) → **TEAK** (dense red fir / subalpine). Sweeping the same parameters across all three exposes structure-dependence, not just density-dependence.

### Oregon / PNW conifer: HJ Andrews

Pair the **field-surveyed ForestGEO/FOREG stem maps** (total-station mapped) with the **2008/2011 Leica ALS50 LiDAR (~8 pts/m² native)** for a genuine tall Douglas-fir / western-hemlock test at the QL1 threshold — structurally the hardest in-domain case and the most relevant to timber-harvest-planning work. Requires manual co-registration of plots to LiDAR; it is the best Oregon airborne pairing.

### Scoring machinery: NeonTreeEvaluation

Its LiDAR is only ~5 pts/m², but the value is the **evaluation package** — it implements detection matching and precision/recall scoring against hand-annotations across the California sites. Use it as the **scoring scaffold** and as the QL2-density data point, with the decimated SOAP ladder supplying the density axis.

### Explicitly excluded for this sweep (with reasons)

- **TreeScope, FOR-instance** — >50 pts/m² ULS/MLS built for 3D instance segmentation, the method family the airborne CHM pipeline explicitly avoids at airborne densities. Including them confounds sensor/platform with density.
- **OFO Emerald Point** — photogrammetry, not LiDAR.
- **FIA** — coordinates unusable for co-location.

### Branch on deployment target

- **If deploying on QL2 national 3DEP coverage:** center the sweep on the **2–4 pts/m²** rungs, and pull matching **3DEP tiles** over the NEON sites for a native-density cross-check alongside the decimated ladder.
- **If deploying on QL1/QL0:** weight toward the **8–20 pts/m²** rungs; the high-density NEON anchor and HJ Andrews are directly representative.

---

<a name="part-3"></a>
## Part 3 — Suggested protocol to proceed

### 1. Assemble the test matrix

| Axis | Levels (suggested) |
|---|---|
| Site / forest structure | SJER (open) · SOAP (mixed conifer) · TEAK (dense fir) · HJA (tall PNW conifer) |
| Density rung (pts/m²) | 20 · 8 · 4 · 2 · 1 (via decimation from the densest available per site) |
| CHM resolution (m) | 0.25 · 0.5 · 1.0 (never finer than pulse spacing) |
| Detection (VWF `ws = a·h + b`, clamped) | slope `a` ∈ {0.05, 0.10, 0.15}; min/max window ∈ {3 m … 5 m}; `hmin = 2 m` |
| Segmentation | `dalponte2016` (`th_seed` ∈ {0.40,0.45,0.50}, `th_cr` ∈ {0.55,0.65}) · marker-controlled watershed · `silva2016` |

**Combinatorics warning.** A full factorial is large. Use a **staged design**: (1) coarse screen on one calibration site to prune dominated settings, (2) refine the survivors across all sites and density rungs. Keep `pit_fill`/`subcircle` fixed early (e.g. `subcircle ≈ 0.2 m`), then sensitivity-test once the rest is settled.

### 2. Evaluation protocol

- **Split calibration from validation plots** per forest type — tune on a subset, report on held-out plots. Otherwise the result is the sweep's training optimum.
- **Stem-level matching by crown class:** 2D distance tolerance + height tolerance for treetop→stem matching (Hungarian/greedy nearest-match within tolerance). Report **recall / precision / F1**, plus crown height and crown-diameter RMSE.
- **Stratify omission/commission** for dominant/co-dominant vs intermediate/suppressed. The suppressed-tree omission is the figure that moves least with parameters and most with density — that contrast is the most useful result the sweep can produce.
- Account for **plot-level positional uncertainty** in the matching tolerance (NEON stems are not survey-grade points).

### 3. Reproducible skeleton (`lasR` / `lidR`)

```r
library(lidR); library(lasR); library(future)
plan(multisession)

# ---- 0. Density ladder from a dense source (e.g. SOAP 2021 ~20 pts/m^2) ----
ctg  <- readLAScatalog("soap_2021/")
opt_chunk_buffer(ctg) <- 25          # >= largest crown radius; avoids seam errors
rungs <- c(20, 8, 4, 2, 1)           # target points/m^2

decimate_to <- function(ctg, target) {
  decimate_points(ctg, homogenize(density = target, res = 5))
}

# ---- 1. Per-rung pipeline ----
run_rung <- function(ctg_d, chm_res, ws_fun, seg) {
  norm <- normalize_height(classify_ground(ctg_d, csf()), tin())
  chm  <- rasterize_canopy(norm, res = chm_res,
            pitfree(thresholds = c(0,10,20), max_edge = c(0,1.5), subcircle = 0.2))
  ttops <- locate_trees(norm, lmf(ws = ws_fun, hmin = 2, shape = "circular"))
  seg_f <- switch(seg,
             dalponte = dalponte2016(chm, ttops, th_seed = 0.45, th_cr = 0.55),
             watershed = watershed(chm),
             silva    = silva2016(chm, ttops))
  list(ttops = ttops, trees = segment_trees(norm, seg_f), chm = chm)
}

# ---- 2. Variable-window allometry (Popescu & Wynne), clamped ----
ws_factory <- function(a, lo = 3, hi = 5) {
  function(h) { y <- a*h + lo; y[h < 2] <- lo; y[h > 20] <- hi; y }
}

# ---- 3. Sweep grid (stage 1: coarse) ----
grid <- expand.grid(rung = rungs,
                    chm_res = c(0.5, 1.0),
                    a = c(0.05, 0.10, 0.15),
                    seg = c("dalponte", "watershed"),
                    stringsAsFactors = FALSE)

results <- lapply(seq_len(nrow(grid)), function(i) {
  g  <- grid[i, ]
  cd <- decimate_to(ctg, g$rung)
  out <- run_rung(cd, g$chm_res, ws_factory(g$a), g$seg)
  score_against_field(out$ttops, field_stems, crown_class = TRUE)  # custom scorer
})
```

> The variable-window-as-a-function and native `summarise()$pulse_density` features live on the `lasR` development branch — pin that build if you depend on them, and verify the multi-threading API in the `lasR3` vignette before scripting the parallel config. For pure header-level density estimates, `lidR::density(ctg)` is cheaper than a streaming pass.

For matching/scoring, the `NeonTreeEvaluation` R package can serve as the harness on the California sites (IoU/centroid matching against hand-annotations); extend it with the crown-class stratification above.

### 4. Suggested repo structure

```
lidar_tree_benchmarks/
├── README.md
├── dataset-research-and-sweep-plan.md      # this document
├── data/                                    # download scripts / DOIs, not raw LiDAR
│   ├── neon_download.R
│   └── hja_plots_notes.md
├── R/
│   ├── density_ladder.R
│   ├── pipeline.R
│   ├── sweep.R
│   └── scoring.R
├── config/sweep_grid.yaml
└── results/                                 # metrics tables + figures (stratified)
```

### 5. Deliverables that make the sweep defensible

- Metric tables **per (site × density × parameter set)**, with recall/precision/F1 and crown/height RMSE **stratified by crown class**.
- Calibration-vs-held-out separation documented.
- A density-sensitivity figure per forest type (the core result), plus a note on decimated-vs-native equivalence at QL2.

---

<a name="references"></a>
## References & data access

**Datasets / portals**
- NEON Data Portal — woody vegetation structure `DP1.10098.001`, discrete-return LiDAR `DP1.30003.001`: https://data.neonscience.org
- NEON LiDAR collection / sensors: https://www.neonscience.org/data-collection/lidar
- NEON D17 (Pacific Southwest, CA): https://nationaldataplatform.org/catalog/dataset/neon-d17-pacific-southwest-california1
- NeonTreeEvaluation: https://github.com/weecology/NeonTreeEvaluation — Zenodo `10.5281/zenodo.5914554`
- IDTReeS 2020: Zenodo `10.5281/zenodo.3934932`; https://idtrees.org
- HJ Andrews data catalog: https://andlter.forestry.oregonstate.edu (LiDAR `GI010`, `GI011`; stem map `TV081`)
- TreeScope: https://treescope.org/data_overview — arXiv `2310.02162`
- FOR-instance: https://zenodo.org/records/8287792 — arXiv `2309.01279`
- Open Forest Observatory: https://openforestobservatory.org/data/ ; photogrammetry workflow: https://openforestobservatory.org/workflows/photogrammetry/
- US Forest Service FIA Spatial Data Services: https://research.fs.usda.gov/programs/fia/sds

**Papers**
- LaRue et al. (2022), *Ecosphere* 13:e4209 — point-density sensitivity / ~8 pts/m² recommendation: https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1002/ecs2.4209
- NEON-SD structural-diversity product (2024), *Nature Scientific Data*: https://www.nature.com/articles/s41597-024-04018-0
- Weinstein et al. (2021), *PLOS Comp. Bio.* `10.1371/journal.pcbi.1009180` (NeonTreeEvaluation)
- Young et al. (2022), *Methods Ecol. Evol.* `10.1111/2041-210X.13860` (OFO / drone individual-tree mapping): https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/2041-210X.13860
- Dubrovin et al. (2024), *Sci. Reports* `10.1038/s41598-024-72669-5` (UAV-LiDAR + field inventory template)
- Coulston et al. — effect of blurred FIA plot coordinates: https://www.srs.fs.usda.gov/pubs/ja/ja_coulston005.pdf
- ForestFormer3D / FOR-instanceV2 (2025), arXiv `2506.16991` (announced, not yet released)

**Method references (from the pipeline document)**
Popescu & Wynne (2004, VWF); Li et al. (2012); Dalponte & Coomes (2016); Silva et al. (2016); Khosravipour et al. (2014, pit-free CHM); Jakubowski et al. (2013); Sparks et al. (2022); Hamraz et al. (2017); Kaartinen et al. (2012); Eysn et al. (2015); Wang et al. (2016). Software: lidR (`r-lidar.github.io/lidRbook`), lasR (`r-lidar.github.io/lasR`), ForestTools, FUSION, PDAL, WhiteboxTools.

---

*Prepared as a planning document for the `lidar_tree_benchmarks` repository.*
