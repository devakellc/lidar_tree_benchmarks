# lidar_tree_benchmarks

Tree-top detection and crown segmentation from airborne LiDAR, with a
density-first workflow and a reproducible **lasR vs lidR** comparison on both a
bundled toy tile and a real USGS 3DEP AOI.

## Documents

Methodology and planning live in [`docs/`](docs/); every results write-up lives
in [`results/`](results/).

- [`docs/treetop-detection-approach.md`](docs/treetop-detection-approach.md)
  — the recommended, density-driven pipeline (pit-free CHM -> variable-window
  local-maximum -> segmentation), tooling, accuracy expectations, and pitfalls.
- [`results/treetop-lasr-vs-lidr-comparison.md`](results/treetop-lasr-vs-lidr-comparison.md)
  — implementation and head-to-head results: toy tile, real 25 ha 3DEP AOI, a
  controlled same-CHM test, a CHM-vs-point-cloud test at high density, crown
  segmentation/metrics (Steps 6-7), and a multi-tile streaming demo.
- [`results/density-ladder-sweep-results.md`](results/density-ladder-sweep-results.md)
  — the NEON density-ladder parameter sweep: detection accuracy vs point
  density, stratified by crown class, on three structure-gradient sites.
- [`results/calibration-validation-results.md`](results/calibration-validation-results.md)
  — held-out calibration/validation split of that sweep: out-of-sample F1 for
  the per-rung best `(chm_res, vwf_a)` (issue #3).
- [`results/native-ql2-crosscheck-results.md`](results/native-ql2-crosscheck-results.md)
  — native USGS 3DEP cross-check of the decimation-as-simulation caveat, by
  crown class (issue #4).
- [`results/temporal-sensitivity-results.md`](results/temporal-sensitivity-results.md)
  — how much the +/-4 yr field-to-LiDAR temporal slack moves recall/precision
  and apex-height bias (exact-2021 re-score; issue #5).
- [`results/pointcloud-detector-results.md`](results/pointcloud-detector-results.md)
  — point-cloud detectors (Li 2012 / lmf) vs CHM-VWF at native density:
  understory-recall deltas and the occlusion floor (issue #6).
- [`results/crown-segmentation-results.md`](results/crown-segmentation-results.md)
  — five crown segmenters scored on crown-diameter RMSE vs NEON field widths,
  by crown class (issue #7).
- [`results/ept-acquisition-sweep-results.md`](results/ept-acquisition-sweep-results.md)
  — lasR remote-EPT acquisition parameter sweep on USGS 3DEP (throughput, not
  detection accuracy).

## Headline findings

- Given the **same CHM**, lasR and lidR local-maximum detectors are effectively
  **identical** (Jaccard 0.95–1.0); the CHM construction drives almost all of
  the end-to-end difference, not the maxima search.
- The CHM resolution, **detection window, and Step-5 smoothing** should all be
  **derived from measured density** (Step 0), not hardcoded. The density-tiered
  3x3 mean smooth before LM (only when density < 8) cuts the toy lasR count by
  ~40% because the underlying `pit_fill` CHM is bumpier than `pitfree`.
- **Crowns (Steps 6-7)** are produced as GPKG polygons: lasR
  `region_growing` + `terra::as.polygons()` for area / apex height; lidR
  `dalponte2016` + `crown_metrics(.stdtreemetrics)` for convex-hull crowns.
- At high density (~14 first-returns/m²) point-cloud **segmentation** (Li 2012)
  recovers ~25% more trees than a CHM — almost all sub-dominant/regen (< 5 m)
  that a 2.5D CHM cannot represent. lasR has point-cloud `local_maximum` but
  no point-cloud segmenter, so that step is lidR/PDAL-only.
- **Multi-tile streaming demo:** lasR auto-buffering reproduces the
  single-file AOI tree count **exactly** across 9 retiled chunks (6,809 =
  6,809), while the lidR LAScatalog path **over-counts by ~13%** at the same
  `opt_chunk_buffer = 20 m`. Prefer lasR streaming for wall-to-wall work.

## Scripts

All under [`scripts/`](scripts/). They expect an environment variable
`CLAUDE_JOB_DIR` pointing at a **working directory** (where `aoi.laz` lives and
outputs are written) — set it to any folder:

```sh
export CLAUDE_JOB_DIR=/path/to/workdir
```

| Script | What it does |
|--------|--------------|
| `detect_lasr.R` / `detect_lidr.R` | Density-first detection on the bundled `MixedConifer.las` (toy). |
| `compare.R` | Spatial matching between two treetop CSVs. |
| `shared_chm.R` | Controlled test: both detectors on one shared CHM. |
| `sweep.R` | Parameter sweep vs the bundled `treeID` reference. |
| `segment_lasr.R` / `segment_lidr.R` | Steps 6-7 on the toy: region-growing / dalponte2016 + crown polygons + metrics. |
| `compare_crowns.R` | Spatial matching + per-pair IoU between two crown GPKGs. |
| `crown_metrics_sweep.R` | Issue #7: seed 5 crown segmenters (dalponte2016, silva2016, marker-free watershed, lasR region_growing, random walker) from shared NEON tree-tops; score crown-diameter RMSE vs field `maxCrownDiameter`/`ninetyCrownDiameter` by crown class. |
| `extract.json` | PDAL pipeline: clip the AOI from the public EPT, reproject 3857 -> UTM 10N, write `aoi.laz`. |
| `detect_lasr_ept_aoi.R` | lasR-native remote EPT AOI pipeline (acquire + process directly in lasR). |
| `detect_lasr_aoi.R` / `detect_lidr_aoi.R` | Full approach on the real 3DEP AOI after PDAL extraction. |
| `shared_chm_aoi.R` | Same-CHM controlled test on the AOI. |
| `segment_lasr_aoi.R` / `segment_lidr_aoi.R` | Steps 6-7 on the AOI: crown polygons + metrics. |
| `pc_vs_chm.R` | CHM-lmf vs point-cloud lmf vs Li 2012 on a sub-clip. |
| `tile_aoi.R` | Retile `aoi.laz` into a tile grid under `tiles/` for the catalog demo. |
| `detect_lasr_catalog.R` / `detect_lidr_catalog.R` | Multi-tile streaming demo (lasR auto-buffering vs `opt_chunk_buffer`). |
| `density_cost.R` | The three detectors vs density: counts + runtime scaling. |
| `extract_big.json` | PDAL: pull a larger ~56 ha block from the EPT (reprojected). |
| `li2012_16core.R` | Retile + 16-core Li 2012 throughput; extrapolates to 10k acres. |
| `calval_split.R` | Calibration/validation split (issue #3): tune `(chm_res, vwf_a)` per density rung on a stratified calibration subset, report held-out F1; multi-seed robustness. |
| `ept_discovery.R` | Find public USGS 3DEP EPT projects covering each NEON site (point-in-polygon vs the entwine boundary index); writes `neon/<SITE>/ql2/ept_candidates.csv`. |
| `native_ql2_crosscheck.R` | Pull the native 3DEP cloud per NEON plot via PDAL (reproject 3857 -> UTM 11N), run the CHM-VWF pipeline, and compare native (+ decimated-to-2) detection to the cached decimated-2 rung by crown class (issue #4). |
| `neon_ground_truth.R` | Build NEON field-stem ground truth (`DP1.10098.001`); writes `ground_truth_stems.csv` with `meas_year`/`dist21` for exact-year filtering. |
| `run_sweep.R` + `sweep_lib.R` | NEON density-ladder sweep: per plot x rung x `chm_res` x `vwf_a`, scored vs stems. `MEAS_YEAR=2021` restricts to exact-year stems (issue #5); use a distinct `OUT=` to keep the +/-4 yr baseline. |
| `analyze_sweep.R` / `compare_sites.R` | Pool the sweep (sum TP / sum n_ref) + figures; cross-site structure gradient SJER -> SOAP -> TEAK. |
| `validate_heights.R` | Apex-vs-field height bias/RMSE at native density; `MEAS_YEAR=2021` writes a distinct `height_pairs_2021.csv` for the temporal cut. |
| `temporal_sensitivity.R` | Pool baseline (+/-4 yr) vs exact-2021 sweeps at modal params per rung; recall/precision/F1/height deltas + height-bias comparison (issue #5). |
| `detect_pc_sweep.R` | Issue #6: point-cloud detectors (lidR lmf-on-points, Li 2012, lasR point `local_maximum`) vs the CHM-VWF baseline at native density, scored on field stems by crown class; quantifies understory-recall deltas and the occlusion floor. Writes `neon/<SITE>/pc_detect_results.csv`. |

## Reproduce

Requirements: R with `lasR` (>= 0.21, dev/`pre-devel` build with EPT parallel
acquisition and variable-window `ws`) and `lidR`; PDAL (>= 2.9) for the EPT
extraction.

```sh
export CLAUDE_JOB_DIR=$(pwd)/work && mkdir -p "$CLAUDE_JOB_DIR"

# Toy tile (no data download needed; uses lasR's bundled MixedConifer.las)
Rscript scripts/detect_lasr.R
Rscript scripts/detect_lidr.R
Rscript scripts/compare.R "$CLAUDE_JOB_DIR/tops_lasr.csv" "$CLAUDE_JOB_DIR/tops_lidr.csv"
Rscript scripts/shared_chm.R                 # same-CHM controlled test
Rscript scripts/sweep.R                      # parameter sweep vs reference
Rscript scripts/segment_lasr.R               # Step 6/7 crowns
Rscript scripts/segment_lidr.R
Rscript scripts/compare_crowns.R "$CLAUDE_JOB_DIR/crowns_lasr.gpkg" "$CLAUDE_JOB_DIR/crowns_lidr.gpkg"

# Real AOI with lasR-only EPT acquisition (pre-devel; auto-partitions into ~32 chunks)
Rscript scripts/detect_lasr_ept_aoi.R

# Real AOI: fetch + reproject ~25 ha with PDAL, then run the full pipeline
(cd "$CLAUDE_JOB_DIR" && pdal pipeline "$OLDPWD/scripts/extract.json")
Rscript scripts/detect_lasr_aoi.R
Rscript scripts/detect_lidr_aoi.R
Rscript scripts/shared_chm_aoi.R
Rscript scripts/pc_vs_chm.R
Rscript scripts/segment_lasr_aoi.R           # AOI crowns
Rscript scripts/segment_lidr_aoi.R
Rscript scripts/compare_crowns.R             # defaults to crowns_*_aoi.gpkg

# Multi-tile streaming demo (approach §3 edge handling)
Rscript scripts/tile_aoi.R                   # retile aoi.laz under work/tiles/
Rscript scripts/detect_lasr_catalog.R
Rscript scripts/detect_lidr_catalog.R
```

Data (`*.laz`, `*.tif`, `*.csv`, `*.gpkg`, and `tiles/`) is gitignored —
regenerate it with the steps above.

## Notes

- The EPT is in EPSG:3857 (Web Mercator); distances there are inflated ~1.32x at
  this latitude. The lasR-only EPT script keeps processing in 3857 for a fully
  native acquisition path; for metric-faithful parameterization, prefer the PDAL
  reproject-to-UTM path before detection.
- Runtimes in the comparison are **not** an engine benchmark (single small
  tiles, and EPT reads from a non-AWS machine are network-bound). lasR's real
  advantage is large-area (>= 100 km²) streaming throughput and low memory.
- The NEON ground-truth builder (`neon_ground_truth.R`) now carries two field
  crown-width columns from `vst_apparentindividual` —
  `maxCrownDiameter` (widest axis) and `ninetyCrownDiameter` (equivalent width)
  — in `ground_truth_stems.csv`. `crown_metrics_sweep.R` scores delineated crown
  diameter against these; see [`crown-segmentation-results.md`](results/crown-segmentation-results.md)
  (issue #7).
