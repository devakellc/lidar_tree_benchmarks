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
- [`docs/neon-lidar-sites.md`](docs/neon-lidar-sites.md)
  — NEON benchmark sites (SJER, SOAP, TEAK): field-ground-truth counts, NEON AOP
  LiDAR acquisition dates, and USGS 3DEP EPT cross-check projects.
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
  — five CHM crown segmenters scored on crown-diameter RMSE vs NEON field widths,
  plus a SOAP-only TreeisoNet `treeOff` crown arm unioned by
  `analyze_crown_metrics.R` (issue #7 / #20).
- [`results/ept-acquisition-sweep-results.md`](results/ept-acquisition-sweep-results.md)
  — lasR remote-EPT acquisition parameter sweep on USGS 3DEP (throughput, not
  detection accuracy).
- [`results/model-benchmark-results.md`](results/model-benchmark-results.md)
  — cross-model density-ladder synthesis (#R10): AMS3D, lmfauto/multichm/ptrees,
  CHM-VWF, TreeisoNet, SegmentAnyTree, native Li 2012, and the native+8
  ForestFormer3D comparison on shared frozen clips, by crown class and height
  band, with head-to-head deltas vs CHM-VWF.
- [`results/crown-allometry-results.md`](./results/crown-allometry-results.md)
  — crown-width + height → DBH / biomass (#S1): how well each segmenter's crowns
  predict field DBH (per class, per rung) + derived AGB. Height dominates DBH;
  crown skill is decoupled from detection F1 (AMS3D best despite worst F1).

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

Most entry points live under [`scripts/`](scripts/); GPU prerequisites and tests
are listed where they belong. The analysis scripts expect an environment
variable `CLAUDE_JOB_DIR` pointing at a **working directory** (where `aoi.laz`
lives and outputs are written) — set it to any folder:

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
| `crown_metrics_sweep.R` | Issue #7: seed 5 crown segmenters (dalponte2016, silva2016, marker-free watershed, lasR region_growing, random walker) from shared NEON tree-tops; score crown-diameter RMSE vs field `maxCrownDiameter`/`ninetyCrownDiameter` by crown class. Issue #33 adds the `RUNGS=native,8,4,2,1` density ladder on the same frozen clips (`rung` column + RMSE/bias-vs-density PNGs). |
| `detect_treeisonet_crowns.R` + `analyze_crown_metrics.R` | TreeisoNet `treeOff` crown arm (#20): the GPU offset net (`gpu/run_treeisonet_crowns.py`) per SOAP plot → per-point instances → `crown_diameter_table` → matched-tree crown-diameter RMSE; `analyze_crown_metrics.R` unions it (SOAP-only) with the #7 CHM segmenters into [`crown-segmentation-results.md`](results/crown-segmentation-results.md). |
| `crown_allometry.R` + `allometry_lib.R` | Crown → DBH / biomass allometry (#S1): joins each matched crown's `d_eq` + field height to field `stemDiameter`/`taxonID`, fits crown-geometry → DBH models per segmenter/class/rung (R²/RMSE/bias), and derives AGB (Jenkins 2003 generic). Pure helpers (`functional_type`, `agb_from_dbh`, `fit_stats`) in `allometry_lib.R`. Writes `neon/<SITE>/crown_allometry.csv` behind [`crown-allometry-results.md`](./results/crown-allometry-results.md). |
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
| `calval_split.R` + `calval_lib.R` | Calibration/validation split (issue #3): tune `(chm_res, vwf_a)` per density rung on a stratified calibration subset, report held-out F1; multi-seed robustness. `calval_lib.R` holds the shared split (`plot_table`, `assign_split`) reused by the multichm arm. |
| `calval_multichm.R` | multichm cal/val (issue #39): on the SAME stratified split as `calval_split.R`, pool multichm vs CHM-VWF (matched-discipline + calib-tuned) on held-out plots per rung; paired complete-case, multi-seed win-fraction; verdict on whether multichm's SOAP/TEAK advantage survives out-of-sample. Reads cached `multichm_sweep_results.csv` + `sweep_results.csv`. |
| `ept_discovery.R` | Find public USGS 3DEP EPT projects covering each NEON site (point-in-polygon vs the entwine boundary index); writes `neon/<SITE>/ql2/ept_candidates.csv`. |
| `neon_download_lidar.R` | Fetch NEON DP1.30003.001 LiDAR tiles overlapping a site/plot set; populates `neon/<SITE>/lidar/` for the density-ladder runs. |
| `verify_geolocation.R` | Audit stem coordinates by re-deriving NEON plot/stem offsets from the API and comparing them with `ground_truth_stems.csv`. |
| `native_ql2_crosscheck.R` | Pull the native 3DEP cloud per NEON plot via PDAL (reproject 3857 -> UTM 11N), run BOTH the CHM-VWF and multichm pipelines (issues #4, #39), and compare each detector's native (+ decimated-to-2) detection to its OWN cached decimated-2 rung by crown class. |
| `bench_lasr_ept_acquisition.R` / `sweep_lasr_ept_params.R` / `sweep_lasr_ept_partitions.R` | EPT acquisition benchmarks for lasR remote reads: throughput, chunking, and partition-parameter sensitivity. |
| `neon_ground_truth.R` | Build NEON field-stem ground truth (`DP1.10098.001`); writes `ground_truth_stems.csv` with `meas_year`/`dist21` for exact-year filtering. |
| `run_sweep.R` + `sweep_lib.R` | NEON density-ladder sweep: per plot x rung x `chm_res` x `vwf_a`, scored vs stems. `MEAS_YEAR=2021` restricts to exact-year stems (issue #5); use a distinct `OUT=` to keep the +/-4 yr baseline. |
| `analyze_sweep.R` / `compare_sites.R` | Pool the sweep (sum TP / sum n_ref) + figures; cross-site structure gradient SJER -> SOAP -> TEAK. |
| `export_geojson.R` | Export benchmark geography as WGS84 GeoJSON: `sites` (convex-hull footprints), `plots` (scoring-box polygons with a `swept` flag + native density, null where unswept), `stems` (field points, `is_tree`); writes the tracked `data/{sites,plots,stems}.geojson`. |
| `validate_heights.R` | Apex-vs-field height bias/RMSE at native density; `MEAS_YEAR=2021` writes a distinct `height_pairs_2021.csv` for the temporal cut. |
| `temporal_sensitivity.R` | Pool baseline (+/-4 yr) vs exact-2021 sweeps at modal params per rung; recall/precision/F1/height deltas + height-bias comparison (issue #5). |
| `detect_pc_sweep.R` + `pc_detect_lib.R` | Issue #6: point-cloud detectors (lidR lmf-on-points, Li 2012, lasR point `local_maximum`) vs the CHM-VWF baseline at native density, scored on field stems by crown class; quantifies understory-recall deltas and the occlusion floor. The three point-cloud apex extractors live in `pc_detect_lib.R` (shared with `detect_pc_ladder.R`). Writes `neon/<SITE>/pc_detect_results.csv`. |
| `detect_pc_ladder.R` | Issue #38: the point-cloud-detector **density ladder** — the same four arms as #6 (`chm_vwf`, `lidr_lmf_pc`, `lidr_li2012`, `lasr_lmax_pc`) at the **native + 8 pts/m²** rungs only on SJER + SOAP + TEAK (4/2/1 out of scope: point segmenters are noise below ~3 first-ret/m²). Uses the seeded `frozen_clip` provider so all arms score identical bytes per (plot, rung); pools by rung with the canonical `pool`/`equal_set_guard`. Writes `neon/<SITE>/pc_detect_ladder_results.csv` + `neon/pc_detect_ladder_pooled.csv` (the pooled CSV carries `cores` + `n_dropped_cells` provenance). Defaults to `CORES=1` for reproducible pooling — lasR `exec` can transiently drop the odd dense-native cell under fork at `CORES>1`, which the equal-set guard then drops for all arms. |
| `model_bench_lib.R` | Shared bridge for the model benchmark (#B2): reduce_instances, crown_diameter_table, seed_for/frozen_clip, pool/equal_set_guard, assert_detection_contract. |
| `detect_ams3d_sweep.R` | AMS3D (crownsegmentr) arm (#B1): adaptive mean-shift crowns over the density ladder, reduced to detections and scored by crown class. |
| `detect_lidrplugins_sweep.R` | lidRplugins competitor arm (#C9): lmfauto/multichm (locate_trees) + ptrees (segment_trees) vs the CHM-VWF baseline over the density ladder. |
| `detect_multichm_sweep.R` | multichm treetop arm on the canonical 3-site density ladder (#37): `lidRplugins::multichm` on the SAME `prepare_clip` lasR path as the cached CHM-VWF `sweep_results.csv` (density-derived `res`, `ws_factory(0.10)`), scored by `score_plot`. Writes `neon/<SITE>/multichm_sweep_results.csv` (one row per plot x rung). Needs only lidR + lidRplugins (CRAN lasR is fine). |
| `analyze_multichm_sweep.R` | Pool the multichm arm and put it head-to-head vs the cached CHM-VWF `sweep_results.csv` on the common (plot, rung) set (same res-rule + `a=0.10`); per-rung + crown-class + height-band tables, Δ of pooled rates, a figure, and the `density-ladder-sweep-results.md` §8 addendum fragment. |
| `detect_li2012_native.R` | Native-only Li 2012 arm (#R10): lidR `li2012` point segmenter on the native frozen clip, reduced to detections via the bridge; the point-segmenter leg of the head-to-head. Writes `neon/<SITE>/li2012_results.csv`. |
| `detect_treeisonet_sweep.R` | TreeisoNet deep-model arm (#M7): runs the headless GPU driver (`gpu/run_treeisonet.py`, cu128/sm_120) on the normalized frozen clip per plot x rung, serially (one GPU), apex-only with a local-canopy-max z-snap, at the calibrated zero-shot `CONF=0.22` default. `VOXEL` accepts either a scalar isotropic override or an anisotropic `x,y,z` vector. Writes `neon/<SITE>/treeisonet_results.csv`. See `docs/superpowers/plans/2026-06-08-gpu-arm-infra-m7-first.md`. |
| `detect_segmentanytree_sweep.R` | SegmentAnyTree deep-model arm (#M6): runs the rebuilt sm_120 Docker image on raw-with-ground frozen clips per plot x rung, reduces `PredInstance` labels to apexes, converts absolute Z to AGL with each clip DTM, and writes checkpointed `neon/<SITE>/segmentanytree_results.csv` rows. `CORES=2` is the tested RTX 5090 throughput setting. See `gpu/segmentanytree-sm120/README.md`. |
| `analyze_model_benchmark.R` | Cross-model synthesis (#R10): unions every arm on the shared frozen clips, equal-set-guards across arms, pools per (detector, rung) by crown class + height band, and writes the density-robustness figures + table fragment behind [`model-benchmark-results.md`](results/model-benchmark-results.md). |
| `compare_model_sites.R` | Cross-site structure gradient for the classical model-benchmark arms (#E11): pools per-site ams3d + lidRplugins results across SJER → SOAP → TEAK, writes `model_cross_site_summary.csv` + per-arm structure-gradient figures. |
| `gpu/setup_treeisonet_env.sh` + `gpu/mirror_weights.sh` | GPU arm prerequisites: create the pinned TreeisoNet venv, mirror weights/configs, and verify the tracked checksum manifest. |
| `tests/run_tests.R` | Unit-test harness for benchmark library code, model runners, extractors, I/O helpers, pooling guards, and synthesis helpers. |

## Reproduce

Requirements: R with `lasR` (>= 0.21, dev/`pre-devel` build with EPT parallel
acquisition and variable-window `ws`) and `lidR`; PDAL (>= 2.9) for the EPT
extraction.

The model-benchmark arms add two more: **crownsegmentr** (CRAN,
`install.packages("crownsegmentr")`) for the AMS3D arm; and **lidRplugins**
for the competitor arm — its CRAN-archived `rgeos`/`rgdal` and Bioconductor
`EBImage` are declared but unused by the detectors we call, so install from a
source clone with those stripped from DESCRIPTION (see
`docs/superpowers/plans/2026-06-07-lidrplugins-competitor-arm.md`, Task 1).

```sh
export CLAUDE_JOB_DIR=$(pwd)/work && mkdir -p "$CLAUDE_JOB_DIR"

# Model benchmark — classical arms per site (SOAP + SJER + TEAK)
Rscript scripts/detect_ams3d_sweep.R       SITE=SOAP PLOTS=ALL CORES=12
Rscript scripts/detect_lidrplugins_sweep.R SITE=SOAP PLOTS=ALL CORES=12
Rscript scripts/detect_ams3d_sweep.R       SITE=SJER PLOTS=ALL CORES=12
Rscript scripts/detect_lidrplugins_sweep.R SITE=SJER PLOTS=ALL CORES=12
Rscript scripts/detect_ams3d_sweep.R       SITE=TEAK PLOTS=ALL CORES=12
Rscript scripts/detect_lidrplugins_sweep.R SITE=TEAK PLOTS=ALL CORES=12
# Synthesize per site, then cross-site
Rscript scripts/analyze_model_benchmark.R  SITE=SOAP
Rscript scripts/analyze_model_benchmark.R  SITE=SJER
Rscript scripts/analyze_model_benchmark.R  SITE=TEAK
Rscript scripts/compare_model_sites.R
# Native-only Li 2012 head-to-head (SOAP)
Rscript scripts/detect_li2012_native.R     SITE=SOAP PLOTS=ALL CORES=12
# Deep GPU arms per site; TreeisoNet CONF defaults to calibrated 0.22
for SITE in SOAP SJER TEAK; do
  Rscript scripts/detect_treeisonet_sweep.R        SITE=$SITE PLOTS=ALL VOXEL=0.8,0.8,2.0
  Rscript scripts/detect_segmentanytree_sweep.R    SITE=$SITE PLOTS=ALL IMAGE=sat-sm120-test CORES=4
  Rscript scripts/detect_forestformer3d_sweep.R    SITE=$SITE PLOTS=ALL REPO=<FF3D repo>
  Rscript scripts/analyze_model_benchmark.R        SITE=$SITE
done

# multichm arm on the canonical 3-site density ladder (#37): same prepare_clip
# lasR path as run_sweep.R; head-to-head vs the cached CHM-VWF sweep_results.csv
Rscript scripts/detect_multichm_sweep.R  SITE=SOAP PLOTS=ALL CORES=12
Rscript scripts/detect_multichm_sweep.R  SITE=SJER PLOTS=ALL CORES=12
Rscript scripts/detect_multichm_sweep.R  SITE=TEAK PLOTS=ALL CORES=12
Rscript scripts/analyze_multichm_sweep.R SITES=SJER,SOAP,TEAK

# Point-cloud-detector density ladder (#38): native + 8 pts/m^2 only, 3 sites.
# CORES=1 keeps the pooling reproducible (no transient lasR fork-drops); the
# script writes the per-site CSVs and prints + writes the pooled rung table.
Rscript scripts/detect_pc_ladder.R       SITES=SJER,SOAP,TEAK CORES=1

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

Run the library tests with `Rscript tests/run_tests.R`. The GPU arms need their
own prerequisites before the first run: TreeisoNet uses
`gpu/setup_treeisonet_env.sh` plus `gpu/mirror_weights.sh`, SegmentAnyTree uses
the Docker image from `gpu/segmentanytree-sm120/`, and ForestFormer3D uses the
ported runtime described under `gpu/forestformer3d-sm120/`.

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
