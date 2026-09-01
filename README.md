# LiDAR Tree Benchmarks

> Reproducible benchmarks for individual-tree detection and crown delineation
> from airborne LiDAR, tested across point densities, forest structures, and
> model families.

This research repository compares density-aware tree-top detection and crown
segmentation workflows. It starts with a small bundled LiDAR tile, extends to a
USGS 3DEP area of interest, and evaluates the principal methods on
field-mapped NEON stems at SJER, SOAP, and TEAK.

The central rule is simple: **measure density first, then derive the canopy
height model (CHM), detection-window, and smoothing parameters from it.** Do
not carry fixed settings from one acquisition into another.

![Overstory and understory recall across the three NEON benchmark
sites.](results/figures/structure_gradient.png)

*CHM-VWF recall across the NEON structure gradient. Solid lines show
overstory recall; dashed lines show understory recall. See the
[density-ladder results](results/density-ladder-sweep-results.md).*

## Start here

Read the [methodology](docs/treetop-detection-approach.md) for the workflow and
its assumptions. Then run the bundled toy comparison; it needs no download:

~~~sh
export CLAUDE_JOB_DIR="$PWD/work"
mkdir -p "$CLAUDE_JOB_DIR"

Rscript scripts/detect_lasr.R
Rscript scripts/detect_lidr.R
Rscript scripts/compare.R \
  "$CLAUDE_JOB_DIR/tops_lasr.csv" \
  "$CLAUDE_JOB_DIR/tops_lidr.csv"
Rscript scripts/shared_chm.R
~~~

This produces density-derived lasR and lidR detections, compares their
locations, and then repeats the local-maximum test on one shared CHM. The
shared-CHM run is the right way to compare the maxima implementations.

| If you want to... | Start with |
| --- | --- |
| Understand the density-first method | [Methodology](docs/treetop-detection-approach.md) |
| Compare lasR and lidR on a known small tile | [Toy workflow](#bundled-toy-tile) |
| Process a real USGS 3DEP AOI | [USGS 3DEP workflow](#usgs-3dep-aoi) |
| Reproduce field-ground-truth accuracy curves | [NEON density ladder](#neon-density-ladder) |
| Compare detector, model, or fusion arms | [Model benchmark](results/model-benchmark-results.md) |
| Find every runnable entry point | [Script reference](#script-reference) |

## What this repository benchmarks

| Area | Methods and outputs |
| --- | --- |
| Tree-top detection | CHM variable-window filtering, point-cloud local maxima, Li 2012, multichm, and learned detector arms |
| Crown delineation | CHM region growing, Dalponte, Silva, watershed, random-walker, and point/instance segmentation approaches |
| Data sources | A bundled lasR tile, a reprojected USGS 3DEP EPT AOI, and 2021 NEON airborne LiDAR with field stems |
| Density robustness | Native data plus 8, 4, 2, and 1 pts/m² density rungs where an arm is meaningful |
| Evaluation | One-to-one apex matching, crown-diameter error, calibrated confidence, uncertainty bands, and point-set IoU, Coverage, and PQ proxy metrics |
| Scale and transfer | Streaming/catalog tests, native-QL2 cross-checks, temporal sensitivity, RGB fusion, and cross-site comparisons |

The model benchmark includes classical CHM and point-cloud methods as well as
AMS3D, Treeiso, TreeisoNet, SegmentAnyTree, ForestFormer3D, DeepForest,
Detectree2, and SAM2Point where their required environments are available.

## Evidence at a glance

These are experiment-specific findings, not universal deployment guarantees.
Follow the linked result documents for datasets, protocols, and caveats.

| Question | Current finding |
| --- | --- |
| Does the engine choose different tops? | Given the same CHM, lasR and lidR local-maximum results nearly agree (Jaccard 0.95); CHM construction drives the larger end-to-end difference. |
| Does density matter? | At SOAP, CHM-VWF F1 remains about 0.35–0.42 over a 15× density range, but recall and precision exchange places and understory detection stays severely occlusion-limited. |
| When do point/instance methods help? | At native SOAP density, SegmentAnyTree improves F1 and understory recall over CHM-VWF, then becomes less competitive at the 2 and 1 pts/m² rungs. |
| Can sparse runs be simulated by decimation? | The native-QL2 cross-check measures that limitation directly rather than assuming decimation is a perfect acquisition surrogate. |
| Which tool streams an AOI cleanly? | In the nine-tile demonstration, lasR reproduced the single-file count exactly; inspect the catalog result before generalizing that result to another tiling scheme. |

- [lasR vs lidR comparison](results/treetop-lasr-vs-lidr-comparison.md)
- [Density-ladder benchmark](results/density-ladder-sweep-results.md)
- [Cross-model benchmark](results/model-benchmark-results.md)
- [Native USGS 3DEP cross-check](results/native-ql2-crosscheck-results.md)
- [EPT acquisition and streaming results](results/ept-acquisition-sweep-results.md)

## Reproduce a workflow

All scripts use CLAUDE_JOB_DIR for their working directory. It defaults to
work below the current directory, but setting it explicitly keeps inputs and
generated outputs separate from the checkout. Scripts accept KEY=VALUE
positional arguments rather than command-line flags.

### Bundled toy tile

The [quick-start commands](#start-here) run the two density-first detection
paths. Run the remaining steps to sweep parameters and create crown products:

~~~sh
Rscript scripts/sweep.R
Rscript scripts/segment_lasr.R
Rscript scripts/segment_lidr.R
Rscript scripts/compare_crowns.R \
  "$CLAUDE_JOB_DIR/crowns_lasr.gpkg" \
  "$CLAUDE_JOB_DIR/crowns_lidr.gpkg"
~~~

The bundled MixedConifer.las file is provided by the lasR installation; it is
not a field-ground-truth benchmark. It is intended for fast, reproducible
pipeline and engine checks.

### USGS 3DEP AOI

The PDAL extraction clips the public EPT and reprojects it from Web Mercator to
UTM before metric-sensitive processing. From the repository root:

~~~sh
REPO_ROOT="$(git rev-parse --show-toplevel)"
(cd "$CLAUDE_JOB_DIR" && pdal pipeline "$REPO_ROOT/scripts/extract.json")

Rscript scripts/detect_lasr_aoi.R
Rscript scripts/detect_lidr_aoi.R
Rscript scripts/shared_chm_aoi.R
Rscript scripts/pc_vs_chm.R
Rscript scripts/segment_lasr_aoi.R
Rscript scripts/segment_lidr_aoi.R
~~~

For a native lasR remote-EPT acquisition path, run
scripts/detect_lasr_ept_aoi.R instead. That path remains in EPSG:3857 and is
best used for acquisition/streaming experiments, not metric-faithful parameter
selection. See the [AOI comparison](results/treetop-lasr-vs-lidr-comparison.md).

### NEON density ladder

This is the primary field-ground-truth benchmark. It downloads the needed NEON
field and LiDAR data, creates density rungs, and scores detections in the
mapped plot cores. A full three-site run needs network access, several GB of
working storage, and meaningful compute time.

~~~sh
for SITE in SJER SOAP TEAK; do
  Rscript scripts/neon_ground_truth.R SITE="$SITE"
  Rscript scripts/neon_download_lidar.R SITE="$SITE" YEAR=2021
  Rscript scripts/run_sweep.R SITE="$SITE" PLOTS=ALL CORES=8
  Rscript scripts/analyze_sweep.R SITE="$SITE"
done
Rscript scripts/compare_sites.R
~~~

Use a distinct OUT value when running an exact-year temporal subset so that it
does not replace the default ±4-year ground-truth baseline:

~~~sh
Rscript scripts/run_sweep.R SITE=SOAP MEAS_YEAR=2021 \
  OUT="$CLAUDE_JOB_DIR/neon/SOAP/sweep_results_2021.csv"
~~~

## Requirements

| Scope | Requirements |
| --- | --- |
| Core toy and AOI workflows | R with lasR, lidR, terra, sf, and data.table |
| Required lasR build | r-lidar/lasR pre-devel; the released 0.21.0 build rejects the variable-window function used by the detection scripts |
| NEON workflows | neonUtilities and jsonlite, plus network access for public NEON products |
| EPT extraction | PDAL 2.9 or later |
| Optional analyses | clue, rpart, crownsegmentr, and lidRplugins as required by the corresponding arm |
| GPU/vision arms | The documented container or conda environment under [gpu](gpu/) for that specific model |
| Tests | testthat |

The [density-ladder setup notes](results/density-ladder-sweep-results.md)
document the pre-development lasR requirement and the feature check behind it.
The model-specific setup instructions live alongside each runtime under
[gpu](gpu/) and in the linked result documents.

There is an important engine distinction: lasR uses a TIN plus post-hoc
pit_fill step, while lidR uses its Khosravipour pitfree method. They should not
be described as the same CHM algorithm.

## Data and reproducibility

- Tracked geographic context lives in [data](data) as GeoJSON site, plot, stem,
  and AOI layers.
- Downloaded LiDAR, rasters, tables, GeoPackages, model weights, and working
  files are intentionally gitignored. Regenerate them in CLAUDE_JOB_DIR.
- The NEON scorer pools counts before calculating rates; it does not average
  plot-level recall or precision. This prevents small plots from dominating a
  site result.
- Field stems support an apex/detection benchmark. The point-set instance
  metrics use a clearly labeled Voronoi-on-stems crown proxy, not hand-drawn
  crown masks.

## Script reference

Run the driver scripts below directly. Supporting libraries are listed last;
they are sourced by drivers and covered by tests rather than run as standalone
workflows.

### Core detection and crown workflows

| Scripts | Purpose |
| --- | --- |
| detect_lasr.R / detect_lidr.R | Density-first detection on the bundled tile |
| compare.R / shared_chm.R / sweep.R | Compare top locations, isolate the same-CHM test, and sweep toy parameters |
| segment_lasr.R / segment_lidr.R / compare_crowns.R | Delineate toy crowns, calculate metrics, and compare polygons |
| detect_lasr_aoi.R / detect_lidr_aoi.R / shared_chm_aoi.R | Run and fairly compare the corresponding USGS AOI paths |
| segment_lasr_aoi.R / segment_lidr_aoi.R / pc_vs_chm.R | AOI crown products and high-density CHM-versus-point-cloud comparison |
| extract.json / extract_big.json | PDAL EPT clip and reproject pipelines |
| detect_lasr_ept_aoi.R | lasR-native remote-EPT acquisition and detection |
| tile_aoi.R / detect_lasr_catalog.R / detect_lidr_catalog.R | Retile an AOI and test multi-tile streaming behavior |
| density_cost.R / li2012_16core.R | Detection density/cost and multi-core Li 2012 throughput studies |
| bench_lasr_ept_acquisition.R / sweep_lasr_ept_params.R / sweep_lasr_ept_partitions.R | EPT throughput, parameter, and partition studies |

### NEON data and density studies

| Scripts | Purpose |
| --- | --- |
| neon_ground_truth.R / verify_geolocation.R | Build and audit field-stem ground truth |
| neon_download_lidar.R / neon_download_aop.R | Download the NEON LiDAR and RGB inputs needed by a selected arm |
| run_sweep.R / analyze_sweep.R / compare_sites.R | Run, pool, and compare the CHM-VWF density ladder |
| calval_split.R / calval_multichm.R | Held-out parameter calibration/validation |
| ept_discovery.R / native_ql2_crosscheck.R | Find covering 3DEP projects and test native-versus-decimated performance |
| validate_heights.R / temporal_sensitivity.R | Height validation and field-to-LiDAR temporal sensitivity |
| detect_pc_sweep.R / detect_pc_ladder.R | Point-cloud detector comparison at native and selected sparse rungs |
| matcher_robustness.R / mc_positional_uncertainty.R | Matching-rule and stem-position-uncertainty sensitivity |

### Model, fusion, and crown analyses

| Scripts | Purpose |
| --- | --- |
| detect_ams3d_sweep.R | Adaptive mean-shift crown segmentation arm |
| detect_lidrplugins_sweep.R / detect_multichm_sweep.R / analyze_multichm_sweep.R | lmfauto, ptrees, and multichm arms and their paired analysis |
| detect_li2012_native.R / detect_treeiso_sweep.R | Native point-cloud segmentation baselines |
| detect_treeisonet_sweep.R / detect_treeisonet_crowns.R | TreeisoNet apex and tree-offset crown arms |
| detect_segmentanytree_sweep.R / detect_forestformer3d_sweep.R | GPU point/instance segmentation arms |
| detect_deepforest_sweep.R / detect_detectree2_sweep.R | RGB detector and crown-width arms |
| detect_sam2point_sweep.R | Promptable seed-to-refine point-cloud arm |
| analyze_model_benchmark.R / compare_model_sites.R | Equal-set-guarded model synthesis and cross-site results |
| score_instances_iou.R / compare_matching_rules.R | Point-set IoU, Coverage, PQ, and metric-ranking sensitivity |
| fuse_detectors.R / calibrate_confidence.R / route_detectors.R | Detector fusion, score calibration, and per-cell routing |
| coverage_gap.R | Re-grade isolated likely-real false positives using cross-family agreement |
| crown_metrics_sweep.R / analyze_crown_metrics.R | Field crown-diameter benchmark and analysis |
| crown_allometry.R | Crown width and height to DBH/biomass analysis |

### Exports and supporting libraries

| Scripts | Purpose |
| --- | --- |
| export_geojson.R / export_stems_ground_truth_geojson.R / export_best_treetops_geojson.R | Export benchmark geography, field stems, and best detections as GeoJSON |
| bootstrap.R / repo_paths.R | Locate the repository and working directory consistently |
| sweep_lib.R / calval_lib.R / pc_detect_lib.R | Shared density-ladder, split, and point-cloud detection helpers |
| model_bench_lib.R / model_runner.R / io_bridge.R | Shared model scoring, runtime, and point-instance I/O helpers |
| route_lib.R / coverage_lib.R / allometry_lib.R | Pure helpers for routing, coverage credit, and allometry |
| crown_metrics_3d.R / crown_metrics_deepmodel.R | Shared crown-metric helpers |

## Documentation and result index

| Read this | For |
| --- | --- |
| [Tree-top detection approach](docs/treetop-detection-approach.md) | Method, parameter rules, tooling, and pitfalls |
| [NEON LiDAR sites](docs/neon-lidar-sites.md) | Site, field-stem, LiDAR, and 3DEP context |
| [Dataset and sweep plan](docs/dataset-research-and-sweep-plan.md) | Benchmark design and evaluation rationale |
| [lasR vs lidR comparison](results/treetop-lasr-vs-lidr-comparison.md) | Toy tile, AOI, same-CHM, crowns, and streaming results |
| [Density-ladder results](results/density-ladder-sweep-results.md) | Cross-density, crown-class, and site results |
| [Model benchmark](results/model-benchmark-results.md) | Classical and deep detector comparison |
| [Crown-segmentation results](results/crown-segmentation-results.md) | Field crown-width error for delineation methods |
| [Point-cloud detector results](results/pointcloud-detector-results.md) | Native-density CHM and point-cloud detector comparison |
| [Instance IoU, Coverage, and PQ](results/instance-iou-pq-results.md) | Mask-aware proxy evaluation |
| [RGB–LiDAR fusion](results/rgb-lidar-fusion-results.md) | DeepForest and Detectree2 results |

Additional targeted analyses:

- [Calibration/validation](results/calibration-validation-results.md),
  [native QL2](results/native-ql2-crosscheck-results.md), and
  [temporal sensitivity](results/temporal-sensitivity-results.md)
- [Detector fusion](results/detector-fusion-results.md),
  [confidence calibration](results/confidence-calibration-results.md), and
  [detector routing](results/detector-routing-results.md)
- [Matcher robustness](results/matcher-robustness-results.md),
  [positional uncertainty](results/positional-uncertainty-results.md), and
  [coverage-gap crediting](results/coverage-gap-results.md)
- [Crown allometry](results/crown-allometry-results.md) and
  [SAM2Point seed-to-refine](results/sam2point-promptable-refine-results.md)

## Tests

Run the library test suite from the repository root:

~~~sh
Rscript tests/run_tests.R
~~~

The tests cover the shared scoring, pooling, detector-extractor, instance-I/O,
model-runner, routing, allometry, and uncertainty helpers. They do not require
the large generated LiDAR working set.
