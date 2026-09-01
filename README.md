# LiDAR Tree Benchmarks

> Reproducible comparisons of individual-tree detection and crown-delineation
> methods from airborne LiDAR.

This research repository evaluates tree-top detectors and crown delineators
from classical CHM methods through point-cloud, deep-learning, and RGB arms.
It starts with a bundled LiDAR tile, extends to a USGS 3DEP area of interest,
and evaluates the methods against field-mapped NEON stems at SJER, SOAP, and
TEAK.

## Methods and results

### Tree-top detection

Detection metrics use one-to-one matching against mapped field stems. Values
below are pooled benchmark results; their site, density, and matched-set scope
are stated so that unlike comparisons are not implied.

| Method family | Implementations in this repository | Headline result | Best fit |
| --- | --- | --- | --- |
| CHM local maxima | lasR and lidR variable-window filtering (CHM-VWF) | With the same CHM, the engines agree closely (Jaccard 0.95). At native density, CHM-VWF reaches F1 0.365 over 699 stems. | Transparent, reproducible baseline and large-area CHM workflows. |
| Multi-layer CHM | multichm and lmfauto | multichm reaches native F1 0.412 over the same three-site population and is the most stable LiDAR arm across sparse rungs. | A robust classical default, especially when point density is modest. |
| Point-cloud detectors | lidR point LMF, Li 2012, and lasR point local maximum | Li 2012 raises understory recall to 0.257 versus 0.200 for CHM-VWF, but pooled F1 is 0.346. | Recovering additional sub-canopy trees when precision trade-off is acceptable. |
| Point/instance models | SegmentAnyTree, TreeisoNet, ForestFormer3D, AMS3D, ptrees, and Treeiso | SegmentAnyTree is the strongest native single arm in the pooled comparison: F1 0.443 and understory recall 0.448. | High-density ALS where compute and model setup are available. |
| RGB detection | DeepForest and Detectree2 | DeepForest achieves F1 0.368 and understory recall 0.348 at SOAP; Detectree2 transfers poorly without domain-matched weights (F1 0.255). | An optical complement or low-density fallback, not the best standalone LiDAR replacement. |
| Multi-detector fusion | Union, majority, layered, and k-of-N consensus | Union lifts pooled recall to 0.764 and understory recall to 0.600, but F1 falls to 0.349 because the stem reference is incomplete for many isolated detections. | Recall-oriented inventories and a calibrated operating-point frontier. |

### Crown delineation

The crown benchmark compares predicted crown diameter with NEON field diameter.
Equivalent-circle diameter (d_eq) is compared with ninetyCrownDiameter, while
max-caliper diameter is compared with maxCrownDiameter. Those are different
geometric targets, so their RMSE values should not be mixed.

| Method | Crown representation | Headline result | Interpretation |
| --- | --- | --- | --- |
| Random walker with per-crown stop rule | CHM regions stopped at a fraction of seed height | Lowest pooled RMSE in the shared-seed classical test: 2.42 m for d_eq and 3.57 m for max-caliper. | Current best classical diameter fit when its Matrix-based workflow is available. |
| lasR region growing | CHM regions converted to polygons | 2.62 m d_eq RMSE and 3.72 m max-caliper RMSE over 225 matched stems. | A strong, compact, engine-native CHM baseline. |
| Dalponte and Silva | Seeded CHM segmenters | Dalponte: 2.70 m d_eq RMSE; Silva: 2.79 m. | Competitive alternatives; seed source and diameter definition matter. |
| AMS3D and ptrees | 3-D point-instance crowns | d_eq RMSE 2.97 m and 3.03 m, respectively; they match 1.4–2.1× more stems than the CHM controls. | Better sub-canopy reach, with different matched populations and a modest RMSE trade-off. |
| Li 2012 | 3-D point-cloud segments | 5.23 m d_eq RMSE and substantial positive diameter bias. | Useful for detection coverage, not the current crown-width choice. |

The [crown-segmentation results](results/crown-segmentation-results.md) include
matched-tree counts, bias, MAE, R², density sensitivity, and the full
three-dimensional comparison. Crown methods should be chosen on that crown
metric, not on the F1 of the detector that supplied their seeds.

### Choosing a method

| Goal | Recommended starting point | Why |
| --- | --- |
| Simple, inspectable tree-top baseline | CHM-VWF | The lasR/lidR same-CHM test shows the peak finder is not the material engine difference. |
| Reliable classical detector at varied density | multichm | It is the strongest stable classical LiDAR arm in the benchmark. |
| Maximum high-density detection F1 | SegmentAnyTree | Best pooled native F1 and understory recall among the evaluated single arms. |
| More understory trees | SegmentAnyTree, Li 2012, or a fusion operating point | These methods raise coverage; choose the precision/recall point explicitly. |
| Optical or sparse-LiDAR complement | DeepForest | It supplies density-independent RGB coverage and a different failure mode. |
| Best classical crown diameter | Random walker with the per-crown stop rule | It has the lowest pooled crown-diameter RMSE in the shared-seed test. |
| Straightforward CHM crown product | lasR region growing or Dalponte | They are close in diameter accuracy and easier to inspect and reproduce. |

Detailed evidence:

- [Cross-model detection benchmark](results/model-benchmark-results.md)
- [Classical and 3-D crown benchmark](results/crown-segmentation-results.md)
- [Point-cloud detector comparison](results/pointcloud-detector-results.md)
- [RGB detector results](results/rgb-lidar-fusion-results.md)
- [Detector-fusion results](results/detector-fusion-results.md)
- [lasR versus lidR implementation comparison](results/treetop-lasr-vs-lidr-comparison.md)

## Benchmark design

The methods are evaluated on a bundled tile, a USGS 3DEP AOI, and 2021 NEON
LiDAR paired with field-mapped stems from three contrasting forest structures.
The NEON experiment holds data preparation and matching rules constant, then
tests methods over native, 8, 4, 2, and 1 pts/m² rungs where a method remains
meaningful. Density is therefore an evaluation condition, not the subject of
the repository.

- CHM resolution, local-maximum window, and smoothing are derived from measured
  density rather than copied as fixed parameters between acquisitions.
- Detection uses one-to-one apex matching; crown delineation uses field crown
  widths; instance metrics use a clearly labelled Voronoi-on-stems proxy rather
  than hand-drawn reference masks.
- Results pool counts or error sums before calculating rates and RMSE, so small
  plots do not dominate a site-level result.

Read the [methodology](docs/treetop-detection-approach.md) for the parameter
rules and [NEON site notes](docs/neon-lidar-sites.md) for the data context.

## Start here

Run the bundled comparison to inspect the two CHM local-maximum implementations.
It needs no data download:

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

| If you want to... | Start with |
| --- | --- |
| Compare detector families and results | [Methods and results](#methods-and-results) |
| Reproduce crown-diameter metrics | [Crown workflow](#bundled-toy-tile) and [crown benchmark](results/crown-segmentation-results.md) |
| Process a real USGS 3DEP AOI | [USGS 3DEP workflow](#usgs-3dep-aoi) |
| Reproduce the field benchmark | [NEON method benchmark](#neon-method-benchmark) |
| Find every runnable entry point | [Script reference](#script-reference) |

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

### NEON method benchmark

This is the primary field-ground-truth benchmark. It downloads the needed NEON
field and LiDAR data, evaluates the CHM-VWF detector over density rungs, and
creates the data products consumed by the other method arms. A full three-site
run needs network access, several GB of working storage, and meaningful
compute time.

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

The [lasR setup notes](results/density-ladder-sweep-results.md) document the
pre-development build requirement and the feature check behind it.
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

### Field data and validation studies

| Scripts | Purpose |
| --- | --- |
| neon_ground_truth.R / verify_geolocation.R | Build and audit field-stem ground truth |
| neon_download_lidar.R / neon_download_aop.R | Download the NEON LiDAR and RGB inputs needed by a selected arm |
| run_sweep.R / analyze_sweep.R / compare_sites.R | Run, pool, and compare the core CHM-VWF field benchmark |
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
