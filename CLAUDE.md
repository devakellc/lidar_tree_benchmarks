# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with
code in this repository.

## What this repo is

A **research/benchmark** repository (not an application): R scripts that detect
tree-tops and delineate crowns from airborne LiDAR, then quantify how detection
accuracy responds to point density. There is no build system, package, or test
suite — every script is a standalone `Rscript` run, and the deliverables are the
markdown result documents under [`results/`](results/) plus CSVs/figures
regenerated under the working directory.

The central methodology lives in [treetop-detection-approach.md](docs/treetop-detection-approach.md):
**measure density first, then derive parameters from it** (Step 0). CHM
resolution, the variable-window-filter (VWF) `ws` function, and the Step-5
pre-LM smoothing are all functions of measured first-return density — never
hardcoded. Read that doc before changing any detection parameter; the result
docs ([treetop-lasr-vs-lidr-comparison.md](results/treetop-lasr-vs-lidr-comparison.md),
[density-ladder-sweep-results.md](results/density-ladder-sweep-results.md)) record the
findings each script produced and must stay consistent with the code.

## Running scripts

Every script reads a **`CLAUDE_JOB_DIR`** environment variable — the working
directory where input `*.laz` lives and all outputs are written. It defaults to
`./work` when unset.

```sh
export CLAUDE_JOB_DIR=$(pwd)/work && mkdir -p "$CLAUDE_JOB_DIR"
Rscript scripts/detect_lasr.R          # toy tile; no data download
```

Scripts that take parameters parse **`KEY=VALUE` positional args** (not flags),
e.g. `Rscript scripts/run_sweep.R SITE=SOAP PLOTS=ALL CORES=8 TOL=4.0`. The
parsing idiom is `strsplit(commandArgs(TRUE), "=")` → a named list `A`.

All generated data (`*.laz`, `*.las`, `*.tif`, `*.csv`, `*.gpkg`, `tiles/`,
`work/`) is **gitignored** — regenerate it; never assume it is present. The
README's script table is the canonical per-script reference for the toy/AOI
pipelines.

## Dependencies (critical)

- **lasR must be the `r-lidar/lasR@pre-devel` build**, not CRAN. The scripts use
  `ws`-as-a-function (variable-window `local_maximum`) and parallel EPT
  acquisition, which CRAN 0.21.0 broke. Installing CRAN lasR will make detection
  scripts fail at the `local_maximum` call. See `memory/lasr-pre-devel-dependency.md`.
- **lidR** (clip/decimate/normalize/segmentation/catalog), **terra**, **sf**,
  **data.table**; **neonUtilities** + **jsonlite** for NEON ground truth;
  **PDAL ≥ 2.9** CLI for EPT extraction (`scripts/extract*.json`).
- Markdown is linted with **rumdl** ([.rumdl.toml](.rumdl.toml)): 80-char prose
  line limit, tables and code blocks exempt. Keep result docs within it.

## Engine split (lasR vs lidR)

This is the architectural backbone — most scripts come in `_lasr` / `_lidr`
pairs that implement the *same* pipeline two ways so the engines can be matched:

- **lasR** builds CHMs (`triangulate` → `rasterize` → `pit_fill`), does
  raster + point-cloud `local_maximum`, region-growing crowns, and streams
  wall-to-wall via auto-buffering. It has **no point-cloud segmenter**.
- **lidR** does clip/decimate/normalize, `dalponte2016` CHM segmentation, **Li
  2012 point-cloud segmentation** (the only path to sub-canopy/regen trees), and
  `LAScatalog` streaming with `opt_chunk_buffer`.

Two findings to preserve when editing: (1) given the *same* CHM the two
local-maximum detectors are near-identical — the CHM construction drives the
difference, so `shared_chm*.R` isolate that; (2) lasR's `pit_fill` is a
TIN+post-hoc fill, **not** the Khosravipour pit-free algorithm in
`lidR::pitfree()` — don't describe them as equivalent.

## The NEON density-ladder sweep (the active workstream)

This is the largest pipeline and spans several files. Flow:

1. [scripts/neon_ground_truth.R](scripts/neon_ground_truth.R) — builds field-stem
   ground truth from NEON woody-veg (DP1.10098.001). Reimplements
   `geoNEON::getLocTOS()` against the public NEON locations API (polar offset
   from named grid points); classifies crown class; writes
   `ground_truth_stems.csv` + `plot_centroids.csv` under `work/neon/<SITE>/`.
2. [scripts/run_sweep.R](scripts/run_sweep.R) + [scripts/sweep_lib.R](scripts/sweep_lib.R)
   — for each plot, decimates to a **density rung** (8/4/2/1 pts/m² + native),
   runs `detect_lasr` over a grid of `chm_res` × `vwf_a`, and scores against
   stems. Output is **long-form**: one row per (plot × rung × chm_res × vwf_a).
3. [scripts/analyze_sweep.R](scripts/analyze_sweep.R) — per-site pooling + figures.
4. [scripts/compare_sites.R](scripts/compare_sites.R) — cross-site structure
   gradient SJER → SOAP → TEAK.

Sweep invariants — get these wrong and the metrics are silently misleading:

- **Pool by summing counts, not averaging rates**: site/rung recall is
  `sum(TP) / sum(n_ref)`, never `mean(per-plot recall)` (small plots would
  dominate). Per-class TP is recovered as `round(rec_class * n_class)`.
- **Plot core tracks plot type**: tower base plots map ±20 m, distributed plots
  only ±10 m (`plot_half()` in sweep_lib). Scoring the full box over a
  distributed plot counts the unmapped ring as false commission.
- **Two density units**: first-return/pulse density (`frdens`) gates CHM
  resolution and the <8 pts/m² smoothing branch; all-return density (`pdens`)
  is the decimation/no-upsampling guard. Don't conflate them.
- **Matching** is global nearest-distance greedy 1:1 with a height-consistency
  gate (`greedy_match`), so a short understory stem can't steal a tall
  neighbour's apex.

## The ML / deep-learning benchmark prep

A separate, **inference-not-yet-run** track that prepares external deep models
(ForestFormer3D, SegmentAnyTree, ForAINet, …) for evaluation. Driven by two YAML
registries in [config/](config/) (`ml_benchmark_methods.yml`,
`ml_benchmark_datasets.yml`), with:

- [scripts/check_ml_method_status.R](scripts/check_ml_method_status.R) — probes
  local deps (docker, nvidia-smi, torch modules) → `work/ml_benchmark/method_status.csv`.
- [scripts/report_ml_benchmark.R](scripts/report_ml_benchmark.R) — renders
  [ml-method-benchmark-results.md](ml-method-benchmark-results.md).
- [scripts/score_point_instance_predictions.R](scripts/score_point_instance_predictions.R)
  — scores any method's exported predictions. **External prediction contract:**
  one row per point with `dataset, plot, point_id, x, y, z, pred_instance_id`
  (+ `ref_instance_id` for per-point instance scoring). Greedy IoU ≥ 0.5
  matching. This contract is how new models plug in — keep it stable.

## Coordinate-system gotcha

The USGS 3DEP EPT is **EPSG:3857** (Web Mercator); distances there are inflated
~1.32× at this latitude. Window sizes, CHM resolution, and density are all
metric-sensitive, so the PDAL extract reprojects to UTM (`extract.json` →
EPSG:32610) before detection. The lasR-only EPT path stays in 3857 for a fully
native acquisition demo — use it for acquisition benchmarking, not for
metric-faithful parameterization.
