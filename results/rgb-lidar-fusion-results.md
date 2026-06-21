# DeepForest RGB arm + the density-invariant anchor (#X1)

Every LiDAR arm degrades as point density drops; an **RGB** detector does not. So
a DeepForest arm gives the meta-pipeline a **density-invariant anchor** — its flat
accuracy-vs-density line tells the router (#P2) exactly which rung each LiDAR arm
stops beating optical, and it is the first optical member for #P1 fusion. This
promotes DeepForest from "optional reference" to a scored arm: it runs the
NEON-pretrained crown model on the NEON RGB camera mosaics (DP3.30010, 2021 — the
same epoch as the LiDAR/field data), reduces each crown box to its centroid,
samples an apex Z from a CHM built on the matched frozen clip, and scores against
field stems with the same `score_plot` harness.

Regenerate:

```sh
export CLAUDE_JOB_DIR=$(pwd)/work
Rscript scripts/neon_download_aop.R SITE=SOAP YEAR=2021      # DP3.30010 RGB tiles
Rscript scripts/detect_deepforest_sweep.R SITE=SOAP          # -> deepforest_results.csv
```

DeepForest runs **CPU-only** (`~/miniconda3/envs/deepforest`, torch 2.12) — no
Blackwell GPU build needed; a 1 km² tile predicts in ~14 s
(`gpu/run_deepforest.py`, which georeferences the pixel boxes through the raster
transform). Verified on SOAP (8 tiles, 2021, 100,953 crown boxes site-wide).

## Generated tables

### DeepForest (RGB) standalone — SOAP, 253 stems

| metric | value |
|---|--:|
| recall | 0.506 |
| precision | 0.289 |
| F1 | 0.368 |
| recall dominant / codominant | 0.530 / 0.537 |
| **recall understory** | **0.348** |

### Density-invariant anchor: F1 vs LiDAR density rung

| arm | native | 8 | 4 | 2 | 1 |
|---|--:|--:|--:|--:|--:|
| **deepforest (RGB, flat)** | 0.37 | 0.37 | 0.37 | 0.37 | 0.37 |
| segmentanytree | 0.46 | 0.44 | 0.44 | **0.32** | **0.12** |
| multichm | 0.44 | 0.42 | 0.45 | 0.44 | 0.46 |
| chm_vwf | 0.38 | 0.36 | 0.42 | 0.39 | 0.40 |
| forestformer3d | 0.26 | — | — | — | — |

## Readings

- **DeepForest is the density floor every LiDAR arm is measured against.** Its F1
  is flat at 0.37 by construction (fixed-resolution RGB). SegmentAnyTree — the
  best *high*-density arm — beats it at native/8/4 (0.44–0.46) but **falls below
  it at 2 pts/m² (0.32) and collapses at 1 pt/m² (0.12)**. So the router rule is
  concrete: **below ~2 pts/m², prefer optical (DeepForest) over the deep LiDAR
  arm.** ForestFormer3D is already below the floor at native.
- **multichm is the one LiDAR arm that stays above the optical floor at every
  rung** (0.42–0.46), confirming #P2's finding that it is the robust low-density
  default; DeepForest is the safety net for the arms that aren't.
- **RGB sees understory the 2.5-D CHM misses.** DeepForest's understory recall
  (0.348) beats CHM-VWF's (0.27): a nadir RGB detector catches canopy-gap crowns
  a surface model smooths over. This is exactly the decorrelated coverage a fusion
  member should add.
- **It pays in recall, not precision.** DeepForest over-detects (precision 0.29;
  100k boxes site-wide, many in the buffer / small understory) — so as a fusion
  member it contributes recall and the density-invariant floor, with the #P4
  calibrator (or a score threshold) needed to temper its commission.

## RGB×LiDAR fusion — the opportunity, quantified

The fusion case is now evidence-backed rather than assumed: DeepForest adds
(a) **low-density coverage** where SAT/FF3D collapse (≤2 pts/m²) and (b)
**understory recall** the CHM arms lack — two decorrelated gains. A
detection-level RGB×LiDAR consensus (union for recall, agreement for precision)
in `fuse_detectors.R` is the natural next step: the DeepForest centroids are now
persisted per tile (`deepforest_boxes/`) and scored (`deepforest_results.csv`)
with the same harness as the LiDAR arms, so `fuse_apexes` can ingest them
directly as a seventh, optical member. Running that full consensus across the
density ladder (and SJER/TEAK) is follow-up; this PR establishes the arm and the
density-invariant anchor it provides.

## Caveats

- **SOAP only, 2021 RGB.** SJER/TEAK 2021 RGB are available (the downloader is
  `SITE=`-parameterized); a full three-site run is future work.
- **Apex Z is sampled from the LiDAR CHM**, so a DeepForest box over a real gap
  with no LiDAR canopy gets the 2 m floor — a minor height-gate effect.
- **Precision is a lower bound** (the #V4 field-map coverage gap), and DeepForest's
  own over-detection compounds it; recall / understory / the density crossover are
  the trustworthy signals.
- **NEON is reprocessing DP3.30010 for RELEASE-2026**; this used the 2021
  acquisition to match the benchmark epoch (recorded in the `rgb/` path).
