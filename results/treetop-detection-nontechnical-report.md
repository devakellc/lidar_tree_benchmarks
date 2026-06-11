# Tree-Top Detection Results: Non-Technical Report

Prepared: 2026-06-11

This report summarizes the tree-top detection experiments in the `results/`
folder. The short version is:

- The system can map visible tree tops, especially dominant and co-dominant
  trees, reasonably consistently.
- It cannot reliably count hidden understory trees from airborne LiDAR. That is
  mostly a visibility problem, not an algorithm problem.
- Forest structure matters more than raw point density. Open forests are easier;
  dense, interlocking conifer crowns are harder.
- A practical production system should choose the detector from measured point
  density, canopy structure, and whether the project values "fewer false tops"
  or "more complete tree counts".

## 1. Dataset

The main benchmark uses NEON 2021 airborne LiDAR and NEON field-surveyed tree
stems. The LiDAR gives a 3D point cloud of the canopy. The field data provides
surveyed tree locations, heights, and canopy classes used as ground truth.

Three California forest sites were used:

| Site | Plain-English forest type | Why it matters | Native first-return density | Field plots / stems used |
| --- | --- | --- | ---: | ---: |
| SJER | Open oak / foothill-pine savanna | Open crowns, broadleaf shapes, easier separation | about 9/m^2 | about 8 plots / 70 stems |
| SOAP | Mixed conifer / deciduous forest | Main anchor case; mixed, multi-layered canopy | about 12/m^2 | 18 plots / 232 stems |
| TEAK | Dense red-fir conifer forest | Closed, packed crowns; hardest structure | about 12/m^2 | about 20 plots / 350-400 stems |

Counts differ slightly between experiments because some comparisons use
"equal-set" filters: a plot is kept only when every algorithm in that comparison
successfully produced a result.

The density ladder was built by starting from native NEON LiDAR and thinning it
to lower point densities: native, 8, 4, 2, and 1 all-return points/m^2. In these
forests, first-return density is roughly half to two-thirds of all-return point
density, so a rung labelled "4" is closer to about 3 first returns/m^2.

The scoring is stem-level:

- Recall: of the field trees, how many were found.
- Precision: of the detected tree tops, how many matched a field tree.
- F1: one combined score balancing recall and precision.

## 2. Algorithms and Ladder-Sweeping Approach

The core ladder sweep used a density-first workflow:

1. Clip each field plot from the LiDAR.
2. Thin the point cloud to the target density rung.
3. Normalize heights so the ground is zero and tree height is above ground.
4. Build a canopy height model, which is like a top-down height map of the
   canopy surface.
5. Detect local peaks in that surface as tree tops.
6. Match detected tops to field trees within a 4 m horizontal tolerance, with a
   height check to avoid obviously wrong matches.

The baseline detector is `CHM-VWF`: a canopy-height-model detector with a
variable-window filter. In simple terms, it searches for high points, using
larger search windows for taller trees and smaller windows for shorter trees.

The sweep varied:

- LiDAR density: native, 8, 4, 2, and 1 all-return points/m^2.
- Canopy height model resolution: 0.25 m, 0.5 m, and 1.0 m where valid.
- Variable-window slope: 0.05, 0.10, and 0.15.
- Sparse-data smoothing: added when first-return density dropped below about
  8/m^2 to reduce noisy false peaks.

Several alternative detectors were also compared:

| Detector family | What it means for non-technical readers | Main finding |
| --- | --- | --- |
| CHM-VWF | Conservative surface peak detector | Stable, cheap, good baseline |
| `multichm` | Looks for peaks across multiple canopy layers | Best classical recall-first option in closed conifer forests |
| Point-cloud local maxima | Searches for peaks directly in the points | Did not beat CHM-VWF; mostly added false tops |
| Li 2012 | 3D point-cloud tree segmentation | Finds a few more understory trees but loses precision and costs more |
| SegmentAnyTree | Deep-learning instance model | Strongest SOAP native F1, but density-sensitive |
| TreeisoNet | Deep-learning / instance model | Stable competitor across SOAP density rungs |
| AMS3D / ptrees / lmfauto | Classical or plugin competitors | Often high recall, but many false detections |

## 3. Results

### Main Density Finding

For the baseline CHM-VWF detector, F1 stayed fairly flat across a large density
range. At SOAP, the main mixed-forest site, CHM-VWF stayed around F1 0.38-0.40
in the model benchmark and about 0.35-0.42 in the density-ladder sweep.

That does not mean density is irrelevant. Density changes the error type:

- At native density, the detector finds more trees but also creates more false
  tops.
- At sparse density, smoothing removes many false tops, so precision improves,
  but some real trees are missed.

The balance stays similar, but the product behavior changes.

### Forest Structure Matters Most

The same baseline detector behaved differently by forest type:

| Site | Structure | Native overstory recall | Sparse-rung overstory recall |
| --- | --- | ---: | ---: |
| SJER | Open oak savanna | about 0.71 | about 0.57 |
| SOAP | Mixed conifer | about 0.57 | about 0.39 |
| TEAK | Dense red fir | about 0.40 | about 0.25 |

The important point: TEAK has conifer crowns that look "ideal" individually, but
the crowns are packed together. The detector often sees several trees as one
merged canopy surface. Open spacing at SJER helps more than crown shape.

### Understory Trees Remain Hard

Understory recall stayed low at SOAP and TEAK. The suppressed class often fell
near zero below the 8 points/m^2 rung. Point-cloud and deep methods recovered
some extra understory trees, but none turned the understory into a reliable
stem census.

The best native point-cloud test, Li 2012, raised pooled understory recall from
0.200 to 0.257 across the three sites. That is only 6 extra understory trees out
of 105, with lower precision and much higher runtime.

### Heights Are Good for Visible Trees

For dominant trees, detected LiDAR apex height matched field height well:

- R^2: 0.85
- Slope: about 0.98
- Bias: about +1 m
- RMSE: about 4.4 m

For sub-canopy trees, height quality drops because the matched LiDAR peak is
often the taller neighboring crown above the field stem. This confirms the same
visibility limit seen in detection.

### Algorithm Winners Depend on the Use Case

At SOAP, the full model benchmark shows:

| Detector | Native recall | Native precision | Native F1 | Understory recall |
| --- | ---: | ---: | ---: | ---: |
| CHM-VWF | 0.48 | 0.32 | 0.38 | 0.27 |
| `multichm` | 0.63 | 0.34 | 0.44 | 0.44 |
| SegmentAnyTree | 0.66 | 0.38 | 0.48 | 0.49 |
| TreeisoNet | 0.60 | 0.29 | 0.39 | 0.33 |
| Li 2012 | 0.59 | 0.26 | 0.36 | 0.33 |

The most useful pattern is:

- SegmentAnyTree is strongest on SOAP at native density and remains competitive
  at the 8 and 4 rungs, but drops badly at 2 and 1.
- `multichm` is the most reliable classical upgrade over CHM-VWF in conifer
  sites. Its advantage survives held-out validation at SOAP and TEAK.
- CHM-VWF remains the conservative, scalable baseline and is safer in open
  savanna-like structure where extra multi-layer detections are often false.
- Point-cloud local maxima and lasR/lidR local maxima are effectively the same
  operator when given the same surface; the canopy surface construction matters
  more than the local-maximum implementation.

### Robustness Checks

The supporting studies do not overturn the main result:

- Calibration/validation: 0.5 m CHM resolution is a strong default for
  decimated SOAP and TEAK rungs, but not universal. Native density and SJER can
  prefer 0.25 m or 1.0 m. The exact VWF slope is second-order.
- Native 3DEP cross-check: decimation predicts sparse-cloud recall within about
  0.07 overall, with a small optimistic bias for overstory trees.
- Temporal sensitivity: the field-to-LiDAR year mismatch has little effect at
  native density and does not change the main conclusions.
- Crown segmentation: after tree tops are found, lasR `region_growing` gives the
  best crown-diameter match, but crown segmentation is a separate downstream
  problem from top detection.

## 4. Outcomes and Suggested Meta-Pipeline

Use a meta-pipeline rather than one fixed detector everywhere. The routing rule
should use four inputs:

1. Measured first-return density, not the project label.
2. Canopy structure: open, mixed, or closed/interlocking.
3. Business priority: fewer false tops or more complete tree counts.
4. Compute budget: CPU-only wall-to-wall processing or GPU/deep-model runs.

Recommended routing:

| Situation | First choice | When to switch |
| --- | --- | --- |
| Low density, below about 4 first returns/m^2 | CHM-VWF at 0.5-1.0 m resolution | Use `multichm` only when recall matters more than false detections |
| Medium density, about 4-8 first returns/m^2 | CHM-VWF or `multichm` | Pick CHM-VWF for open forests; pick `multichm` for closed conifer |
| High density, 8+ first returns/m^2 | CHM-VWF baseline plus `multichm` or SegmentAnyTree | Use SegmentAnyTree where GPU is available and local validation confirms it |
| Open oak / savanna | CHM-VWF, coarser CHM, larger windows | Avoid recall-first methods unless false tops are acceptable |
| Mixed conifer | `multichm` or SegmentAnyTree | Use CHM-VWF for conservative production maps |
| Dense red fir / closed conifer | `multichm` | Expect lower recall; report understory as low-confidence |
| Understory inventory | Do not rely on airborne LiDAR tops alone | Use field plots, terrestrial/mobile LiDAR, or treat results as partial |

Operational recommendation:

1. Run a cheap preflight over each tile: density, canopy cover, height
   variability, and crown spacing.
2. Route each tile or stand to CHM-VWF, `multichm`, or a deep/GPU arm.
3. Keep CHM-VWF as the reference baseline everywhere.
4. Validate each forest type with held-out plots before changing defaults.
5. Report visible-overstory results separately from understory estimates.

For large tiled jobs, use the lasR streaming path for the baseline production
run because it reproduced single-file tiled results exactly in the comparison.
Use lidR and GPU model pipelines for targeted algorithm arms, smaller jobs, or
validation runs.

## Source Result Files

This report summarizes:

- `density-ladder-sweep-results.md`
- `model-benchmark-results.md`
- `pointcloud-detector-results.md`
- `treetop-lasr-vs-lidr-comparison.md`
- `calibration-validation-results.md`
- `native-ql2-crosscheck-results.md`
- `temporal-sensitivity-results.md`
- `cross_site_summary.csv`

Adjacent, non-primary evidence:

- `crown-segmentation-results.md` for downstream crown delineation.
- `ept-acquisition-sweep-results.md` for large EPT acquisition/runtime behavior.
