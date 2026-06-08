# Cross-model density-ladder benchmark (SOAP)

Cross-model synthesis (#R10) of every tree detector currently runnable on the
NEON SOAP density ladder, scored on the **same frozen clips** by the **same**
field-stem harness, so the arms are directly comparable. It unifies the AMS3D
instance segmenter (`crownsegmentr`), the lidRplugins competitors (`lmfauto`,
`multichm`, `ptrees`), the CHM-VWF baseline (`detect_lasr`), and a native-only
Li 2012 point segmenter. The deep GPU models from the deep-research report
(SegmentAnyTree, TreeisoNet-ALS, ForestFormer3D) are deferred to the GPU track
(#M6–#M8) and will slot into these tables as further arms.

Regenerate:

```sh
export CLAUDE_JOB_DIR=$(pwd)/work
Rscript scripts/detect_ams3d_sweep.R       SITE=SOAP PLOTS=ALL CORES=12
Rscript scripts/detect_lidrplugins_sweep.R SITE=SOAP PLOTS=ALL CORES=12
Rscript scripts/detect_li2012_native.R     SITE=SOAP PLOTS=ALL CORES=12
Rscript scripts/analyze_model_benchmark.R  SITE=SOAP
```

## What this is

- **Population.** 18 SOAP plots, 232 field stems pooled, five density rungs
  (native ≈ 11.8, 8, 4, 2, 1 first-return pulses/m²). Every arm runs on the
  byte-identical frozen normalized clip per (plot, rung), so differences are the
  detector, not the input.
- **Scoring.** Each detector is reduced to apex detections `(x, y, z)` and
  scored by the unchanged `score_plot`/`greedy_match` (global nearest-distance
  1:1 with a height-consistency gate) against field stems in the plot core.
  Pooling sums counts (recall = ΣTP / Σn_ref), never averages per-plot rates.
- **Equal-set guard.** The density-ladder comparison keeps only (plot, rung)
  cells scored by **all five** ladder arms; here that dropped **0** cells (every
  arm returned a row for every cell, a 0-row table where it found nothing), so
  all five arms are compared on the identical 18-plot population at each rung.
- **Why Li 2012 is a separate table.** It is run native-only (a point segmenter
  is meaningless at 1–2 pts/m²), so including it in the ladder guard would drop
  every decimated cell. It appears only in the native point-segmenter
  head-to-head below.

## Read recall and precision together

Reducing an instance/point segmenter to apex detections **flatters it on
recall**: AMS3D and ptrees post the highest recall at native density (0.79 and
0.85) precisely because they emit many more apexes (precision 0.11 and 0.15)
than the CHM detectors. Recall alone therefore rewards over-segmentation; F1 is
the honest single-number summary, and per-class/precision columns are reported
alongside throughout.

## Per-crown-class recall, density ladder

| detector | rung | frdens | n_plots | n_ref | recall | precision | F1 | rec_dominant | rec_codominant | rec_understory |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ams3d | native | 11.78 | 18 | 232 | 0.79 | 0.11 | 0.19 | 0.79 | 0.82 | 0.69 |
| ams3d | 8 | 5.79 | 18 | 232 | 0.81 | 0.18 | 0.29 | 0.83 | 0.82 | 0.71 |
| ams3d | 4 | 3.00 | 18 | 232 | 0.74 | 0.24 | 0.37 | 0.79 | 0.75 | 0.62 |
| ams3d | 2 | 1.53 | 18 | 232 | 0.68 | 0.35 | 0.46 | 0.76 | 0.67 | 0.53 |
| ams3d | 1 | 0.78 | 18 | 232 | 0.44 | 0.42 | 0.43 | 0.52 | 0.43 | 0.33 |
| chm_vwf | native | 11.78 | 18 | 232 | 0.48 | 0.32 | 0.38 | 0.54 | 0.51 | 0.27 |
| chm_vwf | 8 | 5.79 | 18 | 232 | 0.33 | 0.48 | 0.39 | 0.44 | 0.33 | 0.13 |
| chm_vwf | 4 | 3.00 | 18 | 232 | 0.35 | 0.48 | 0.40 | 0.45 | 0.35 | 0.16 |
| chm_vwf | 2 | 1.53 | 18 | 232 | 0.34 | 0.46 | 0.39 | 0.46 | 0.33 | 0.13 |
| chm_vwf | 1 | 0.78 | 18 | 232 | 0.34 | 0.45 | 0.39 | 0.42 | 0.33 | 0.18 |
| lmfauto | native | 11.78 | 18 | 232 | 0.51 | 0.25 | 0.34 | 0.56 | 0.55 | 0.29 |
| lmfauto | 8 | 5.79 | 18 | 232 | 0.57 | 0.26 | 0.36 | 0.58 | 0.64 | 0.36 |
| lmfauto | 4 | 3.00 | 18 | 232 | 0.66 | 0.25 | 0.37 | 0.76 | 0.70 | 0.40 |
| lmfauto | 2 | 1.53 | 18 | 232 | 0.76 | 0.23 | 0.35 | 0.87 | 0.79 | 0.51 |
| lmfauto | 1 | 0.78 | 18 | 232 | 0.81 | 0.19 | 0.31 | 0.89 | 0.82 | 0.62 |
| multichm | native | 11.78 | 18 | 232 | 0.63 | 0.34 | 0.44 | 0.68 | 0.67 | 0.44 |
| multichm | 8 | 5.79 | 18 | 232 | 0.61 | 0.32 | 0.42 | 0.62 | 0.67 | 0.44 |
| multichm | 4 | 3.00 | 18 | 232 | 0.64 | 0.35 | 0.45 | 0.66 | 0.71 | 0.42 |
| multichm | 2 | 1.53 | 18 | 232 | 0.65 | 0.37 | 0.47 | 0.65 | 0.66 | 0.58 |
| multichm | 1 | 0.78 | 18 | 232 | 0.60 | 0.36 | 0.45 | 0.62 | 0.64 | 0.47 |
| ptrees | native | 11.78 | 18 | 232 | 0.85 | 0.15 | 0.26 | 0.83 | 0.88 | 0.78 |
| ptrees | 8 | 5.79 | 18 | 232 | 0.61 | 0.26 | 0.36 | 0.65 | 0.65 | 0.42 |
| ptrees | 4 | 3.00 | 18 | 232 | 0.50 | 0.34 | 0.41 | 0.58 | 0.54 | 0.27 |
| ptrees | 2 | 1.53 | 18 | 232 | 0.31 | 0.37 | 0.34 | 0.45 | 0.28 | 0.11 |
| ptrees | 1 | 0.78 | 18 | 232 | 0.21 | 0.38 | 0.27 | 0.28 | 0.21 | 0.07 |

Reading it:

- **multichm is the standout classical competitor** — the only arm that beats
  CHM-VWF on pooled F1 at **every** rung, while also lifting understory recall.
  Its recall is the most density-stable of all (0.60–0.65 across the whole
  ladder).
- **AMS3D** trades precision for recall: top-tier recall and understory recall
  (0.69–0.71 at the dense rungs) but precision 0.11–0.18, so its F1 trails
  CHM-VWF until the sparse rungs, where the field thins and precision recovers.
- **ptrees** is the native-density recall leader (0.85) but the least
  density-robust — recall collapses to 0.21 by 1 pt/m² as the point segmenter is
  starved of returns.
- **lmfauto inverts** the usual trend: recall *rises* as density falls
  (0.51 → 0.81) because its parameter-free window grows coarser and locks onto
  dominant crowns; precision falls in step.
- **CHM-VWF** holds the best precision (0.32–0.48) and the most even
  precision/recall trade, the reference every other arm is measured against.

## Per-height-band recall

Bands: short < 8 m, mid 8–15 m, tall ≥ 15 m.

| detector | rung | rec_h_tall | n_h_tall | rec_h_mid | n_h_mid | rec_h_short | n_h_short |
| --- | --- | --- | --- | --- | --- | --- | --- |
| ams3d | native | 0.58 | 60 | 0.79 | 90 | 0.94 | 79 |
| ams3d | 8 | 0.63 | 60 | 0.80 | 90 | 0.94 | 79 |
| ams3d | 4 | 0.63 | 60 | 0.78 | 90 | 0.77 | 79 |
| ams3d | 2 | 0.70 | 60 | 0.74 | 90 | 0.57 | 79 |
| ams3d | 1 | 0.67 | 60 | 0.47 | 90 | 0.24 | 79 |
| chm_vwf | native | 0.63 | 60 | 0.38 | 90 | 0.46 | 79 |
| chm_vwf | 8 | 0.43 | 60 | 0.29 | 90 | 0.28 | 79 |
| chm_vwf | 4 | 0.48 | 60 | 0.30 | 90 | 0.29 | 79 |
| chm_vwf | 2 | 0.53 | 60 | 0.23 | 90 | 0.29 | 79 |
| chm_vwf | 1 | 0.50 | 60 | 0.30 | 90 | 0.23 | 79 |
| lmfauto | native | 0.62 | 60 | 0.44 | 90 | 0.48 | 79 |
| lmfauto | 8 | 0.67 | 60 | 0.53 | 90 | 0.52 | 79 |
| lmfauto | 4 | 0.80 | 60 | 0.54 | 90 | 0.68 | 79 |
| lmfauto | 2 | 0.83 | 60 | 0.73 | 90 | 0.73 | 79 |
| lmfauto | 1 | 0.92 | 60 | 0.74 | 90 | 0.78 | 79 |
| multichm | native | 0.73 | 60 | 0.60 | 90 | 0.58 | 79 |
| multichm | 8 | 0.77 | 60 | 0.57 | 90 | 0.54 | 79 |
| multichm | 4 | 0.82 | 60 | 0.57 | 90 | 0.58 | 79 |
| multichm | 2 | 0.82 | 60 | 0.64 | 90 | 0.51 | 79 |
| multichm | 1 | 0.83 | 60 | 0.52 | 90 | 0.51 | 79 |
| ptrees | native | 0.87 | 60 | 0.83 | 90 | 0.85 | 79 |
| ptrees | 8 | 0.73 | 60 | 0.54 | 90 | 0.58 | 79 |
| ptrees | 4 | 0.62 | 60 | 0.38 | 90 | 0.54 | 79 |
| ptrees | 2 | 0.45 | 60 | 0.24 | 90 | 0.25 | 79 |
| ptrees | 1 | 0.37 | 60 | 0.13 | 90 | 0.16 | 79 |

The short-tree band is where density bites hardest. AMS3D and ptrees own short
trees at native density (0.94 and 0.85) but both collapse on the sparse rungs
(AMS3D short 0.94 → 0.24; ptrees 0.85 → 0.16). CHM-VWF short recall is flat and
mediocre (~0.23–0.46) at all densities. lmfauto's tall-tree recall climbs to
0.92 at 1 pt/m² — the coarse-window effect again, now visibly concentrated in
the dominant canopy.

## Density-robustness curves

Two figures under `work/neon/SOAP/figs/` (regenerated, gitignored):

- `model_recall_vs_density.png` — overall recall vs first-return density, one
  line per ladder arm. The shapes summarise the table: multichm flat, CHM-VWF
  flat-and-low, AMS3D a gentle decline, ptrees a steep collapse, lmfauto rising
  into low density.
- `model_understory_vs_density.png` — understory recall (intermediate +
  suppressed) vs density. AMS3D and ptrees lead at the top of the ladder and
  converge down toward CHM-VWF at the bottom.

## Native point-segmenter head-to-head

The report's central claim is that point/instance methods reach sub-canopy stems
the CHM cannot. At native density, on the identical 18-plot / 45-understory-stem
population:

| detector | n_plots | n_ref | recall | precision | F1 | rec_understory | n_understory |
| --- | --- | --- | --- | --- | --- | --- | --- |
| ams3d | 18 | 232 | 0.79 | 0.11 | 0.19 | 0.69 | 45 |
| chm_vwf | 18 | 232 | 0.48 | 0.32 | 0.38 | 0.27 | 45 |
| li2012 | 18 | 232 | 0.59 | 0.26 | 0.36 | 0.33 | 45 |
| ptrees | 18 | 232 | 0.85 | 0.15 | 0.26 | 0.78 | 45 |

Every point/instance segmenter beats CHM-VWF on understory recall — ptrees 0.78
and AMS3D 0.69 emphatically, Li 2012 0.33 modestly — confirming the thesis, but
all of them pay for it in precision, and only Li 2012 stays near CHM-VWF on F1.
The sub-canopy gain is real; whether it is worth the commission cost depends on
the downstream use.

## Head-to-head deltas vs CHM-VWF

Each arm minus the CHM-VWF baseline at the same rung (positive = the arm beats
the baseline):

| detector | rung | d_recall | d_F1 | d_understory |
| --- | --- | --- | --- | --- |
| ams3d | native | 0.31 | -0.19 | 0.42 |
| ams3d | 8 | 0.47 | -0.10 | 0.58 |
| ams3d | 4 | 0.39 | -0.04 | 0.47 |
| ams3d | 2 | 0.34 | 0.08 | 0.40 |
| ams3d | 1 | 0.11 | 0.04 | 0.16 |
| lmfauto | native | 0.03 | -0.04 | 0.02 |
| lmfauto | 8 | 0.24 | -0.03 | 0.22 |
| lmfauto | 4 | 0.31 | -0.04 | 0.24 |
| lmfauto | 2 | 0.43 | -0.04 | 0.38 |
| lmfauto | 1 | 0.47 | -0.07 | 0.44 |
| multichm | native | 0.15 | 0.06 | 0.18 |
| multichm | 8 | 0.28 | 0.03 | 0.31 |
| multichm | 4 | 0.29 | 0.05 | 0.27 |
| multichm | 2 | 0.31 | 0.08 | 0.44 |
| multichm | 1 | 0.27 | 0.07 | 0.29 |
| ptrees | native | 0.37 | -0.12 | 0.51 |
| ptrees | 8 | 0.28 | -0.03 | 0.29 |
| ptrees | 4 | 0.16 | 0.00 | 0.11 |
| ptrees | 2 | -0.03 | -0.05 | -0.02 |
| ptrees | 1 | -0.12 | -0.12 | -0.11 |

Everyone beats CHM-VWF on recall and understory recall at the denser rungs, but
**only multichm beats it on F1 at every rung** (+0.03 to +0.08): the others buy
recall with precision the CHM keeps. ptrees crosses below the baseline on every
metric by 1–2 pts/m², the clearest "point segmenter needs points" signal in the
study.

## Appendix: zero-shot ledger

Every arm here is **classical / parameter-derived**, applied zero-shot with no
NEON-specific fitting; the only knobs are the density-first parameters the repo
already derives (CHM resolution and the VWF window from measured first-return
density — see [treetop-detection-approach.md](../docs/treetop-detection-approach.md)),
plus literature-default crown allometry for AMS3D (Ferraz 2016). The triage,
zero-shot protocol, and weights-mirror policy are in
[model-benchmark-plan.md](../docs/model-benchmark-plan.md) (#A0). The pretrained
deep models (SegmentAnyTree, TreeisoNet-ALS, ForestFormer3D) are deferred to the
GPU track (#M6–#M8) and will be added as further arms in this same synthesis.

## Caveats

- **Decimation is not native low density.** The sparse rungs are decimated from
  the native cloud, which removes returns but preserves the native scan geometry
  and pulse pattern; a genuinely low-density acquisition differs. Carry the
  native USGS 3DEP cross-check
  ([native-ql2-crosscheck-results.md](native-ql2-crosscheck-results.md), #4) as
  the precedent for how much this matters.
- **Discrete-return ALS vs ULS/UAS.** The instance segmenters in the source
  report were developed/trained largely on dense ULS/UAS or TLS point clouds;
  NEON is discrete-return airborne LiDAR at far lower density. These results
  characterise zero-shot transfer to ALS, not the methods at their design
  density.
- **Field-stem ground truth reduces every model to detections.** Scoring is
  apex recall/precision against mapped stems, not point-level instance IoU; an
  over-segmenter is penalised only through precision, and a model that recovers
  crown shape well but misses the stem location is not credited for the shape.
- **Single site, partial Li 2012.** SOAP only; the cross-site structure gradient
  (SJER, TEAK) and a full-ladder Li 2012 are deferred (#E11). The native-only Li
  2012 here answers the dense-input sub-canopy question, not its density
  response.
