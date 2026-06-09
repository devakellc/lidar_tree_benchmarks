# Cross-model density-ladder benchmark (SOAP)

Cross-model synthesis (#R10) of every tree detector currently runnable on the
NEON SOAP density ladder, scored on the **same frozen clips** by the **same**
field-stem harness, so the arms are directly comparable. It unifies the AMS3D
instance segmenter (`crownsegmentr`), the lidRplugins competitors (`lmfauto`,
`multichm`, `ptrees`), the CHM-VWF baseline (`detect_lasr`), a native-only
Li 2012 point segmenter, and two **deep GPU models** run zero-shot on the
RTX 5090 — the TreeisoNet-ALS instance segmenter (#M7) and ForestFormer3D
(#M8, native + 8 only; see its section below). SegmentAnyTree (#M6) stays
deferred behind a remaining sm_120 inference blocker and will slot in as a
further arm.

Regenerate:

```sh
export CLAUDE_JOB_DIR=$(pwd)/work
Rscript scripts/detect_ams3d_sweep.R       SITE=SOAP PLOTS=ALL CORES=12
Rscript scripts/detect_lidrplugins_sweep.R SITE=SOAP PLOTS=ALL CORES=12
Rscript scripts/detect_li2012_native.R     SITE=SOAP PLOTS=ALL CORES=12
Rscript scripts/detect_treeisonet_sweep.R                       # GPU, serial; CONF=0.22 default
Rscript scripts/detect_forestformer3d_sweep.R SITE=SOAP REPO=<FF3D repo>  # GPU, serial; native+8 (~8 min)
Rscript scripts/analyze_model_benchmark.R  SITE=SOAP
```

## What this is

- **Population.** 18 SOAP plots, 232 field stems pooled, five density rungs —
  the undecimated native cloud plus all-return decimation targets of 8/4/2/1
  pts/m² (the `lidR::homogenize` unit), i.e. first-return ≈ 11.8 / 5.8 / 3.0 /
  1.5 / 0.8 pulses/m² as reported in the `frdens` column and used as the
  density-curve x-axis. Every arm runs on the byte-identical frozen normalized
  clip per (plot, rung), so differences are the detector, not the input.
- **Scoring.** Each detector is reduced to apex detections `(x, y, z)` and
  scored by the unchanged `score_plot`/`greedy_match` (global nearest-distance
  1:1 with a height-consistency gate) against field stems in the plot core.
  Pooling sums counts (recall = ΣTP / Σn_ref), never averages per-plot rates.
- **Equal-set guard.** The density-ladder comparison keeps only (plot, rung)
  cells scored by **all six** ladder arms; here that dropped **0** cells (every
  arm returned a row for every cell, a 0-row table where it found nothing —
  TreeisoNet included), so all six arms are compared on the identical 18-plot
  population at each rung, and the five classical arms' numbers are unchanged by
  adding the deep model.
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
| treeisonet | native | 11.78 | 18 | 232 | 0.27 | 0.08 | 0.12 | 0.25 | 0.28 | 0.27 |
| treeisonet | 8 | 5.79 | 18 | 232 | 0.06 | 0.10 | 0.07 | 0.07 | 0.04 | 0.07 |
| treeisonet | 4 | 3.00 | 18 | 232 | 0.02 | 0.21 | 0.04 | 0.03 | 0.02 | 0.00 |
| treeisonet | 2 | 1.53 | 18 | 232 | 0.00 | 0.00 | NA | 0.00 | 0.00 | 0.00 |
| treeisonet | 1 | 0.78 | 18 | 232 | 0.00 | NA | NA | 0.00 | 0.00 | 0.00 |

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
- **TreeisoNet (the deep model)** is the **weakest arm**: native recall 0.27 /
  F1 0.12, and it collapses to ~0 by rung 2 — see the dedicated section below.

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
| treeisonet | native | 0.42 | 60 | 0.29 | 90 | 0.14 | 79 |
| treeisonet | 8 | 0.07 | 60 | 0.09 | 90 | 0.00 | 79 |
| treeisonet | 4 | 0.05 | 60 | 0.01 | 90 | 0.00 | 79 |
| treeisonet | 2 | 0.00 | 60 | 0.00 | 90 | 0.00 | 79 |
| treeisonet | 1 | 0.00 | 60 | 0.00 | 90 | 0.00 | 79 |

The short-tree band is where density bites hardest. AMS3D and ptrees own short
trees at native density (0.94 and 0.85) but both collapse on the sparse rungs
(AMS3D short 0.94 → 0.24; ptrees 0.85 → 0.16). CHM-VWF short recall is flat and
mediocre (~0.23–0.46) at all densities. lmfauto's tall-tree recall climbs to
0.92 at 1 pt/m² — the coarse-window effect again, now visibly concentrated in
the dominant canopy. TreeisoNet is worst in every band even at native (tall 0.42
vs CHM-VWF 0.63, short 0.14 vs 0.46) and zero below rung 8.

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
| treeisonet | 18 | 232 | 0.27 | 0.08 | 0.12 | 0.27 | 45 |

The *classical* point/instance segmenters beat CHM-VWF on understory recall —
ptrees 0.78 and AMS3D 0.69 emphatically, Li 2012 0.33 modestly — confirming the
thesis, all paying for it in precision (only Li 2012 stays near CHM-VWF on F1).
The lone **deep** arm, TreeisoNet, is the exception: understory 0.27 only *ties*
CHM-VWF and its overall recall (0.27) is the lowest of any native arm. The
sub-canopy gain is real for the classical point methods; the zero-shot deep
model does not deliver it here.

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
| treeisonet | native | -0.21 | -0.26 | 0.00 |
| treeisonet | 8 | -0.28 | -0.32 | -0.07 |
| treeisonet | 4 | -0.33 | -0.37 | -0.16 |
| treeisonet | 2 | -0.34 | NA | -0.13 |
| treeisonet | 1 | -0.34 | NA | -0.18 |

The *classical* arms beat CHM-VWF on recall and understory recall at the denser
rungs, but **only multichm beats it on F1 at every rung** (+0.03 to +0.08): the
others buy recall with precision the CHM keeps. ptrees crosses below the
baseline by 1–2 pts/m², the clearest "point segmenter needs points" signal among
them. TreeisoNet is **below CHM-VWF on every metric at every rung** (Δrecall
−0.21 → −0.34) — the deep model does not clear even the classical baseline.

## Deep model: TreeisoNet zero-shot (#M7)

TreeisoNet-ALS (`treeisonet_als_reclamation`, esegformer3D, NRCan/TreeAIBox) is
the first deep instance segmenter in the benchmark, run zero-shot on the RTX
5090 (Blackwell sm_120, torch cu128) via `gpu/run_treeisonet.py`: TreeLoc →
`postPeakExtraction` tops → snap each apex z to the local canopy max (the same
canopy-surface height the CHM arms use). One fixed confidence threshold
(`conf = 0.22`) was calibrated once for pooled F1 on a 5-plot subset, then
applied unchanged across the ladder — no per-plot tuning, postPeakExtraction
defaults otherwise.

The result is a clear negative. TreeisoNet is the weakest arm at every rung and
below the CHM-VWF baseline throughout (native recall 0.27 / precision 0.08 / F1
0.12; understory 0.27 merely ties CHM-VWF). It **collapses with density** far
faster than any classical arm — recall 0.27 → 0.06 → ~0 from native to rung 4,
and exactly zero at rungs 2 and 1 (precision then undefined, shown `NA`). Lower
thresholds do not help: at `conf = 0.10` it floods the core with hundreds of
spurious tops (precision ~0.01), so no threshold yields a usable F1 (peak ≈ 0.06
in calibration).

This is the zero-shot domain-shift failure made concrete: a model trained on
dense ULS/UAV/TLS (its published ALS auto-mIoU is ~0.59) does not transfer to
sparse discrete-return NEON ALS without fine-tuning, and at NEON densities it is
beaten by a tuned classical CHM detector. Read it as a *floor* the heavier deep
arms must clear to justify their cost — ForestFormer3D (#M8, next section) clears
it comfortably but still trails the classical baseline; SegmentAnyTree (#M6)
remains to be run.

## ForestFormer3D (#M8) — zero-shot, native + 8 only

ForestFormer3D (ICCV 2025), ported to the RTX 5090 (sm_120; torch 2.7 / cu128 —
see [`gpu/forestformer3d-sm120/`](../gpu/forestformer3d-sm120/README.md)), run
zero-shot on the **top of the ladder only** (native + 8 pts/m²; it is the
heaviest arm, ~8 min for a full SOAP run on one GPU). Each plot core is tiled into
16 m-radius cylinders (`SPACING = 24` m → 8 m overlap; 9 cylinders per tower
plot, 4 per distributed), with one model pass per plot; cross-block instances are
then merged by apex-cluster union-find (`MERGE_TOL = 2` m, **different blocks
only**, so a model's within-cylinder over-segmentation is scored honestly) before
the shared apex reducer. It is reported as its **own** native+8 pool against the
CHM-VWF baseline and the other deep arm over the identical 18-plot / 232-stem
set, and is **not** part of the five-rung ladder above (an additive comparison,
so it cannot shrink the ladder's equal-set population).

| detector | rung | recall | precision | F1 | rec_dominant | rec_understory |
| --- | --- | --- | --- | --- | --- | --- |
| forestformer3d | native | 0.40 | 0.19 | 0.25 | 0.46 | 0.20 |
| forestformer3d | 8 | 0.45 | 0.25 | 0.32 | 0.55 | 0.22 |
| chm_vwf | native | 0.48 | 0.32 | 0.38 | 0.54 | 0.27 |
| chm_vwf | 8 | 0.33 | 0.48 | 0.39 | 0.44 | 0.13 |
| treeisonet | native | 0.27 | 0.08 | 0.12 | 0.25 | 0.27 |
| treeisonet | 8 | 0.06 | 0.10 | 0.07 | 0.07 | 0.07 |

Two readings. (1) **Zero-shot, FF3D trails the tuned classical CHM-VWF baseline**
(F1 0.25–0.32 vs 0.38–0.39) — the expected dense-ULS/TLS → sparse-ALS domain gap,
the same story as TreeisoNet but on a query-based transformer. (2) But FF3D
**clears the TreeisoNet floor comfortably** (F1 0.32 vs 0.07 at rung 8; ~4× the
dominant-tree recall), so among pretrained deep arms it transfers far better to
NEON ALS. Notably FF3D's F1 *rises* from native to rung 8 (0.25 → 0.32): the
denser native cloud yields more spurious instances (precision 0.19), and
decimating toward its training density sharpens both precision and recall — a
mild signal its sweet spot sits below NEON-native density. Understory recall
(~0.20) is non-trivial but well under the canopy-class recall, consistent with a
canopy-trained model.

## Appendix: zero-shot ledger

Every arm except TreeisoNet is **classical / parameter-derived**, applied
zero-shot with no NEON-specific fitting; the only knobs are the density-first
parameters the repo already derives (CHM resolution and the VWF window from
measured first-return density — see
[treetop-detection-approach.md](../docs/treetop-detection-approach.md)), plus
literature-default crown allometry for AMS3D (Ferraz 2016). TreeisoNet and
ForestFormer3D are the **pretrained deep** arms, also zero-shot: published
weights, no fine-tuning. The triage, zero-shot protocol, and weights-mirror
policy are in [model-benchmark-plan.md](../docs/model-benchmark-plan.md) (#A0).
SegmentAnyTree (#M6) remains deferred behind a remaining sm_120 inference blocker
and will be added as a further arm in this same synthesis.

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
  density — TreeisoNet's collapse and ForestFormer3D's sub-baseline F1 (both
  above) are the concrete evidence, and #M6 is expected to share the risk.
- **Field-stem ground truth reduces every model to detections.** Scoring is
  apex recall/precision against mapped stems, not point-level instance IoU; an
  over-segmenter is penalised only through precision, and a model that recovers
  crown shape well but misses the stem location is not credited for the shape.
- **Single site, partial Li 2012.** SOAP only; the cross-site structure gradient
  (SJER, TEAK) and a full-ladder Li 2012 are deferred (#E11). The native-only Li
  2012 here answers the dense-input sub-canopy question, not its density
  response.
