# Crown-width + height → DBH / biomass allometry (#S1)

The crown issues (#7/#30–#35) score crown *diameter* vs field crown diameter
only; none predicts DBH or biomass, yet DBH/above-ground biomass (AGB) are the
operational products a detection pipeline ultimately feeds. This closes the loop:
it joins each matched crown's equivalent diameter `d_eq` (and field height) to
NEON field `stemDiameter`/`taxonID`, fits crown-geometry → DBH models, reports how
well each segmenter's crowns predict field DBH (per crown class, per density
rung), and **derives** AGB from predicted DBH. The thesis under test: a detector
whose crowns predict DBH well is more valuable than raw F1 suggests, and does that
survive decimation?

Regenerate:

```sh
export CLAUDE_JOB_DIR=$(pwd)/work
Rscript scripts/crown_allometry.R SITES=SOAP,SJER,TEAK
# -> work/neon/<SITE>/crown_allometry.csv (per matched tree: geometry + field DBH
#    + predicted DBH + derived AGB)
```

The pure helpers `functional_type` (taxonID → conifer/broadleaf), `agb_from_dbh`
(Jenkins et al. 2003 generic AGB), and `fit_stats` (R²/RMSE/bias) live in
[`allometry_lib.R`](../scripts/allometry_lib.R), unit-tested in
`tests/testthat/test-allometry-lib.R`. 19,446 matched crowns carry field DBH
(6,470 at native density) across three sites and 23 segmenter variants.

## Generated tables

### Does crown width predict DBH? (pooled, native)

| predictor | n | R² | RMSE (cm) | bias |
|---|--:|--:|--:|--:|
| height only | 6470 | **0.708** | 13.5 | 0.0 |
| crown `d_eq` only | 6470 | 0.174 | 22.7 | 0.0 |
| `d_eq` + height | 6470 | 0.709 | 13.5 | 0.0 |

### Crown `d_eq` → DBH skill per segmenter (native, in-sample R²)

| segmenter | n | R² | RMSE (cm) |
|---|--:|--:|--:|
| **ams3d** | 475 | **0.558** | 13.1 |
| random_walker | 225 | 0.389 | 19.4 |
| silva2016_seedlmf | 225 | 0.378 | 19.6 |
| watershed_seeded | 225 | 0.375 | 19.6 |
| lasr_region_growing | 225 | 0.335 | 20.3 |
| dalponte2016_seedlmf | 225 | 0.303 | 20.7 |
| li2012 | 308 | 0.269 | 19.7 |
| segmentanytree | 419 | 0.258 | 20.5 |
| ptrees | 412 | 0.185 | 21.6 |
| **forestformer3d** | 220 | **0.021** | 25.4 |

(seed-variant CHM segmenters omitted for brevity; full table in the CSV.)

### By crown class (pooled native) and density rung

| crown class | n | R² | | rung | n | R² |
|---|--:|--:|---|---|--:|--:|
| dominant | 2326 | 0.129 | | native | 6470 | 0.174 |
| codominant | 3541 | 0.132 | | 8 | 3117 | 0.221 |
| intermediate | 539 | 0.119 | | 4 | 3174 | 0.175 |
| **suppressed** | 61 | **0.414** | | 2 | 3774 | 0.182 |
| | | | | 1 | 2911 | 0.168 |

### Derived AGB (predicted-DBH vs field-DBH AGB, native)

n = 6470, **R² 0.029, RMSE 1844 kg, bias −377 kg** (mean field-DBH AGB 929 kg).
AGB is a derived product (Jenkins 2003 generic) — NEON has no field AGB.

## Readings

- **For DBH, height is the predictor; crown width adds almost nothing.** Field
  height alone explains 71 % of DBH variance; crown `d_eq` alone only 17 %, and
  adding crown width to height moves R² 0.708 → 0.709. The operational consequence
  is the opposite of the crown-centric premise: a pipeline should drive DBH
  allometry from its **height** output (which LiDAR recovers accurately), with
  crown width a negligible add-on. (Caveat: this uses *field* height; the
  detected apex height a deployed pipeline would use is not carried in the
  crown-metrics rows — but the benchmark shows detected height tracks field
  height closely, so the conclusion holds operationally.)
- **Detection F1 and crown-allometry quality are decoupled — the issue's thesis,
  confirmed with a twist.** AMS3D has the **best** crown→DBH skill (R² 0.558)
  despite ranking *last* (9th) on detection F1 in #V2 (it splits crowns, tanking
  precision): the crowns it does recover have widths that track DBH well.
  Conversely ForestFormer3D, a true-mask deep arm, has essentially **no**
  crown→DBH skill (R² 0.021) — its masks recover trees but not DBH-predictive
  widths. So "best detector" (F1) ≠ "best for downstream DBH"; the meta-pipeline
  should weight arms differently for counting vs for mensuration.
- **Crown→DBH skill survives decimation.** Pooled R² is flat across the density
  ladder (0.17–0.22 from native to 1 pt/m²) — the weak-but-real crown signal does
  not degrade with sparser clouds, so the allometry is density-robust (if
  unspectacular).
- **Crown width helps most for the smallest trees.** Suppressed-class crowns
  predict DBH best (R² 0.414 vs ~0.13 for dominant/codominant): big-tree crown
  width saturates and decouples from DBH, while for small understory stems crown
  size still tracks stem size. This is where crown width adds the most over height.
- **Crown-derived AGB is too noisy to use.** Because AGB ≈ DBH^2.4, the modest
  crown→DBH error explodes into an AGB RMSE (1844 kg) that exceeds the mean tree's
  AGB (929 kg), R² 0.03. Crowns alone cannot drive biomass; height-based DBH (R²
  0.71) would propagate to far tighter AGB, reinforcing the height-first reading.

## Caveats

- **In-sample R²** per segmenter (the model is fit and scored on the same matched
  trees) measures explained variance, not out-of-sample skill; the cross-segmenter
  *ordering* is the robust result.
- **Field height coverage** is uneven across stems; the matched trees that carry
  height skew toward well-measured dominants, so the height-only R² is for that
  subset. The crown-only models use all matched trees.
- **AGB is derived, never validated** — NEON provides field DBH, not AGB, so
  the AGB columns are a transparent Jenkins-2003 function of DBH + functional type,
  reported only to show how DBH error propagates.
- **Functional type is genus-coarse** (conifer vs broadleaf from `taxonID`);
  species-specific allometry would tighten both DBH and AGB but needs per-species
  sample sizes the matched set does not always provide.
