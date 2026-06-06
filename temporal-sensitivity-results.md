# Ground-Truth Temporal-Sensitivity — Results

*Bounding how much the **+/-4 yr field-to-LiDAR temporal slack** in the
density-ladder ground truth moves detection recall/precision and apex-height
bias, by re-scoring on the **exact-2021** stem subset and comparing to the
[`density-ladder-sweep-results.md`](density-ladder-sweep-results.md) baseline.
Implements GitHub **issue #5**. Last run: 2026-06-05.*

---

## Why this check exists

The 2021 NEON airborne LiDAR is scored against NEON Woody Plant Vegetation
Structure (`DP1.10098.001`) field stems. NEON does **not** re-survey every plot
every year, so the ground-truth builder
([`neon_ground_truth.R`](scripts/neon_ground_truth.R)) pairs each stem with the
`apparentindividual` measurement **nearest the 2021 acquisition within +/-4 yr**
and records two columns: `meas_year` (calendar year of that measurement) and
`dist21` (`|meas_year - 2021|`). A stem measured in 2018 or 2023 has up to four
growing seasons of height growth, mortality, or stem turnover between the tape
and the laser — a mismatch that can both move a tree's apex and change whether
it is even alive/present.

**Exact-2021 coverage of the live-tree ground truth** (the
`live & is_tree & mapped` subset the sweep scores):

| Site | live-tree stems | measured in 2021 | exact-2021 share |
|---|---|---|---|
| TEAK | 483 | 233 | 48% |
| SOAP | 268 | 52 | 19% |
| **SJER** | 113 | **0** | **0%** |

**SJER is excluded** from this analysis: zero of its mapped live stems were
measured in 2021, so an exact-year re-score has no data. The SJER baseline in
the density-ladder results therefore carries the **full** +/-4 yr slack with no
way to bound it from within-site exact-year data; treat its numbers as the
upper end of temporal exposure.

## Method

Re-ran the sweep and the height validation on the exact-2021 stem subset only,
holding everything else fixed, then pooled at the **modal detection
parameters** (`chm_res = 0.5`, `vwf_a = 0.10`) per density rung and differenced
against the +/-4 yr baseline.

- **Filter.** `run_sweep.R` / `validate_heights.R` gain a `MEAS_YEAR` arg; when
  set, after the existing `live & is_tree & mapped` filter they additionally
  keep only `meas_year == MEAS_YEAR`. Unset = byte-identical to the baseline.
- **Pooling** reuses `analyze_sweep.R`'s rule: pooled
  `recall = sum(TP)/sum(n_ref)`, `precision = sum(tp_core)/sum(n_det)` with
  `tp_core = round(precision * n_det)`, and TP-weighted height RMSE — never a
  mean of per-plot rates (small plots would dominate).
- **Plot floor.** The `>= 6`-stems-per-plot rule is re-applied *after* the
  exact-year cut, so plots that fall below six 2021 stems drop out. This shrinks
  the comparison set sharply (see `n_plots` below) and is the main caveat.
- **Height bias** = `apex_z - field_h` (positive = LiDAR apex taller than the
  field tape), position-only match within 3 m at native density.

Commands (baseline CSVs are **never** overwritten — when `MEAS_YEAR` is set both
`run_sweep.R` and `validate_heights.R` auto-default their output to a distinct
`*_<YEAR>.csv` filename, so the explicit `OUT=` below is optional):

```sh
export CLAUDE_JOB_DIR=/path/to/work
# OUT= omitted -> auto-defaults to neon/<SITE>/sweep_results_2021.csv
Rscript scripts/run_sweep.R SITE=TEAK PLOTS=ALL CORES=6 TOL=4 MEAS_YEAR=2021
Rscript scripts/run_sweep.R SITE=SOAP PLOTS=ALL CORES=6 TOL=4 MEAS_YEAR=2021
Rscript scripts/validate_heights.R SITES=TEAK,SOAP TOL=3 MEAS_YEAR=2021
Rscript scripts/temporal_sensitivity.R SITES=TEAK,SOAP
```

---

## TL;DR — the temporal slack is a minor confound, not a metric-breaker

- **Recall.** At **native** density the exact-2021 re-score moves pooled recall
  by essentially **zero** (TEAK -0.0004, SOAP +0.000). On decimated rungs the
  exact-2021 subset scores **modestly lower** recall — up to **-0.05** at TEAK
  and **-0.13** at SOAP — but SOAP's exact-2021 cut is only **2 plots / 48
  stems**, too thin to read as a real density interaction.
- **Precision** tracks recall: within **+/-0.01** of baseline at native, down a
  few points on decimated rungs (TEAK worst -0.09 at rung 8; SOAP -0.15 at
  rung 4).
- **Height bias** barely moves where it can be measured: TEAK signed bias falls
  from **+5.52 m to +4.94 m (delta -0.58 m)** and height RMSE *improves* by
  ~0.2-0.7 m across rungs. SOAP's bias rises **+1.33 m** but on only **20**
  matched pairs.
- **Bound.** Across the well-populated site (TEAK), the +/-4 yr temporal slack
  accounts for **<= ~5 recall points, <= ~9 precision points, and <= ~0.6 m of
  height bias** at the modal parameters. The headline density-ladder
  conclusions (flat F1 across density; crown class as the dominant axis;
  overstory apex near-1:1) are **robust** to the temporal mismatch. SJER cannot
  be bounded this way (0% exact-2021).

---

## TEAK — +/-4 yr baseline vs exact-2021 (chm_res=0.5, vwf_a=0.10)

`n_ref` = mapped live stems in plot cores; `delta` = exact-2021 minus baseline.

| rung | cut | n_plots | n_ref | n_det | TP | recall | precision | F1 | h_rmse |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| native | baseline +/-4yr | 19 | 353 | 293 | 129 | 0.365 | 0.396 | 0.380 | 3.23 |
| native | exact-2021 | 5 | 189 | 167 | 69 | 0.365 | 0.389 | 0.377 | 2.68 |
| native | **delta** | -14 | -164 | -126 | -60 | **-0.000** | **-0.007** | -0.003 | **-0.56** |
| 8 | baseline +/-4yr | 20 | 396 | 203 | 102 | 0.258 | 0.453 | 0.328 | 3.92 |
| 8 | exact-2021 | 5 | 189 | 97 | 40 | 0.212 | 0.361 | 0.267 | 3.71 |
| 8 | **delta** | -15 | -207 | -106 | -62 | **-0.046** | **-0.092** | -0.062 | -0.22 |
| 4 | baseline +/-4yr | 20 | 396 | 201 | 101 | 0.255 | 0.453 | 0.326 | 4.18 |
| 4 | exact-2021 | 5 | 189 | 93 | 41 | 0.217 | 0.398 | 0.281 | 3.52 |
| 4 | **delta** | -15 | -207 | -108 | -60 | **-0.038** | **-0.055** | -0.046 | -0.66 |
| 2 | baseline +/-4yr | 20 | 396 | 196 | 90 | 0.227 | 0.413 | 0.293 | 4.35 |
| 2 | exact-2021 | 5 | 189 | 87 | 40 | 0.212 | 0.425 | 0.283 | 3.89 |
| 2 | **delta** | -15 | -207 | -109 | -50 | **-0.016** | **+0.012** | -0.011 | -0.46 |
| 1 | baseline +/-4yr | 20 | 396 | 184 | 92 | 0.232 | 0.451 | 0.307 | 4.58 |
| 1 | exact-2021 | 5 | 189 | 90 | 37 | 0.196 | 0.411 | 0.265 | 3.90 |
| 1 | **delta** | -15 | -207 | -94 | -55 | **-0.037** | **-0.040** | -0.041 | -0.68 |

TEAK keeps 5 plots / 189 stems after the exact-year + 6-stem floor. Native
recall is unchanged to four decimals; the largest recall drop (-0.05) is at
rung 8 and is within the noise of a 5-plot subset. Height RMSE is uniformly
**lower** on exact-2021 — consistent with removing trees whose field height is
years out of date.

## SOAP — +/-4 yr baseline vs exact-2021 (chm_res=0.5, vwf_a=0.10)

| rung | cut | n_plots | n_ref | n_det | TP | recall | precision | F1 | h_rmse |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| native | baseline +/-4yr | 18 | 232 | 318 | 116 | 0.500 | 0.330 | 0.398 | 3.60 |
| native | exact-2021 | 2 | 48 | 64 | 24 | 0.500 | 0.344 | 0.407 | 2.98 |
| native | **delta** | -16 | -184 | -254 | -92 | **+0.000** | **+0.014** | +0.010 | **-0.62** |
| 8 | baseline +/-4yr | 18 | 232 | 155 | 72 | 0.310 | 0.432 | 0.361 | 4.34 |
| 8 | exact-2021 | 2 | 48 | 30 | 9 | 0.188 | 0.300 | 0.231 | 3.64 |
| 8 | **delta** | -16 | -184 | -125 | -63 | **-0.123** | **-0.132** | -0.131 | -0.70 |
| 4 | baseline +/-4yr | 18 | 232 | 160 | 84 | 0.362 | 0.513 | 0.424 | 4.41 |
| 4 | exact-2021 | 2 | 48 | 30 | 11 | 0.229 | 0.367 | 0.282 | 4.19 |
| 4 | **delta** | -16 | -184 | -130 | -73 | **-0.133** | **-0.146** | -0.142 | -0.22 |
| 2 | baseline +/-4yr | 18 | 232 | 155 | 80 | 0.345 | 0.458 | 0.393 | 4.04 |
| 2 | exact-2021 | 2 | 48 | 29 | 12 | 0.250 | 0.345 | 0.290 | 4.04 |
| 2 | **delta** | -16 | -184 | -126 | -68 | **-0.095** | **-0.113** | -0.104 | +0.00 |
| 1 | baseline +/-4yr | 18 | 232 | 152 | 78 | 0.336 | 0.487 | 0.398 | 4.88 |
| 1 | exact-2021 | 2 | 48 | 31 | 15 | 0.313 | 0.452 | 0.369 | 4.11 |
| 1 | **delta** | -16 | -184 | -121 | -63 | **-0.024** | **-0.035** | -0.028 | -0.77 |

SOAP collapses to **2 plots / 48 stems** at exact-2021 (only 19% of its stems
are 2021-measured, and few plots clear the 6-stem floor once cut). Native recall
is again unchanged (0.500 both cuts); the larger decimated-rung swings (up to
-0.13) ride on those 2 plots and should **not** be read as a density
interaction. The native and rung-1 rows — where the subset is largest — show
deltas of <= ~0.03, consistent with TEAK.

## Apex-vs-field height bias (native, res=0.5, a=0.10, match <= 3 m)

Signed bias = `apex_z - field_h` (positive = LiDAR taller than the field tape).
Baseline = +/-4 yr pairs restricted to TEAK,SOAP; exact-2021 = `meas_year==2021`.

| Site | cut | n pairs | bias (m) | RMSE (m) |
|---|---|---:|---:|---:|
| TEAK | baseline +/-4yr | 144 | +5.52 | 11.26 |
| TEAK | exact-2021 | 65 | +4.94 | 10.20 |
| TEAK | **delta** | -79 | **-0.58** | **-1.06** |
| SOAP | baseline +/-4yr | 107 | +2.72 | 6.88 |
| SOAP | exact-2021 | 20 | +4.05 | 7.06 |
| SOAP | **delta** | -87 | **+1.33** | **+0.18** |
| TEAK+SOAP | baseline +/-4yr | 251 | +4.32 | 9.64 |
| TEAK+SOAP | exact-2021 | 85 | +4.73 | 9.55 |
| TEAK+SOAP | **delta** | -166 | **+0.41** | **-0.09** |

The pooled bias moves only **+0.41 m** and pooled RMSE is flat. TEAK (the larger
exact-2021 sample) shows the *expected* sign — removing stale field heights
nudges bias **down** -0.58 m and RMSE down -1.06 m. SOAP's +1.33 m rise is on
20 pairs and is dominated by sampling, not a real temporal trend.

---

## Bounding statement (issue #5)

At the modal detection parameters, the +/-4 yr field-to-LiDAR temporal slack in
the ground truth moves the headline metrics by **at most**, on the
well-populated site (TEAK, exact-2021 = 5 plots / 189 stems):

- **recall:** <= ~0.05 (0.000 at native);
- **precision:** <= ~0.09 (within 0.01 at native);
- **apex-height bias:** <= ~0.6 m, in the direction that *improves* agreement.

SOAP corroborates the near-zero native effect (recall delta 0.000) but its
exact-2021 cut (2 plots / 48 stems) is too small to bound the decimated rungs.
**SJER cannot be bounded at all** (0% exact-2021 coverage). The temporal
mismatch is therefore a **minor, recall-conservative confound** that does not
overturn any density-ladder conclusion; if anything the baseline slightly
*over-states* height bias on TEAK by carrying years-old field heights.

### Caveats

- The 6-stem plot floor, re-applied after the exact-year cut, is what shrinks
  the comparison (TEAK 20 -> 5 plots, SOAP 18 -> 2). The decimated-rung deltas
  ride on these small subsets; only the native and largest-subset rows are
  well-supported.
- `meas_year` is the year of the *nearest* `apparentindividual` record, not a
  guarantee of a 2021 stem-mapping; position still comes from the (stable)
  mapped polar offset.
- SJER (0% exact-2021) and any 0%-coverage acquisition year remain exposed to
  the full slack — its density-ladder numbers should be read as the upper bound
  of temporal exposure.
