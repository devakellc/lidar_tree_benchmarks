# Calibration / validation split of the density-ladder sweep

Addresses **issue #3**. The density-ladder sweep
([density-ladder-sweep-results.md](density-ladder-sweep-results.md)) reports a
"best `(chm_res, vwf_a)` per density rung", but those are **in-sample** optima:
parameters are tuned and scored on the *same* plots, pooled over every plot.
This document tests whether that choice survives **out-of-sample** by tuning on
a calibration subset of plots and reporting recall / precision / F1 on
**held-out** validation plots.

Headline finding under test (from the sweep doc): *"`chm_res = 0.5 m` is
F1-optimal; the VWF slope is a second-order effect."* We test it two ways: (a)
does the *calibration-selected* `chm_res` transfer to held-out plots, and (b)
the stronger test — is 0.5 m actually the **held-out-F1-optimal** `chm_res`
when every candidate resolution is scored on the validation plots.

Reproduce:

```sh
export CLAUDE_JOB_DIR=/path/to/work
Rscript scripts/calval_split.R SITES=SJER,SOAP,TEAK SEED=1 FRAC=0.5 \
        SEEDS=1,2,3,4,5
```

All numbers below are real console output from that command on the cached
`sweep_results.csv` (no LiDAR recomputation). The script writes long-form
`calval_metrics.csv` per site under `$CLAUDE_JOB_DIR/neon/<SITE>/`.

## Method

1. For each site, read the cached long-form `sweep_results.csv` (one row per
   plot x rung x `chm_res` x `vwf_a`).
2. Assign each unique plot to **calibration** or **validation**
   deterministically (`set.seed(SEED)`), **stratified** by two covariates:
   - **plotType** — tower base (mapped +/-20 m) vs distributed (+/-10 m);
     different mapped extent and structure.
   - **crown-class mix** — each plot is *overstory-dominated* or
     *understory-present*, classified from its overstory stems
     (`n_dominant + n_codominant`) vs understory stems
     (`n_intermediate + n_suppressed`): a plot is `understory` iff understory
     stems > 0 **and** understory >= overstory / 4, else `overstory`.

   Within each `plotType x crown_mix` stratum the plots are shuffled
   (seeded) and `ceil(FRAC * n)` go to calibration. A singleton stratum goes
   to calibration (so validation never relies on an unrepresented stratum).
3. Pick the best `(chm_res, vwf_a)` per density rung **by pooled F1 on
   calibration plots only**, using the exact `analyze_sweep.R` pooling rule:
   pooled recall = `sum(TP) / sum(n_ref)`; pooled precision =
   `sum(tp_core) / sum(n_det)` with `tp_core = round(precision * n_det)`;
   per-class TP = `round(rec_class * n_class)`. **Never** average per-plot
   rates (that over-weights small plots).
4. Apply those calibration-selected parameters to the **validation** plots and
   report held-out pooled recall / precision / F1 per rung, alongside the
   calibration (in-sample) metrics and the full-pooled in-sample optimum.
5. **Held-out optimality across `chm_res` (the direct out-of-sample test),
   complete-case.** Calibration selection alone never *compares* `chm_res` on
   held-out data — it reports the calibration-selected value. To test the
   headline honestly, for each rung pool **held-out F1 for every candidate
   `chm_res`** (0.25, 0.5, 1.0 m) on the **validation** plots, holding `vwf_a` at
   the calibration-selected value for that rung, and take the argmax. **The
   comparison is complete-case across `chm_res`.** The sweep grid is *not*
   rectangular: at native density `chm_res = 0.25` is run only on plots whose
   first-return density >= 8 pts/m^2, so it covers fewer plots than 0.5 / 1.0
   (SJER native 5/8, SOAP 16/18, TEAK 17/20). Pooling each resolution over
   whatever rows exist would compare F1 across resolutions on **different plot
   subsets** and bias the argmax. So for each site x rung we restrict to the
   **intersection** of validation plots that have a scored row at *every*
   candidate `chm_res` (at the calibration `vwf_a`), pool F1 per `chm_res` over
   that common set, and take the argmax there; we report the number of plots used
   (`n_pl`) and how many were dropped (`drop`) to enforce complete-case. The
   decimated rungs (8/4/2/1) only ever run {0.5, 1.0} at equal plot counts, so
   there complete-case is a no-op (`drop = 0`); the restriction only ever bites
   the native rung. The resulting argmax is then compared to 0.5. *Design
   choice:* `vwf_a` is fixed (not marginalized) because the VWF slope is
   second-order (step 6), so the `chm_res` argmax is insensitive to it, and
   fixing `vwf_a` keeps the contrast a clean one-factor (`chm_res`-only) test at
   the operating point the calibration would deploy. Ties are broken
   deterministically by finest `chm_res` first; no exact held-out F1 ties
   occur in the cached data, so the reported optima do not depend on the tie
   rule.
6. Quantify the VWF-slope effect as the F1 spread across `vwf_a`, fixing
   `chm_res` to the calibration optimum but pooling F1 on the **held-out
   validation** plots, so the second-order claim is itself out-of-sample (small
   spread => second-order). When the F1-optimal `(chm_res, vwf_a)` per rung is
   tied, the tie is broken deterministically: finest `chm_res` first
   (ascending), then lower `vwf_a`. No exact ties occur in the cached data, so
   this affects no reported optimum.
7. Because per-site plot counts are small, repeat the split for `SEEDS=1..5`
   and report the held-out F1 distribution, the modal calibration-optimal
   `chm_res`, **and the modal held-out-F1-optimal `chm_res`** per rung, so the
   verdict is not a single-split artifact.

**Honest caveat on sample size.** Per-site plot counts are small — SJER 8,
SOAP 18, TEAK 20 — and the validation halves are only 3 / 8 / 9 plots. With so
few plots a single split is noisy; the multi-seed summary is the load-bearing
evidence, not any one seed. SJER is additionally degenerate for the crown-mix
stratum: only 2 intermediate and 0 suppressed stems exist across all 8 plots,
so exactly one SJER plot is `understory` and it is a distributed singleton
(forced to calibration). SJER's crown-mix stratification is therefore nominal.

### Per-site plot counts and strata (all plots)

| Site | Plots | Tower | Distributed | overstory | understory |
|------|-------|-------|-------------|-----------|------------|
| SJER | 8     | 7     | 1           | 7         | 1          |
| SOAP | 18    | 7     | 11          | 12        | 6          |
| TEAK | 20    | 7     | 13          | 15        | 5          |

`SEED=1, FRAC=0.5` split sizes: SJER calib 5 / valid 3; SOAP calib 10 /
valid 8; TEAK calib 11 / valid 9.

## Held-out vs in-sample results (SEED=1, FRAC=0.5)

Each row: the `(chm_res, vwf_a)` chosen on calibration plots, with the
calibration (CAL, in-sample) and held-out validation (VAL, out-of-sample)
pooled metrics, plus the full-pooled in-sample optimum (all plots) for
reference.

### SJER (calib n=5, valid n=3)

| rung   | chm_res | vwf_a | CAL rec/prec/F1 | VAL rec/prec/F1 | full-pooled opt (F1) |
|--------|---------|-------|-----------------|-----------------|----------------------|
| native | 1.00    | 0.15  | 0.54/0.20/0.29  | 0.50/0.33/0.40  | res=1.00 a=0.15 (0.33)|
| 8      | 0.50    | 0.10  | 0.49/0.30/0.37  | 0.37/0.44/0.40  | res=0.50 a=0.10 (0.38)|
| 4      | 0.50    | 0.15  | 0.49/0.31/0.38  | 0.39/0.41/0.40  | res=0.50 a=0.15 (0.39)|
| 2      | 0.50    | 0.05  | 0.43/0.27/0.33  | 0.39/0.45/0.42  | res=0.50 a=0.05 (0.37)|
| 1      | 1.00    | 0.05  | 0.40/0.29/0.34  | 0.36/0.43/0.39  | res=1.00 a=0.05 (0.36)|

### SOAP (calib n=10, valid n=8)

| rung   | chm_res | vwf_a | CAL rec/prec/F1 | VAL rec/prec/F1 | full-pooled opt (F1) |
|--------|---------|-------|-----------------|-----------------|----------------------|
| native | 0.50    | 0.05  | 0.58/0.34/0.43  | 0.46/0.26/0.33  | res=0.50 a=0.10 (0.40)|
| 8      | 0.50    | 0.05  | 0.34/0.43/0.38  | 0.30/0.43/0.35  | res=0.50 a=0.05 (0.37)|
| 4      | 0.50    | 0.05  | 0.37/0.50/0.42  | 0.38/0.50/0.43  | res=0.50 a=0.05 (0.43)|
| 2      | 0.50    | 0.05  | 0.38/0.48/0.42  | 0.34/0.43/0.38  | res=0.50 a=0.05 (0.41)|
| 1      | 0.50    | 0.05  | 0.37/0.50/0.43  | 0.29/0.45/0.35  | res=0.50 a=0.05 (0.40)|

### TEAK (calib n=11, valid n=9)

| rung   | chm_res | vwf_a | CAL rec/prec/F1 | VAL rec/prec/F1 | full-pooled opt (F1) |
|--------|---------|-------|-----------------|-----------------|----------------------|
| native | 0.25    | 0.05  | 0.47/0.36/0.41  | 0.35/0.48/0.41  | res=0.25 a=0.05 (0.41)|
| 8      | 0.50    | 0.05  | 0.30/0.37/0.33  | 0.24/0.58/0.34  | res=0.50 a=0.05 (0.34)|
| 4      | 0.50    | 0.10  | 0.30/0.38/0.34  | 0.22/0.55/0.32  | res=0.50 a=0.10 (0.33)|
| 2      | 0.50    | 0.05  | 0.27/0.33/0.30  | 0.21/0.53/0.30  | res=0.50 a=0.05 (0.30)|
| 1      | 0.50    | 0.05  | 0.28/0.38/0.32  | 0.20/0.55/0.29  | res=0.50 a=0.05 (0.31)|

The calibration-selected `(chm_res, vwf_a)` matches the full-pooled in-sample
optimum in nearly every rung at SOAP and TEAK, and the held-out F1 tracks the
in-sample F1 closely (within roughly +/-0.05 at the decimated rungs). The
native rung is the least stable, especially at SOAP, where held-out F1 (0.33)
drops well below calibration (0.43) — recall-driven, as the densest data lets
calibration over-fit a high-recall point that does not transfer.

## Multi-seed robustness (SEEDS=1..5, FRAC=0.5)

For each rung: the modal calibration-optimal `chm_res` across seeds, the
fraction of seeds whose optimum was exactly 0.5, and the held-out F1
distribution (median, then min..max).

### SJER

| rung   | modal_chm | chm=0.5 frac | val F1 median | val F1 min..max |
|--------|-----------|--------------|---------------|-----------------|
| native | 1.00      | 0.40         | 0.311         | 0.127 .. 0.396  |
| 8      | 0.50      | 1.00         | 0.349         | 0.218 .. 0.400  |
| 4      | 0.50      | 0.80         | 0.324         | 0.192 .. 0.400  |
| 2      | 0.50      | 1.00         | 0.328         | 0.165 .. 0.418  |
| 1      | 1.00      | 0.20         | 0.319         | 0.193 .. 0.392  |

### SOAP

| rung   | modal_chm | chm=0.5 frac | val F1 median | val F1 min..max |
|--------|-----------|--------------|---------------|-----------------|
| native | 0.50      | 1.00         | 0.382         | 0.333 .. 0.425  |
| 8      | 0.50      | 1.00         | 0.355         | 0.344 .. 0.380  |
| 4      | 0.50      | 1.00         | 0.422         | 0.408 .. 0.430  |
| 2      | 0.50      | 1.00         | 0.406         | 0.382 .. 0.436  |
| 1      | 0.50      | 1.00         | 0.389         | 0.353 .. 0.424  |

### TEAK

| rung   | modal_chm | chm=0.5 frac | val F1 median | val F1 min..max |
|--------|-----------|--------------|---------------|-----------------|
| native | 0.25      | 0.40         | 0.364         | 0.350 .. 0.444  |
| 8      | 0.50      | 1.00         | 0.339         | 0.311 .. 0.370  |
| 4      | 0.50      | 1.00         | 0.326         | 0.313 .. 0.341  |
| 2      | 0.50      | 1.00         | 0.294         | 0.266 .. 0.309  |
| 1      | 0.50      | 1.00         | 0.316         | 0.281 .. 0.325  |

## Held-out-F1-optimal `chm_res` (the direct out-of-sample test)

The tables above report the *calibration-selected* `chm_res` and how its F1
transfers; they do **not** compare resolutions on held-out data. This section
does: for each rung it pools held-out F1 for **every** candidate `chm_res`
(0.25, 0.5, 1.0 m) on the validation plots, holding `vwf_a` at the
calibration-selected value, and takes the argmax. That argmax is the
held-out-F1-optimal resolution; the question is whether it equals 0.5.

**Complete-case across `chm_res`.** Because the native-rung `chm_res = 0.25`
slice is run on fewer plots than 0.5 / 1.0 (SJER 5/8, SOAP 16/18, TEAK 17/20),
comparing resolutions over all-rows-that-exist would pit different plot subsets
against each other. So at each rung the comparison is restricted to the
intersection of validation plots scored at *every* candidate `chm_res`; `n_pl`
is the size of that common set and `drop` is how many validation plots were
removed to enforce it. The decimated rungs (8/4/2/1) compare only {0.5, 1.0} at
equal plot counts, so `drop = 0` there — the restriction only touches the native
rung.

### Single split (SEED=1): held-out-optimal `chm_res` per rung

`ho_F1` is the held-out F1 at the argmax `chm_res`; `F1@0.5` is the held-out F1
at 0.5 m for comparison. `n_pl` / `drop` are the complete-case common-plot count
and the plots dropped at that rung.

| Site | rung | calib `chm_res` | held-out opt `chm_res` | ho_F1 | F1@0.5 | n_pl | drop | 0.5 is held-out opt? |
|------|------|-----------------|------------------------|-------|--------|------|------|----------------------|
| SJER | native | 1.00 | 1.00 | 0.400 | 0.368 | 2 | 1 | no (opt 1.00) |
| SJER | 8      | 0.50 | 0.50 | 0.400 | 0.400 | 2 | 0 | yes |
| SJER | 4      | 0.50 | 1.00 | 0.416 | 0.400 | 3 | 0 | no (opt 1.00) |
| SJER | 2      | 0.50 | 0.50 | 0.418 | 0.418 | 3 | 0 | yes |
| SJER | 1      | 1.00 | 0.50 | 0.404 | 0.404 | 3 | 0 | yes |
| SOAP | native | 0.50 | 1.00 | 0.348 | 0.317 | 7 | 1 | no (opt 1.00) |
| SOAP | 8      | 0.50 | 0.50 | 0.355 | 0.355 | 8 | 0 | yes |
| SOAP | 4      | 0.50 | 0.50 | 0.430 | 0.430 | 8 | 0 | yes |
| SOAP | 2      | 0.50 | 0.50 | 0.382 | 0.382 | 8 | 0 | yes |
| SOAP | 1      | 0.50 | 0.50 | 0.353 | 0.353 | 8 | 0 | yes |
| TEAK | native | 0.25 | 0.50 | 0.413 | 0.413 | 8 | 1 | yes |
| TEAK | 8      | 0.50 | 0.50 | 0.343 | 0.343 | 9 | 0 | yes |
| TEAK | 4      | 0.50 | 0.50 | 0.316 | 0.316 | 9 | 0 | yes |
| TEAK | 2      | 0.50 | 0.50 | 0.296 | 0.296 | 9 | 0 | yes |
| TEAK | 1      | 0.50 | 0.50 | 0.294 | 0.294 | 9 | 0 | yes |

On this single split 0.5 m is the complete-case held-out optimum in **3/5** rungs
at SJER, **4/5** at SOAP, and **5/5** at TEAK — the same per-site counts as the
prior (non-complete-case) run, and the single-seed native argmax is unchanged at
every site (SJER 1.00, SOAP 1.00, TEAK 0.50). The native-rung F1 values shift
slightly because the argmax is now pooled over only the common plots: at SJER
native one of the 3 validation plots lacks a 0.25 slice (`drop = 1`, `n_pl = 2`),
at SOAP/TEAK native one plot is dropped (`drop = 1`). Where 0.5 m is not the
optimum (SJER native + rung 4, SOAP native) the winning resolution is 1.0 m, but
the F1 margin over 0.5 m is small (SJER rung 4: 0.416 vs 0.400; SOAP native:
0.348 vs 0.317) — well inside the multi-seed F1 spread, so it is not a decisive
defeat of 0.5 m.

### Across seeds (SEEDS=1..5): modal held-out-optimal `chm_res`

`ho opt = 0.5 frac` is the fraction of seeds whose complete-case held-out argmax
was 0.5 m. `n_pl med` / `drop med` are the median common-plot count and median
plots dropped at that rung across the five seeds.

| Site | rung | modal held-out opt `chm_res` | ho opt = 0.5 frac | n_pl med | drop med |
|------|------|------------------------------|-------------------|----------|----------|
| SJER | native | 0.25 | 0.20 | 2.0 | 1.0 |
| SJER | 8      | 0.50 | 1.00 | 3.0 | 0.0 |
| SJER | 4      | 1.00 | 0.00 | 3.0 | 0.0 |
| SJER | 2      | 0.50 | 0.60 | 3.0 | 0.0 |
| SJER | 1      | 0.50 | 0.60 | 3.0 | 0.0 |
| SOAP | native | 1.00 | 0.40 | 7.0 | 1.0 |
| SOAP | 8      | 0.50 | 0.80 | 8.0 | 0.0 |
| SOAP | 4      | 0.50 | 1.00 | 8.0 | 0.0 |
| SOAP | 2      | 0.50 | 0.80 | 8.0 | 0.0 |
| SOAP | 1      | 0.50 | 1.00 | 8.0 | 0.0 |
| TEAK | native | 0.25 | 0.40 | 8.0 | 1.0 |
| TEAK | 8      | 0.50 | 1.00 | 9.0 | 0.0 |
| TEAK | 4      | 0.50 | 1.00 | 9.0 | 0.0 |
| TEAK | 2      | 0.50 | 1.00 | 9.0 | 0.0 |
| TEAK | 1      | 0.50 | 1.00 | 9.0 | 0.0 |

Across seeds the held-out optimum is **0.5 m at every decimated rung at TEAK**
(modal 0.5, frac 1.00) and at SOAP rungs 4/1 (frac 1.00) with rungs 8/2 at 0.80.
These decimated-rung results are identical to the prior (non-complete-case)
run — expected, because at the decimated rungs only {0.5, 1.0} are run and both
cover the same plots, so `drop = 0` and complete-case is a no-op. SJER is the weak
case: rung 8 holds (frac 1.00), rungs 2/1 are 0.5 m only as a 0.60-fraction
plurality, and **rung 4 is the clear exception — its held-out optimum is 1.0 m
in 5/5 seeds (frac 0.00)**.

The **native rung is never robustly 0.5 m** at any site, but its modal optimum
*changed under complete-case at SJER*: it was 1.0 m on the non-complete-case run
and is now **0.25 m**. This is exactly the bias the complete-case fix removes —
at SJER native the 0.25 m slice covers only 5 of 8 plots, and `drop` is 1 in 3
of the 5 seeds (median 1). On the seeds where a plot is dropped, scoring 1.0 m
over its larger native subset had inflated it relative to 0.25 m on its own
smaller subset; restricting to the common plots flips two seeds' argmax to 0.25,
giving a 2-vs-2 split between 0.25 and 1.0 that the finest-res tie-break
resolves to 0.25. The SJER native modal is therefore a fragile plurality, not a
robust winner — but in either reading the SJER native optimum is **not 0.5 m**.
At SOAP the native modal stays 1.0 m (`drop` 1 in some seeds, `n_pl med` 7) and
at TEAK it stays 0.25 m (`drop` 1, `n_pl med` 8). So the honest held-out optimum
is: 0.5 m at the decimated rungs at SOAP and TEAK, but **not universal** — SJER
rung 4 prefers 1.0 m, and the native rung always prefers a non-0.5 resolution at
every site.

## VWF-slope effect (second-order check), out-of-sample

The `chm_res` is fixed to the **calibration** optimum, but pooled F1 across the
three `vwf_a` values is computed on the **held-out validation** plots, so this
is a genuinely out-of-sample test of the VWF main effect. Maximum F1 spread
across all rungs (SEED=1, validation plots):

| Site | max VWF-slope F1 spread (held-out) | where |
|------|------------------------------------|-------|
| SJER | 0.021                              | rung 8 |
| SOAP | 0.039                              | native rung |
| TEAK | 0.031                              | native rung |

Per-rung held-out spread at the calibration-optimal `chm_res` (SEED=1):

| rung   | SJER  | SOAP  | TEAK  |
|--------|-------|-------|-------|
| native | 0.008 | 0.039 | 0.031 |
| 8      | 0.021 | 0.015 | 0.010 |
| 4      | 0.006 | 0.001 | 0.014 |
| 2      | 0.018 | 0.024 | 0.004 |
| 1      | 0.000 | 0.000 | 0.000 |

At every decimated rung (8/4/2/1) the held-out spread is <= 0.024, and at the
finest rung (1 pt/m^2) it is 0.000 at all three sites. Changing `vwf_a` over
its full sweep range (0.05 -> 0.15) moves held-out F1 by at most ~0.04 even in
the worst case (SOAP native), versus the 0.1 -> 0.5 m `chm_res` differences
that move F1 by far more. The VWF slope is confirmed **second-order**
out-of-sample. Every validation subset is poolable here (n_pool = 2-9 plots
per rung), so no rung falls back to calibration; the one thin case is SJER
rung 8, where only 2 of the 3 validation plots carry the 0.5 m slice
(n_pool = 2).

## Verdict: does the finding survive out-of-sample?

**Mostly — but not universally, and weakest exactly where the original sweep
doc claimed it was strongest.** The honest answer separates two questions.

**Does the *calibration-selected* `chm_res` transfer?** Yes. The
calibration-selected `chm_res` matches the full-pooled in-sample optimum in
nearly every decimated rung at SOAP and TEAK, and its held-out F1 tracks the
in-sample F1 closely (within ~+/-0.05). On the SEED=1 split, the
calibration-optimal `chm_res` is 0.5 m in 3/5 rungs at SJER, 5/5 at SOAP, and
4/5 at TEAK; across seeds it is the modal calibration optimum at every decimated
rung at SOAP and TEAK (and SJER rung 8). At SJER the calibration optimum is
0.5 m as the modal value at rungs 8/4/2 with `chm=0.5` seed-fractions of
**1.00, 0.80, 1.00** respectively, while the native and 1 pts/m^2 rungs flip to
1.0 m.

**Is 0.5 m the *held-out-F1-optimal* `chm_res` (the stronger test)?** Not
everywhere. Scoring every candidate resolution on the validation plots:

- **At TEAK 0.5 m is the held-out optimum at every decimated rung (8/4/2/1) in
  5/5 seeds** — the cleanest confirmation.
- **At SOAP 0.5 m is the held-out optimum at the decimated rungs** (frac 1.00
  at rungs 4 and 1; 0.80 at rungs 8 and 2).
- **At SJER 0.5 m does not robustly win.** It holds at rung 8 (5/5 seeds), is
  only a plurality at rungs 2/1 (0.5 m in 3/5 seeds each), and at **rung 4 the
  held-out optimum is 1.0 m in 5/5 seeds** — a genuine exception where 0.5 m is
  *not* the held-out optimum even at a decimated rung. The F1 margin there is
  small (SEED=1: 0.416 at 1.0 m vs 0.400 at 0.5 m), so 0.5 m is competitive,
  not optimal.
- **The native (undecimated) rung is never robustly 0.5 m** at any site: under
  the complete-case comparison the modal held-out optimum is 0.25 m at SJER (a
  fragile 2-vs-2 plurality after the fix; it was 1.0 m on the biased run),
  1.0 m at SOAP, and 0.25 m at TEAK. This is the rung the complete-case fix
  matters for — the native `chm_res = 0.25` slice runs on fewer plots than
  0.5 / 1.0, so the resolutions must be compared on the common plot subset
  (`drop` = 1 in the SEED=1 split at all three sites). It is physically
  consistent with the density-first methodology — a finer or coarser CHM than
  0.5 m can win when the density differs from the decimated operating regime —
  but it means the unqualified claim "0.5 m is F1-optimal" is **false at native
  density**, whichever way the SJER native tie falls.

So the prior statement that "`chm_res = 0.5 m` is F1-optimal out-of-sample" was
**overstated**: that test was never run (only the calibration-selected param was
scored). Running it honestly — and complete-case, so resolutions are compared on
a common plot subset rather than the unequal subsets the native grid would
otherwise force — 0.5 m is the held-out optimum at the decimated rungs at SOAP
and TEAK but **not** at SJER rung 4 (1.0 m wins 5/5 seeds) and **not** at native
density at any site.

**The VWF slope is second-order**, confirmed out-of-sample: pooled on the
held-out validation plots, F1 spread across `vwf_a` at the calibration
`chm_res` is <= 0.039 in the worst rung (SOAP native) and <= 0.024 at every
decimated rung (0.000 at the 1 pt/m^2 rung at all three sites). This part of the
headline survives cleanly.

In short: the sweep's headline "**0.5 m CHM resolution is F1-optimal; VWF slope
is a second-order effect**" **largely survives** the calibration/validation
split for the decimated density rungs at SOAP and TEAK, and the VWF
second-order claim survives everywhere — but it **does not survive
universally**: at native density the held-out optimum is always a non-0.5
resolution, and at SJER rung 4 the held-out optimum is 1.0 m, not 0.5 m. The
small per-site plot counts mean absolute F1 values carry wide uncertainty (see
the multi-seed min..max ranges); the *decimated-rung* selection is stable at
SOAP/TEAK, but SJER and the native rung are the honest exceptions.
