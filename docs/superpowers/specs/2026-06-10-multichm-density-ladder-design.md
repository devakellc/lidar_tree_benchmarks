# multichm arm on the 3-site density ladder — design (issue #37)

*Status: approved 2026-06-10. Adds a `multichm` treetop detector arm to the
canonical NEON density-ladder pipeline (`run_sweep.R` →
`density-ladder-sweep-results.md`), which today scores **CHM-VWF only**.*

## Problem

`multichm` (Eysn-style multi-layer CHM local maxima, `lidRplugins::multichm`) is
the strongest classical arm on the SOAP **model** benchmark (F1 0.44 vs CHM-VWF
0.38 at native; beats CHM-VWF at every rung). But the canonical **density-ladder**
doc — the one operators read for "how does detection respond to point density" —
scores only CHM-VWF (lasR `local_maximum_raster` on a `pit_fill` CHM). The model
benchmark and the density ladder also differ in clip provenance (frozen-clip
`pitfree` vs `prepare_clip` `pit_fill`), so the model-benchmark multichm numbers
cannot simply be quoted into the density-ladder story. We need a multichm arm
built on the **same `prepare_clip` lasR-based path** as the cached CHM-VWF
`sweep_results.csv`, so the comparison is internally consistent.

## Scope

In scope: a multichm detector arm over the 3-site ladder (SJER, SOAP, TEAK) with
the existing plots, rungs, scoring, and pooling; a head-to-head vs the cached
CHM-VWF `sweep_results.csv`; pooled tables, a figure, and a results-doc addendum.

Out of scope (per the issue): a multichm internal-parameter sweep
(`layer_thickness`, `dist_2d`/`dist_3d`). Fixed literature defaults unless a clear
main effect appears.

## Architecture

Two new scripts, mirroring the `run_sweep.R` / `analyze_sweep.R` split, plus one
unit test. No existing script is modified except the README table and the results
doc.

### 1. `scripts/detect_multichm_sweep.R` (detection)

Parallel to `run_sweep.R`, **not** the frozen-clip model-benchmark path.

- Reuses `sweep_lib.R::prepare_clip` — the *same* clip → `decimate_points` →
  `normalize_height` the cached CHM-VWF `sweep_results.csv` was built on. This is
  the load-bearing choice: identical clip provenance keeps the density-ladder
  comparison honest.
- Per plot × rung {native, 8, 4, 2, 1 pts/m²}, with the existing **no-upsampling
  guard** (rung target vs the plot's native all-return density) and
  **plot-type-aware core** (`plot_half`: tower ±20 m, distributed ±10 m).
- Detector: `det_multichm_run(las, res, a)` — a small, unit-testable wrapper that
  runs `lidR::locate_trees(las, lidRplugins::multichm(res = res, ws =
  ws_factory(a)))` and returns the lowercase `x, y, z` detection contract (the
  same shape `score_plot` consumes; mirrors `detect_lidrplugins_sweep.R::det_multichm`
  - `.tops_to_det`). `res = 0.25` if first-return density ≥ 8 else `0.5`; `a =
  0.10`. multichm geometry is 2-D, so the apex `z` is read from the tops' `Z`
  attribute (the same fallback `.tops_to_det` uses).
- Scoring: `score_plot(stems, det, tol_xy = TOL=4, core_*)` — unchanged greedy
  height-gated 1:1 matcher, crown-class + height-band recall.
- Output: `work/neon/<SITE>/multichm_sweep_results.csv`, **one row per
  (plot × rung)** (multichm has no `chm_res`/`vwf_a` axes here — single
  density-derived `res`, fixed `a`). Columns: identifiers (`site, plot, plotType,
  rung, pdens, frdens, chm_res, n_apex`) + the full `score_plot` block, plus
  `tp_core = round(precision * n_det)` for poolers.
- Args: `SITE`, `PLOTS`, `CORES`, `TOL`, `A` (KEY=VALUE), matching the sibling
  scripts. `run_main()` guarded by `if (sys.nframe() == 0L)` so the file can be
  `source()`d by the test without executing.

Dependency note: this arm needs only `lidR` + `lidRplugins`; the CHM-VWF baseline
comes from the **cached** `sweep_results.csv`, so lasR `pre-devel` is *not*
required to run the multichm arm (CRAN lasR 0.21.0 is fine — `prepare_clip` uses
only lidR, and we never call `local_maximum_raster(ws = f)`).

### 2. `scripts/analyze_multichm_sweep.R` (pooling + comparison)

- Loads `multichm_sweep_results.csv` for the site(s) and the cached CHM-VWF
  `sweep_results.csv`.
- CHM-VWF comparison cell: for each cached `(plot, rung)`, select the row whose
  `chm_res` equals that row's own density-derived value (`0.25` if `frdens ≥ 8`
  else `0.5`) and `vwf_a == 0.10` — the **same discipline** multichm follows, so
  the head-to-head is apples-to-apples (both pick `res` by the `frdens ≥ 8` rule
  and use `a = 0.10`).
- Restricts both arms to the **common (plot, rung) set** before pooling (a 2-arm
  equal-set guard).
- Pools per rung with the canonical `model_bench_lib.R::pool` (sum TP / sum n_ref;
  per-class TP recovered as `round(rec_<cls> · n_<cls>)`; understory =
  intermediate + suppressed; height bands when present). Deltas are computed
  **between pooled rates**, never as a mean of per-row deltas.
- Emits: per-site pooled multichm table by rung (recall/precision/F1 + crown class
  - height band); head-to-head Δ table (multichm − CHM-VWF) per rung; a cross-site
  roll-up; a figure (multichm vs CHM-VWF recall/F1 across the measured
  first-return density); a markdown fragment for the doc addendum. Writes
  `work/neon/<SITE>/multichm_summary_by_rung.csv` and
  `work/neon/multichm_vs_chmvwf.csv`.
- Args: `SITES` (comma list, default `SJER,SOAP,TEAK`).

### 3. `tests/testthat/test-multichm-sweep.R`

Mirrors `test-lidrplugins-extractors.R`: `source()` the detect script, assert
`det_multichm_run(synth_las_normalized(), res = 0.5, a = 0.10)` returns a base
`data.frame` (not `sf`) with exactly numeric `x, y, z`; and that a no-tree flat
cloud yields a 0-row frame (or NULL) honoring the contract. No NEON data required.

## Data flow

```text
NEON catalog ─prepare_clip(rung)→ normalized .laz ─readLAS→ LAS
   → locate_trees(multichm(res=f(frdens), ws=a·h+3 clamp[3,5]))
   → x,y,z tops → score_plot vs field stems (core, TOL=4)
   → multichm_sweep_results.csv (one row / plot×rung)
analyze: + cached sweep_results.csv (CHM-VWF, same res-rule, a=0.10)
   → equal-set (plot,rung) → pool per arm per rung → Δ tables + figure + md
```

## Validation / acceptance

- `Rscript tests/run_tests.R` green (includes the new extractor test).
- 3-site `multichm_sweep_results.csv` produced from the cached NEON data
  (`CLAUDE_JOB_DIR` = main checkout `work`).
- Pooled multichm-vs-CHM-VWF table + figure rendered; the addendum in
  `results/density-ladder-sweep-results.md` reports the **real** numbers and is
  rumdl-clean (80-char prose).
- README script table + reproduce block list both new scripts.

## Risks / mitigations

- *Decimation realization differs from the cached run* (`prepare_clip` is
  unseeded, exactly like `run_sweep.R`). The comparison pools over plots, so it
  is robust to the realization; each arm picks `res` from its **own** measured
  `frdens`. Documented as a caveat.
- *multichm `z` is 2-D geometry.* Read apex height from the `Z` attribute, asserted
  by the contract; height RMSE is reported but the detection metrics (which drive
  the story) do not depend on `z` precision.
- *Writing into the shared main `work` dir.* Only the new
  `multichm_sweep_results.csv` / summary CSVs are written (gitignored, additive);
  nothing existing is clobbered.
