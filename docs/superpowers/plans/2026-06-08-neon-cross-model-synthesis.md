# NEON Cross-Model Density-Ladder Synthesis (#R10) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> superpowers:subagent-driven-development (recommended) or
> superpowers:executing-plans to implement this plan task-by-task. Steps use
> checkbox (`- [ ]`) syntax for tracking.

**Goal:** Produce the cross-model density-ladder result doc
(`results/model-benchmark-results.md`) and its generator
(`scripts/analyze_model_benchmark.R`), unifying every arm scored on the shared
frozen SOAP clips — AMS3D, lmfauto, multichm, ptrees, CHM-VWF, and a native-only
Li 2012 — into per-crown-class and per-height-band recall tables,
density-robustness curves, and head-to-head deltas vs CHM-VWF.

**Architecture:** The synthesis is pure analysis: it reads the per-arm long-form
CSVs each detector arm already writes under `work/neon/SOAP/`, unions them
(schema-harmonised), applies the existing `equal_set_guard` so every cross-arm
comparison rests on an identical (plot,rung) population, and pools per
(detector,rung) with the canonical `pool()` (sum counts, never average rates).
Two data gaps are closed first: AMS3D is regenerated over all 18 plots (it had
only `SOAP_001`), and a thin native-only Li 2012 arm is added so the spec's
"head-to-head vs CHM-VWF + Li 2012" is real and on identical inputs. `pool()` is
extended additively to also pool the three height bands, which it does not do
today, because the spec makes per-height-band recall a primary table.

**Tech Stack:** R (`Rscript`), lidR 4.3.2 (`li2012`, `segment_trees`),
data.table (`rbindlist(fill=TRUE)`), base-R `png` graphics, testthat 3.x. No new
dependencies.

---

## Data reality this plan is built around

- `work/neon/SOAP/lidrplugins_results.csv` — 360 rows: 18 plots x 5 rungs x
  {lmfauto, multichm, ptrees, chm_vwf}. Has extra `chm_res` + `tp_core` columns.
  **Done.**
- `work/neon/SOAP/ams3d_results.csv` — currently only `SOAP_001` (5 rows).
  **Task 0 regenerates it over all 18 plots** (running in background at plan
  time). Lacks `chm_res`/`tp_core`.
- Li 2012 — **does not exist** on the frozen clips. `scripts/detect_pc_sweep.R`
  runs it native-only on its *own* clips, so it is not comparable. **Task 2 adds
  a native-only frozen-clip Li 2012 arm.**
- All 90 frozen clips (18 plots x 5 rungs) already exist under
  `work/neon/SOAP/frozen/`, so Tasks 0 and 2 are cache-hits on clips plus
  segmentation compute — no re-clipping.

**Arm-set design (because Li 2012 is native-only):**

- **Density ladder** comparison uses `LADDER_ARMS =
  c("ams3d","lmfauto","multichm","ptrees","chm_vwf")` — all present at all 5
  rungs; `equal_set_guard` on this set.
- **Native point-segmenter head-to-head** uses `NATIVE_ARMS =
  c("chm_vwf","ams3d","ptrees","li2012")` at the native rung only; its own
  `equal_set_guard`. Putting Li 2012 in the ladder guard would drop every
  8/4/2/1 cell (Li 2012 absent there), so it is deliberately a separate table.

---

## Task 0: Regenerate the AMS3D arm over all 18 plots (data prep)

**Files:** none created/modified — runs the existing
`scripts/detect_ams3d_sweep.R`.

This is launched in the background before Task 1 and only needs to be *complete*
before Task 4 (the synthesis driver) consumes it. No code change.

- [ ] **Step 1: Launch (already running at plan time)**

```bash
export CLAUDE_JOB_DIR=$(pwd)/work
nohup nice -n 10 Rscript scripts/detect_ams3d_sweep.R SITE=SOAP PLOTS=ALL CORES=12 \
  > work/neon/SOAP/ams3d_fullrun.log 2>&1 &
```

- [ ] **Step 2: Verify completion before Task 4**

```bash
tail -3 work/neon/SOAP/ams3d_fullrun.log
Rscript -e 'd<-read.csv("work/neon/SOAP/ams3d_results.csv"); cat("plots:",length(unique(d$plot)),"rows:",nrow(d),"\n")'
```

Expected: `DONE` line in the log; `plots: 18 rows: 90` (18 plots x 5 rungs;
cells where the clip was unusable may be fewer — that is fine, the guard handles
it).

---

## Task 1: Extend `pool()` to pool height bands (the bridge)

**Files:**

- Modify: `scripts/model_bench_lib.R` (the `pool()` function, ~line 116-142)
- Test: `tests/testthat/test-pool-guard.R`

**Why a bridge change:** `pool()` is *our* shared library (not the frozen
`sweep_lib.R` scorer). The spec makes per-height-band recall a primary table,
and `pool()` currently pools only the 4 crown classes + understory. The
extension is **additive and guarded** (only acts when `n_h_short/mid/tall` are
present), so the two existing callers (`detect_lidrplugins_sweep.R`, and any
analysis) only gain columns. This keeps a single canonical "sum counts, never
average rates" pooler rather than forking a second one in the analysis script.

- [ ] **Step 1: Write the failing test**

Add to `tests/testthat/test-pool-guard.R`:

```r
test_that("pool pools height bands by summed counts when present", {
  hrow <- function(plot, n_ref, TP, tp_core, n_tall, rec_tall) {
    data.frame(site = "SOAP", plot = plot, rung = "native", detector = "x",
               n_ref = n_ref, TP = TP, n_det = TP, precision = tp_core / TP,
               recall = TP / n_ref,
               n_dominant = 0, rec_dominant = NA_real_,
               n_codominant = 0, rec_codominant = NA_real_,
               n_intermediate = 0, rec_intermediate = NA_real_,
               n_suppressed = 0, rec_suppressed = NA_real_,
               n_h_short = 0, rec_h_short = NA_real_,
               n_h_mid = 0, rec_h_mid = NA_real_,
               n_h_tall = n_tall, rec_h_tall = rec_tall)
  }
  df <- rbind(hrow("p1", 10, 5, 5, 4, 0.50),   # tall: 2 of 4
              hrow("p2", 10, 6, 6, 6, 0.50))    # tall: 3 of 6
  p <- pool(df)
  expect_equal(p$n_h_tall, 10)                  # summed n
  expect_equal(p$rec_h_tall, 0.5)               # (2+3)/(4+6), NOT mean(0.5,0.5)
  expect_true(is.na(p$rec_h_short))             # band absent -> NA, n=0
  expect_equal(p$n_h_short, 0)
})

test_that("pool omits height-band columns when n_h_* absent", {
  df <- rbind(mk_row("SOAP","p1","8","ams3d",10,5,6,4),
              mk_row("SOAP","p2","8","ams3d",90,45,50,40))
  p <- pool(df)
  expect_false("rec_h_tall" %in% names(p))      # additive: no bands in, none out
})
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-pool-guard.R")'`
Expected: FAIL — `pool()` returns no `n_h_tall`/`rec_h_tall` columns
(`p$n_h_tall` is NULL).

- [ ] **Step 3: Implement the additive height-band pooling**

In `scripts/model_bench_lib.R`, in `pool()`, insert the height-band block
immediately before the final `out` return (after the `out$n_understory <-
nref_u` line):

```r
  hbands <- c("short", "mid", "tall")
  if (all(paste0("n_h_", hbands) %in% names(df))) for (b in hbands) {
    nref <- sum(df[[paste0("n_h_", b)]], na.rm = TRUE)
    tp   <- sum(ifelse(df[[paste0("n_h_", b)]] > 0,
                       round(df[[paste0("rec_h_", b)]] * df[[paste0("n_h_", b)]]), 0),
                na.rm = TRUE)
    out[[paste0("rec_h_", b)]] <- if (nref) tp / nref else NA_real_
    out[[paste0("n_h_", b)]]   <- nref
  }
```

Also update the function's header comment (line ~110-115) to note it pools
height bands when `n_h_*` are present.

- [ ] **Step 4: Run the test to verify it passes**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-pool-guard.R")'`
Expected: PASS (all tests, including the two pre-existing ones).

- [ ] **Step 5: Commit**

```bash
git add scripts/model_bench_lib.R tests/testthat/test-pool-guard.R
git commit -m "feat(#R10): pool() pools height bands additively when present"
```

---

## Task 2: Native-only Li 2012 frozen-clip arm

**Files:**

- Create: `scripts/detect_li2012_native.R`
- Test: `tests/testthat/test-li2012-extractor.R`

**Why:** The spec wants a head-to-head vs Li 2012, the canonical lidR
point-cloud segmenter and the repo's only sub-canopy path.
`detect_pc_sweep.R`'s Li 2012 is native-only on different clips, so it is not
comparable. This arm runs Li 2012 on the **same native frozen clip** per plot,
collapses per-point `treeID` through the bridge's `reduce_instances`, and scores
with the existing harness — making the comparison faithful. Native-only by
design: Li 2012 is the dense-input understory test; scoring it at 1 pt/m² is
meaningless for a sub-canopy segmenter (and the report's interest in it is
exactly the native question).

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-li2012-extractor.R`:

```r
source(file.path("..", "..", "scripts", "detect_li2012_native.R"), local = TRUE)

test_that("det_li2012 returns the detection contract on a 2-tree clip", {
  las <- synth_las_normalized()
  det <- det_li2012(las, hmin = 2)
  expect_s3_class(det, "data.frame")
  expect_false(inherits(det, "sf"))
  expect_identical(names(det), c("x", "y", "z"))
  expect_true(all(vapply(det, is.numeric, logical(1))))
  expect_gte(nrow(det), 1L)                 # finds at least one apex
})

test_that("det_li2012 returns a 0-row frame (not NULL) on a no-canopy clip", {
  las <- synth_las_normalized()
  las@data$Z <- pmin(las@data$Z, 0.5)        # crush everything below hmin
  det <- det_li2012(las, hmin = 2)
  expect_s3_class(det, "data.frame")
  expect_identical(names(det), c("x", "y", "z"))
  expect_equal(nrow(det), 0L)                # ran-but-empty -> legit recall 0
})
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-li2012-extractor.R")'`
Expected: FAIL — `det_li2012` not found / file not sourceable.

- [ ] **Step 3: Write the arm**

Create `scripts/detect_li2012_native.R` (mirror `detect_lidrplugins_sweep.R`'s
locator + driver; native rung only):

```r
#!/usr/bin/env Rscript
# Native-only Li 2012 arm (#R10) of the NEON model benchmark. Runs lidR's
# li2012 point-cloud segmenter on the SAME native frozen clip per plot,
# collapses per-point treeID through the bridge's reduce_instances(), and
# scores against field stems with the existing harness. Native-only by design:
# li2012 is the dense-input sub-canopy test; decimated rungs are meaningless for
# a point segmenter, and CHM-VWF/ptrees/ams3d already carry the full ladder.
#
# Usage:  Rscript scripts/detect_li2012_native.R [SITE=SOAP] [PLOTS=ALL]
#             [CORES=6] [TOL=4]
# Output: $CLAUDE_JOB_DIR/neon/<SITE>/li2012_results.csv (one row per plot,
#         detector "li2012", rung "native").
suppressMessages({ library(lidR); library(data.table); library(parallel) })
options(lidR.progress = FALSE)
d <- Sys.getenv("CLAUDE_JOB_DIR", file.path(getwd(), "work"))
.find <- function(rel) Find(file.exists, c(file.path("scripts", rel),
                                           file.path("..", "..", "scripts", rel),
                                           file.path(getwd(), "scripts", rel)))
source(.find("sweep_lib.R"))
source(.find("model_bench_lib.R"))

# li2012 on a normalized clip -> per-point treeID -> reduce_instances apex set.
# Guard the no-canopy case (parity with ptrees): too few returns above hmin
# means no trees anyway, so return a 0-row frame without entering the segmenter.
det_li2012 <- function(las, hmin = 2, dt1 = 1.5, dt2 = 2, R = 2) {
  empty <- data.frame(x = numeric(), y = numeric(), z = numeric())
  if (sum(las$Z >= hmin) < 1) { assert_detection_contract(empty); return(empty) }
  seg <- tryCatch(lidR::segment_trees(las,
                    lidR::li2012(dt1 = dt1, dt2 = dt2, R = R, hmin = hmin)),
                  error = function(e) NULL)
  if (is.null(seg)) return(NULL)               # crash -> skip (guard drops cell)
  if (!"treeID" %in% names(seg@data)) return(NULL)
  det <- reduce_instances(seg@data, id_col = "treeID", x = "X", y = "Y", z = "Z")
  assert_detection_contract(det)
  det
}

args  <- strsplit(commandArgs(TRUE), "=")
A     <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
SITE  <- if (is.null(A$SITE))  "SOAP" else A$SITE
PLOTS <- if (is.null(A$PLOTS) || A$PLOTS == "ALL") NULL else strsplit(A$PLOTS, ",")[[1]]
CORES <- as.integer(if (is.null(A$CORES)) 6 else A$CORES)
TOL   <- as.numeric(if (is.null(A$TOL)) 4.0 else A$TOL)
MINTREES <- 6

run_main <- function() {
  nd  <- file.path(d, "neon", SITE)
  gt  <- read.csv(file.path(nd, "ground_truth_stems.csv"), stringsAsFactors = FALSE)
  pc  <- read.csv(file.path(nd, "plot_centroids.csv"),     stringsAsFactors = FALSE)
  gt  <- gt[gt$live & gt$is_tree & !is.na(gt$E), ]
  laz <- list.files(file.path(nd, "lidar"), pattern = "\\.laz$",
                    recursive = TRUE, full.names = TRUE)
  ctg <- readLAScatalog(laz, progress = FALSE)
  counts <- table(gt$plotID)
  keep   <- names(counts)[counts >= MINTREES]
  if (!is.null(PLOTS)) keep <- intersect(keep, PLOTS)
  keep   <- intersect(keep, pc$plotID)
  cat(sprintf("[%s] li2012 plots: %d (%s)\n", SITE, length(keep),
              paste(keep, collapse = ",")))

  run_plot <- function(pid) {
    ci <- pc[pc$plotID == pid, ][1, ]
    cx <- ci$easting; cy <- ci$northing
    ph <- plot_half(ci$plotType)
    stems <- gt[gt$plotID == pid & abs(gt$E - cx) <= ph & abs(gt$N - cy) <= ph, ]
    if (nrow(stems) < 1) return(NULL)
    prep <- tryCatch(frozen_clip(ctg, SITE, pid, NA, cx, cy, ph,
                                 out_root = file.path(nd, "frozen")),
                     error = function(e) NULL)
    if (is.null(prep)) return(NULL)
    las <- tryCatch(readLAS(prep$normalized), error = function(e) NULL)
    if (is.null(las) || is.empty(las)) return(NULL)
    det <- det_li2012(las, hmin = 2)
    if (is.null(det)) return(NULL)
    sc <- tryCatch(score_plot(stems, det, tol_xy = TOL, core_cx = cx,
                              core_cy = cy, core_half = ph),
                   error = function(e) NULL)
    if (is.null(sc)) return(NULL)
    cbind(data.frame(site = SITE, plot = pid, plotType = ci$plotType,
                     detector = "li2012", rung = "native",
                     pdens = round(prep$pdens, 2), frdens = round(prep$frdens, 2),
                     n_apex = nrow(det)), sc)
  }

  res_list <- mclapply(keep, function(p)
                tryCatch(run_plot(p), error = function(e) {
                  message("plot ", p, " failed: ", conditionMessage(e)); NULL }),
                mc.cores = CORES, mc.preschedule = FALSE)
  results <- do.call(rbind, Filter(Negate(is.null), res_list))
  if (is.null(results) || !nrow(results)) { cat("no li2012 results\n"); return(invisible()) }
  results$tp_core <- round(results$precision * results$n_det)
  write.csv(results, file.path(nd, "li2012_results.csv"), row.names = FALSE)
  cat(sprintf("[%s] li2012 DONE: %d rows -> %s\n", SITE, nrow(results),
              file.path(nd, "li2012_results.csv")))
}

if (sys.nframe() == 0L) run_main()
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-li2012-extractor.R")'`
Expected: PASS (both tests).

- [ ] **Step 5: Run the arm over all 18 plots (data prep for Task 4)**

```bash
export CLAUDE_JOB_DIR=$(pwd)/work
Rscript scripts/detect_li2012_native.R SITE=SOAP PLOTS=ALL CORES=12
Rscript -e 'd<-read.csv("work/neon/SOAP/li2012_results.csv"); cat("plots:",length(unique(d$plot)),"rows:",nrow(d),"\n")'
```

Expected: `li2012 DONE` line; `plots: 18 rows: 18` (one native row per plot;
treeless plots may be fewer).

- [ ] **Step 6: Commit**

```bash
git add scripts/detect_li2012_native.R tests/testthat/test-li2012-extractor.R
git commit -m "feat(#R10): native-only Li 2012 frozen-clip arm for the head-to-head"
```

---

## Task 3: Synthesis core functions (`analyze_model_benchmark.R`)

**Files:**

- Create: `scripts/analyze_model_benchmark.R` (functions only this task; driver
  in Task 4)
- Test: `tests/testthat/test-model-benchmark-analysis.R`

**Responsibility:** Pure, testable functions — union/harmonise the arm CSVs,
pool per (detector,rung) over the guarded common cells, and compute deltas vs
the CHM-VWF baseline. No file I/O or plotting in these functions (that is Task
4's `run_main`).

The union must reconcile schema (`chm_res`/`tp_core` present only on the
lidRplugins arm) and **recompute `tp_core` for every row** from `round(precision

- n_det)` so pooling is identical across arms regardless of which CSV carried it
(an incoming `tp_core` with NAs for the ams3d/li2012 rows would silently
undercount their precision).

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-model-benchmark-analysis.R`:

```r
source(file.path("..", "..", "scripts", "analyze_model_benchmark.R"), local = TRUE)

# Two arms with DIFFERENT schemas (one has chm_res+tp_core, one does not),
# mimicking lidrplugins_results.csv vs ams3d_results.csv.
mk_arm <- function(det, rung, plot, n_ref, TP, n_det, with_extra) {
  base <- data.frame(site = "SOAP", plot = plot, plotType = "distributed",
                     detector = det, rung = rung, pdens = 14, frdens = 7,
                     n_apex = n_det, n_ref = n_ref, n_det = n_det, TP = TP,
                     recall = TP / n_ref, precision = TP / n_det, F1 = NA_real_,
                     height_rmse = 1,
                     rec_dominant = TP / n_ref, n_dominant = n_ref,
                     rec_codominant = NA_real_, n_codominant = 0,
                     rec_intermediate = NA_real_, n_intermediate = 0,
                     rec_suppressed = NA_real_, n_suppressed = 0,
                     rec_h_short = NA_real_, n_h_short = 0,
                     rec_h_mid = NA_real_, n_h_mid = 0,
                     rec_h_tall = TP / n_ref, n_h_tall = n_ref)
  if (with_extra) { base$chm_res <- 0.5; base$tp_core <- TP }
  base
}

test_that("harmonize_union fills missing cols and recomputes tp_core for all rows", {
  a <- mk_arm("chm_vwf", "8", "p1", 10, 6, 8, TRUE)
  b <- mk_arm("ams3d",   "8", "p1", 10, 5, 20, FALSE)   # no chm_res/tp_core
  u <- harmonize_union(list(a, b))
  expect_true(all(c("chm_res", "tp_core") %in% names(u)))
  expect_equal(nrow(u), 2L)
  # tp_core recomputed = round(precision * n_det) for BOTH rows
  expect_equal(u$tp_core[u$detector == "ams3d"], round((5 / 20) * 20))   # 5
  expect_true(is.na(u$chm_res[u$detector == "ams3d"]))                   # filled NA
})

test_that("pool_arms guards to common cells and pools per (detector,rung)", {
  u <- harmonize_union(list(
    mk_arm("chm_vwf", "8", "p1", 10, 6, 8, TRUE),
    mk_arm("ams3d",   "8", "p1", 10, 5, 20, FALSE),
    mk_arm("chm_vwf", "8", "p2", 10, 4, 9, TRUE)))   # p2 has no ams3d -> dropped
  pooled <- pool_arms(u, arms = c("chm_vwf", "ams3d"), rungs = "8")
  expect_setequal(unique(pooled$detector), c("chm_vwf", "ams3d"))
  # only p1/8 survives the guard, so chm_vwf recall = 6/10 (NOT pooling p2)
  expect_equal(pooled$recall[pooled$detector == "chm_vwf"], 0.6)
  expect_equal(pooled$n_plots[pooled$detector == "chm_vwf"], 1L)
  expect_true("rec_h_tall" %in% names(pooled))         # height bands flow through
})

test_that("deltas_vs_baseline subtracts the baseline arm per rung", {
  pooled <- data.frame(detector = c("chm_vwf", "ams3d"), rung = c("8", "8"),
                       recall = c(0.6, 0.5), F1 = c(0.5, 0.4),
                       rec_understory = c(0.2, 0.4))
  dl <- deltas_vs_baseline(pooled, baseline = "chm_vwf")
  expect_equal(dl$d_recall[dl$detector == "ams3d"], -0.1)
  expect_equal(dl$d_understory[dl$detector == "ams3d"], 0.2)
  expect_false("chm_vwf" %in% dl$detector)            # baseline itself excluded
})
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-model-benchmark-analysis.R")'`
Expected: FAIL — functions not found / file not sourceable.

- [ ] **Step 3: Write the functions**

Create `scripts/analyze_model_benchmark.R` (functions + a robust locator; the
driver is added in Task 4):

```r
#!/usr/bin/env Rscript
# Cross-model density-ladder synthesis (#R10). Unifies every arm scored on the
# shared frozen SOAP clips (AMS3D, lmfauto, multichm, ptrees, CHM-VWF, and the
# native-only Li 2012) into per-class + per-height-band recall tables, density-
# robustness curves, and head-to-head deltas vs CHM-VWF. Pools with the canonical
# pool() (sum counts, never average rates) over equal_set_guard common cells.
#
# Usage:  Rscript scripts/analyze_model_benchmark.R [SITE=SOAP]
# Reads:  $CLAUDE_JOB_DIR/neon/<SITE>/{ams3d,lidrplugins,li2012}_results.csv
# Writes: summary CSVs + figs/ + a generated markdown table fragment under the
#         same dir; the narrative results/model-benchmark-results.md is authored.
suppressMessages({ library(data.table) })
.find <- function(rel) Find(file.exists, c(file.path("scripts", rel),
                                           file.path("..", "..", "scripts", rel),
                                           file.path(getwd(), "scripts", rel)))
source(.find("model_bench_lib.R"))

RUNG_LEVELS <- c("native", "8", "4", "2", "1")

# Union arm data.frames with differing schemas; recompute tp_core uniformly so
# pooling is identical across arms regardless of which CSV carried the column.
harmonize_union <- function(dfs) {
  u <- as.data.frame(data.table::rbindlist(dfs, fill = TRUE, use.names = TRUE))
  u$tp_core <- round(u$precision * u$n_det)
  u$rung <- factor(u$rung, levels = RUNG_LEVELS)
  u
}

# Restrict to `arms` and `rungs`, apply equal_set_guard per (site,plot,rung), then
# pool() each (detector,rung). Carries frdens (mean over the cell's plots) as the
# density-curve x. Returns a long table; attr "dropped" = guard's dropped cells.
pool_arms <- function(u, arms, rungs = RUNG_LEVELS) {
  sub <- u[u$detector %in% arms & as.character(u$rung) %in% rungs, ]
  sub$rung <- as.character(sub$rung)
  g <- equal_set_guard(sub, arms = arms)
  out <- list()
  for (a in arms) for (rl in rungs) {
    s <- g[g$detector == a & g$rung == rl, ]
    if (!nrow(s)) next
    out[[length(out) + 1]] <- cbind(
      detector = a, rung = rl, frdens = round(mean(s$frdens), 2), pool(s))
  }
  res <- do.call(rbind, out)
  if (!is.null(res)) res$rung <- factor(res$rung, levels = rungs)
  attr(res, "dropped") <- attr(g, "dropped")
  res
}

# Per rung, each non-baseline arm's recall/F1/understory minus the baseline arm's.
deltas_vs_baseline <- function(pooled, baseline = "chm_vwf") {
  base <- pooled[pooled$detector == baseline, ]
  others <- pooled[pooled$detector != baseline, ]
  m <- merge(others, base[, c("rung", "recall", "F1", "rec_understory")],
             by = "rung", suffixes = c("", "_base"))
  data.frame(detector = m$detector, rung = m$rung,
             d_recall = m$recall - m$recall_base,
             d_F1 = m$F1 - m$F1_base,
             d_understory = m$rec_understory - m$rec_understory_base)
}
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-model-benchmark-analysis.R")'`
Expected: PASS (all three tests).

- [ ] **Step 5: Commit**

```bash
git add scripts/analyze_model_benchmark.R tests/testthat/test-model-benchmark-analysis.R
git commit -m "feat(#R10): cross-model synthesis core (union, pool_arms, deltas)"
```

---

## Task 4: Synthesis driver + outputs (CSVs, figures, generated tables)

**Files:**

- Modify: `scripts/analyze_model_benchmark.R` (append `run_main()` + `if
  (sys.nframe()==0L)`)

**Responsibility:** Wire the pure functions to the real data: load whichever arm
CSVs exist, build the ladder table (LADDER_ARMS) and the native point-segmenter
table (NATIVE_ARMS), compute deltas, write summary CSVs and base-R figures, and
emit a machine-generated markdown table fragment that Task 5's narrative
includes.

- [ ] **Step 1: Append the driver to `scripts/analyze_model_benchmark.R`**

```r
LADDER_ARMS <- c("ams3d", "lmfauto", "multichm", "ptrees", "chm_vwf")
NATIVE_ARMS <- c("chm_vwf", "ptrees", "ams3d", "li2012")

# Render a long pooled table to a GitHub-markdown block (selected columns).
.md_table <- function(df, cols, digits = 2) {
  hdr <- paste0("| ", paste(cols, collapse = " | "), " |")
  sep <- paste0("|", paste(rep(" --- ", length(cols)), collapse = "|"), "|")
  body <- apply(df[, cols, drop = FALSE], 1, function(r)
    paste0("| ", paste(vapply(r, function(v) {
      if (is.numeric(v) && !is.na(suppressWarnings(as.numeric(v))))
        formatC(as.numeric(v), format = "f", digits = digits) else as.character(v)
    }, character(1)), collapse = " | "), " |"))
  paste(c(hdr, sep, body), collapse = "\n")
}

run_main <- function() {
  args <- strsplit(commandArgs(TRUE), "=")
  A    <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
  SITE <- if (is.null(A$SITE)) "SOAP" else A$SITE
  d    <- Sys.getenv("CLAUDE_JOB_DIR", file.path(getwd(), "work"))
  nd   <- file.path(d, "neon", SITE)

  arm_files <- c(ams3d = "ams3d_results.csv", lidrplugins = "lidrplugins_results.csv",
                 li2012 = "li2012_results.csv")
  dfs <- lapply(file.path(nd, arm_files), function(f)
    if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE) else NULL)
  dfs <- Filter(Negate(is.null), dfs)
  if (!length(dfs)) stop("no arm result CSVs found under ", nd)
  u <- harmonize_union(dfs)
  cat("arms present:", paste(sort(unique(u$detector)), collapse = ", "), "\n")

  ladder <- pool_arms(u, arms = intersect(LADDER_ARMS, unique(u$detector)))
  native <- pool_arms(u, arms = intersect(NATIVE_ARMS, unique(u$detector)),
                      rungs = "native")
  dl <- deltas_vs_baseline(ladder, baseline = "chm_vwf")

  for (nm in c("ladder", "native", "dl"))
    write.csv(get(nm), file.path(nd, paste0("model_bench_", nm, ".csv")),
              row.names = FALSE)
  if (length(attr(ladder, "dropped")))
    cat(sprintf("ladder guard dropped %d (plot,rung) cells\n",
                length(attr(ladder, "dropped"))))

  ## figures: recall and understory recall vs first-return density, per arm
  fig <- file.path(nd, "figs"); dir.create(fig, showWarnings = FALSE)
  arms <- intersect(LADDER_ARMS, unique(ladder$detector))
  pal  <- setNames(c("#1b9e77","#d95f02","#7570b3","#e7298a","#666666")[seq_along(arms)], arms)
  draw <- function(file, ycol, ylab, main) {
    png(file.path(fig, file), 1000, 750, res = 130); on.exit(dev.off())
    par(mar = c(4.2, 4.2, 2.5, 1)); first <- TRUE
    for (a in arms) {
      s <- ladder[ladder$detector == a, ]; s <- s[order(s$frdens), ]
      if (!nrow(s)) next
      if (first) { plot(s$frdens, s[[ycol]], type = "o", pch = 19, ylim = c(0, 1),
                        log = "x", xlab = "first-return density (pulses/m^2)",
                        ylab = ylab, main = main, col = pal[a], lwd = 2); first <- FALSE
      } else lines(s$frdens, s[[ycol]], type = "o", pch = 19, col = pal[a], lwd = 2)
    }
    legend("topleft", arms, col = pal[arms], lwd = 2, pch = 19, bty = "n", cex = 0.85)
  }
  draw("model_recall_vs_density.png", "recall", "recall (overall)",
       "SOAP cross-model recall vs density")
  draw("model_understory_vs_density.png", "rec_understory", "understory recall",
       "SOAP understory recall vs density")

  ## generated markdown fragment (tables only; narrative is authored separately)
  cols_l <- c("detector","rung","frdens","n_plots","n_ref","recall","precision",
              "F1","rec_dominant","rec_codominant","rec_understory")
  cols_h <- c("detector","rung","rec_h_tall","n_h_tall","rec_h_mid","n_h_mid",
              "rec_h_short","n_h_short")
  cols_n <- c("detector","n_plots","n_ref","recall","precision","F1",
              "rec_understory","n_understory")
  cols_d <- c("detector","rung","d_recall","d_F1","d_understory")
  frag <- c("<!-- generated by scripts/analyze_model_benchmark.R; do not edit by hand -->",
            "", "#### Density ladder, per crown class (pooled, equal-set)", "",
            .md_table(ladder[order(ladder$detector, ladder$rung), ], cols_l), "",
            "#### Density ladder, per height band", "",
            .md_table(ladder[order(ladder$detector, ladder$rung), ], cols_h), "",
            "#### Native point-segmenter head-to-head", "",
            .md_table(native[order(native$detector), ], cols_n), "",
            "#### Head-to-head deltas vs CHM-VWF (recall/F1/understory)", "",
            .md_table(dl[order(dl$detector, dl$rung), ], cols_d))
  writeLines(frag, file.path(nd, "model_bench_tables.md"))
  cat(sprintf("synthesis -> %s {ladder,native,dl}.csv, figs/, model_bench_tables.md\n", nd))
}

if (sys.nframe() == 0L) run_main()
```

- [ ] **Step 2: Run the full test suite (sourcing under `sys.nframe` is safe)**

Run: `Rscript tests/run_tests.R`
Expected: all test files PASS (the `if (sys.nframe()==0L)` guard means sourcing
the script under testthat does not invoke `run_main`).

- [ ] **Step 3: Run the synthesis on the real data (Tasks 0 & 2 must be complete)**

```bash
export CLAUDE_JOB_DIR=$(pwd)/work
Rscript scripts/analyze_model_benchmark.R SITE=SOAP
ls work/neon/SOAP/model_bench_*.csv work/neon/SOAP/model_bench_tables.md work/neon/SOAP/figs/model_*vs_density.png
cat work/neon/SOAP/model_bench_tables.md
```

Expected: prints `arms present: ams3d, chm_vwf, li2012, lmfauto, multichm,
ptrees`, a non-zero "ladder guard dropped" count is fine, and writes the three
CSVs, two PNGs, and the markdown fragment. Sanity-check the printed tables
(recall in [0,1], chm_vwf present at every rung).

- [ ] **Step 4: Commit**

```bash
git add scripts/analyze_model_benchmark.R
git commit -m "feat(#R10): synthesis driver — pooled tables, density figures, md fragment"
```

---

## Task 5: Author the result doc + README row

**Files:**

- Create: `results/model-benchmark-results.md`
- Modify: `README.md` (script table: add the new scripts; results section: link
  the doc)

**Responsibility:** The deliverable. Weave a hand-written narrative around the
generated tables (paste the verified blocks from `model_bench_tables.md`), with
the structure the spec requires. Keep prose within rumdl's 80-char limit
(tables/code exempt).

- [ ] **Step 1: Write `results/model-benchmark-results.md`**

Required sections (per spec #R10):

1. **What this is** — cross-model density-ladder synthesis on shared frozen SOAP
   clips; the arm set and why Li 2012 is a separate native-only table; equal-set
   guard note (every cross-arm comparison on an identical (plot,rung)
   population; report the dropped-cell count).
2. **Per-crown-class recall, density ladder** — the `ladder` table; narrate
   AMS3D-vs-CHM density robustness and where each arm wins per class.
3. **Per-height-band recall** — the height-band table; tall/mid/short trends.
4. **Density-robustness curves** — reference the two PNGs; describe the slopes.
5. **Native point-segmenter head-to-head** — the `native` table; CHM-VWF vs
   ptrees/AMS3D/Li 2012 for understory at native density.
6. **Head-to-head deltas vs CHM-VWF** — the `dl` table; one-paragraph reading.
7. **Zero-shot ledger appendix** — these are classical/zero-shot detectors (no
   NEON-specific tuning); cite `model-benchmark-plan.md` (#A0). Note the GPU
   deep models (SegmentAnyTree/TreeisoNet/ForestFormer3D) are deferred
   (#M6-#M8).
8. **Caveats** — (a) decimation != native low density (carry the existing repo
   caveat; cite the native-QL2 cross-check #4); (b) discrete-return ALS vs the
   ULS/UAS the instance models were trained on; (c) field-stem ground truth
   reduces every model to detections — apex recall, not point-IoU; (d) Li 2012
   full-ladder and other sites (SJER/TEAK) are deferred (#E11).

Pull the actual numbers from `work/neon/SOAP/model_bench_tables.md` (do not
invent them). Lint: `rumdl check results/model-benchmark-results.md`.

- [ ] **Step 2: Add the scripts to the README script table and link the doc**

Add rows for `detect_li2012_native.R`, `analyze_model_benchmark.R` (and confirm
`detect_ams3d_sweep.R`/`detect_lidrplugins_sweep.R` rows exist), and link
`results/model-benchmark-results.md` from the results section. Match the
existing table's column format.

- [ ] **Step 3: Lint and verify**

Run: `rumdl check results/model-benchmark-results.md README.md`
Expected: no errors (fix any 80-char prose overflows).

- [ ] **Step 4: Commit**

```bash
git add results/model-benchmark-results.md README.md
git commit -m "docs(#R10): cross-model density-ladder result doc + README rows"
```

---

## Self-review notes (gaps to watch during execution)

- **`pool()` change is the one bridge edit** — keep it additive; the two
  existing callers must still pass. Re-run `test-pool-guard.R` and the
  lidRplugins smoke.
- **`tp_core` recompute on union** is load-bearing: without it the ams3d/li2012
  rows (no incoming `tp_core`) would pool precision as if `tp_core = NA -> 0`.
- **Arm-set split** (ladder excludes li2012; native table includes it) is what
  keeps the guard from nuking the 8/4/2/1 rungs. Do not put li2012 in
  `LADDER_ARMS`.
- **Numbers in the doc come from the generated fragment**, never hand-typed.
- If Task 0 / Task 2 leave a plot with no AMS3D or no Li 2012 cell, the guard
  drops it — report the count, do not silently proceed.
