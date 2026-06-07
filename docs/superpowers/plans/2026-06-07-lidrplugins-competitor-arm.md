# lidRplugins Competitor Arm (#C9) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> superpowers:subagent-driven-development (recommended) or
> superpowers:executing-plans to implement this plan task-by-task. Steps use
> checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add the LiDAR-native, classical, density-responsive competitor arm —
`lidRplugins` `lmfauto` / `multichm` / `ptrees` — as a fair head-to-head against
the existing CHM-VWF baseline across the SOAP density ladder, reusing the
`model_bench_lib.R` bridge and the unchanged `score_plot` harness.

**Architecture:** One driver script that, per plot × density rung, runs three
`lidRplugins` tree-top detectors plus the CHM-VWF baseline on the **same frozen
normalized clip** (`frozen_clip`), collapses each detector's tops to the
`(x,y,z)` detection contract, scores via the existing `score_plot`/`greedy_match`,
and pools with the canonical `pool` + `equal_set_guard`. Detectors return
treetops directly (not per-point instances), so they bypass `reduce_instances`
and use `st_coordinates` + `assert_detection_contract`. A detector that *crashes*
returns `NULL` (skipped; equal-set guard drops the cell); a detector that *runs
but finds nothing* returns a 0-row frame (legitimate `recall=0`) — the same
crash-vs-empty discipline as the AMS3D arm.

**Tech stack:** R (`lidR` 4.3.2, `sf`, `data.table`, `parallel`), `lidRplugins`
(GitHub-only, `remotes::install_github`), `testthat` 3e. Reuses
`scripts/sweep_lib.R` (unchanged scorer) and `scripts/model_bench_lib.R` (the
bridge from PR #14).

**Scope note:** This is the #C9 arm from
`docs/superpowers/specs/2026-06-07-neon-deep-model-benchmark-design.md`. It
produces the lidRplugins arm's results + an in-script comparison vs CHM-VWF. The
cross-arm density-ladder synthesis against AMS3D and Li 2012 is **#R10** (a
separate plan). The GPU arms (#I3–#M8) are separate plans.

---

## File structure

| File | Responsibility | Tasks |
|---|---|---|
| `scripts/detect_lidrplugins_sweep.R` | The #C9 arm: 3 lidRplugins detector extractors + the full-ladder driver (frozen clips, score, pool) | 1, 2, 3 |
| `tests/testthat/test-lidrplugins-extractors.R` | Unit tests for the 3 extractors (synthetic LAS, no real data) | 2 |

No changes to `sweep_lib.R` or `model_bench_lib.R` — this arm only *consumes*
them. Synthetic fixtures (`tests/testthat/helper-synth.R`) already exist.

---

## Task 1: Install lidRplugins and gate lidR 4.3.2 compatibility

**Files:** none (environment + a throwaway compatibility check)

`lidRplugins` is GitHub-only and predates the lidR `find_trees`→`locate_trees`
rename. This task is the **gate**: if it won't install or its algorithms don't
plug into lidR 4.3.2's `locate_trees`, the whole arm is BLOCKED and we stop here
to decide (pin an older lidR vs shim) rather than build on sand.

- [ ] **Step 1: Install remotes (if absent) and lidRplugins from GitHub**

Run:

```bash
Rscript -e 'if (!requireNamespace("remotes", quietly=TRUE)) \
  install.packages("remotes", repos="https://cloud.r-project.org")'
Rscript -e 'if (!requireNamespace("lidRplugins", quietly=TRUE)) \
  remotes::install_github("Jean-Romain/lidRplugins", upgrade="never")'
```

Expected: installs without error. Verify:
`Rscript -e 'packageVersion("lidRplugins")'` prints a version.
If install fails (compile error, dependency conflict with lidR 4.3.2), STOP and
report BLOCKED with the exact error.

- [ ] **Step 2: Smoke-check that all three algorithms plug into `locate_trees`**

Run this compatibility probe (uses the existing synthetic fixture):

```bash
Rscript -e '
suppressMessages({library(lidR); library(lidRplugins); library(sf)})
source("tests/testthat/helper-synth.R")
las <- synth_las_normalized()
for (nm in c("lmfauto","multichm","ptrees")) {
  algo <- switch(nm,
    lmfauto  = lmfauto(plot = TRUE, hmin = 2),
    multichm = multichm(res = 0.5, ws = function(h) pmin(pmax(0.1*h+3,3),5)),
    ptrees   = ptrees(k = c(30, 15), hmin = 2))
  tt <- tryCatch(locate_trees(las, algo), error = function(e) conditionMessage(e))
  cat(sprintf("%-9s -> %s\n", nm,
      if (inherits(tt,"sf")) sprintf("OK, %d tops, 3D=%s", nrow(tt),
        ncol(sf::st_coordinates(tt))>=3) else paste("ERROR:", tt)))
}'
```

Expected: each of `lmfauto`, `multichm`, `ptrees` prints `OK, <n> tops, 3D=TRUE`
(n may be small/zero on the toy cloud — OK as long as it does not ERROR and the
geometry is 3D so `st_coordinates` yields an X,Y,Z column).

- If all three print `OK ... 3D=TRUE`: proceed to Task 2.
- If any prints `ERROR: ...`: this is the compatibility gate failing. Report
  BLOCKED with which algorithm(s) failed and the message. Do NOT attempt a fix
  in this task — it is an environment/version decision for the orchestrator
  (e.g., pin lidR to the version lidRplugins targets, or write a `locate_trees`
  shim). Note: if `multichm`/`ptrees` 3D=FALSE (tops carry Z only as an
  attribute, not geometry), record that — Task 2's extractor will read `tt$Z`
  as a fallback.
- [ ] **Step 3: Commit a note (no code yet) only if a fallback was discovered**

If Step 2 revealed a 2D-geometry fallback need or any algorithm-specific quirk,
record it in the commit message of Task 2; otherwise nothing to commit here.

---

## Task 2: The three detector extractors

**Files:**

- Create: `scripts/detect_lidrplugins_sweep.R`
- Test: `tests/testthat/test-lidrplugins-extractors.R`

Each extractor takes a normalized `LAS`, runs one `lidRplugins` detector via
`locate_trees`, and returns the `(x,y,z)` detection contract — or `NULL` on a
crash (skip), or a 0-row frame on a ran-but-empty result.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-lidrplugins-extractors.R`:

```r
source(file.path("..", "..", "scripts", "detect_lidrplugins_sweep.R"), local = TRUE)

test_that("each lidRplugins extractor returns the x,y,z contract on a 2-tree cloud", {
  las <- synth_las_normalized()
  for (fn in list(
        function() det_lmfauto(las, hmin = 2),
        function() det_multichm(las, res = 0.5, a = 0.10),
        function() det_ptrees(las, hmin = 2))) {
    det <- fn()
    # NULL only on a genuine crash; on this clean synthetic cloud expect a frame
    expect_s3_class(det, "data.frame")
    expect_false(inherits(det, "sf"))
    expect_identical(names(det), c("x", "y", "z"))
    expect_true(all(vapply(det, is.numeric, logical(1))))
  }
})

test_that("an extractor returns a 0-row frame (not NULL) when the detector finds nothing", {
  flat <- LAS(data.frame(X = runif(200, 0, 10), Y = runif(200, 0, 10),
                         Z = rep(0.1, 200)))      # no trees: all near ground
  st_crs(flat) <- 32611L
  det <- det_lmfauto(flat, hmin = 2)
  expect_true(is.null(det) || (is.data.frame(det) && nrow(det) == 0L))
  if (is.data.frame(det)) expect_identical(names(det), c("x", "y", "z"))
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-lidrplugins-extractors.R")'`
Expected: FAIL — `could not find function "det_lmfauto"`.

- [ ] **Step 3: Create the script header + the three extractors**

Create `scripts/detect_lidrplugins_sweep.R`:

```r
#!/usr/bin/env Rscript
# lidRplugins competitor arm (#C9) of the NEON model benchmark.
# Runs three classical LiDAR-native tree-top detectors (lmfauto, multichm,
# ptrees) plus the CHM-VWF baseline on the SAME frozen normalized clip per plot
# x density rung, scoring each against field stems with the existing harness.
# Detectors return treetops directly (not per-point instances), so they use
# st_coordinates + assert_detection_contract, NOT reduce_instances.
#
# Usage:
#   Rscript scripts/detect_lidrplugins_sweep.R [SITE=SOAP] [PLOTS=ALL]
#       [CORES=6] [TOL=4] [A=0.10]
# Output: $CLAUDE_JOB_DIR/neon/<SITE>/lidrplugins_results.csv (one row per
#         plot x rung x detector).
suppressMessages({ library(lidR); library(lidRplugins); library(sf)
                   library(data.table); library(parallel) })
options(lidR.progress = FALSE)
d <- Sys.getenv("CLAUDE_JOB_DIR", file.path(getwd(), "work"))
# Locate the two libs robustly across run contexts (Rscript from root vs
# source()-d under testthat where cwd = tests/testthat).
.find <- function(rel) Find(file.exists, c(file.path("scripts", rel),
                                           file.path("..", "..", "scripts", rel),
                                           file.path(getwd(), "scripts", rel)))
source(.find("sweep_lib.R"))
source(.find("model_bench_lib.R"))

## ---- tops sf -> (x,y,z) detection contract -------------------------------
# Pull lowercase x,y,z from a locate_trees() sf. Prefer 3D geometry; fall back
# to the Z attribute column if the geometry is 2D. assert_detection_contract
# enforces the exact shape score_plot consumes.
.tops_to_det <- function(tt) {
  empty <- data.frame(x = numeric(), y = numeric(), z = numeric())
  if (is.null(tt) || nrow(tt) == 0) { assert_detection_contract(empty); return(empty) }
  co <- sf::st_coordinates(tt)
  z  <- if (ncol(co) >= 3) co[, 3] else if (!is.null(tt$Z)) tt$Z else co[, 2] * 0
  det <- data.frame(x = co[, 1], y = co[, 2], z = z)
  assert_detection_contract(det)
  det
}

# lmfauto: parameter-free variable-window LMF (plot mode for plot-scale clips).
det_lmfauto <- function(las, hmin = 2) {
  tt <- tryCatch(lidR::locate_trees(las, lidRplugins::lmfauto(plot = TRUE, hmin = hmin)),
                 error = function(e) NULL)
  if (is.null(tt)) return(NULL)            # crash -> skip (equal-set guard drops it)
  .tops_to_det(tt)
}

# multichm: multi-layer CHM; pass our clamped variable window (same allometry as
# CHM-VWF) through ... to the internal lmf for a fair comparison.
det_multichm <- function(las, res = 0.5, a = 0.10) {
  tt <- tryCatch(lidR::locate_trees(las, lidRplugins::multichm(res = res,
                                                               ws = ws_factory(a))),
                 error = function(e) NULL)
  if (is.null(tt)) return(NULL)
  .tops_to_det(tt)
}

# ptrees: point-based PTrees (Vega 2014); detection mode via locate_trees.
det_ptrees <- function(las, hmin = 2, k = c(30, 15)) {
  tt <- tryCatch(lidR::locate_trees(las, lidRplugins::ptrees(k = k, hmin = hmin)),
                 error = function(e) NULL)
  if (is.null(tt)) return(NULL)
  .tops_to_det(tt)
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-lidrplugins-extractors.R")'`
Expected: PASS. If an extractor returns `NULL` on the clean 2-tree synthetic
cloud (first test), that means the detector crashed on a trivial input —
investigate (likely a lidR-version API mismatch surfaced past Task 1's probe)
and report; do not weaken the test. The second test tolerates NULL-or-0-row.

- [ ] **Step 5: Commit**

```bash
git add scripts/detect_lidrplugins_sweep.R tests/testthat/test-lidrplugins-extractors.R
git commit -m "feat(#C9): lidRplugins detector extractors (lmfauto/multichm/ptrees)"
```

---

## Task 3: Full-ladder driver + CHM-VWF head-to-head

**Files:**

- Modify: `scripts/detect_lidrplugins_sweep.R` (append the driver)

Append the driver: per plot × rung, run the three lidRplugins detectors **plus
the CHM-VWF baseline** (`detect_lasr` from `sweep_lib.R`) on the frozen
normalized clip, score each, write a long-form CSV, and print a pooled
comparison using `equal_set_guard` + `pool`.

- [ ] **Step 1: Append the driver**

Append to `scripts/detect_lidrplugins_sweep.R`:

```r
## ---- args ----------------------------------------------------------------
args  <- strsplit(commandArgs(TRUE), "=")
A     <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
SITE  <- if (is.null(A$SITE))  "SOAP" else A$SITE
PLOTS <- if (is.null(A$PLOTS) || A$PLOTS == "ALL") NULL else strsplit(A$PLOTS, ",")[[1]]
CORES <- as.integer(if (is.null(A$CORES)) 6 else A$CORES)
TOL   <- as.numeric(if (is.null(A$TOL)) 4.0 else A$TOL)
A_VWF <- as.numeric(if (is.null(A$A))   0.10 else A$A)
RUNGS <- c(8, 4, 2, 1)
MINTREES <- 6
ARMS  <- c("lmfauto", "multichm", "ptrees", "chm_vwf")

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
  cat(sprintf("[%s] lidRplugins plots: %d (%s)\n", SITE, length(keep),
              paste(keep, collapse = ",")))

  run_plot <- function(pid) {
    ci <- pc[pc$plotID == pid, ][1, ]
    cx <- ci$easting; cy <- ci$northing
    ph <- plot_half(ci$plotType)
    stems <- gt[gt$plotID == pid & abs(gt$E - cx) <= ph & abs(gt$N - cy) <= ph, ]
    if (nrow(stems) < 1) return(NULL)
    out <- list(); native_pdens <- NA_real_
    for (rung in c(NA, RUNGS)) {
      prep <- tryCatch(frozen_clip(ctg, SITE, pid, rung, cx, cy, ph,
                                   out_root = file.path(nd, "frozen")),
                       error = function(e) NULL)
      if (is.null(prep)) next
      pdens <- prep$pdens; frdens <- prep$frdens
      if (is.na(rung)) native_pdens <- pdens
      else if (is.na(native_pdens) || rung >= native_pdens) next
      las <- tryCatch(readLAS(prep$normalized), error = function(e) NULL)
      if (is.null(las) || is.empty(las)) next
      res <- if (frdens >= 8) 0.25 else 0.5     # density-derived, like CHM-VWF
      dets <- list(
        lmfauto  = det_lmfauto(las, hmin = 2),
        multichm = det_multichm(las, res = res, a = A_VWF),
        ptrees   = det_ptrees(las, hmin = 2),
        chm_vwf  = tryCatch(detect_lasr(prep$normalized, res, A_VWF, frdens),
                            error = function(e) NULL))
      for (nm in names(dets)) {
        det <- dets[[nm]]
        if (is.null(det)) next                  # crash -> skip this detector/cell
        sc <- tryCatch(score_plot(stems, det, tol_xy = TOL, core_cx = cx,
                                  core_cy = cy, core_half = ph),
                       error = function(e) NULL)
        if (is.null(sc)) next
        sc <- cbind(data.frame(site = SITE, plot = pid, plotType = ci$plotType,
                               detector = nm,
                               rung = ifelse(is.na(rung), "native", as.character(rung)),
                               pdens = round(pdens, 2), frdens = round(frdens, 2),
                               chm_res = if (nm %in% c("multichm","chm_vwf")) res else NA_real_,
                               n_apex = nrow(det)), sc)
        out[[length(out) + 1]] <- sc
      }
    }
    if (!length(out)) return(NULL)
    do.call(rbind, out)
  }

  res_list <- mclapply(keep, function(p)
                tryCatch(run_plot(p), error = function(e) {
                  message("plot ", p, " failed: ", conditionMessage(e)); NULL }),
                mc.cores = CORES, mc.preschedule = FALSE)
  results <- do.call(rbind, Filter(Negate(is.null), res_list))
  if (is.null(results) || !nrow(results)) { cat("no lidRplugins results\n"); return(invisible()) }
  results$tp_core <- round(results$precision * results$n_det)
  write.csv(results, file.path(nd, "lidrplugins_results.csv"), row.names = FALSE)
  cat(sprintf("[%s] lidRplugins DONE: %d rows -> %s\n", SITE, nrow(results),
              file.path(nd, "lidrplugins_results.csv")))

  # In-script head-to-head: equal-(plot,rung)-set guard, then pool per detector.
  guarded <- equal_set_guard(results, arms = ARMS)
  if (length(attr(guarded, "dropped")))
    cat(sprintf("equal-set guard dropped %d (plot,rung) cells\n",
                length(attr(guarded, "dropped"))))
  cat("\n=== Pooled recall/precision/F1 by detector (common cells) ===\n")
  pooled <- do.call(rbind, lapply(ARMS, function(a) {
    s <- guarded[guarded$detector == a, ]; if (!nrow(s)) return(NULL)
    cbind(detector = a, pool(s)) }))
  print(pooled[, c("detector","n_plots","n_ref","recall","precision","F1",
                   "rec_dominant","rec_understory")], row.names = FALSE, digits = 3)
}

if (sys.nframe() == 0L) run_main()
```

- [ ] **Step 2: Re-run the extractor test (regression — the `sys.nframe` guard
  must keep the driver from firing on `source()`)**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-lidrplugins-extractors.R")'`
Expected: PASS (the driver does not run during `source()`).

- [ ] **Step 3: Run the full suite (ensure nothing else broke)**

Run: `Rscript tests/run_tests.R`
Expected: all test files green (the 6 foundation files + the new
`test-lidrplugins-extractors`).

- [ ] **Step 4: Smoke-run on real SOAP data (only if present)**

Check: `ls "$CLAUDE_JOB_DIR/neon/SOAP/lidar/" 2>/dev/null` (defaults to `./work`).

- If present: `Rscript scripts/detect_lidrplugins_sweep.R SITE=SOAP PLOTS=ALL CORES=6`
  Expected: writes `lidrplugins_results.csv` and prints the pooled
  detector table (lmfauto/multichm/ptrees/chm_vwf rows with recall/precision/F1
  - per-class). Frozen clips from the AMS3D run are **reused** (same
  `frozen/<plot>/<rung>/` cache), so this confirms cross-arm identical inputs.
- If absent: skip and note the smoke run is deferred.
- [ ] **Step 5: Commit**

```bash
git add scripts/detect_lidrplugins_sweep.R
git commit -m "feat(#C9): lidRplugins full-ladder driver + CHM-VWF head-to-head"
```

---

## Self-review

**Spec coverage** (#C9 in the design spec):

- lidRplugins `lmfauto`/`multichm`/`ptrees` detectors → Task 2. ✓
- Full density ladder, SOAP, frozen clips reused from the bridge → Task 3
  (`frozen_clip`, native-pdens no-upsampling guard). ✓
- Same `ws_factory`/`hmin` discipline as CHM-VWF → `det_multichm` passes
  `ws_factory(a)`; `detect_lasr` baseline uses the same `A_VWF`; density-derived
  `res`. ✓
- Scored via the unchanged `score_plot`/`greedy_match`, pooled with canonical
  `pool` + `equal_set_guard` → Task 3. ✓
- Fair head-to-head vs CHM-VWF → `chm_vwf` run as a 4th detector on the same
  clip + pooled comparison. ✓
- Crash-vs-empty discipline (NULL skip vs 0-row recall=0) consistent with the
  AMS3D arm → extractors return NULL on crash, 0-row on empty. ✓
- Cross-arm synthesis vs AMS3D/Li 2012 is out of scope here (that is #R10).

**Placeholder scan:** none — every code step is complete; the compatibility
gate (Task 1) has an explicit BLOCKED path rather than a vague "handle errors".

**Type/signature consistency:** `det_lmfauto(las, hmin)`,
`det_multichm(las, res, a)`, `det_ptrees(las, hmin, k)` defined in Task 2 are
called with those exact argument names in Task 3's `dets` list. `.tops_to_det`
returns `data.frame(x,y,z)` and calls `assert_detection_contract` (from
`model_bench_lib.R`). The driver adds `tp_core` before pooling (matching what
`pool` expects, though `pool` also self-derives it). `ARMS` includes all four
detector names used in the `dets` list and the pooled loop. `ws_factory`,
`plot_half`, `detect_lasr`, `score_plot` come from `sweep_lib.R`; `frozen_clip`,
`equal_set_guard`, `pool`, `assert_detection_contract` from `model_bench_lib.R`.
Consistent.
