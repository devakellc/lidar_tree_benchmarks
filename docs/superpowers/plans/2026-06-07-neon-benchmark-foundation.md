# NEON Deep-Model Benchmark — Foundation Plan (#A0 + #B1 + #B2)

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> superpowers:subagent-driven-development (recommended) or
> superpowers:executing-plans to implement this plan task-by-task. Steps use
> checkbox (`- [ ]`) syntax for tracking.

**Goal:** Stand up the R-native foundation of the NEON model benchmark — a
deterministic frozen-clip provider, a universal predictions→detections reducer,
a canonical pooler with a per-(plot,rung) equal-set guard, a scorer-contract
conformance harness, and the first model arm (AMS3D) — all feeding the existing
`score_plot`/`greedy_match` harness unchanged.

**Architecture:** Walking skeleton first (B): build the AMS3D arm end-to-end on
the existing `sweep_lib.R` to prove the bridge in-language with zero GPU. Then
extract the reusable, hardened bridge (A) into `scripts/model_bench_lib.R` and
refactor the arm onto it. Every later arm (lidRplugins, the GPU models) consumes
this library.

**Tech stack:** R (`lidR`, `lasR`, `sf`, `terra`, `data.table`),
`crownsegmentr` (CRAN, AMS3D), `testthat` 3e for unit tests, `jsonlite` for the
clip manifest. Detections are scored by the existing `sweep_lib.R` functions.

**Scope note:** This is one of several plans from
`docs/superpowers/specs/2026-06-07-neon-deep-model-benchmark-design.md`. It
covers **#A0, #B1, #B2** only. `#C9` (lidRplugins) and the GPU fan-out
(`#I3`–`#M8`) are separate plans that reuse the library built here.

---

## File structure

| File | Responsibility | Tasks |
|---|---|---|
| `docs/model-benchmark-plan.md` | #A0 decision record: run/defer/drop set, zero-shot protocol, weights-mirror policy | 1 |
| `tests/run_tests.R` | Rscript entry point: runs `testthat::test_dir("tests/testthat")` | 2 |
| `tests/testthat/helper-synth.R` | Synthetic-LAS / labelled-points fixtures shared by tests | 2 |
| `scripts/detect_ams3d_sweep.R` | #B1→#B2 AMS3D arm: full ladder, scores via the harness, long-form CSV | 3, 4, 11 |
| `scripts/model_bench_lib.R` | #B2 shared bridge: `reduce_instances`, `crown_diameter_table`, `seed_for`, `frozen_clip`, `pool`, `equal_set_guard`, `assert_detection_contract` | 5–10 |
| `tests/testthat/test-*.R` | Unit tests per library function (synthetic inputs, no real data) | 3, 5–10 |

The library is the durable deliverable; the AMS3D script is its first consumer
and regression check. Tests use **synthetic in-memory inputs only** — real NEON
tiles are gitignored and may be absent, so no automated test depends on them.

---

## Task 1: #A0 — Triage + standing-protocol decision record

**Files:**

- Create: `docs/model-benchmark-plan.md`

This is a documentation task (no code, no test). It pins the decisions so later
arms don't relitigate them.

- [ ] **Step 1: Write the decision record**

Create `docs/model-benchmark-plan.md` with exactly this content:

````markdown
# NEON Model Benchmark — Triage & Standing Protocol (#A0)

Decision record for the deep-model benchmark. Source of truth for the runnable
set and the rules every arm must follow. Derived from
[the design spec](superpowers/specs/2026-06-07-neon-deep-model-benchmark-design.md).

## Runnable set

| Model | Disposition | Rung scope | Input variant |
|---|---|---|---|
| AMS3D (`crownsegmentr`) | run — walking skeleton | full ladder | normalized |
| lidRplugins (`lmfauto`/`multichm`/`ptree`) | run — competitor | full ladder | normalized |
| SegmentAnyTree | run | full ladder | raw-with-ground |
| TreeisoNet-ALS | run | full ladder | raw-with-ground |
| ForestFormer3D | run | native + 8 only | raw-with-ground |
| ForAINet | defer | — | — |
| HFC | defer (agreement) | — | — |
| TreeLearn | drop | — | — |
| Dersch graph-cut | drop | — | — |
| DeepForest (RGB) | optional reference | n/a (density-invariant) | RGB + CHM |

Baselines: existing CHM-VWF and Li 2012 arms.

## Zero-shot protocol

- No fine-tuning. Published weights / classical defaults only.
- Per-arm knobs MAY be set from prior literature/allometry; they MUST NOT be
  tuned on SOAP scoring outcomes.
- Every non-default knob is recorded in the run's ledger (the
  `*_ledger.csv` written next to each arm's results), with the value and its
  literature source.
- AMS3D knobs to record: `crown_diameter_to_tree_height`,
  `crown_length_to_tree_height`, `segment_crowns_only_above`.

## Weights-mirror policy

- Mirror every model weight/image to a project store before first use; record
  SHA256 + source URL. Upstream links (TreeisoNet personal server, ForAINet
  Dropbox) are single points of failure.
- Tracked in the GPU-infra plan (#I5), not here.

## Density framing

Rungs are all-return pts/m². SOAP native ≈ 20 all-return / ≈ 12 first-return.
Deep models clear their floor (10) at native; feasibility risk is the 4/2/1
rungs and the discrete-return-vs-ULS structural gap, not native sparsity.
````

- [ ] **Step 2: Lint and commit**

Run: `rumdl check docs/model-benchmark-plan.md`
Expected: `Success: No issues found`

```bash
git add docs/model-benchmark-plan.md
git commit -m "docs(#A0): triage + zero-shot protocol + mirror policy"
```

---

## Task 2: Test harness scaffolding

**Files:**

- Create: `tests/run_tests.R`
- Create: `tests/testthat/helper-synth.R`

The repo has no test suite; `testthat` 3.3.2 is installed. This adds a
minimal Rscript-runnable harness and the synthetic fixtures all unit tests use.

- [ ] **Step 1: Write the test runner**

Create `tests/run_tests.R`:

```r
#!/usr/bin/env Rscript
# Runs the benchmark library unit tests. Usage: Rscript tests/run_tests.R
suppressMessages(library(testthat))
res <- test_dir("tests/testthat", reporter = "summary", stop_on_failure = TRUE)
```

- [ ] **Step 2: Write the synthetic fixtures helper**

Create `tests/testthat/helper-synth.R`:

```r
# Synthetic fixtures shared across tests. No real NEON data is ever required.
suppressMessages({ library(lidR); library(data.table) })

# A labelled point table: two well-separated "trees", each a small vertical
# cluster. Columns X,Y,Z (metres) + a known instance id. Tree A apex (10,10,18),
# tree B apex (40,12,12).
synth_labelled_points <- function() {
  a <- data.table(X = c(10, 10.2, 9.8, 10.1), Y = c(10, 10.1, 9.9, 10.0),
                  Z = c(18, 14, 9, 4),  id = 1L)
  b <- data.table(X = c(40, 40.1, 39.9),       Y = c(12, 12.1, 11.9),
                  Z = c(12, 8, 3),      id = 2L)
  noise <- data.table(X = 25, Y = 25, Z = 1.0,  id = NA_integer_)  # unassigned
  rbind(a, b, noise)
}

# A minimal NORMALIZED LAS (ground at 0) with two clusters, for arm smoke tests.
synth_las_normalized <- function() {
  set.seed(1)
  mk <- function(cx, cy, top, n = 60) {
    data.frame(X = rnorm(n, cx, 0.7), Y = rnorm(n, cy, 0.7),
               Z = runif(n, 0, top))
  }
  df <- rbind(mk(10, 10, 18), mk(40, 12, 12))
  df$Z[1]  <- 18; df$Y[1]  <- 10; df$X[1]  <- 10   # guarantee a tree-A apex
  df$Z[61] <- 12; df$Y[61] <- 12; df$X[61] <- 40   # guarantee a tree-B apex
  las <- LAS(df)
  projection(las) <- 32611
  las
}
```

- [ ] **Step 3: Verify the harness runs (with zero tests yet)**

Run: `Rscript tests/run_tests.R`
Expected: runs without error (reports 0 failures; "No tests found" or an empty
summary is acceptable at this point).

- [ ] **Step 4: Commit**

```bash
git add tests/run_tests.R tests/testthat/helper-synth.R
git commit -m "test: add testthat runner + synthetic fixtures"
```

---

## Task 3: #B1 — AMS3D apex extractor (walking skeleton, part 1)

**Files:**

- Create: `scripts/detect_ams3d_sweep.R`
- Test: `tests/testthat/test-ams3d-extractor.R`

Build the AMS3D detector function first, test-first. It runs `crownsegmentr`'s
`segment_tree_crowns()` on a normalized clip and collapses each `crown_id` to a
max-Z apex `data.frame(x,y,z)` — the same apex rule the Li 2012 arm uses.

**Install dependency first.**

- [ ] **Step 1: Install crownsegmentr**

Run:

```bash
Rscript -e 'if (!requireNamespace("crownsegmentr", quietly=TRUE)) \
  install.packages("crownsegmentr", repos="https://cloud.r-project.org")'
```

Expected: installs `crownsegmentr` (and `dbscan`) without error.
Verify: `Rscript -e 'packageVersion("crownsegmentr")'` prints a version.

- [ ] **Step 2: Write the failing test**

Create `tests/testthat/test-ams3d-extractor.R`:

```r
source(file.path("scripts", "detect_ams3d_sweep.R"), local = TRUE)

test_that("det_ams3d returns the x,y,z detection contract on a 2-tree cloud", {
  las <- synth_las_normalized()
  det <- det_ams3d(las, cd_ratio = 0.4, cl_ratio = 0.8, min_above = 2)
  expect_s3_class(det, "data.frame")
  expect_false(inherits(det, "sf"))
  expect_identical(names(det), c("x", "y", "z"))
  expect_true(all(vapply(det, is.numeric, logical(1))))
  expect_gte(nrow(det), 1L)          # at least one crown found
  expect_true(max(det$z) > 8)        # apex is a real height, normalized
})

test_that("det_ams3d returns an empty 0-row frame, never NULL, when no crown", {
  empty <- LAS(data.frame(X = 0, Y = 0, Z = 0.1))
  projection(empty) <- 32611
  det <- det_ams3d(empty, cd_ratio = 0.4, cl_ratio = 0.8, min_above = 2)
  expect_s3_class(det, "data.frame")
  expect_identical(names(det), c("x", "y", "z"))
  expect_equal(nrow(det), 0L)
})
```

- [ ] **Step 3: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-ams3d-extractor.R")'`
Expected: FAIL — `could not find function "det_ams3d"` (the script doesn't
define it yet).

- [ ] **Step 4: Write the minimal script with det_ams3d**

Create `scripts/detect_ams3d_sweep.R`:

```r
#!/usr/bin/env Rscript
# AMS3D (adaptive mean shift, crownsegmentr) arm of the NEON model benchmark.
# Walking skeleton (#B1): runs segment_tree_crowns on the normalized clip per
# plot x density rung, collapses each crown_id to its max-Z apex (x,y,z), and
# scores against field stems with the existing sweep_lib harness.
#
# Usage:
#   Rscript scripts/detect_ams3d_sweep.R [SITE=SOAP] [PLOTS=ALL] [CORES=6]
#       [TOL=4] [CD_RATIO=0.4] [CL_RATIO=0.8]
# Output: $CLAUDE_JOB_DIR/neon/<SITE>/ams3d_results.csv (one row per
#         plot x rung) + ams3d_ledger.csv (recorded zero-shot knobs).
suppressMessages({ library(lidR); library(sf); library(data.table)
                   library(parallel); library(crownsegmentr) })
options(lidR.progress = FALSE)
d <- Sys.getenv("CLAUDE_JOB_DIR", file.path(getwd(), "work"))
source(file.path("scripts", "sweep_lib.R"))

## ---- AMS3D apex extractor ------------------------------------------------
# las: a NORMALIZED lidR::LAS (ground at 0). Returns data.frame(x,y,z) of
# per-crown apexes (max-Z point of each crown_id). cd_ratio/cl_ratio are the
# crown-diameter/length-to-tree-height allometric ratios (zero-shot knobs).
det_ams3d <- function(las, cd_ratio = 0.4, cl_ratio = 0.8, min_above = 2) {
  empty <- data.frame(x = numeric(), y = numeric(), z = numeric())
  if (is.empty(las) || npoints(las) < 5) return(empty)
  seg <- tryCatch(
    crownsegmentr::segment_tree_crowns(
      point_cloud = las,
      crown_diameter_to_tree_height = cd_ratio,
      crown_length_to_tree_height   = cl_ratio,
      segment_crowns_only_above     = min_above,
      ground_height                 = NULL,        # input already normalized
      crown_id_column_name          = "crown_id"),
    error = function(e) NULL)
  if (is.null(seg)) return(empty)
  dt <- as.data.table(seg@data)[!is.na(crown_id), .(X, Y, Z, crown_id)]
  if (!nrow(dt)) return(empty)
  ap <- dt[, .(x = X[which.max(Z)], y = Y[which.max(Z)], z = max(Z)),
           by = crown_id]
  data.frame(x = ap$x, y = ap$y, z = ap$z)
}
```

- [ ] **Step 5: Run test to verify it passes**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-ams3d-extractor.R")'`
Expected: PASS (2 tests). If AMS3D finds >2 crowns on the synthetic cloud that
is fine — the test only asserts `>= 1` and the contract shape.

- [ ] **Step 6: Commit**

```bash
git add scripts/detect_ams3d_sweep.R tests/testthat/test-ams3d-extractor.R
git commit -m "feat(#B1): AMS3D apex extractor + contract test"
```

---

## Task 4: #B1 — AMS3D arm driver over the full ladder (walking skeleton, part 2)

**Files:**

- Modify: `scripts/detect_ams3d_sweep.R` (append the driver)

Wire the extractor into a full-ladder per-plot sweep using the **existing**
`prepare_clip`/`score_plot` from `sweep_lib.R`, mirroring `run_sweep.R`'s ladder
and no-upsampling guard. This is the end-to-end spike; it is integration code
verified by a smoke run (real data permitting), not a unit test.

- [ ] **Step 1: Append the driver to the script**

Append to `scripts/detect_ams3d_sweep.R`:

```r
## ---- args ----------------------------------------------------------------
args  <- strsplit(commandArgs(TRUE), "=")
A     <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
SITE  <- if (is.null(A$SITE))  "SOAP" else A$SITE
PLOTS <- if (is.null(A$PLOTS) || A$PLOTS == "ALL") NULL else strsplit(A$PLOTS, ",")[[1]]
CORES <- as.integer(if (is.null(A$CORES)) 6 else A$CORES)
TOL   <- as.numeric(if (is.null(A$TOL)) 4.0 else A$TOL)
CD    <- as.numeric(if (is.null(A$CD_RATIO)) 0.4 else A$CD_RATIO)
CL    <- as.numeric(if (is.null(A$CL_RATIO)) 0.8 else A$CL_RATIO)
RUNGS <- c(8, 4, 2, 1)
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
  cat(sprintf("[%s] AMS3D plots: %d (%s)\n", SITE, length(keep),
              paste(keep, collapse = ",")))
  tmpdir <- file.path(tempdir(), paste0("ams3d_", SITE))
  dir.create(tmpdir, showWarnings = FALSE, recursive = TRUE)

  run_plot <- function(pid) {
    ci <- pc[pc$plotID == pid, ][1, ]
    cx <- ci$easting; cy <- ci$northing
    ph <- plot_half(ci$plotType)
    stems <- gt[gt$plotID == pid & abs(gt$E - cx) <= ph & abs(gt$N - cy) <= ph, ]
    if (nrow(stems) < 1) return(NULL)
    out <- list(); native_pdens <- NA_real_
    for (rung in c(NA, RUNGS)) {
      prep <- tryCatch(prepare_clip(ctg, cx, cy, rung, tmpdir, core_half = ph),
                       error = function(e) NULL)
      if (is.null(prep)) next
      pdens <- prep$pdens; frdens <- prep$frdens
      if (is.na(rung)) native_pdens <- pdens
      else if (is.na(native_pdens) || rung >= native_pdens) { unlink(prep$file); next }
      las <- tryCatch(readLAS(prep$file), error = function(e) NULL)
      if (!is.null(las) && !is.empty(las)) {
        det <- det_ams3d(las, cd_ratio = CD, cl_ratio = CL, min_above = 2)
        sc  <- score_plot(stems, det, tol_xy = TOL, core_cx = cx, core_cy = cy,
                          core_half = ph)
        sc  <- cbind(data.frame(site = SITE, plot = pid, plotType = ci$plotType,
                                detector = "ams3d",
                                rung = ifelse(is.na(rung), "native", as.character(rung)),
                                pdens = round(pdens, 2), frdens = round(frdens, 2),
                                n_apex = nrow(det)), sc)
        out[[length(out) + 1]] <- sc
      }
      unlink(prep$file)
    }
    if (!length(out)) return(NULL)
    do.call(rbind, out)
  }

  res_list <- mclapply(keep, function(p)
                tryCatch(run_plot(p), error = function(e) {
                  message("plot ", p, " failed: ", conditionMessage(e)); NULL }),
                mc.cores = CORES, mc.preschedule = FALSE)
  results <- do.call(rbind, Filter(Negate(is.null), res_list))
  if (is.null(results)) { cat("no AMS3D results\n"); return(invisible()) }
  write.csv(results, file.path(nd, "ams3d_results.csv"), row.names = FALSE)
  write.csv(data.frame(knob = c("crown_diameter_to_tree_height",
                                "crown_length_to_tree_height",
                                "segment_crowns_only_above"),
                       value = c(CD, CL, 2),
                       source = "literature default (Ferraz 2016 allometry)"),
            file.path(nd, "ams3d_ledger.csv"), row.names = FALSE)
  cat(sprintf("[%s] AMS3D DONE: %d rows -> %s\n", SITE, nrow(results),
              file.path(nd, "ams3d_results.csv")))
}

if (sys.nframe() == 0L) run_main()
```

Note the `if (sys.nframe() == 0L) run_main()` guard: it lets the test
`source()` the script to get `det_ams3d` without triggering the data-dependent
driver.

- [ ] **Step 2: Re-run the extractor test to confirm sourcing still works**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-ams3d-extractor.R")'`
Expected: PASS (2 tests) — the `sys.nframe()` guard prevents the driver from
running during `source()`.

- [ ] **Step 3: Smoke-run on real SOAP data (only if present)**

Check:

```bash
ls "$CLAUDE_JOB_DIR/neon/SOAP/ground_truth_stems.csv" \
   "$CLAUDE_JOB_DIR/neon/SOAP/lidar/" 2>/dev/null
```

- If present:
  Run: `Rscript scripts/detect_ams3d_sweep.R SITE=SOAP PLOTS=ALL CORES=6`
  Expected: prints `AMS3D DONE: N rows` and writes `ams3d_results.csv` with
  columns `site,plot,plotType,detector,rung,pdens,frdens,n_apex,n_ref,n_det,TP,
  recall,precision,F1,height_rmse,rec_dominant,...`.
- If absent: skip; note in the commit that the smoke run is deferred until SOAP
  data is fetched (`scripts/neon_ground_truth.R` + `scripts/neon_download_lidar.R`).
- [ ] **Step 4: Commit**

```bash
git add scripts/detect_ams3d_sweep.R
git commit -m "feat(#B1): AMS3D full-ladder driver (walking skeleton end-to-end)"
```

---

## Task 5: #B2 — `reduce_instances()` universal reducer

**Files:**

- Create: `scripts/model_bench_lib.R`
- Test: `tests/testthat/test-reduce-instances.R`

Generalize the per-instance max-Z apex collapse (currently duplicated inline in
`det_ams3d` and `det_lidr_li2012`) into one reducer that every arm uses.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-reduce-instances.R`:

```r
source(file.path("scripts", "model_bench_lib.R"), local = TRUE)

test_that("reduce_instances collapses each id to its max-Z apex as x,y,z", {
  pts <- synth_labelled_points()          # tree 1 apex (10,10,18); tree 2 (40,12,12)
  det <- reduce_instances(pts, id_col = "id")
  expect_s3_class(det, "data.frame")
  expect_false(inherits(det, "sf"))
  expect_identical(names(det), c("x", "y", "z"))
  expect_equal(nrow(det), 2L)             # NA-id noise point dropped
  det <- det[order(det$x), ]
  expect_equal(det$z, c(18, 12))
  expect_equal(det$x, c(10, 40))
  expect_equal(det$y, c(10, 12))
})

test_that("reduce_instances returns a 0-row x,y,z frame on all-NA ids", {
  pts <- data.table::data.table(X = 1, Y = 1, Z = 1, id = NA_integer_)
  det <- reduce_instances(pts, id_col = "id")
  expect_identical(names(det), c("x", "y", "z"))
  expect_equal(nrow(det), 0L)
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-reduce-instances.R")'`
Expected: FAIL — `could not find function "reduce_instances"`.

- [ ] **Step 3: Write the minimal library with reduce_instances**

Create `scripts/model_bench_lib.R`:

```r
#!/usr/bin/env Rscript
# Shared bridge for the NEON model benchmark (#B2). Sourced by every model arm.
# Provides: reduce_instances (universal predictions->detections collapse),
# crown_diameter_table (optional side metric), seed_for + frozen_clip
# (deterministic two-variant clip provider), pool + equal_set_guard (canonical
# pooling), assert_detection_contract (scorer-contract conformance harness).
suppressMessages({ library(data.table); library(lidR); library(terra)
                   library(jsonlite) })

## ---- universal reducer: labelled points -> data.frame(x,y,z) -------------
# pts: data.frame/data.table with coordinate columns (default X,Y,Z) and an
# instance-id column (default "crown_id"; NA = unassigned, dropped). Apex per
# instance = its max-Z point. Returns a base data.frame with exactly x,y,z.
reduce_instances <- function(pts, id_col = "crown_id",
                             x = "X", y = "Y", z = "Z") {
  empty <- data.frame(x = numeric(), y = numeric(), z = numeric())
  dt <- as.data.table(pts)
  if (!nrow(dt) || !id_col %in% names(dt)) return(empty)
  dt <- dt[!is.na(dt[[id_col]]), c(x, y, z, id_col), with = FALSE]
  if (!nrow(dt)) return(empty)
  setnames(dt, c(x, y, z, id_col), c("X", "Y", "Z", "ID"))
  ap <- dt[, .(x = X[which.max(Z)], y = Y[which.max(Z)], z = max(Z)), by = ID]
  data.frame(x = ap$x, y = ap$y, z = ap$z)
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-reduce-instances.R")'`
Expected: PASS (2 tests).

- [ ] **Step 5: Commit**

```bash
git add scripts/model_bench_lib.R tests/testthat/test-reduce-instances.R
git commit -m "feat(#B2): universal predictions->detections reducer"
```

---

## Task 6: #B2 — `crown_diameter_table()` side metric with point-count floor

**Files:**

- Modify: `scripts/model_bench_lib.R`
- Test: `tests/testthat/test-crown-diameter.R`

Crown diameter is **not** in the detection contract (the scorer ignores it). It
is a separate optional add-on keyed off the same instance groups, with a
point-count floor so a caliper over a handful of points is never reported.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-crown-diameter.R`:

```r
source(file.path("scripts", "model_bench_lib.R"), local = TRUE)

test_that("crown_diameter_table reports d_eq, d_caliper, n_pts per instance", {
  # one wide instance (id=1, ~4 m caliper) and one sparse instance (id=2, 2 pts)
  pts <- data.table::data.table(
    X = c(0, 4, 0, 4,  100, 100.5),
    Y = c(0, 0, 4, 4,  100, 100),
    Z = c(5, 5, 5, 5,  3,   3),
    crown_id = c(1, 1, 1, 1,  2, 2))
  tab <- crown_diameter_table(pts, id_col = "crown_id", min_pts = 4)
  expect_identical(sort(names(tab)),
                   sort(c("id", "n_pts", "d_eq", "d_caliper")))
  r1 <- tab[tab$id == 1, ]
  expect_equal(r1$n_pts, 4L)
  expect_equal(round(r1$d_caliper, 2), round(sqrt(32), 2))  # diagonal of 4x4
  # sparse instance below the floor: diameters NA, but row + n_pts still present
  r2 <- tab[tab$id == 2, ]
  expect_equal(r2$n_pts, 2L)
  expect_true(is.na(r2$d_eq) && is.na(r2$d_caliper))
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-crown-diameter.R")'`
Expected: FAIL — `could not find function "crown_diameter_table"`.

- [ ] **Step 3: Add crown_diameter_table to the library**

Append to `scripts/model_bench_lib.R`:

```r
## ---- optional crown-diameter side metric (NOT in the detection contract) -
# Per instance: d_eq = 2*sqrt(area/pi) from the 2-D convex hull area;
# d_caliper = max pairwise distance among the instance's points. Both NA when
# the instance has fewer than `min_pts` points (a diameter over a few points is
# noise). Returns data.frame(id, n_pts, d_eq, d_caliper).
crown_diameter_table <- function(pts, id_col = "crown_id", min_pts = 5) {
  dt <- as.data.table(pts)
  dt <- dt[!is.na(dt[[id_col]])]
  if (!nrow(dt)) return(data.frame(id = integer(), n_pts = integer(),
                                   d_eq = numeric(), d_caliper = numeric()))
  setnames(dt, id_col, "ID")
  one <- function(s) {
    n <- nrow(s)
    if (n < min_pts) return(list(n_pts = n, d_eq = NA_real_, d_caliper = NA_real_))
    h  <- grDevices::chull(s$X, s$Y)
    hx <- s$X[h]; hy <- s$Y[h]
    area <- abs(sum(hx * c(hy[-1], hy[1]) - c(hx[-1], hx[1]) * hy)) / 2
    cal  <- max(dist(cbind(s$X, s$Y)))
    list(n_pts = n, d_eq = 2 * sqrt(area / pi), d_caliper = cal)
  }
  res <- dt[, one(.SD), by = ID, .SDcols = c("X", "Y")]
  data.frame(id = res$ID, n_pts = res$n_pts, d_eq = res$d_eq,
             d_caliper = res$d_caliper)
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-crown-diameter.R")'`
Expected: PASS (1 test, 5 expectations).

- [ ] **Step 5: Commit**

```bash
git add scripts/model_bench_lib.R tests/testthat/test-crown-diameter.R
git commit -m "feat(#B2): crown-diameter side metric with point-count floor"
```

---

## Task 7: #B2 — `seed_for()` + `frozen_clip()` deterministic provider

**Files:**

- Modify: `scripts/model_bench_lib.R`
- Test: `tests/testthat/test-frozen-clip.R`

The provider decimates the raw clip **exactly once** with a seed keyed by
`site/plot/rung`, then derives both the raw-with-ground and normalized variants
from that *same* decimated set, persists the DTM (raw→height transform) and a
manifest, and reuses cached outputs on re-call. This is what guarantees every
arm sees an identical point subset.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-frozen-clip.R`:

```r
source(file.path("scripts", "model_bench_lib.R"), local = TRUE)
suppressMessages(library(lidR))

test_that("seed_for is deterministic and varies by key", {
  expect_equal(seed_for("SOAP", "SOAP_001", 8), seed_for("SOAP", "SOAP_001", 8))
  expect_false(seed_for("SOAP", "SOAP_001", 8) == seed_for("SOAP", "SOAP_001", 4))
  expect_false(seed_for("SOAP", "SOAP_001", 8) == seed_for("SOAP", "SOAP_002", 8))
  expect_type(seed_for("SOAP", "SOAP_001", 8), "integer")
})

test_that("seeded homogenize decimation is reproducible", {
  set.seed(123)
  big <- LAS(data.frame(X = runif(5000, 0, 50), Y = runif(5000, 0, 50),
                        Z = runif(5000, 0, 30)))
  s <- seed_for("SOAP", "SOAP_001", 8)
  set.seed(s); a <- decimate_points(big, homogenize(density = 8, res = 5))
  set.seed(s); b <- decimate_points(big, homogenize(density = 8, res = 5))
  expect_equal(npoints(a), npoints(b))
  expect_equal(a@data$X, b@data$X)     # identical point subset, not just count
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-frozen-clip.R")'`
Expected: FAIL — `could not find function "seed_for"`.

- [ ] **Step 3: Add seed_for and frozen_clip to the library**

Append to `scripts/model_bench_lib.R`:

```r
## ---- deterministic seed from a (site, plot, rung) key --------------------
# Maps the key string to a stable non-negative 31-bit integer (FNV-1a), so
# decimation is reproducible across arms and runs without external deps.
seed_for <- function(site, plot, rung) {
  key <- paste(site, plot, ifelse(is.na(rung), "native", rung), sep = "|")
  h <- 2166136261
  for (b in utf8ToInt(key)) h <- bitwAnd((bitwXor(h, b) * 16777619), 0x7FFFFFFF)
  as.integer(h)
}

## ---- frozen clip provider: one decimated set -> two variants + DTM -------
# Decimates the plot clip ONCE (seeded by site/plot/rung), then derives both the
# raw-with-ground and the normalized variant from that SAME decimated set, plus
# a TIN DTM and a JSON manifest. Caches under out_root; reuses on re-call.
# rung = NA means native (no decimation). Returns a list with file paths +
# densities, or NULL if the clip is unusable.
frozen_clip <- function(ctg, site, plot, rung, cx, cy, core_half, out_root,
                        buffer = 25) {
  rdir <- file.path(out_root, site, plot, ifelse(is.na(rung), "native", rung))
  fp <- list(rawground   = file.path(rdir, "clip_rawground.laz"),
             normalized  = file.path(rdir, "clip_normalized.laz"),
             dtm         = file.path(rdir, "ground_dtm.tif"),
             manifest    = file.path(rdir, "manifest.json"))
  if (file.exists(fp$manifest)) {                    # cached -> reuse verbatim
    mf <- jsonlite::read_json(fp$manifest, simplifyVector = TRUE)
    return(c(fp, list(pdens = mf$pdens, frdens = mf$frdens, seed = mf$seed)))
  }
  dir.create(rdir, showWarnings = FALSE, recursive = TRUE)
  half <- core_half + buffer
  las  <- lidR::clip_rectangle(ctg, cx - half, cy - half, cx + half, cy + half)
  if (lidR::is.empty(las) || lidR::npoints(las) < 100) return(NULL)
  seed <- seed_for(site, plot, rung)
  if (!is.na(rung)) { set.seed(seed)
    las <- lidR::decimate_points(las, lidR::homogenize(density = rung, res = 5)) }
  if (sum(las$Classification == 2L) < 10) return(NULL)   # need ground for DTM
  dtm <- lidR::rasterize_terrain(las, res = 1, algorithm = lidR::tin())
  nrm <- lidR::normalize_height(las, lidR::tin(), na.rm = TRUE)
  nrm <- lidR::filter_poi(nrm, Z >= -1, Z < 80)
  area   <- (2 * half)^2
  pdens  <- lidR::npoints(nrm) / area
  frdens <- sum(nrm$ReturnNumber == 1L) / area
  lidR::writeLAS(las, fp$rawground)
  lidR::writeLAS(nrm, fp$normalized)
  terra::writeRaster(dtm, fp$dtm, overwrite = TRUE)
  jsonlite::write_json(list(site = site, plot = plot,
                            rung = ifelse(is.na(rung), "native", rung),
                            seed = seed, n_raw = lidR::npoints(las),
                            n_norm = lidR::npoints(nrm),
                            pdens = round(pdens, 3), frdens = round(frdens, 3),
                            buffer = buffer, core_half = core_half),
                       fp$manifest, auto_unbox = TRUE, pretty = TRUE)
  c(fp, list(pdens = pdens, frdens = frdens, seed = seed))
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-frozen-clip.R")'`
Expected: PASS (2 tests). (The test exercises `seed_for` and seeded `homogenize`
reproducibility directly; `frozen_clip`'s file I/O is covered by the smoke run
in Task 11.)

- [ ] **Step 5: Commit**

```bash
git add scripts/model_bench_lib.R tests/testthat/test-frozen-clip.R
git commit -m "feat(#B2): deterministic seed + frozen-clip two-variant provider"
```

---

## Task 8: #B2 — canonical `pool()` + `equal_set_guard()`

**Files:**

- Modify: `scripts/model_bench_lib.R`
- Test: `tests/testthat/test-pool-guard.R`

Extract the one canonical pooler (today copy-pasted in four scripts) and the
per-(plot,rung) equal-set guard. Pooling sums counts (never averages rates); the
guard drops any `(site,plot,rung)` where some arm failed, so a rung's cross-arm
table compares identical plot populations.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-pool-guard.R`:

```r
source(file.path("scripts", "model_bench_lib.R"), local = TRUE)

mk_row <- function(site, plot, rung, det, n_ref, TP, n_det, tp_core) {
  data.frame(site = site, plot = plot, rung = rung, detector = det,
             n_ref = n_ref, TP = TP, n_det = n_det, tp_core = tp_core,
             n_dominant = n_ref, rec_dominant = TP / n_ref,
             n_codominant = 0, rec_codominant = NA_real_,
             n_intermediate = 0, rec_intermediate = NA_real_,
             n_suppressed = 0, rec_suppressed = NA_real_,
             frdens = 10, secs = 1)
}

test_that("pool sums counts: recall = sum(TP)/sum(n_ref)", {
  df <- rbind(mk_row("SOAP","p1","8","ams3d",10,5,6,4),
              mk_row("SOAP","p2","8","ams3d",90,45,50,40))
  p <- pool(df)
  expect_equal(p$n_ref, 100); expect_equal(p$TP, 50)
  expect_equal(p$recall, 0.5)            # NOT mean(0.5, 0.5)-by-rate weighting
  expect_equal(p$precision, 44 / 56)
})

test_that("equal_set_guard drops (site,plot,rung) missing any arm", {
  df <- rbind(mk_row("SOAP","p1","8","ams3d",10,5,6,4),
              mk_row("SOAP","p1","8","li2012",10,6,7,5),
              mk_row("SOAP","p2","8","ams3d",10,5,6,4))  # li2012 missing on p2/8
  g <- equal_set_guard(df, arms = c("ams3d","li2012"))
  expect_equal(nrow(g), 2L)                       # only p1/8 survives
  expect_true(all(paste(g$site,g$plot,g$rung) == "SOAP p1 8"))
  expect_equal(attr(g, "dropped"), "SOAP::p2::8")
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-pool-guard.R")'`
Expected: FAIL — `could not find function "pool"`.

- [ ] **Step 3: Add pool and equal_set_guard to the library**

Append to `scripts/model_bench_lib.R`:

```r
## ---- canonical pooler: sum counts, never average rates -------------------
# df: long-form scored rows (one per site x plot x rung x detector subset).
# Pools to a single row: recall = sum(TP)/sum(n_ref), precision =
# sum(tp_core)/sum(n_det), F1 from the pooled rates; per-class recall recovered
# from round(rec_<cls> * n_<cls>); understory = intermediate + suppressed.
POOL_CLASSES <- c("dominant", "codominant", "intermediate", "suppressed")
pool <- function(df, classes = POOL_CLASSES) {
  out <- data.frame(
    n_plots = length(unique(paste(df$site, df$plot, df$rung))),
    n_ref = sum(df$n_ref), n_det = sum(df$n_det), TP = sum(df$TP),
    recall = sum(df$TP) / sum(df$n_ref),
    precision = sum(df$tp_core, na.rm = TRUE) / sum(df$n_det))
  out$F1 <- if (!is.na(out$recall) && !is.na(out$precision) &&
                (out$recall + out$precision) > 0)
    2 * out$recall * out$precision / (out$recall + out$precision) else NA_real_
  for (cl in classes) {
    nref <- sum(df[[paste0("n_", cl)]], na.rm = TRUE)
    tp   <- sum(ifelse(df[[paste0("n_", cl)]] > 0,
                       round(df[[paste0("rec_", cl)]] * df[[paste0("n_", cl)]]), 0),
                na.rm = TRUE)
    out[[paste0("rec_", cl)]] <- if (nref) tp / nref else NA_real_
    out[[paste0("n_", cl)]]   <- nref
  }
  nref_u <- sum(df$n_intermediate, na.rm = TRUE) + sum(df$n_suppressed, na.rm = TRUE)
  tp_u <- sum(ifelse(df$n_intermediate > 0,
                     round(df$rec_intermediate * df$n_intermediate), 0), na.rm = TRUE) +
          sum(ifelse(df$n_suppressed > 0,
                     round(df$rec_suppressed * df$n_suppressed), 0), na.rm = TRUE)
  out$rec_understory <- if (nref_u) tp_u / nref_u else NA_real_
  out$n_understory   <- nref_u
  out
}

## ---- per-(plot,rung) equal-set guard -------------------------------------
# Keep only (site,plot,rung) cells scored by EVERY arm; drop the rest so a rung's
# cross-arm comparison uses an identical plot population. Returns the filtered
# df with the dropped cell keys in attr(.,"dropped").
equal_set_guard <- function(df, arms, key_cols = c("site", "plot", "rung")) {
  k <- do.call(paste, c(df[key_cols], sep = "::"))
  df$.k <- k
  n_arms <- tapply(df$detector, k, function(v) length(unique(v)))
  common <- names(n_arms)[n_arms == length(arms)]
  dropped <- setdiff(unique(k), common)
  out <- df[df$.k %in% common, setdiff(names(df), ".k"), drop = FALSE]
  attr(out, "dropped") <- dropped
  out
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-pool-guard.R")'`
Expected: PASS (2 tests).

- [ ] **Step 5: Commit**

```bash
git add scripts/model_bench_lib.R tests/testthat/test-pool-guard.R
git commit -m "feat(#B2): canonical pool() + per-(plot,rung) equal-set guard"
```

---

## Task 9: #B2 — `assert_detection_contract()` conformance harness

**Files:**

- Modify: `scripts/model_bench_lib.R`
- Test: `tests/testthat/test-contract.R`

A guard every arm calls before scoring: the detection table must be a base
`data.frame` with exactly numeric `x,y,z`, not an `sf` object, and empty must be
a 0-row frame (never NULL). Catches the common `sf::st_coordinates` capital-XYZ
mistake.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-contract.R`:

```r
source(file.path("scripts", "model_bench_lib.R"), local = TRUE)

test_that("assert_detection_contract passes a valid x,y,z frame", {
  ok <- data.frame(x = 1, y = 2, z = 3)
  expect_true(assert_detection_contract(ok))
  expect_true(assert_detection_contract(data.frame(x = numeric(), y = numeric(),
                                                    z = numeric())))
})

test_that("assert_detection_contract rejects contract violations", {
  expect_error(assert_detection_contract(NULL), "NULL")
  expect_error(assert_detection_contract(data.frame(X = 1, Y = 2, Z = 3)),
               "x, y, z")                                   # capital cols
  expect_error(assert_detection_contract(list(x = 1, y = 2, z = 3)),
               "data.frame")                                # not a data.frame
  expect_error(assert_detection_contract(
    data.frame(x = "a", y = 2, z = 3)), "numeric")          # non-numeric
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-contract.R")'`
Expected: FAIL — `could not find function "assert_detection_contract"`.

- [ ] **Step 3: Add assert_detection_contract to the library**

Append to `scripts/model_bench_lib.R`:

```r
## ---- scorer-contract conformance harness ---------------------------------
# Asserts `det` satisfies what score_plot consumes: a base data.frame (not sf),
# with exactly numeric columns x, y, z (lowercase). Empty must be a 0-row frame,
# never NULL. Returns TRUE invisibly or stops with a specific message.
assert_detection_contract <- function(det) {
  if (is.null(det)) stop("detection table is NULL; return a 0-row data.frame instead")
  if (!is.data.frame(det)) stop("detection table must be a base data.frame")
  if (inherits(det, "sf")) stop("detection table must not be an sf object")
  if (!identical(names(det), c("x", "y", "z")))
    stop("detection table must have exactly columns x, y, z (lowercase)")
  if (!all(vapply(det, is.numeric, logical(1))))
    stop("detection columns x, y, z must be numeric")
  invisible(TRUE)
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-contract.R")'`
Expected: PASS (2 tests, 6 expectations).

- [ ] **Step 5: Commit**

```bash
git add scripts/model_bench_lib.R tests/testthat/test-contract.R
git commit -m "feat(#B2): detection-contract conformance harness"
```

---

## Task 10: #B2 — full library test pass

**Files:** none (verification only)

- [ ] **Step 1: Run the whole suite**

Run: `Rscript tests/run_tests.R`
Expected: all tests across `test-ams3d-extractor`, `test-reduce-instances`,
`test-crown-diameter`, `test-frozen-clip`, `test-pool-guard`, `test-contract`
PASS; `stop_on_failure = TRUE` means a non-zero exit on any failure.

- [ ] **Step 2: Commit (only if any fixup was needed)**

If a cross-test inconsistency surfaced and you fixed it:

```bash
git add -A
git commit -m "test(#B2): green full library suite"
```

---

## Task 11: #B1↔#B2 — refactor the AMS3D arm onto the shared library

**Files:**

- Modify: `scripts/detect_ams3d_sweep.R`

Now that the bridge exists, refactor the walking skeleton to consume it:
`det_ams3d` uses `reduce_instances`; the driver uses `frozen_clip` (deterministic
clips) and `assert_detection_contract`. This is the regression check that the
extracted library reproduces the skeleton.

- [ ] **Step 1: Source the library and use reduce_instances in det_ams3d**

In `scripts/detect_ams3d_sweep.R`, after the existing
`source(file.path("scripts", "sweep_lib.R"))` line, add:

```r
source(file.path("scripts", "model_bench_lib.R"))
```

Replace the body of `det_ams3d` (the `dt`/`ap`/`data.frame(...)` tail after the
`if (is.null(seg)) return(empty)` line) with:

```r
  det <- reduce_instances(seg@data, id_col = "crown_id", x = "X", y = "Y", z = "Z")
  assert_detection_contract(det)
  det
}
```

- [ ] **Step 2: Use frozen_clip in the driver**

In `run_plot` inside `run_main`, replace the `prepare_clip(...)` call and the
subsequent `readLAS(prep$file)` with the frozen provider's normalized variant:

```r
    prep <- tryCatch(frozen_clip(ctg, SITE, pid, rung, cx, cy, ph,
                                 out_root = file.path(nd, "frozen")),
                     error = function(e) NULL)
    if (is.null(prep)) next
    pdens <- prep$pdens; frdens <- prep$frdens
    if (is.na(rung)) native_pdens <- pdens
    else if (is.na(native_pdens) || rung >= native_pdens) next
    las <- tryCatch(readLAS(prep$normalized), error = function(e) NULL)
```

Remove the now-unneeded `unlink(prep$file)` lines (frozen clips are cached
deliberately, not deleted).

- [ ] **Step 3: Re-run the extractor test (regression)**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-ams3d-extractor.R")'`
Expected: PASS (2 tests) — `det_ams3d` still satisfies the contract via
`reduce_instances`.

- [ ] **Step 4: Re-run the full suite**

Run: `Rscript tests/run_tests.R`
Expected: all PASS.

- [ ] **Step 5: Smoke-run on real SOAP data (only if present)**

If `"$CLAUDE_JOB_DIR/neon/SOAP/lidar/"` exists:
Run: `Rscript scripts/detect_ams3d_sweep.R SITE=SOAP PLOTS=ALL CORES=6`
Expected: writes `ams3d_results.csv` (full ladder rows) and reuses
`frozen/SOAP/<plot>/<rung>/` clips on a second run (idempotent — check the
manifests exist and the second run is faster).

- [ ] **Step 6: Commit**

```bash
git add scripts/detect_ams3d_sweep.R
git commit -m "refactor(#B1): AMS3D arm consumes model_bench_lib (frozen clips + reducer)"
```

---

## Self-review

**Spec coverage** (against §3–§4 of the design spec):

- #A0 triage + zero-shot protocol + mirror policy → Task 1. ✓
- #B1 AMS3D arm end-to-end, full ladder → Tasks 3, 4, 11. ✓
- #B2 frozen-clip provider (one decimation, two variants, DTM, manifest, seed)
  → Task 7. ✓
- #B2 universal reducer (x,y,z, standardized apex) → Task 5. ✓
- #B2 crown-diameter side metric with point-count floor → Task 6. ✓
- #B2 canonical `pool()` + per-(plot,rung) equal-set guard → Task 8. ✓
- #B2 conformance harness → Task 9. ✓
- Raw-ground→height DTM transform → persisted in Task 7 (`ground_dtm.tif`);
  *consumed* by the raw-ground GPU arms in their own plan (not this one, which
  only has the normalized-input AMS3D arm). Noted, not a gap.

**Out of this plan's scope (separate plans, by design):** #C9 lidRplugins,

## I3–#I5 GPU infra, #M6–#M8 GPU arms, #R10 result doc, #E11–#E14 extensions

**Placeholder scan:** none — every code step has complete code; every run step
has an exact command + expected result.

**Type/signature consistency:** `reduce_instances(pts, id_col, x, y, z)` defined
in Task 5 is called with `id_col="crown_id", x="X", y="Y", z="Z"` in Task 11;
`pool()` columns (`n_ref,TP,n_det,tp_core,n_<cls>,rec_<cls>`) match the row
schema `score_plot` emits plus the `tp_core` the driver adds; `frozen_clip`
returns `$normalized`/`$pdens`/`$frdens` as consumed in Task 11;
`assert_detection_contract` enforces the exact `names == c("x","y","z")` the
reducer produces. Consistent.
