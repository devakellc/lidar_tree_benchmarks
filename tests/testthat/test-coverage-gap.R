# Unit tests for the #V5 coverage-gap crediting helpers: co_detect_credit +
# fp_points (sweep_lib.R), the pool() corrected-precision extension
# (model_bench_lib.R), and the coverage_lib.R loaders (arm_family,
# read_arm_cache, boxes_to_dets). All synthetic; no NEON data.
source(file.path("..", "..", "scripts", "sweep_lib.R"), local = TRUE)
source(file.path("..", "..", "scripts", "model_bench_lib.R"), local = TRUE)
source(file.path("..", "..", "scripts", "coverage_lib.R"), local = TRUE)

## ---- co_detect_credit ------------------------------------------------------
# An isolated FP is credited as probable-real only when witness detections from
# >= min_fam DISTINCT families sit within r of it.
test_that("co_detect_credit requires min_fam distinct witness families within r", {
  fpx <- c(0, 10, 20); fpy <- c(0, 0, 0)
  wx <- c(0.5, -0.5, 10.3, 10.8); wy <- c(0, 0, 0, 0)
  wfam <- c("chm", "pc", "rgb", "rgb")
  # FP1: chm + pc within 2 m -> 2 families -> credited.
  # FP2: two rgb witnesses -> 1 family -> not credited.
  # FP3: no witnesses -> not credited.
  expect_equal(co_detect_credit(fpx, fpy, wx, wy, wfam, r = 2, min_fam = 2),
               c(TRUE, FALSE, FALSE))
  # min_fam = 1: a single confirming family is enough (FP2 flips).
  expect_equal(co_detect_credit(fpx, fpy, wx, wy, wfam, r = 2, min_fam = 1),
               c(TRUE, TRUE, FALSE))
})

test_that("co_detect_credit honours the radius", {
  # Witnesses at 1.9 m and 2.1 m from the FP: only the first is inside r = 2.
  cr <- co_detect_credit(0, 0, c(1.9, 2.1), c(0, 0), c("chm", "pc"),
                         r = 2, min_fam = 2)
  expect_false(cr)                       # only one family inside r
  expect_true(co_detect_credit(0, 0, c(1.9, 2.1), c(0, 0), c("chm", "pc"),
                               r = 2.5, min_fam = 2))
})

test_that("co_detect_credit handles empty inputs", {
  expect_equal(co_detect_credit(numeric(0), numeric(0),
                                c(1, 2), c(0, 0), c("chm", "pc")), logical(0))
  expect_equal(co_detect_credit(c(0, 5), c(0, 0),
                                numeric(0), numeric(0), character(0)),
               c(FALSE, FALSE))
})

## ---- fp_points -------------------------------------------------------------
# Core false positives with coordinates + the near/isolated split, consistent
# with score_plot's fp_near/fp_isolated counts on the same input.
test_that("fp_points returns core FPs with the same near/isolated split as score_plot", {
  # One stem at the core centre, matched by det1. det2 sits 3 m from the stem
  # (near, over-seg); det3 sits isolated at (15, 15); det4 is outside the core.
  stems <- data.frame(E = 100, N = 100, height = 20, crown_class = "dominant")
  det <- data.frame(x = c(100.5, 103, 115, 130), y = c(100, 100, 115, 100),
                    z = c(19, 18, 12, 15))
  sc <- score_plot(stems, det, tol_xy = 4, core_cx = 100, core_cy = 100,
                   core_half = 20)
  fp <- fp_points(stems, det, tol_xy = 4, core_cx = 100, core_cy = 100,
                  core_half = 20)
  expect_equal(nrow(fp), sc$fp_near + sc$fp_isolated)
  expect_equal(sum(!fp$isolated), sc$fp_near)
  expect_equal(sum(fp$isolated), sc$fp_isolated)
  expect_setequal(fp$x, c(103, 115))     # matched det1 and out-of-core det4 excluded
})

test_that("fp_points returns a 0-row frame when there is nothing to flag", {
  stems <- data.frame(E = 100, N = 100, height = 20, crown_class = "dominant")
  none <- fp_points(stems, data.frame(x = numeric(), y = numeric(), z = numeric()),
                    tol_xy = 4, core_cx = 100, core_cy = 100, core_half = 20)
  expect_equal(nrow(none), 0L)
  expect_true(all(c("x", "y", "z", "isolated") %in% names(none)))
  # a single matched detection -> no FPs
  one <- fp_points(stems, data.frame(x = 100.5, y = 100, z = 19),
                   tol_xy = 4, core_cx = 100, core_cy = 100, core_half = 20)
  expect_equal(nrow(one), 0L)
})

## ---- pool() corrected precision --------------------------------------------
test_that("pool derives corrected precision from summed fp_credited counts", {
  df <- data.frame(site = "S", plot = c("a", "b"), rung = "native",
                   n_ref = c(10, 10), n_det = c(10, 5), TP = c(6, 3),
                   tp_core = c(6, 3), recall = c(0.6, 0.3),
                   precision = c(0.6, 0.6), fp_credited = c(2, 1))
  p <- pool(df)
  expect_equal(p$precision, 9 / 15)                 # raw unchanged
  expect_equal(p$fp_credited, 3L)
  expect_equal(p$precision_cred, 9 / 12)            # denominator loses credits
  expect_equal(p$F1_cred, 2 * 0.45 * 0.75 / (0.45 + 0.75))
})

test_that("pool omits corrected columns when fp_credited is absent", {
  df <- data.frame(site = "S", plot = "a", rung = "native",
                   n_ref = 10, n_det = 10, TP = 6, tp_core = 6,
                   recall = 0.6, precision = 0.6)
  p <- pool(df)
  expect_null(p$precision_cred)
  expect_null(p$F1_cred)
})

test_that("pool corrected precision never exceeds 1 (denominator floored at TP)", {
  # Degenerate: crediting claims more FPs than exist beyond the TPs.
  df <- data.frame(site = "S", plot = "a", rung = "native",
                   n_ref = 10, n_det = 5, TP = 4, tp_core = 4,
                   recall = 0.4, precision = 0.8, fp_credited = 3)
  p <- pool(df)
  expect_equal(p$precision_cred, 1.0)               # 4 / max(5-3, 4)
})

## ---- arm_family ------------------------------------------------------------
test_that("arm_family maps every benchmark arm to its modality family", {
  expect_equal(arm_family(c("chm_vwf", "lmfauto", "multichm")),
               rep("chm", 3))
  expect_equal(arm_family(c("ptrees", "ams3d", "li2012", "lidr_li2012",
                            "lidr_lmf_pc", "lasr_lmax_pc")), rep("pc", 6))
  expect_equal(arm_family(c("treeisonet", "segmentanytree", "forestformer3d")),
               rep("deep", 3))
  expect_equal(arm_family(c("deepforest", "detectree2")), rep("rgb", 2))
  expect_true(is.na(arm_family("nonesuch")))
})

## ---- read_arm_cache --------------------------------------------------------
test_that("read_arm_cache reads plain and param-suffixed cache files", {
  tmp <- file.path(tempdir(), "btc_test"); dir.create(tmp, showWarnings = FALSE)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  det <- data.frame(x = c(1, 2), y = c(3, 4), z = c(10, 20))
  write.csv(det, file.path(tmp, "ams3d__SOAP__SOAP_001__2.csv"), row.names = FALSE)
  write.csv(det, file.path(tmp, "chm_vwf__SOAP__SOAP_001__4__res0.5__a0.05.csv"),
            row.names = FALSE)
  expect_equal(read_arm_cache(tmp, "ams3d", "SOAP", "SOAP_001", "2")$z, c(10, 20))
  expect_equal(nrow(read_arm_cache(tmp, "chm_vwf", "SOAP", "SOAP_001", "4")), 2L)
  expect_null(read_arm_cache(tmp, "ams3d", "SOAP", "SOAP_001", "native"))
  # a "li2012" query must not swallow "lidr_li2012" files
  write.csv(det, file.path(tmp, "lidr_li2012__SOAP__SOAP_002__8.csv"),
            row.names = FALSE)
  expect_null(read_arm_cache(tmp, "li2012", "SOAP", "SOAP_002", "8"))
})

## ---- credit eligibility vs ANY mapped stem (#96 review) ----------------------
test_that("fp_points flags FPs near any mapped stem as credit-ineligible", {
  # Stem A is matched; stem B is unmatched (its only nearby detection fails the
  # height gate). An FP 2.5 m from unmapped... from UNMATCHED stem B is
  # `isolated` under the #V4 split (not near a MATCHED stem) but must NOT be
  # credit-eligible: the field map has a stem right there.
  stems <- data.frame(E = c(100, 115), N = c(100, 100), height = c(20, 15),
                      crown_class = c("dominant", "codominant"))
  det <- data.frame(x = c(100.5, 112.5, 140), y = c(100, 100, 130),
                    z = c(19, 2.0, 10))      # det2 height-gated off stem B
  fp <- fp_points(stems, det, tol_xy = 4, core_cx = 120, core_cy = 115,
                  core_half = 50)
  expect_equal(nrow(fp), 2L)                 # det2 + det3 are core FPs
  expect_true(all(fp$isolated))              # neither is near a MATCHED stem
  expect_equal(fp$credit_eligible[order(fp$x)], c(FALSE, TRUE))
})

test_that("fp_points eligibility sees mapped stems outside the scored core", {
  # NEON maps stems whose reconstructed position lands just outside the plot
  # core (15-47 per site here). They never enter n_ref, but a core FP sitting on
  # one is still explained by the field map, so it may not be credited.
  stems <- data.frame(E = 100, N = 100, height = 20, crown_class = "dominant")
  outside <- data.frame(E = 121, N = 100, height = 12, crown_class = "codominant")
  det <- data.frame(x = c(100.5, 119), y = c(100, 100), z = c(19, 11))
  args <- list(det = det, tol_xy = 4, core_cx = 100, core_cy = 100, core_half = 20)
  core_only <- do.call(fp_points, c(list(stems = stems), args))
  expect_equal(nrow(core_only), 1L)                    # the x=119 detection
  expect_true(core_only$credit_eligible)               # blind to the outside stem
  with_ring <- do.call(fp_points, c(list(stems = stems), args,
                                    list(elig_stems = rbind(stems, outside))))
  expect_true(with_ring$isolated)                      # score_plot split unchanged
  expect_false(with_ring$credit_eligible)              # ... but never creditable
})

test_that("credit_isolated honours credit_eligible on both target and witness sides", {
  fps <- data.frame(x = c(0, 20), y = 0, z = 10, isolated = TRUE,
                    credit_eligible = c(TRUE, FALSE))
  wit <- data.frame(x = c(0.5, 0.5, 20.2, 20.4), y = 0,
                    fam = c("pc", "rgb", "pc", "rgb"), isolated = TRUE,
                    credit_eligible = c(TRUE, TRUE, TRUE, TRUE),
                    stringsAsFactors = FALSE)
  # target FP at x=20 is ineligible even though witnesses sit on it
  expect_equal(credit_isolated(fps, wit, "chm", r = 2, min_fam = 2), 1L)
  # ineligible witnesses never testify: strike the rgb witness at the eligible FP
  wit2 <- wit; wit2$credit_eligible[2] <- FALSE
  expect_equal(credit_isolated(fps, wit2, "chm", r = 2, min_fam = 2), 0L)
})

## ---- param-pinned cache reads (#96 review) -----------------------------------
test_that("read_selection carries the chm_vwf parameter columns when present", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  write.csv(data.frame(site = "SOAP", method = c("chm_vwf", "ams3d"),
                       rung = c("4", "2"), chm_res = c(0.5, NA),
                       vwf_a = c(0.05, NA)), tmp, row.names = FALSE)
  sel <- read_selection(tmp, "SOAP")
  expect_true(all(c("chm_res", "vwf_a") %in% names(sel)))
  expect_equal(sel$chm_res[sel$method == "chm_vwf"], 0.5)
})

test_that("read_arm_cache prefers the exact param-pinned variant over the glob", {
  tmp <- file.path(tempdir(), "btc_pin"); dir.create(tmp, showWarnings = FALSE)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  write.csv(data.frame(x = 1, y = 2, z = 1),
            file.path(tmp, "chm_vwf__SOAP__SOAP_001__4__res0.25__a0.1.csv"),
            row.names = FALSE)
  write.csv(data.frame(x = 1, y = 2, z = 99),
            file.path(tmp, "chm_vwf__SOAP__SOAP_001__4__res0.5__a0.05.csv"),
            row.names = FALSE)
  det <- read_arm_cache(tmp, "chm_vwf", "SOAP", "SOAP_001", "4",
                        params = c("res0.5", "a0.05"))
  expect_equal(det$z, 99)                    # pinned, not glob-first (res0.25)
  # pinned variant absent -> warning fallback to the glob
  expect_warning(
    d2 <- read_arm_cache(tmp, "chm_vwf", "SOAP", "SOAP_001", "4",
                         params = c("res9", "a9")),
    "variant")
  expect_equal(nrow(d2), 1L)
})

## ---- credit_isolated --------------------------------------------------------
test_that("credit_isolated credits only isolated FPs, never from the target's own family", {
  fps <- data.frame(x = c(0, 10), y = c(0, 0), z = c(12, 9),
                    isolated = c(TRUE, FALSE))          # the near FP never credits
  wit <- data.frame(x = c(0.4, 0.6, 10.2), y = c(0, 0, 0),
                    fam = c("pc", "rgb", "pc"), isolated = TRUE,
                    stringsAsFactors = FALSE)
  # target family "chm": pc + rgb witnesses at the isolated FP -> 1 credit; the
  # near FP at x=10 is not eligible even though a witness sits on it.
  expect_equal(credit_isolated(fps, wit, target_fam = "chm",
                               r = 2, min_fam = 2), 1L)
  # target family "pc": its own family's witnesses are struck, leaving only rgb
  # (1 family) -> nothing credits.
  expect_equal(credit_isolated(fps, wit, target_fam = "pc",
                               r = 2, min_fam = 2), 0L)
  # non-isolated witnesses (over-seg of a mapped tree) never testify
  wit2 <- within(wit, isolated <- c(TRUE, FALSE, TRUE))
  expect_equal(credit_isolated(fps, wit2, target_fam = "chm",
                               r = 2, min_fam = 2), 0L)
  # empty inputs
  expect_equal(credit_isolated(fps[0, ], wit, "chm"), 0L)
  expect_equal(credit_isolated(fps, wit[0, ], "chm"), 0L)
})

test_that("credit_isolated counts one credit per co-detected tree, not per duplicate FP", {
  # Two target FPs 1 m apart sit over the same witness cluster: crediting both
  # would launder over-segmentation of one unmapped tree. Exactly one credit;
  # the far FP at x=20 has its own witnesses and credits separately.
  fps <- data.frame(x = c(0, 1, 20), y = 0, z = 10, isolated = TRUE)
  wit <- data.frame(x = c(0.5, 0.5, 20.2, 20.4), y = 0,
                    fam = c("pc", "rgb", "pc", "rgb"), isolated = TRUE,
                    stringsAsFactors = FALSE)
  expect_equal(credit_isolated(fps, wit, "chm", r = 2, min_fam = 2), 2L)
})

## ---- read_selection ---------------------------------------------------------
test_that("read_selection maps each arm to its selected best rung for a site", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  write.csv(data.frame(site = c("SOAP", "SOAP", "SJER"),
                       method = c("ams3d", "chm_vwf", "ams3d"),
                       rung = c("2", "4", "1")), tmp, row.names = FALSE)
  sel <- read_selection(tmp, "SOAP")
  expect_equal(nrow(sel), 2L)
  expect_equal(sel$rung[sel$method == "ams3d"], "2")
  expect_null(read_selection(file.path(tempdir(), "nonexistent.csv"), "SOAP"))
  expect_null(read_selection(tmp, "TEAK"))
})

test_that("read_arm_cache warns when multiple parameter variants match one rung", {
  tmp <- file.path(tempdir(), "btc_multi"); dir.create(tmp, showWarnings = FALSE)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  det <- data.frame(x = 1, y = 2, z = 3)
  write.csv(det, file.path(tmp, "chm_vwf__SOAP__SOAP_001__4__res0.5__a0.05.csv"),
            row.names = FALSE)
  write.csv(det, file.path(tmp, "chm_vwf__SOAP__SOAP_001__4__res0.25__a0.1.csv"),
            row.names = FALSE)
  expect_warning(r <- read_arm_cache(tmp, "chm_vwf", "SOAP", "SOAP_001", "4"),
                 "variant")
  expect_equal(nrow(r), 1L)
})

## ---- cached_rungs ----------------------------------------------------------
test_that("cached_rungs discovers which rungs an arm cached for a plot", {
  tmp <- file.path(tempdir(), "btc_rungs"); dir.create(tmp, showWarnings = FALSE)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  det <- data.frame(x = 1, y = 2, z = 3)
  write.csv(det, file.path(tmp, "ams3d__SOAP__SOAP_001__2.csv"), row.names = FALSE)
  write.csv(det, file.path(tmp, "chm_vwf__SOAP__SOAP_001__4__res0.5__a0.05.csv"),
            row.names = FALSE)
  write.csv(det, file.path(tmp, "chm_vwf__SOAP__SOAP_001__native__res0.25__a0.1.csv"),
            row.names = FALSE)
  expect_equal(cached_rungs(tmp, "ams3d", "SOAP", "SOAP_001")$rung, "2")
  expect_setequal(cached_rungs(tmp, "chm_vwf", "SOAP", "SOAP_001")$rung,
                  c("4", "native"))
  expect_equal(nrow(cached_rungs(tmp, "ptrees", "SOAP", "SOAP_001")), 0L)
  # plot anchoring: SOAP_001 must not swallow SOAP_0011
  write.csv(det, file.path(tmp, "ams3d__SOAP__SOAP_0011__8.csv"), row.names = FALSE)
  expect_equal(cached_rungs(tmp, "ams3d", "SOAP", "SOAP_001")$rung, "2")
})

## ---- boxes_to_dets ---------------------------------------------------------
test_that("boxes_to_dets reads apex z from the CHM and floors off-CHM boxes", {
  skip_if_not_installed("terra")
  chm <- terra::rast(xmin = 0, xmax = 10, ymin = 0, ymax = 10,
                     resolution = 1, vals = 7.5)
  boxes <- data.frame(x = c(2.5, 50), y = c(2.5, 50), score = c(0.9, 0.8))
  det <- boxes_to_dets(boxes, chm)
  expect_equal(det$z, c(7.5, 2.0))                  # on-CHM height; off-CHM floor
  expect_equal(det$x, boxes$x)
  # score filter drops low-confidence boxes
  det2 <- boxes_to_dets(boxes, chm, score_min = 0.85)
  expect_equal(nrow(det2), 1L)
})

test_that("boxes_to_dets handles empty box sets", {
  skip_if_not_installed("terra")
  chm <- terra::rast(xmin = 0, xmax = 10, ymin = 0, ymax = 10,
                     resolution = 1, vals = 7.5)
  det <- boxes_to_dets(data.frame(x = numeric(), y = numeric(),
                                  score = numeric()), chm)
  expect_equal(nrow(det), 0L)
  expect_true(all(c("x", "y", "z") %in% names(det)))
})
