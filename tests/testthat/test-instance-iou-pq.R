# Unit tests for the #V1 point-set IoU / Coverage / Panoptic-Quality scorer
# helpers in model_bench_lib.R. All fixtures are synthetic integer label vectors
# over a shared point substrate (NA = unassigned) and small coordinate sets, so
# no real NEON data is ever required. Every expected value is hand-computed in
# the comments.
source(file.path("..", "..", "scripts", "model_bench_lib.R"), local = TRUE)

## ---- point_set_iou -------------------------------------------------------
# 10-point substrate. pred id1 = pts 1-4 (size 4), id2 = pts 5-6 (size 2).
# ref  id1 = pts 1-3 (size 3), id2 = pts 5-7 (size 3). Overlaps:
#   (pred1,ref1): pts 1-3 -> inter 3, union 4+3-3=4, iou 0.75
#   (pred2,ref2): pts 5-6 -> inter 2, union 2+3-2=3, iou 2/3
test_that("point_set_iou computes pairwise intersection/union/iou over the substrate", {
  pred <- c(1L, 1L, 1L, 1L, 2L, 2L, NA, NA, NA, NA)
  ref  <- c(1L, 1L, 1L, NA, 2L, 2L, 2L, NA, NA, NA)
  io <- point_set_iou(pred, ref)
  expect_identical(names(io),
                   c("pred_id", "ref_id", "inter", "pred_size", "ref_size",
                     "union", "iou"))
  expect_equal(nrow(io), 2L)                       # only the two overlapping pairs
  r11 <- io[io$pred_id == 1 & io$ref_id == 1, ]
  expect_equal(r11$inter, 3); expect_equal(r11$union, 4); expect_equal(r11$iou, 0.75)
  r22 <- io[io$pred_id == 2 & io$ref_id == 2, ]
  expect_equal(r22$inter, 2); expect_equal(r22$iou, 2 / 3)
})

test_that("point_set_iou is 0-row safe on all-unassigned and length-0 input", {
  z <- point_set_iou(c(NA, NA), c(NA, NA))
  expect_equal(nrow(z), 0L)
  expect_identical(names(z),
                   c("pred_id", "ref_id", "inter", "pred_size", "ref_size",
                     "union", "iou"))
  expect_equal(nrow(point_set_iou(integer(0), integer(0))), 0L)
})

test_that("point_set_iou rejects mismatched lengths", {
  expect_error(point_set_iou(c(1L, 1L), c(1L)))
})

## ---- iou_match (greedy max-IoU, IoU>=0.5 gate) ----------------------------
# Three candidate pairs; pred2 overlaps BOTH ref1 (iou 0.9) and ref2 (iou 0.6).
# Greedy best-first must take (pred2,ref1,0.9) first, blocking (pred2,ref2);
# (pred1,ref2,0.55) survives; (pred3,ref3,0.4) is below the gate -> dropped.
test_that("iou_match is greedy best-first, 1:1, and gated at IoU>=0.5", {
  io <- data.frame(
    pred_id = c(2L, 2L, 1L, 3L),
    ref_id  = c(1L, 2L, 2L, 3L),
    iou     = c(0.9, 0.6, 0.55, 0.40))
  m <- iou_match(io, gate = 0.5)
  expect_identical(names(m), c("pred_id", "ref_id", "iou"))
  expect_equal(nrow(m), 2L)
  expect_true(any(m$pred_id == 2 & m$ref_id == 1))   # best pair won
  expect_false(any(m$pred_id == 2 & m$ref_id == 2))  # pred2 already used
  expect_true(any(m$pred_id == 1 & m$ref_id == 2))   # ref2 free for pred1
  expect_false(any(m$ref_id == 3))                   # 0.40 below the 0.5 gate
})

test_that("iou_match is 0-row safe", {
  m <- iou_match(point_set_iou(c(NA_integer_), c(NA_integer_)))
  expect_equal(nrow(m), 0L)
})

## ---- coverage_table (mean max-IoU over references) -----------------------
# pred1 = pts 1-2, pred2 = pts 3-4. ref1 = pts 1-2 (iou 1.0 vs pred1),
# ref2 = pt 3 (pred2 covers 1 of 2 pred-pts -> inter1/union2 = 0.5),
# ref3 = pt 5 (pred NA there -> no overlap -> 0). mean = (1.0+0.5+0)/3.
test_that("coverage_table reports per-reference max IoU including zero-overlap refs", {
  pred <- c(1L, 1L, 2L, 2L, NA)
  ref  <- c(1L, 1L, 2L, NA, 3L)
  cv <- coverage_table(pred, ref)
  expect_identical(names(cv), c("ref_id", "max_iou"))
  expect_setequal(cv$ref_id, c(1, 2, 3))
  expect_equal(cv$max_iou[cv$ref_id == 1], 1.0)
  expect_equal(cv$max_iou[cv$ref_id == 2], 0.5)
  expect_equal(cv$max_iou[cv$ref_id == 3], 0.0)     # zero-overlap ref still listed
  expect_equal(mean(cv$max_iou), 0.5)
})

test_that("coverage_table is 0-row safe when there are no references", {
  expect_equal(nrow(coverage_table(c(1L, 1L), c(NA, NA))), 0L)
})

## ---- panoptic_quality (PQ = SQ*RQ) ---------------------------------------
# 2 matched pairs iou 0.75 & 2/3 (sum 17/12). n_pred=2, n_ref=3 -> 1 ref unmatched.
# TP=2 FP=0 FN=1. SQ = (17/12)/2. RQ = 2/(2+0.5*0+0.5*1)=0.8. PQ=SQ*RQ.
# RQ is exactly the F1 over matched pairs.
test_that("panoptic_quality returns summable accumulators and SQ/RQ/PQ", {
  matched <- data.frame(pred_id = c(1L, 2L), ref_id = c(1L, 2L),
                        iou = c(0.75, 2 / 3))
  pq <- panoptic_quality(matched, n_pred = 2L, n_ref = 3L)
  expect_identical(names(pq),
                   c("TP", "FP", "FN", "sum_iou", "n_pred", "n_ref",
                     "SQ", "RQ", "PQ", "precision", "recall", "F1"))
  expect_equal(pq$TP, 2); expect_equal(pq$FP, 0); expect_equal(pq$FN, 1)
  expect_equal(pq$sum_iou, 0.75 + 2 / 3)
  expect_equal(pq$SQ, (0.75 + 2 / 3) / 2)
  expect_equal(pq$RQ, 0.8)
  expect_equal(pq$PQ, pq$SQ * pq$RQ)
  expect_equal(pq$precision, 1.0)
  expect_equal(pq$recall, 2 / 3)
  expect_equal(pq$F1, pq$RQ)                         # RQ == F1 over matched pairs
})

test_that("panoptic_quality is well-defined with zero true positives", {
  pq <- panoptic_quality(
    data.frame(pred_id = integer(), ref_id = integer(), iou = numeric()),
    n_pred = 3L, n_ref = 4L)
  expect_equal(pq$TP, 0); expect_equal(pq$FP, 3); expect_equal(pq$FN, 4)
  expect_equal(pq$sum_iou, 0)
  expect_true(is.na(pq$SQ))                          # mean over no matches
  expect_equal(pq$RQ, 0); expect_equal(pq$PQ, 0)
})

## ---- assign_points_to_stems (Voronoi-on-stems within per-stem radius) -----
# Stem A (0,0) r=2, Stem B (10,0) r=1. A point is assigned to the NEAREST stem
# only if it falls within THAT stem's radius, else NA (proxy crown extent).
test_that("assign_points_to_stems assigns nearest stem gated by its crown radius", {
  px <- c(0.5, 9.5, 5.0, 10.8, 11.5)
  py <- c(0.0, 0.0, 0.0,  0.0,  0.0)
  a <- assign_points_to_stems(px, py, sx = c(0, 10), sy = c(0, 0),
                              radius = c(2, 1))
  expect_equal(a, c(1L, 2L, NA, 2L, NA))
})

test_that("assign_points_to_stems treats a non-finite crown radius as no capture", {
  # A stem with NA/Inf radius (no measured crown width and no fallback) must
  # capture no points -- never error on an NA subscript.
  a <- assign_points_to_stems(px = c(0.1, 10.1), py = c(0, 0),
                              sx = c(0, 10), sy = c(0, 0),
                              radius = c(NA, 1))
  expect_equal(a, c(NA_integer_, 2L))
})

test_that("assign_points_to_stems handles no stems and no points", {
  expect_equal(assign_points_to_stems(c(1, 2), c(1, 2), numeric(0), numeric(0), numeric(0)),
               c(NA_integer_, NA_integer_))
  expect_equal(length(assign_points_to_stems(numeric(0), numeric(0), 0, 0, 2)), 0L)
})

## ---- transfer_labels (grid-hash nearest-neighbour, dependency-free) -------
# Source labelled points (0,0)=id5, (10,10)=id7. Each destination point takes the
# label of the nearest source point within tol (XY); beyond tol -> NA.
test_that("transfer_labels assigns nearest source label within tolerance", {
  out <- transfer_labels(sx = c(0, 10), sy = c(0, 10), sid = c(5L, 7L),
                         dx = c(0.1, 10.0, 5.0, 0.4, 0.6),
                         dy = c(0.0, 10.0, 5.0, 0.0, 0.0),
                         tol = 0.5)
  expect_equal(out, c(5L, 7L, NA, 5L, NA))
})

test_that("transfer_labels is empty-safe on both sides", {
  expect_equal(transfer_labels(numeric(0), numeric(0), integer(0),
                              c(1, 2), c(1, 2), tol = 0.5),
               c(NA_integer_, NA_integer_))
  expect_equal(length(transfer_labels(0, 0, 9L, numeric(0), numeric(0), 0.5)), 0L)
})

## ---- score_instance_cell (one plot x rung x model accumulator row) --------
# 10-pt substrate. pred1=pts1-4, pred2=pts5-6 (n_pred 2). ref1=pts1-3,
# ref2=pts5-7, ref3=pts9-10 (n_ref 3). IoU(1,1)=0.75, IoU(2,2)=2/3; ref3 has no
# overlap. classes: ref1 & ref3 dominant, ref2 suppressed.
#   TP=2 FP=0 FN=1; sum_iou=0.75+2/3; sum_maxiou=same (ref3 contributes 0).
#   dominant: n=2 tp=1(ref1) fn=1(ref3) sumiou=0.75 sumcov=0.75
#   suppressed: n=1 tp=1(ref2) fn=0 sumiou=2/3 sumcov=2/3
test_that("score_instance_cell builds summable, crown-class-stratified accumulators", {
  pred <- c(1L, 1L, 1L, 1L, 2L, 2L, NA, NA, NA, NA)
  ref  <- c(1L, 1L, 1L, NA, 2L, 2L, 2L, NA, 3L, 3L)
  rcls <- c("1" = "dominant", "2" = "suppressed", "3" = "dominant")
  row <- score_instance_cell(pred, ref, ref_class = rcls, gate = 0.5)
  expect_equal(row$n_pred, 2); expect_equal(row$n_ref, 3)
  expect_equal(row$TP, 2); expect_equal(row$FP, 0); expect_equal(row$FN, 1)
  expect_equal(row$sum_iou, 0.75 + 2 / 3)
  expect_equal(row$sum_maxiou, 0.75 + 2 / 3)         # ref3 max-iou 0
  expect_equal(row$coverage, (0.75 + 2 / 3) / 3)
  expect_equal(row$PQ, row$SQ * row$RQ)
  # per-class accumulators
  expect_equal(row$n_dominant, 2);   expect_equal(row$tp_dominant, 1)
  expect_equal(row$fn_dominant, 1);  expect_equal(row$sumiou_dominant, 0.75)
  expect_equal(row$sumcov_dominant, 0.75)
  expect_equal(row$n_suppressed, 1); expect_equal(row$tp_suppressed, 1)
  expect_equal(row$sumiou_suppressed, 2 / 3)
  expect_equal(row$n_codominant, 0); expect_equal(row$tp_codominant, 0)
})

test_that("score_instance_cell is well-defined on an empty reference", {
  row <- score_instance_cell(c(1L, 1L), c(NA, NA), ref_class = character(0))
  expect_equal(row$n_ref, 0); expect_equal(row$TP, 0)
  expect_true(is.na(row$recall) || row$recall == 0 || is.nan(row$recall))
})

## ---- pool_pq (sum accumulators, recompute rates) -------------------------
# Pooling two identical cells (the fixture above) must SUM accumulators and
# recompute pooled rates from the sums (NEVER average the per-cell rates):
#   TP 4, FN 2, FP 0, sum_iou 2*(0.75+2/3); SQ = sum_iou/TP (unchanged 0.7083),
#   RQ = 4/(4+1)=0.8; recall = 4/6; coverage = sum_maxiou/6.
#   dominant pooled: n 4, tp 2 -> rec_dominant 0.5; SQ_dominant = 1.5/2 = 0.75.
test_that("pool_pq sums panoptic accumulators and recomputes rates", {
  pred <- c(1L, 1L, 1L, 1L, 2L, 2L, NA, NA, NA, NA)
  ref  <- c(1L, 1L, 1L, NA, 2L, 2L, 2L, NA, 3L, 3L)
  rcls <- c("1" = "dominant", "2" = "suppressed", "3" = "dominant")
  cell <- score_instance_cell(pred, ref, ref_class = rcls)
  cell$site <- "X"; cell$plot <- "P"; cell$rung <- "native"; cell$model <- "m"
  pooled <- pool_pq(rbind(cell, cell))
  expect_equal(pooled$TP, 4); expect_equal(pooled$FN, 2); expect_equal(pooled$FP, 0)
  expect_equal(pooled$sum_iou, 2 * (0.75 + 2 / 3))
  expect_equal(pooled$SQ, (0.75 + 2 / 3) / 2)        # SQ invariant to duplication
  expect_equal(pooled$RQ, 0.8)
  expect_equal(pooled$recall, 4 / 6)
  expect_equal(pooled$precision, 1.0)
  expect_equal(pooled$coverage, 2 * (0.75 + 2 / 3) / 6)
  expect_equal(pooled$rec_dominant, 0.5)             # 2 TP / 4 dominant refs
  expect_equal(pooled$SQ_dominant, 0.75)             # 1.5 matched-iou / 2 TP
})

## ---- #V6 generic instance-LAZ loader ----------------------------------------
test_that("make_laz_loader round-trips a persisted classical instance cloud", {
  src <- file.path("..", "..", "scripts")
  source(file.path(src, "io_bridge.R"), local = TRUE)
  source(file.path(src, "score_instances_iou.R"), local = TRUE)
  las <- synth_las_normalized()
  las@data$treeID <- rep(c(1L, 2L), each = 60)
  las@data$treeID[c(5, 65)] <- NA_integer_          # unassigned points
  f <- tempfile(fileext = ".laz")
  on.exit(unlink(f), add = TRUE)
  write_instances_laz(las, f, id_col = "treeID")
  loader <- make_laz_loader("treeID")
  pts <- loader(f)
  expect_identical(names(pts), c("X", "Y", "id"))
  expect_equal(nrow(pts), 118L)                     # unassigned dropped
  expect_setequal(unique(pts$id), c(1L, 2L))
  # wrong id field -> NULL (schema failure -> cell skipped)
  expect_null(make_laz_loader("nonesuch")(f))
})
