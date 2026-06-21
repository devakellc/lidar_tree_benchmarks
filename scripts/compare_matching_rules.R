#!/usr/bin/env Rscript
.bs_ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
.bs_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
bs <- Find(file.exists, c(
  if (!is.null(.bs_ofile) && length(.bs_ofile) && nzchar(.bs_ofile))
    file.path(dirname(.bs_ofile), "bootstrap.R"),
  if (length(.bs_file)) file.path(dirname(sub("^--file=", "", .bs_file[1])),
                                  "bootstrap.R"),
  file.path("scripts", "bootstrap.R"),
  file.path("..", "..", "scripts", "bootstrap.R"),
  file.path(getwd(), "scripts", "bootstrap.R")))
if (!length(bs)) stop("bootstrap.R not found", call. = FALSE)
source(bs[1]); rm(bs, .bs_ofile, .bs_file)

# #V2 arm-ranking sensitivity: distance vs IoU0.5/Coverage leaderboard flips (#69).
#
# If the router (#P2) selects 'best arm per rung', that decision is only
# trustworthy if the ranking is stable across reasonable metrics. Building on the
# #V1 scorer, this recomputes the SOAP leaderboard under three rules and flags
# where arms reorder, testing the hypothesis that apex-distance matching
# over-credits high-recall/low-precision point segmenters (ptrees, AMS3D -- split
# crowns) that point-set IoU/RQ should demote, while the true-mask arms
# (SegmentAnyTree, ForestFormer3D) hold or rise.
#
# HONEST DATA LIMIT: only the deep arms persist per-point instance labels, so the
# #V1 instance_iou_pq.csv scores ONLY segmentanytree + forestformer3d. The
# classical/point arms (chm_vwf, multichm, lmfauto, ptrees, ams3d, li2012,
# treeisonet) are marked `n/a (no masks)` for the IoU/Coverage columns -- exactly
# the equal_set_guard honesty the issue requires. The over-crediting hypothesis
# for ptrees/AMS3D therefore CANNOT be tested until their per-point labels are
# materialized; what IS testable is whether the two mask arms reorder between
# distance and IoU, and how far distance over-credits them (the #V1 collapse).
#
# Usage:  Rscript scripts/compare_matching_rules.R SITE=SOAP
# Reads: the per-arm distance leaderboard CSVs (the union analyze_model_benchmark
#   builds: {lidrplugins,ams3d,li2012,treeisonet,segmentanytree,forestformer3d}_
#   results.csv) + the #V1 instance_iou_pq.csv. Writes work/neon/<SITE>/
#   matching_rule_ranks.csv and prints the leaderboards + Kendall tau.
suppressMessages({ library(data.table) })
d <- .job_dir()
source(.find("model_bench_lib.R"))   # pool

args <- strsplit(commandArgs(TRUE), "=")
A    <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
SITE <- if (is.null(A$SITE)) "SOAP" else A$SITE
RUNG <- if (is.null(A$RUNG)) "native" else A$RUNG
nd   <- file.path(d, "neon", SITE)

## ---- distance leaderboard: pool each detector at the rung ------------------
DIST_FILES <- c("lidrplugins_results.csv", "ams3d_results.csv", "li2012_results.csv",
                "treeisonet_results.csv", "segmentanytree_results.csv",
                "forestformer3d_results.csv")
read_dist <- function() {
  rows <- list()
  for (f in DIST_FILES) {
    p <- file.path(nd, f); if (!file.exists(p)) next
    df <- read.csv(p, stringsAsFactors = FALSE)
    det_col <- intersect(c("detector", "algo"), names(df))[1]
    if (is.na(det_col)) next
    df$rung <- as.character(df$rung)
    df <- df[df$rung == RUNG, , drop = FALSE]; if (!nrow(df)) next
    for (a in unique(df[[det_col]])) {
      sub <- df[df[[det_col]] == a, , drop = FALSE]
      if (is.null(sub$site)) sub$site <- SITE
      p2 <- tryCatch(pool(sub), error = function(e) NULL)
      if (!is.null(p2)) rows[[length(rows) + 1]] <-
        data.frame(arm = a, dist_recall = p2$recall, dist_precision = p2$precision,
                   dist_F1 = p2$F1, dist_rec_understory = p2$rec_understory,
                   n_ref = p2$n_ref, stringsAsFactors = FALSE)
    }
  }
  d <- do.call(rbind, rows)
  d[!duplicated(d$arm), , drop = FALSE]              # one row per arm
}

## ---- IoU leaderboard from instance_iou_pq.csv (mask arms only) -------------
read_iou <- function() {
  p <- file.path(nd, "instance_iou_pq.csv")
  if (!file.exists(p)) return(NULL)
  iq <- as.data.table(read.csv(p, stringsAsFactors = FALSE))
  iq <- iq[as.character(rung) == RUNG]
  if (!nrow(iq)) return(NULL)
  # instance_iou_pq.csv carries the score_instance_cell accumulators (#V1):
  # n_ref, TP, sum_iou, sum_maxiou -- pooled here by SUM (the #V1 rule).
  iq[, .(iou_recall = sum(TP) / sum(n_ref),
         coverage   = sum(sum_maxiou) / sum(n_ref),
         sq         = sum(sum_iou) / sum(TP)), by = .(arm = model)]
}

run_main <- function() {
  dist <- read_dist(); iou <- read_iou()
  if (is.null(dist) || !nrow(dist)) { cat("no distance leaderboard\n"); return(invisible()) }
  lb <- merge(dist, if (is.null(iou)) data.frame(arm = character()) else as.data.frame(iou),
              by = "arm", all.x = TRUE)
  lb$has_mask <- !is.na(lb$iou_recall)
  lb <- lb[order(-lb$dist_F1), , drop = FALSE]
  lb$rank_dist <- rank(-lb$dist_F1, ties.method = "min")
  # IoU rank only among mask arms
  lb$rank_iou <- NA_integer_
  mk <- which(lb$has_mask)
  if (length(mk)) lb$rank_iou[mk] <- rank(-lb$iou_recall[mk], ties.method = "min")
  lb$rank_cov <- NA_integer_
  if (length(mk)) lb$rank_cov[mk] <- rank(-lb$coverage[mk], ties.method = "min")

  cat(sprintf("\n===== MATCHING-RULE LEADERBOARD (%s, %s) =====\n", SITE, RUNG))
  cat(sprintf("%-16s %8s %8s %8s | %9s %9s %7s\n",
              "arm", "dist_R", "dist_F1", "rank", "iou_R@.5", "Coverage", "rank"))
  for (i in seq_len(nrow(lb))) {
    r <- lb[i, ]
    iou_s <- if (r$has_mask) sprintf("%9.3f %9.3f %7d", r$iou_recall, r$coverage, r$rank_iou)
             else sprintf("%9s %9s %7s", "n/a", "n/a", "-")
    cat(sprintf("%-16s %8.3f %8.3f %8d | %s\n",
                r$arm, r$dist_recall, r$dist_F1, r$rank_dist, iou_s))
  }
  # rank correlation over the COMMON (mask) arms
  cm <- lb[lb$has_mask, , drop = FALSE]
  cat(sprintf("\nComparable (mask) arms: %d (%s)\n", nrow(cm), paste(cm$arm, collapse = ", ")))
  if (nrow(cm) >= 3) {
    tau <- suppressWarnings(cor(cm$dist_F1, cm$iou_recall, method = "kendall"))
    rho <- suppressWarnings(cor(cm$dist_F1, cm$iou_recall, method = "spearman"))
    cat(sprintf("Kendall tau-b (dist_F1 vs iou_recall) = %.3f ; Spearman rho = %.3f\n", tau, rho))
  } else {
    cat("Kendall tau / Spearman rho: UNDEFINED (need >=3 comparable arms; only the\n")
    cat("two deep arms persist masks). Reporting the pairwise order instead.\n")
    if (nrow(cm) == 2) {
      o_d <- cm$arm[order(-cm$dist_F1)]; o_i <- cm$arm[order(-cm$iou_recall)]
      cat(sprintf("  distance order: %s ; IoU order: %s -> %s\n",
                  paste(o_d, collapse = ">"), paste(o_i, collapse = ">"),
                  if (identical(o_d, o_i)) "NO FLIP" else "FLIP"))
      # how far distance over-credits each mask arm vs IoU
      for (i in seq_len(nrow(cm)))
        cat(sprintf("  %s: dist_recall %.3f vs iou_recall@0.5 %.3f (ratio %.1fx)\n",
                    cm$arm[i], cm$dist_recall[i], cm$iou_recall[i],
                    cm$dist_recall[i] / cm$iou_recall[i]))
    }
  }
  o <- file.path(nd, "matching_rule_ranks.csv")
  write.csv(lb, o, row.names = FALSE)
  cat(sprintf("\nwrote %s\n", o))
}

if (sys.nframe() == 0L) run_main()
