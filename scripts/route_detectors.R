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

# #P2 routing study (GitHub issue #65).
#
# The repo's only routing knob today is measured density (the frdens/pdens guards
# in sweep_lib), and the SJER->SOAP->TEAK structure gradient is described only
# descriptively. But the ladder already shows the crossovers a router would key
# on: SegmentAnyTree beats CHM-VWF at native/8/4 then falls below it at 2/1
# pts/m2; multichm is the flattest arm and the strongest low-density baseline.
# This study quantifies the achievable meta-pipeline by SELECTING a detector per
# cell from cheap deploy-time features, instead of shipping one fixed arm.
#
# Pipeline:
#   1. Assemble the per-cell x per-arm scored ladder from the per-arm *_results.csv
#      (CHM-VWF + multichm at the canonical density-derived chm_res / vwf_a=0.10;
#      SegmentAnyTree, native Li 2012, native+8 ForestFormer3D as-scored), keyed by
#      (site, plot, rung). Restrict to the equal set where the three laddered arms
#      (chm_vwf, multichm, segmentanytree) are all present.
#   2. Compute deploy-time structure features over each cell's clip_normalized.laz
#      with lidR: rumple_index (canopy rugosity), canopy cover, height CV, gap
#      fraction, mean canopy height -- plus measured frdens/pdens. NO field data.
#   3. Label each cell with its argmax-F1 arm (oracle_pick) and fit an
#      interpretable rpart CART, evaluated leave-one-PLOT-out (a plot's rungs are
#      not independent).
#   4. Score four policies by the canonical pool() (sum counts, never average
#      rates): every single arm, the fixed best-single arm, the ORACLE per-cell
#      argmax, and the LEARNED rpart router (held-out predictions). Report pooled
#      F1 + understory recall, and emit the routing table the workflow ships with.
#
# Usage:
#   Rscript scripts/route_detectors.R SITES=SOAP,SJER,TEAK
# Reads: work/neon/<SITE>/{sweep_results,multichm_sweep_results,
#   segmentanytree_results,li2012_results,forestformer3d_results}.csv and the
#   frozen normalized clips. Writes: work/neon/<SITE>/router_policy.csv
#   (per-cell features + oracle label + held-out router prediction) and prints the
#   fitted CART + the policy comparison behind results/detector-routing-results.md.
suppressMessages({ library(lidR); library(data.table); library(rpart) })
options(lidR.progress = FALSE, lidR.verbose = FALSE)
d <- .job_dir()
source(.find("sweep_lib.R"))
source(.find("model_bench_lib.R"))   # pool, equal_set_guard
source(.find("route_lib.R"))

args  <- strsplit(commandArgs(TRUE), "=")
A     <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
SITES <- if (!is.null(A$SITES)) strsplit(A$SITES, ",")[[1]] else
  if (!is.null(A$SITE)) A$SITE else c("SOAP", "SJER", "TEAK")
VWF_A   <- as.numeric(if (is.null(A$VWF_A)) 0.10 else A$VWF_A)
CORE_ARMS <- c("chm_vwf", "multichm", "segmentanytree")  # the equal-set ladder
FEATURES  <- c("frdens", "pdens", "rumple", "cover", "height_cv", "gap", "mean_ht")
# canonical CHM resolution rule (matches the benchmark): fine when dense.
canon_res <- function(frdens) ifelse(frdens >= 8, 0.25, 0.5)

POOL_NUM <- c("n_ref", "n_det", "TP", "recall", "precision",
              "rec_dominant", "n_dominant", "rec_codominant", "n_codominant",
              "rec_intermediate", "n_intermediate", "rec_suppressed", "n_suppressed")

## ---- per-arm laddered rows -> normalized (site,plot,rung,arm,F1,counts) ----
norm_cols <- function(df, site, arm) {
  df$site <- if ("site" %in% names(df)) df$site else site
  df$arm  <- arm
  for (c in c("F1", POOL_NUM, "frdens", "pdens"))
    if (is.null(df[[c]])) df[[c]] <- NA_real_
  df$rung <- as.character(df$rung)
  df[, c("site", "plot", "rung", "arm", "F1", "frdens", "pdens", POOL_NUM),
     drop = FALSE]
}
load_arm <- function(nd, site, arm) {
  f <- switch(arm,
    chm_vwf        = "sweep_results.csv",
    multichm       = "multichm_sweep_results.csv",
    segmentanytree = "segmentanytree_results.csv",
    li2012         = "li2012_results.csv",
    forestformer3d = "forestformer3d_results.csv")
  p <- file.path(nd, f); if (!file.exists(p)) return(NULL)
  df <- read.csv(p, stringsAsFactors = FALSE)
  if (arm == "chm_vwf")
    df <- df[df$vwf_a == VWF_A & df$chm_res == canon_res(df$frdens), , drop = FALSE]
  if (arm == "multichm" && "chm_res" %in% names(df))
    df <- df[df$chm_res == canon_res(df$frdens), , drop = FALSE]
  if (!nrow(df)) return(NULL)
  norm_cols(df, site, arm)
}

## ---- deploy-time structure features over one normalized clip --------------
clip_features <- function(clip) {
  las <- tryCatch(suppressWarnings(lidR::readLAS(clip)), error = function(e) NULL)
  if (is.null(las) || lidR::is.empty(las)) return(NULL)
  z <- las$Z; cz <- z[z > 2]                          # canopy returns
  if (length(cz) < 20) return(NULL)
  chm <- tryCatch(suppressWarnings(lidR::rasterize_canopy(las, res = 0.5,
                    algorithm = lidR::p2r())), error = function(e) NULL)
  rump <- if (!is.null(chm)) tryCatch(lidR::rumple_index(chm), error = function(e) NA_real_) else NA_real_
  cover <- mean(z > 2)                                # fraction of returns in canopy
  data.frame(rumple = as.numeric(rump), cover = cover,
             height_cv = stats::sd(cz) / mean(cz), gap = 1 - cover,
             mean_ht = mean(cz))
}

## ---- per-site assembly: long ladder + per-cell features -------------------
build_site <- function(site) {
  nd <- file.path(d, "neon", site)
  long <- rbindlist(Filter(Negate(is.null),
    lapply(c(CORE_ARMS, "li2012", "forestformer3d"), function(a) load_arm(nd, site, a))),
    fill = TRUE)
  if (!nrow(long)) return(NULL)
  long <- as.data.frame(long)
  # equal set: cells where all three laddered arms are present
  long$detector <- long$arm
  guarded <- equal_set_guard(long, arms = CORE_ARMS)
  if (!nrow(guarded)) return(NULL)
  cells <- unique(guarded[, CELL_KEYS])
  # deploy-time features per surviving cell (frdens/pdens carried on the rows)
  meta <- unique(guarded[guarded$arm == "chm_vwf",
                         c(CELL_KEYS, "frdens", "pdens")])
  feats <- rbindlist(lapply(seq_len(nrow(cells)), function(i) {
    cl <- cells[i, ]
    clip <- file.path(nd, "frozen", site, cl$plot, cl$rung, "clip_normalized.laz")
    ff <- if (file.exists(clip)) clip_features(clip) else NULL
    if (is.null(ff)) return(NULL)
    cbind(cl, ff)
  }), fill = TRUE)
  if (!nrow(feats)) return(NULL)
  feats <- merge(as.data.frame(feats), meta, by = CELL_KEYS, all.x = TRUE)
  list(long = guarded, feats = feats)
}

## ---- driver ---------------------------------------------------------------
run_main <- function() {
  t0 <- Sys.time()
  longs <- list(); featl <- list()
  for (site in SITES) {
    b <- tryCatch(build_site(site), error = function(e) {
      message("site ", site, ": ", conditionMessage(e)); NULL })
    if (!is.null(b)) { longs[[site]] <- b$long; featl[[site]] <- b$feats }
  }
  long <- do.call(rbind, longs); feats <- do.call(rbind, featl)
  if (is.null(long) || !nrow(long) || is.null(feats) || !nrow(feats)) {
    cat("No routable cells assembled.\n"); return(invisible()) }
  # keep only cells that have BOTH a scored ladder and computed features
  fk <- do.call(paste, c(feats[CELL_KEYS], sep = "\r"))
  lk <- do.call(paste, c(long[CELL_KEYS], sep = "\r"))
  long <- long[lk %in% fk, , drop = FALSE]
  # label: per-cell argmax-F1 arm, joined onto the feature table
  lab <- oracle_pick(long, value = "F1")
  feats <- merge(feats, lab, by = CELL_KEYS)
  feats$arm <- factor(feats$arm)
  feats <- feats[stats::complete.cases(feats[, FEATURES]), , drop = FALSE]
  cat(sprintf("Routable cells: %d across %d sites; candidate arms: %s\n",
              nrow(feats), length(longs), paste(levels(feats$arm), collapse = ",")))
  cat("\nOracle label distribution (per-cell best arm):\n"); print(table(feats$arm))

  ## leave-one-PLOT-out CV: predict each plot's cells from a CART on the others
  feats$plotkey <- paste(feats$site, feats$plot, sep = "::")
  fold <- lopo_folds(feats$plotkey)
  feats$pred <- NA_character_
  form <- stats::as.formula(paste("arm ~", paste(FEATURES, collapse = " + ")))
  for (fo in sort(unique(fold))) {
    tr <- feats[fold != fo, , drop = FALSE]; te <- which(fold == fo)
    tr$arm <- droplevels(tr$arm)
    if (nlevels(tr$arm) < 2) { feats$pred[te] <- as.character(tr$arm[1]); next }
    fit <- rpart::rpart(form, data = tr, method = "class",
                        control = rpart::rpart.control(minbucket = 8, cp = 0.01))
    feats$pred[te] <- as.character(predict(fit, feats[te, , drop = FALSE], type = "class"))
  }

  ## score the policies with pool() over the SAME routable cells
  pol_pool <- function(picks) {
    rows <- select_policy_rows(long, picks)
    if (!nrow(rows)) return(NULL)
    rows$rung <- as.character(rows$rung); pool(rows)
  }
  cand <- levels(feats$arm)
  single <- setNames(lapply(cand, function(a)
    pol_pool(data.frame(feats[CELL_KEYS], arm = a, stringsAsFactors = FALSE))), cand)
  best_arm <- cand[which.max(vapply(cand, function(a)
    { p <- single[[a]]; if (is.null(p)) -1 else p$F1 }, numeric(1)))]
  oracle  <- pol_pool(data.frame(feats[CELL_KEYS], arm = as.character(feats$arm),
                                 stringsAsFactors = FALSE))
  learned <- pol_pool(data.frame(feats[CELL_KEYS], arm = feats$pred,
                                 stringsAsFactors = FALSE))

  cat("\n===== ROUTING POLICIES (pooled, equal-set ladder) =====\n")
  cat(sprintf("%-22s %6s %7s %7s %9s\n", "policy", "n_ref", "F1", "recall", "rec_und"))
  row <- function(nm, p) if (!is.null(p))
    cat(sprintf("%-22s %6d %7.3f %7.3f %9.3f\n", nm, p$n_ref, p$F1, p$recall, p$rec_understory))
  for (a in cand) row(paste0("single:", a), single[[a]])
  row(paste0("FIXED-BEST (", best_arm, ")"), single[[best_arm]])
  row("LEARNED router (LOPO)", learned)
  row("ORACLE (per-cell best)", oracle)
  if (!is.null(learned) && !is.null(single[[best_arm]]))
    cat(sprintf("\nlearned vs fixed-best dF1 = %+.3f ; oracle headroom dF1 = %+.3f\n",
                learned$F1 - single[[best_arm]]$F1, oracle$F1 - single[[best_arm]]$F1))

  ## the deployable CART (fit on all data) + its rules
  feats$arm <- droplevels(feats$arm)
  full <- rpart::rpart(form, data = feats, method = "class",
                       control = rpart::rpart.control(minbucket = 8, cp = 0.01))
  cat("\n===== LEARNED ROUTING TABLE (CART on all cells) =====\n")
  print(full)
  for (site in names(longs))
    write.csv(feats[feats$site == site, ],
              file.path(d, "neon", site, "router_policy.csv"), row.names = FALSE)
  cat(sprintf("\nDONE: %d cells in %.1f min; wrote router_policy.csv per site\n",
              nrow(feats), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

if (sys.nframe() == 0L) run_main()
