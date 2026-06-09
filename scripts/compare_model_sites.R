#!/usr/bin/env Rscript
# Script: cross-site structure-gradient for classical model-benchmark arms (#E11)
# Mirrors compare_sites.R but reads model-benchmark arm CSVs
# Usage: Rscript scripts/compare_model_sites.R [SITES=SJER,SOAP,TEAK] [ARMS=chm_vwf,lmfauto,multichm,ptrees,ams3d]
# Reads:  $CLAUDE_JOB_DIR/neon/<SITE>/{ams3d,lidrplugins}_results.csv
# Writes: $CLAUDE_JOB_DIR/neon/model_cross_site_summary.csv
#         $CLAUDE_JOB_DIR/neon/figs/model_structure_gradient.png  (all arms, per-site lines)
#         $CLAUDE_JOB_DIR/neon/figs/model_structure_gradient_<arm>.png  (per arm)
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

source(.find("analyze_model_benchmark.R"))

args  <- strsplit(commandArgs(TRUE), "=")
A     <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
SITES <- strsplit(if (is.null(A$SITES)) "SJER,SOAP,TEAK" else A$SITES, ",")[[1]]
ARMS  <- strsplit(if (is.null(A$ARMS))  "chm_vwf,lmfauto,multichm,ptrees,ams3d"
                  else A$ARMS, ",")[[1]]

run_main <- function() {
  d <- .job_dir()
  site_pooled <- list()

  for (st in SITES) {
    nd <- file.path(d, "neon", st)
    arm_files <- c(file.path(nd, "ams3d_results.csv"),
                   file.path(nd, "lidrplugins_results.csv"))
    dfs <- lapply(arm_files, function(f)
      if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE) else NULL)
    dfs <- Filter(Negate(is.null), dfs)
    if (!length(dfs)) { cat("skip", st, "(no arm CSVs)\n"); next }

    u <- harmonize_union(dfs)
    arms_here <- intersect(ARMS, unique(u$detector))
    if (!length(arms_here)) { cat("skip", st, "(no matching arms)\n"); next }

    ladder <- pool_arms(u, arms = arms_here)
    if (is.null(ladder) || !nrow(ladder)) { cat("skip", st, "(pool returned empty)\n"); next }

    site_pooled[[st]] <- cbind(site = st, ladder)
    cat(sprintf("\n=== %s (arms=%d; rungs=%d; recall %.2f-%.2f) ===\n",
                st, length(arms_here), length(unique(ladder$rung)),
                min(ladder$recall, na.rm = TRUE), max(ladder$recall, na.rm = TRUE)))
    print(ladder[, intersect(c("detector", "rung", "frdens", "recall",
                               "rec_understory"), names(ladder))],
          row.names = FALSE, digits = 2)
  }

  if (!length(site_pooled))
    stop("no sites produced results; run arm sweeps first", call. = FALSE)

  comb <- do.call(rbind, Filter(Negate(is.null), site_pooled))
  write.csv(comb, file.path(d, "neon", "model_cross_site_summary.csv"), row.names = FALSE)

  fig <- file.path(d, "neon", "figs")
  dir.create(fig, showWarnings = FALSE, recursive = TRUE)

  site_cols <- c(SJER = "#1b9e77", SOAP = "#d95f02", TEAK = "#7570b3")
  sites_present <- names(site_pooled)

  .draw_sites <- function(sub, main_label, cex = 1.0) {
    xlim <- range(sub$frdens, na.rm = TRUE)
    plot(NA, xlim = xlim, ylim = c(0, 1), log = "x",
         xlab = "first-return density (pulses/m^2)", ylab = "recall",
         main = main_label)
    for (st in sites_present) {
      s <- sub[sub$site == st, ]
      if (!nrow(s)) next
      o <- order(s$frdens)
      lines(s$frdens[o], s$recall[o], type = "o", pch = 19,
            col = site_cols[st], lwd = 2)
      uu <- s$rec_understory[o]
      uu[is.na(s$n_understory[o]) | s$n_understory[o] < 5] <- NA
      lines(s$frdens[o], uu, type = "o", pch = 1,
            col = site_cols[st], lwd = 2, lty = 2)
    }
    legend("topright", sites_present, col = site_cols[sites_present],
           lwd = 2, pch = 19, bty = "n", cex = cex)
  }

  ## per-arm PNGs
  arms_in_comb <- intersect(ARMS, unique(comb$detector))
  for (arm in arms_in_comb) {
    local({
      png(file.path(fig, paste0("model_structure_gradient_", arm, ".png")),
          1050, 760, res = 130); on.exit(dev.off())
      par(mar = c(4.2, 4.2, 2.5, 1))
      .draw_sites(comb[comb$detector == arm, ],
                  paste0(arm, ": overall (solid) vs understory (dashed) by site"))
    })
  }

  ## overall summary PNG: one panel per arm in a grid
  n_arms <- length(arms_in_comb)
  nr <- if (n_arms <= 3) 1L else if (n_arms <= 6) 2L else 3L
  nc <- ceiling(n_arms / nr)
  local({
    png(file.path(fig, "model_structure_gradient.png"),
        width = 500 * nc, height = 440 * nr, res = 130); on.exit(dev.off())
    par(mfrow = c(nr, nc), mar = c(4.0, 4.0, 2.2, 0.8))
    for (arm in arms_in_comb)
      .draw_sites(comb[comb$detector == arm, ], main_label = arm, cex = 0.75)
  })

  cat(sprintf("\ncombined -> %s\nfigs     -> %s\n",
              file.path(d, "neon", "model_cross_site_summary.csv"), fig))
  cat(sprintf("per-arm PNGs: %s\noverall PNG:  %s\n",
              paste(paste0("model_structure_gradient_", arms_in_comb, ".png"),
                    collapse = ", "),
              "model_structure_gradient.png"))
}

if (sys.nframe() == 0L) run_main()
