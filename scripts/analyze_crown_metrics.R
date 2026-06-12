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

# Crown-diameter comparison (#20, extended #30/#34): union the five CHM
# segmenters from #7 (crown_metrics_results.csv), the 3-D point-instance
# segmenters from #30 (crown_metrics_3d_results.csv: li2012, ptrees, ams3d), the
# SOAP-only TreeisoNet treeOff arm (treeisonet_crown_metrics.csv), and the #34
# deep instance-segmentation arms SegmentAnyTree + ForestFormer3D
# (segmentanytree_crown_metrics.csv, forestformer3d_crown_metrics.csv), then
# score pooled crown-diameter RMSE/MAE/bias/R2 per algorithm. All sources share
# the canonical crown-metrics columns and a per-algo `algo` label, so the union
# is a straight rbind; pooling is by SUMMED squared errors over matched trees,
# never a mean of per-plot rates.
#
# Multi-site: pass SITES=SJER,SOAP,TEAK to pool each algorithm across all the
# requested sites (the 3-D arms run on all three; the TreeisoNet arm is SOAP-only
# and contributes only when SOAP is requested). The legacy single-site form
# SITE=SOAP is preserved -- it pools that one site only, reproducing the original
# SOAP-only behaviour. Sites whose metric CSVs are absent are skipped.
#
# Usage:
#   Rscript scripts/analyze_crown_metrics.R SITES=SJER,SOAP,TEAK
#   Rscript scripts/analyze_crown_metrics.R SITE=SOAP        # legacy single-site
# Output: prints + $CLAUDE_JOB_DIR/neon/crown_compare_tables.md (multi-site) or
#   $CLAUDE_JOB_DIR/neon/<SITE>/crown_compare_tables.md (single legacy site).
args <- strsplit(commandArgs(TRUE), "=")
A    <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
# SITES wins when given; else fall back to the legacy single SITE; else default
# to the cross-site set. The single-SITE path keeps the original output layout.
if (!is.null(A$SITES)) {
  SITES <- strsplit(A$SITES, ",")[[1]]; single <- FALSE
} else if (!is.null(A$SITE)) {
  SITES <- A$SITE; single <- TRUE
} else {
  SITES <- c("SJER", "SOAP", "TEAK"); single <- FALSE
}
d <- .job_dir()

# Mirrors crown_metrics_sweep.R::err_stats (pooled SSE, never a mean of per-plot
# rates). det = detected diameter, fld = field diameter.
err_stats <- function(det, fld) {
  ok <- is.finite(det) & is.finite(fld); det <- det[ok]; fld <- fld[ok]
  n <- length(det)
  if (n < 2) return(data.frame(n = n, rmse = NA, mae = NA, bias = NA, r2 = NA))
  e <- det - fld; ss_res <- sum((fld - det)^2); ss_tot <- sum((fld - mean(fld))^2)
  data.frame(n = n, rmse = sqrt(mean(e^2)), mae = mean(abs(e)), bias = mean(e),
             r2 = if (ss_tot > 0) 1 - ss_res / ss_tot else NA_real_)
}

cols <- c("site", "plot", "algo", "crown_class", "individualID",
          "d_eq", "d_caliper", "area", "field_maxCD", "field_ninetyCD")
read_metric_csv <- function(path, required = cols, optional = FALSE) {
  if (!file.exists(path)) {
    if (optional) return(NULL)
    stop(sprintf("missing required metric file: %s", path))
  }
  x <- tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(x)) {
    if (optional) return(NULL)
    stop(sprintf("could not read required metric file: %s", path))
  }
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    msg <- sprintf("%s missing columns: %s", path, paste(missing, collapse = ", "))
    if (optional) { warning(msg); return(NULL) }
    stop(msg)
  }
  x[, required, drop = FALSE]
}

# Union every crown-metric source available for one site: the #7 CHM arms
# (required when single-site for back-compat; optional when pooling several
# sites so a missing site is skipped, not fatal), the #30 3-D arms, and the
# SOAP-only TreeisoNet arm. Each is filtered to its own site so a stale cross-
# site row cannot leak in. Returns a canonical-cols frame (0 rows if none).
collect_site <- function(site, chm_required) {
  nd <- file.path(d, "neon", site)
  parts <- list()
  chm <- read_metric_csv(file.path(nd, "crown_metrics_results.csv"),
                         optional = !chm_required)
  if (!is.null(chm)) parts[[length(parts) + 1]] <- chm[chm$site == site, cols,
                                                        drop = FALSE]
  td <- read_metric_csv(file.path(nd, "crown_metrics_3d_results.csv"),
                        optional = TRUE)
  if (!is.null(td)) parts[[length(parts) + 1]] <- td[td$site == site, cols,
                                                      drop = FALSE]
  ti <- read_metric_csv(file.path(nd, "treeisonet_crown_metrics.csv"),
                        optional = TRUE)
  if (!is.null(ti)) parts[[length(parts) + 1]] <- ti[ti$site == site, cols,
                                                      drop = FALSE]
  # #34 deep-model arms (SegmentAnyTree, ForestFormer3D): one optional per-model
  # CSV each, written by crown_metrics_deepmodel.R. Absent until the GPU instance
  # clouds exist for the site, so optional like the 3-D / TreeisoNet arms.
  for (dm in c("segmentanytree_crown_metrics.csv",
               "forestformer3d_crown_metrics.csv")) {
    x <- read_metric_csv(file.path(nd, dm), optional = TRUE)
    if (!is.null(x)) parts[[length(parts) + 1]] <- x[x$site == site, cols,
                                                     drop = FALSE]
  }
  if (!length(parts)) return(NULL)
  do.call(rbind, parts)
}

# Single-site keeps the original behaviour (CHM file is required); multi-site
# pooling treats each site's CHM file as optional so the run does not abort on a
# site that has not been crown-scored yet.
res <- do.call(rbind, Filter(Negate(is.null),
  lapply(SITES, function(s) collect_site(s, chm_required = single))))
if (is.null(res) || !nrow(res))
  stop(sprintf("no crown metric rows for SITES=%s", paste(SITES, collapse = ",")))
algos    <- unique(res$algo)
site_lbl <- paste(SITES, collapse = "+")
cat(sprintf("[%s] crown-diameter comparison: %d matched-tree rows, %d algos\n",
            site_lbl, nrow(res), length(algos)))

md <- c("<!-- generated by scripts/analyze_crown_metrics.R; do not edit by hand -->",
        "", sprintf("#### Crown-diameter accuracy on %s (pooled matched trees)",
                    site_lbl), "")
for (defn in list(c("d_eq", "field_ninetyCD", "Equivalent-circle d_eq vs ninetyCrownDiameter"),
                  c("d_caliper", "field_maxCD", "Max-caliper d_caliper vs maxCrownDiameter"))) {
  cat(sprintf("\n--- %s ---\n", defn[3]))
  md <- c(md, sprintf("**%s**", defn[3]), "",
          "| algo | n | rmse | mae | bias | r2 |",
          "| --- | --- | --- | --- | --- | --- |")
  tab <- do.call(rbind, lapply(algos, function(a)
    cbind(algo = a, err_stats(res[[defn[1]]][res$algo == a],
                              res[[defn[2]]][res$algo == a]))))
  tab <- tab[order(tab$rmse), ]
  for (i in seq_len(nrow(tab))) {
    cat(sprintf("  %-22s n=%d rmse=%.2f mae=%.2f bias=%+.2f r2=%+.3f\n",
                tab$algo[i], tab$n[i], tab$rmse[i], tab$mae[i], tab$bias[i], tab$r2[i]))
    md <- c(md, sprintf("| %s | %d | %.2f | %.2f | %+.2f | %+.3f |",
                tab$algo[i], tab$n[i], tab$rmse[i], tab$mae[i], tab$bias[i], tab$r2[i]))
  }
  md <- c(md, "")
}
# Single legacy site writes per-site (unchanged); multi-site writes one pooled
# table under neon/.
out <- if (single) file.path(d, "neon", SITES, "crown_compare_tables.md") else
  file.path(d, "neon", "crown_compare_tables.md")
writeLines(md, out)
cat(sprintf("\n-> %s\n", out))
