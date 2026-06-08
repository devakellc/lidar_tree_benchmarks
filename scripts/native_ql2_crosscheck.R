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

# Native USGS 3DEP QL2 cross-check for the density-ladder sweep (issue #4).
#
# The density-ladder sweep DECIMATES a dense NEON cloud to simulate sparse
# acquisitions ("decimation-as-simulation", a stated caveat). This script pulls
# the *native* USGS 3DEP cloud over the same NEON plots, runs the SAME CHM-VWF
# detection pipeline (detect_lasr + score_plot from sweep_lib), and compares it
# to the cached DECIMATED 2-pts/m^2 rung by crown class. Deliverable: does
# decimation faithfully predict native-sparse (QL2 ~2 pulses/m^2) behaviour?
#
# HONESTY ON COVERAGE. The entwine/USGS public catalog does NOT have a QL2
# acquisition over every NEON site. Each site's covering EPT is verified by its
# *measured* first-return (pulse) density over a probe AOI, not by its name:
#   * "native_ql2" -- native pulse density falls in the QL2 band ~[1.5, 4]:
#     this is the real native-vs-decimated equivalence test. `native_full` is
#     the native-QL2 detection.
#   * "dense"       -- native pulse density is well above QL2 (a high-density
#     survey, e.g. CA_SierraNevada_*_B22 at ~30-40 pulses/m^2). There is NO
#     public native QL2 here; we report that honestly. `native_full` is the
#     dense detection (context only) and `native_dec2` (the dense cloud
#     decimated to 2) is a SECONDARY cross-source diagnostic, NOT native QL2.
#
# Pipeline per plot:
#   centroid (UTM 11N) -> +/-PAD m box -> reproject box to EPT native SRS ->
#   PDAL readers.ept(bounds) -> drop noise (class 7/18) + withheld points ->
#   filters.reprojection out_srs=OUTCRS -> writers.las -> per-plot laz.
#   Then Step-0 density -> detect_lasr(res-by-density, a=A) -> score_plot
#   (stems reprojected 32611 -> OUTCRS so they share the detection CRS).
# Pool per site exactly like the sweep: recall = sum(TP)/sum(n_ref); per-class
# TP recovered as round(rec_class * n_class).
#
# Usage:
#   export CLAUDE_JOB_DIR=.../work
#   Rscript scripts/native_ql2_crosscheck.R SITES=SOAP,SJER,TEAK \
#       [EPT_SOAP=<url>] [EPT_SJER=<url>] [EPT_TEAK=<url>] [OUTCRS=EPSG:32611] \
#       [TOL=4] [A=0.10] [PAD=45] [MAXPLOTS=ALL] [PLOTS=SJER_050,...] \
#       [QL2_LO=1.5] [QL2_HI=4]
#
# If EPT_<SITE> is not given, the script auto-discovers via ept_discovery.R
# (writes work/neon/<SITE>/ql2/ept_candidates.csv) and uses the preferred one.
suppressMessages({
  library(lidR); library(lasR); library(terra); library(sf); library(jsonlite)
})
options(lidR.progress = FALSE)

src_dir <- dirname(sub("^--file=", "",
                       grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(src_dir) || !nzchar(src_dir)) src_dir <- "scripts"
source(file.path(src_dir, "sweep_lib.R"))     # detect_lasr, score_plot, plot_half
source(file.path(src_dir, "ept_discovery.R")) # discover_ept (for auto-resolution)

JOB <- .job_dir()
NEON_EPSG <- 32611                            # SOAP/SJER/TEAK = UTM 11N

## ---- args ----------------------------------------------------------------
args <- strsplit(commandArgs(TRUE), "=")
A <- setNames(lapply(args, function(z) paste(z[-1], collapse = "=")),
              sapply(args, `[`, 1))
SITES <- if (is.null(A$SITES)) c("SOAP", "SJER", "TEAK") else
  strsplit(A$SITES, ",")[[1]]
OUTCRS <- if (is.null(A$OUTCRS)) "EPSG:32611" else A$OUTCRS
OUT_EPSG <- as.integer(sub("EPSG:", "", OUTCRS))
TOL  <- as.numeric(if (is.null(A$TOL)) 4.0 else A$TOL)
AVWF <- as.numeric(if (is.null(A$A)) 0.10 else A$A)
PAD  <- as.numeric(if (is.null(A$PAD)) 45 else A$PAD)
QL2_LO <- as.numeric(if (is.null(A$QL2_LO)) 1.5 else A$QL2_LO)
QL2_HI <- as.numeric(if (is.null(A$QL2_HI)) 4.0 else A$QL2_HI)
MAXPLOTS <- if (is.null(A$MAXPLOTS) || A$MAXPLOTS == "ALL") NA_integer_ else
  as.integer(A$MAXPLOTS)
ONLY_PLOTS <- if (is.null(A$PLOTS)) NULL else strsplit(A$PLOTS, ",")[[1]]

# Per-site EPT override argument: EPT_<SITE>=<url>.
ept_override <- function(site) {
  k <- paste0("EPT_", site)
  if (!is.null(A[[k]]) && nzchar(A[[k]])) return(A[[k]]) else NA_character_
}

# Resolve the EPT URL for a site: explicit override, else cached candidate, else
# run discovery now. Returns the preferred EPT URL or NA.
resolve_ept <- function(site) {
  ov <- ept_override(site)
  if (!is.na(ov)) return(ov)
  cand_path <- file.path(JOB, "neon", site, "ql2", "ept_candidates.csv")
  if (!file.exists(cand_path)) {
    discover_ept(site, JOB, site_epsg = NEON_EPSG)   # writes candidates csv
  }
  if (!file.exists(cand_path)) return(NA_character_)
  cand <- read.csv(cand_path, stringsAsFactors = FALSE)
  pick <- cand$url[which(cand$preferred %in% TRUE)]
  if (length(pick) && !is.na(pick[1]) && nzchar(pick[1])) return(pick[1])
  ok <- cand$url[!is.na(cand$url) & nzchar(cand$url)]
  if (length(ok)) ok[1] else NA_character_
}

# Read the native horizontal SRS (EPSG code) from an ept.json.
ept_srs <- function(url) {
  j <- tryCatch(jsonlite::fromJSON(url), error = function(e) NULL)
  if (is.null(j)) return(NA_integer_)
  s <- j$srs
  if (!is.null(s$horizontal) && nzchar(as.character(s$horizontal)))
    return(as.integer(s$horizontal))
  NA_integer_
}

# Provenance manifest for a cached per-plot laz. The cache key is NOT just the
# plot id: a laz pulled under a different EPT/PAD/OUTCRS is a different cloud and
# must be re-pulled. We write a sidecar `<laz>.json` recording (ept_url, pad,
# outcrs) and only reuse the cached laz when all three match. This prevents a
# stale laz (e.g. from a prior EPT or a different reprojection target) from being
# silently reused.
manifest_path <- function(out_laz) paste0(out_laz, ".json")

write_manifest <- function(out_laz, ept_url, pad, outcrs) {
  writeLines(jsonlite::toJSON(
    list(ept_url = ept_url, pad = pad, outcrs = outcrs),
    auto_unbox = TRUE, pretty = TRUE), manifest_path(out_laz))
}

# TRUE iff a cached laz exists AND its sidecar manifest matches the provenance
# we are about to use (ept_url + pad + outcrs). Missing/mismatched manifest -> not
# a valid cache hit -> caller re-pulls.
cache_valid <- function(out_laz, ept_url, pad, outcrs) {
  if (!file.exists(out_laz) || file.info(out_laz)$size < 500) return(FALSE)
  mp <- manifest_path(out_laz)
  if (!file.exists(mp)) return(FALSE)
  m <- tryCatch(jsonlite::fromJSON(mp), error = function(e) NULL)
  if (is.null(m)) return(FALSE)
  identical(as.character(m$ept_url), as.character(ept_url)) &&
    isTRUE(as.numeric(m$pad) == as.numeric(pad)) &&
    identical(as.character(m$outcrs), as.character(outcrs))
}

# Reproject a set of NEON-native points (E,N in EPSG:32611, the woody-veg stem
# CRS) into OUTCRS for scoring (fix 2). Detections come out of detect_lasr in
# OUTCRS (the laz is reprojected by filters.reprojection out_srs=OUTCRS), so the
# field stems MUST be moved into the SAME CRS before matching -- otherwise any
# non-default OUTCRS silently mismatches stems against detections. When OUTCRS is
# the native EPSG:32611 this is an identity transform. Returns a 2-col matrix.
to_outcrs <- function(e, n) {
  if (OUT_EPSG == NEON_EPSG) return(cbind(x = e, y = n))   # identity fast-path
  p <- st_as_sf(data.frame(x = e, y = n), coords = c("x", "y"), crs = NEON_EPSG)
  xy <- st_coordinates(st_transform(p, OUT_EPSG))
  cbind(x = xy[, 1], y = xy[, 2])
}

# Pull a single plot's native cloud via PDAL readers.ept, reproject to OUTCRS.
# box (UTM) -> EPT SRS bounds -> ept reader -> drop noise -> reproject -> laz.
pull_plot <- function(cx, cy, pad, ept_url, ept_epsg, out_laz) {
  box <- st_as_sf(data.frame(id = 1:4,
                             x = c(cx - pad, cx + pad, cx + pad, cx - pad),
                             y = c(cy - pad, cy - pad, cy + pad, cy + pad)),
                  coords = c("x", "y"), crs = NEON_EPSG)
  box_ept <- st_transform(box, ept_epsg)
  bb <- st_bbox(box_ept)
  bounds <- sprintf("([%f, %f], [%f, %f])",
                    bb["xmin"], bb["xmax"], bb["ymin"], bb["ymax"])
  # Drop noise classes (7 low noise, 18 high noise) AND any point flagged
  # Withheld (fix 3): a genuine withheld exclusion so the code matches the
  # documented "drop class 7/18 + withheld". `Withheld![1:1]` keeps points whose
  # Withheld flag is NOT 1 (i.e. removes withheld). Multiple range conditions in
  # one limits string are ANDed by filters.range.
  pipe <- list(
    list(type = "readers.ept", filename = ept_url, bounds = bounds),
    list(type = "filters.range",                       # drop noise + withheld
         limits = paste0("Classification![7:7],",
                         "Classification![18:18],Withheld![1:1]")),
    list(type = "filters.reprojection", out_srs = OUTCRS),
    list(type = "writers.las", filename = out_laz, compression = "laszip")
  )
  pj <- tempfile(fileext = ".json")
  writeLines(jsonlite::toJSON(pipe, auto_unbox = TRUE, pretty = TRUE), pj)
  rc <- suppressWarnings(system2("pdal", c("pipeline", shQuote(pj)),
                                 stdout = TRUE, stderr = TRUE))
  st <- attr(rc, "status")
  unlink(pj)
  if (!is.null(st) && st != 0)
    return(list(ok = FALSE, msg = paste(tail(rc, 4), collapse = " | ")))
  if (!file.exists(out_laz) || file.info(out_laz)$size < 500)
    return(list(ok = FALSE, msg = "empty/no output laz"))
  write_manifest(out_laz, ept_url, pad, OUTCRS)   # record provenance for cache
  list(ok = TRUE, msg = "")
}

# Run detect_lasr on a (possibly decimated) LAS; returns det + measured density.
# res_override: when non-NA, pin the CHM resolution instead of deriving it from
# the (decimated) first-return density. Used by the decimated-2 arm so its CHM
# resolution MATCHES the NEON dec2 arm it is compared against (chm_res=0.5).
# Without this, a native cloud decimated to ~2 ppsm would fall into the res=1.0
# branch (frdens<4) while NEON dec2 is pooled at res=0.5 -- a resolution confound
# that contaminates the "sensor difference" delta. With it, the residual delta is
# genuinely sensor/return-structure, not CHM resolution.
detect_on <- function(src_las, rung, pad, res_override = NA_real_) {
  l <- src_las
  area0 <- (2 * pad)^2
  pdens_native <- npoints(src_las) / area0
  if (!is.na(rung)) {
    if (pdens_native <= rung) return(NULL)             # never upsample
    l <- decimate_points(l, homogenize(density = rung, res = 5))
  }
  if (sum(l$Classification == 2L) < 10) return(NULL)
  nrm <- tryCatch(normalize_height(l, tin(), na.rm = TRUE),
                  error = function(e) NULL)
  if (is.null(nrm)) return(NULL)
  nrm <- filter_poi(nrm, Z >= -1, Z < 80)
  a2   <- (2 * pad)^2
  pdn  <- npoints(nrm) / a2
  frdn <- sum(nrm$ReturnNumber == 1L) / a2
  fn   <- tempfile(fileext = ".laz"); writeLAS(nrm, fn)
  res  <- if (!is.na(res_override)) res_override else
          if (frdn >= 8) 0.25 else if (frdn >= 4) 0.5 else 1.0
  # NB the <8 ppsm pre-LM smoothing branch in detect_lasr keys on dens=frdn, not
  # on res, so pinning res does not change which smoothing path the dec2 arm
  # takes -- only the CHM grid is matched to the NEON dec2 arm.
  det  <- tryCatch(detect_lasr(fn, res = res, a = AVWF, dens = frdn),
                   error = function(e) { message(conditionMessage(e)); NULL })
  unlink(fn)
  if (is.null(det)) return(NULL)
  list(det = det, res = res, pdens = pdn, frdens = frdn)
}

## ---- per-site run ---------------------------------------------------------
run_site <- function(site) {
  ept_url <- resolve_ept(site)
  cat(sprintf("\n########## %s  EPT=%s\n", site, ept_url))
  if (is.na(ept_url)) { cat("  no EPT URL resolved; skipping\n"); return(NULL) }
  ept_epsg <- ept_srs(ept_url)
  if (is.na(ept_epsg)) { cat("  could not read SRS from ept.json; skipping\n")
    return(NULL) }
  cat(sprintf("  native EPT SRS = EPSG:%d ; reproject -> %s\n", ept_epsg, OUTCRS))

  nd  <- file.path(JOB, "neon", site)
  pc  <- read.csv(file.path(nd, "plot_centroids.csv"))
  gt  <- read.csv(file.path(nd, "ground_truth_stems.csv"))
  gt  <- gt[gt$live & gt$is_tree & !is.na(gt$E), ]
  MINTREES <- 6
  ql2 <- file.path(nd, "ql2")
  dir.create(ql2, showWarnings = FALSE, recursive = TRUE)

  rows <- list()
  pids <- pc$plotID
  if (!is.null(ONLY_PLOTS)) pids <- intersect(pids, ONLY_PLOTS)
  done <- 0L
  for (pid in pids) {
    if (!is.na(MAXPLOTS) && done >= MAXPLOTS) break
    ci <- pc[pc$plotID == pid, ][1, ]
    cx <- ci$easting; cy <- ci$northing
    ph <- plot_half(ci$plotType)                       # tower +/-20, dist +/-10
    stems <- gt[gt$plotID == pid &
                abs(gt$E - cx) <= ph & abs(gt$N - cy) <= ph, ]
    if (nrow(stems) < MINTREES) next

    out_laz <- file.path(ql2, sprintf("%s.laz", pid))
    # Cache reuse is provenance-gated (fix 1): a pre-existing laz is a cache HIT
    # ONLY when its sidecar manifest exists AND matches (ept_url, pad, outcrs).
    # Provenance is NEVER inferred from file size -- a manifest-less laz (or one
    # whose manifest records a different EPT/PAD/OUTCRS) is treated as a cache
    # MISS and re-pulled, which writes a fresh manifest. This closes the stale-
    # data hole: a laz from an earlier EPT/pad/reprojection target can never be
    # silently reused just because it is non-empty on disk.
    if (!cache_valid(out_laz, ept_url, PAD, OUTCRS)) {
      if (file.exists(out_laz))
        cat(sprintf("  %s: cache miss (no/invalid manifest) -> re-pull\n", pid))
      pr <- pull_plot(cx, cy, PAD, ept_url, ept_epsg, out_laz)
      if (!pr$ok) { cat(sprintf("  %s: PULL FAILED (%s)\n", pid, pr$msg)); next }
    } else {
      cat(sprintf("  %s: cache HIT (provenance match)\n", pid))
    }
    las <- tryCatch(readLAS(out_laz), error = function(e) NULL)
    if (is.null(las) || is.empty(las) || npoints(las) < 200) {
      cat(sprintf("  %s: too few native points; skip\n", pid)); next }
    if (sum(las$Classification == 2L) < 10) {
      cat(sprintf("  %s: <10 ground pts; skip\n", pid)); next }

    # native densities over the padded box (all-return + first-return/pulse)
    area     <- (2 * PAD)^2
    pdens_n  <- npoints(las) / area
    frdens_n <- sum(las$ReturnNumber == 1L) / area

    # Move stems + core centre from native EPSG:32611 into OUTCRS so they share
    # the detection CRS (fix 2). Stem-core SELECTION above stayed native (matches
    # the AOI pull); SCORING below is done entirely in OUTCRS.
    stems_o <- stems
    s_xy <- to_outcrs(stems$E, stems$N)
    stems_o$E <- s_xy[, "x"]; stems_o$N <- s_xy[, "y"]
    c_xy <- to_outcrs(cx, cy)
    cx_o <- c_xy[1, "x"]; cy_o <- c_xy[1, "y"]

    # native_full = detection at the survey's native density (its true sparsity);
    #   CHM resolution derived from its (high) native first-return density.
    # native_dec2 = native cloud decimated to ALL-RETURN density 2 (pdens, the
    #   sweep's homogenize unit), CHM resolution PINNED to 0.5 to match the NEON
    #   dec2 arm (round-1 fix: resolution-matched equivalence). Each arm carries a
    #   (rung, res_override) pair.
    variants <- list(native_full = list(rung = NA_real_, res = NA_real_),
                     native_dec2 = list(rung = 2,        res = 0.5))
    got_any <- FALSE
    for (vn in names(variants)) {
      v <- detect_on(las, variants[[vn]]$rung, PAD,
                     res_override = variants[[vn]]$res)
      if (is.null(v)) { cat(sprintf("  %s [%s]: skip\n", pid, vn)); next }
      sc <- score_plot(stems_o, v$det, tol_xy = TOL, core_cx = cx_o,
                       core_cy = cy_o, core_half = ph)
      sc <- cbind(data.frame(site = site, plot = pid, plotType = ci$plotType,
                             variant = vn,
                             native_pdens = round(pdens_n, 2),
                             native_frdens = round(frdens_n, 2),
                             pdens = round(v$pdens, 2), frdens = round(v$frdens, 2),
                             chm_res = v$res, vwf_a = AVWF), sc)
      rows[[length(rows) + 1]] <- sc
      got_any <- TRUE
      cat(sprintf("  %s [%s]: native_fr=%.2f -> fr=%.2f res=%.2f n_ref=%d TP=%d recall=%.2f prec=%.2f\n",
                  pid, vn, frdens_n, v$frdens, v$res, sc$n_ref, sc$TP,
                  sc$recall, sc$precision))
    }
    if (got_any) done <- done + 1L
  }
  if (!length(rows)) { cat(sprintf("  %s: NO plots scored\n", site)); return(NULL) }
  res_df <- do.call(rbind, rows)
  res_df$ept_url <- ept_url
  res_df$ept_srs <- ept_epsg
  # classify the site by MEDIAN native first-return density at the probe AOIs
  med_fr <- median(unique(res_df[, c("plot", "native_frdens")])$native_frdens)
  cls <- if (med_fr >= QL2_LO && med_fr <= QL2_HI) "native_ql2" else
    if (med_fr > QL2_HI) "dense" else "sub_ql2"
  res_df$site_class <- cls
  res_df$native_frdens_median <- round(med_fr, 2)
  write.csv(res_df, file.path(ql2, "ql2_detect_results.csv"), row.names = FALSE)
  cat(sprintf("  --> %s native median first-return density = %.2f pulses/m^2 (class=%s)\n",
              site, med_fr, cls))
  res_df
}

## ---- pooling (matches the sweep's count-summing rule) ---------------------
# Pooling rules (matched to analyze_sweep.R):
#  * recall    = sum(TP) / sum(n_ref)  -- TP is the global match count (incl.
#                boundary matches), n_ref is all core stems. This is correct.
#  * precision = sum(tp_core) / sum(n_det). The PER-ROW precision from score_plot
#                is already tp_core/n_det (core-only numerator). Pooling
#                sum(TP)/sum(n_det) would mix a boundary-inclusive numerator with
#                a core-only denominator and INFLATE precision (round-1 fix). So
#                we recover tp_core = round(precision * n_det) per row -- exactly as
#                analyze_sweep.R:17 does -- and pool that.
pool_native <- function(df) {
  classes <- c("dominant", "codominant", "intermediate", "suppressed")
  hbands  <- c("short", "mid", "tall")
  tp_core <- round(df$precision * df$n_det); tp_core[is.na(tp_core)] <- 0
  out <- list()
  out[["overall"]] <- data.frame(
    scope = "overall",
    n_ref = sum(df$n_ref), TP = sum(df$TP),
    recall = sum(df$TP) / sum(df$n_ref),
    n_det = sum(df$n_det),
    precision = if (sum(df$n_det)) sum(tp_core) / sum(df$n_det) else NA_real_)
  for (cls in classes) {
    nc <- df[[paste0("n_", cls)]]; rc <- df[[paste0("rec_", cls)]]
    tp <- round(rc * nc); tp[is.na(tp)] <- 0
    nref <- sum(nc, na.rm = TRUE)
    out[[cls]] <- data.frame(scope = cls, n_ref = nref, TP = sum(tp),
      recall = if (nref) sum(tp) / nref else NA_real_,
      n_det = NA_integer_, precision = NA_real_)
  }
  for (b in hbands) {
    nc <- df[[paste0("n_h_", b)]]; rc <- df[[paste0("rec_h_", b)]]
    if (is.null(nc)) next
    tp <- round(rc * nc); tp[is.na(tp)] <- 0
    nref <- sum(nc, na.rm = TRUE)
    out[[paste0("h_", b)]] <- data.frame(scope = paste0("h_", b),
      n_ref = nref, TP = sum(tp),
      recall = if (nref) sum(tp) / nref else NA_real_,
      n_det = NA_integer_, precision = NA_real_)
  }
  do.call(rbind, out)
}

# Pool the cached DECIMATED 2-rung the same way (rung==2,res==0.5,a==0.10),
# restricted to a plot set so the comparison is paired.
pool_decimated <- function(site, plots) {
  f <- file.path(JOB, "neon", site, "sweep_results.csv")
  if (!file.exists(f)) return(NULL)
  d <- read.csv(f, stringsAsFactors = FALSE)
  d <- d[d$rung == "2" & d$chm_res == 0.5 & d$vwf_a == 0.10, ]
  if (!is.null(plots)) d <- d[d$plot %in% plots, ]
  if (!nrow(d)) return(NULL)
  pool_native(d)
}

## ---- main -----------------------------------------------------------------
all_native <- list(); all_cmp <- list()
for (site in SITES) {
  nv <- tryCatch(run_site(site), error = function(e) {
    cat("  ", site, " errored: ", conditionMessage(e), "\n"); NULL })
  if (is.null(nv) || !nrow(nv)) next
  all_native[[site]] <- nv
  scls <- nv$site_class[1]

  # The equivalence test depends on the site class:
  #  native_ql2 -> compare native_full (== native QL2) vs neon_dec2
  #  dense      -> there is NO native QL2; compare the dense survey decimated to
  #                2 (native_dec2) vs neon_dec2 as a cross-source diagnostic.
  test_variant <- if (scls == "native_ql2") "native_full" else "native_dec2"
  test_plots <- unique(nv$plot[nv$variant == test_variant])

  pools <- list()
  for (vn in unique(nv$variant)) {
    sub <- nv[nv$variant == vn, ]
    pv  <- pool_native(sub); pv$engine <- vn
    pv$n_plots <- length(unique(sub$plot))
    pools[[vn]] <- pv
  }
  pd <- pool_decimated(site, test_plots)
  if (!is.null(pd)) { pd$engine <- "neon_dec2"; pd$n_plots <- length(test_plots) }

  site_tbl <- do.call(rbind, c(pools, list(neon_dec2 = pd)))
  site_tbl$site <- site
  site_tbl$site_class <- scls
  site_tbl$native_frdens_median <- nv$native_frdens_median[1]
  site_tbl$test_variant <- test_variant
  all_cmp[[site]] <- site_tbl

  cat(sprintf("\n=== %s pooled  (class=%s, native median fr=%.2f, n=%d plots) ===\n",
              site, scls, nv$native_frdens_median[1], length(unique(nv$plot))))
  print(site_tbl[, c("engine", "scope", "n_plots", "n_ref", "TP", "recall",
                     "precision")], row.names = FALSE)

  if (!is.null(pd) && test_variant %in% names(pools)) {
    nt <- pools[[test_variant]]
    eqv <- merge(nt[, c("scope", "recall", "n_ref", "TP")],
                 pd[, c("scope", "recall", "n_ref", "TP")],
                 by = "scope", suffixes = c("_native", "_neonDec2"))
    eqv$delta_recall <- eqv$recall_native - eqv$recall_neonDec2
    lab <- if (scls == "native_ql2") "native_QL2 vs neon_dec2"
           else "native_dec2(cross-source) vs neon_dec2"
    cat(sprintf("\n--- %s EQUIVALENCE: %s ---\n", site, lab))
    print(eqv, row.names = FALSE, digits = 3)
  }
}

if (length(all_cmp)) {
  comb <- do.call(rbind, all_cmp)
  outf <- file.path(JOB, "neon", "native_ql2_vs_decimated.csv")
  write.csv(comb, outf, row.names = FALSE)
  cat("\nWROTE", outf, "\n")
}
cat("\nDONE native_ql2_crosscheck\n")
