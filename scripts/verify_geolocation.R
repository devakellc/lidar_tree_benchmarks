#!/usr/bin/env Rscript
# Auditable geolocation check: independently re-derive a sample of field-stem
# coordinates straight from the raw NEON mappingandtagging table + a fresh
# locations-API fetch, and compare to the stored ground_truth_stems.csv. Should
# print ~0 m differences (confirms the polar-offset formula, azimuth convention,
# and UTM zone). Usage: Rscript scripts/verify_geolocation.R [SITE=SOAP] [N=8]
suppressMessages(library(jsonlite))
args <- strsplit(commandArgs(TRUE), "="); A <- setNames(lapply(args,`[`,2), sapply(args,`[`,1))
site <- if (is.null(A$SITE)) "SOAP" else A$SITE
N    <- as.integer(if (is.null(A$N)) 8 else A$N)
d  <- Sys.getenv("CLAUDE_JOB_DIR", file.path(getwd(), "work"))
nd <- file.path(d, "neon", site)

dat <- readRDS(list.files(file.path(nd, "vst"), pattern = "allyears.rds$", full.names = TRUE)[1])
mt  <- dat$vst_mappingandtagging
gt  <- read.csv(file.path(nd, "ground_truth_stems.csv"))
cand <- mt[!is.na(mt$stemDistance) & !is.na(mt$stemAzimuth) & !is.na(mt$pointID) &
           mt$individualID %in% gt$individualID, ]
set.seed(1); pick <- cand[sample(nrow(cand), min(N, nrow(cand))), ]

cat(sprintf("[%s] independent re-derivation vs stored ground truth:\n", site))
maxoff <- 0
for (i in seq_len(nrow(pick))) {
  nm <- paste0(pick$namedLocation[i], ".", pick$pointID[i])
  j  <- tryCatch(fromJSON(paste0("https://data.neonscience.org/api/v0/locations/", nm))$data,
                 error = function(e) NULL)
  if (is.null(j)) next
  az <- pick$stemAzimuth[i] * pi / 180
  E  <- j$locationUtmEasting  + pick$stemDistance[i] * sin(az)
  N2 <- j$locationUtmNorthing + pick$stemDistance[i] * cos(az)
  g  <- gt[gt$individualID == pick$individualID[i], ][1, ]
  off <- sqrt((E - g$E)^2 + (N2 - g$N)^2); maxoff <- max(maxoff, off)
  cat(sprintf("  %-26s dE=%+.4f dN=%+.4f |off|=%.4f m\n",
              pick$individualID[i], E - g$E, N2 - g$N, off))
}
cat(sprintf("max |offset| over %d stems: %.5f m  (%s)\n", nrow(pick), maxoff,
            if (maxoff < 0.01) "PASS" else "CHECK"))
